/*
solve_bpm.cu (CUDA)

BPM ソルバー (Douglas-Gunn ADI 法, bpm/FDBPMpropagator.cu の CUDA カーネルを使用)
sol/solve_bpm.cpp の CUDA 版
*/

#include "obpm.h"
#include "obpm_cuda.h"
#include "obpm_prototype.h"

#include "hdf5.h"
#define FILE_NAME "time_series_data.h5"
// 光活性化関数曲線 (powersweep 指定時に出力)
#define FN_activation "activation_curve.csv"

#include "bpm/bpm_prototype.h"
#include "bpm/wabpm.h"

static void cuda_check(cudaError_t code, const char *file, int line)
{
	if (code != cudaSuccess) {
		fprintf(stderr, "*** CUDA error : %s (%s:%d)\n", cudaGetErrorString(code), file, line);
		exit(1);
	}
}
#define CUDA_CHECK(ans) cuda_check((ans), __FILE__, __LINE__)

void solve_bpm(int io, double *tdft, FILE *fp)
{
	// HDF5ファイルの作成
	hid_t file_id;
	// local
	hid_t dataset_id, dataspace_id;
	herr_t status;

	char str[BUFSIZ];

	if (!GPU) {
		fprintf(stderr, "*** BPM solver (CUDA) requires GPU.\n");
		exit(1);
	}

	// セルの幅（空間ステップ）を計算 (Xn/Yn/Zn は節点座標で要素数は Nx+1/Ny+1/Nz+1)
	double Dx = (Xn[Nx] - Xn[0]) / Nx;
	double Dy = (Yn[Ny] - Yn[0]) / Ny;
	double Dz = (Zn[Nz] - Zn[0]) / Nz;
	sprintf(str, "cell step : %.6e %.6e %.6e", Dx, Dy, Dz);
	if (io) fprintf(fp, "%s\n", str);
	fprintf(stdout, "%s\n", str);

	// 不均一メッシュの警告 : BPM カーネルは等間隔グリッドを仮定するため、
	// 不均一の場合は平均セル幅で計算される
	{
		const char *axisname[] = {"x", "y", "z"};
		double *node[] = {Xn, Yn, Zn};
		const int ncell[] = {Nx, Ny, Nz};
		for (int a = 0; a < 3; a++) {
			double dmin = node[a][1] - node[a][0];
			double dmax = dmin;
			for (int i = 1; i < ncell[a]; i++) {
				double d = node[a][i + 1] - node[a][i];
				if (d < dmin) dmin = d;
				if (d > dmax) dmax = d;
			}
			if ((dmax - dmin) > (1e-6 * dmax)) {
				sprintf(str, "*** warning : %smesh is not uniform (min=%.6e, max=%.6e). BPM uses the average cell width.",
						axisname[a], dmin, dmax);
				if (io) fprintf(fp, "%s\n", str);
				fprintf(stdout, "%s\n", str);
			}
		}
	}

	// 波長 : OpenFDTD の周波数データ (frequency2 優先) から決定する
	double lambda;
	if      (NFreq2 > 0) lambda = SPEED_OF_LIGHT / Freq2[0];
	else if (NFreq1 > 0) lambda = SPEED_OF_LIGHT / Freq1[0];
	else                 lambda = 1.55e-6;  // 周波数未指定時の既定値
	const double k0 = (2 * PI) / lambda;

	// 材料毎の複素屈折率テーブル : n = sqrt(epsr - i*sigma/(omega*eps0))
	// 損失を imag(n) > 0 として格納する (applyMultiplier の exp(d*imag(n)), d < 0 の規約)
	floatcomplex *n_mat = (floatcomplex *)malloc(NMaterial * sizeof(floatcomplex));
	{
		const double omega = 2 * PI * SPEED_OF_LIGHT / lambda;
		for (int64_t m = 0; m < NMaterial; m++) {
			double epsr, sgm;
			if (Material[m].type == 2) {
				// 分散性材料は einf で近似 (BPM の単一波長計算)
				epsr = Material[m].einf;
				sgm  = 0;
			}
			else {
				epsr = Material[m].epsr;
				sgm  = Material[m].esgm;
			}
			if (epsr <= 0) epsr = 1;  // PEC 等は BPM では真空扱い
			const double s = sgm / (omega * EPS0);
			const double mag = sqrt((epsr * epsr) + (s * s));
			floatcomplex v = {(float)sqrt((mag + epsr) / 2), (float)sqrt((mag - epsr) / 2)};
			n_mat[m] = v;
		}
	}

	// 材料毎の二光子吸収 (TPA) 係数テーブル [m/W]
	// 入力 (tpa = <material id> <beta[cm/GW]>) から単位変換して格納する
	//   1 cm/GW = 1e-11 m/W (例 : 424 cm/GW = 4.24e-9 m/W)
	double *tpa_mat = (double *)malloc(NMaterial * sizeof(double));
	for (int64_t m = 0; m < NMaterial; m++) {
		tpa_mat[m] = 0;
	}
	for (int n = 0; n < NTpaB; n++) {
		tpa_mat[TpaB[n].m] = TpaB[n].beta * 1e-11;
	}
	const int haveTpa = (NTpaB > 0);

	// BPM パラメータ (BPM-MATLAB の FDBPM に対応)
	struct parameters P_var;
	struct parameters *P = &P_var;
	P->Nx = Nx;
	P->Ny = Ny;
	P->dx = (float)Dx;
	P->dy = (float)Dy;
	P->dz = (float)Dz;
	P->iz_start = 0;
	P->iz_end = Nz;
	// 対称境界 (symmetry = <x|y|xy> [sym|anti]) :
	// 鏡像面は指定軸のメッシュ始端。カーネルの命名は「対称面の軸」基準のため、
	// x 方向の鏡像 (x = xmin) は ySymmetry、y 方向の鏡像 (y = ymin) は xSymmetry。
	P->ySymmetry = (unsigned char)BPM.symx;   // x = xmin の鏡像面
	P->xSymmetry = (unsigned char)BPM.symy;   // y = ymin の鏡像面

	// 参照屈折率 : 入力 (refindex =) を優先し、
	// 未指定なら計算領域中心 (z=先頭スライス) の材料から取得する
	if (BPM.n0 > 0) {
		P->n_0 = (float)BPM.n0;
	}
	else {
		int64_t nn0 = NA(Nx / 2, Ny / 2, 0);
		P->n_0 = CREALF(n_mat[iEx[nn0]]);
	}

	// 位相・損失係数 d = -dz*k0
	P->d = (float)(-P->dz * k0);

	// 複素屈折率分布 : 3D 配列として一括でデバイスへ転送する
	// (カーネル側は z 方向に区分一定で参照する)
	P->Nx_n = Nx;
	P->Ny_n = Ny;
	P->Nz_n = Nz;
	P->dz_n = P->dz;
	P->n_in = (floatcomplex *)malloc((size_t)P->Nx_n * P->Ny_n * P->Nz_n * sizeof(floatcomplex));
	for (int iz = 0; iz < Nz; iz++) {
		for (int iy = 0; iy < Ny; iy++) {
			for (int ix = 0; ix < Nx; ix++) {
				size_t index = (size_t)ix + ((size_t)P->Nx_n * iy) + ((size_t)P->Nx_n * P->Ny_n * iz);
				P->n_in[index] = n_mat[iEx[NA(ix, iy, iz)]];
			}
		}
	}

	// テーパ (taper = 出口/入口の横方向スケール比) とツイスト (twist = [deg/m]) :
	// 屈折率分布を座標変換 (相似縮小 + 回転) で z に沿って変化させる。
	// 変換時は入力断面 (n_in の先頭スライス) を 2D RIP として参照する
	// (カーネル bpm/FDBPMpropagator.cu の変換branch と同一)。
	const long nzstep = P->iz_end - P->iz_start;
	if (BPM.taper <= 0) {
		sprintf(str, "*** warning : taper must be > 0 (given %g), ignored.", BPM.taper);
		if (io) fprintf(fp, "%s\n", str);
		fprintf(stdout, "%s\n", str);
		BPM.taper = 1;
	}
	P->taperPerStep = (nzstep > 0) ? (float)((1.0 - BPM.taper) / nzstep) : 0.0f;
	P->twistPerStep = (float)(BPM.twist * DTOR * P->dz);
	const int xform = (P->taperPerStep != 0.0f) || (P->twistPerStep != 0.0f);
	if (xform) {
		sprintf(str, "taper : ratio = %.4f (per step %.6e), twist : %.4f [deg/m] (total %.2f [deg])",
		        BPM.taper, P->taperPerStep, BPM.twist, BPM.twist * (Zn[Nz] - Zn[0]));
		if (io) fprintf(fp, "%s\n", str);
		fprintf(stdout, "%s\n", str);
	}
	// 曲げ : 入力 (bend = RoC [dir] [rho_e]) があれば等価屈折率法で反映する
	if (BPM.RoC > 0) {
		P->rho_e = (float)BPM.rho_e;
		P->RoC = (float)BPM.RoC;
		P->sinBendDirection = (float)sin(BPM.bendDir * DTOR);
		P->cosBendDirection = (float)cos(BPM.bendDir * DTOR);
	}
	else {
		P->rho_e = 0.0;
		P->RoC = INFINITY;  // 直線導波路 (曲げ半径無限大)
		P->sinBendDirection = 0.0;
		P->cosBendDirection = 0.0;
	}

	// ADI 法の係数 ax = dz/(4i*dx^2*k0*n0), ay = dz/(4i*dy^2*k0*n0)
	P->ax = -I * (float)(P->dz / (4 * P->dx * P->dx * k0 * P->n_0));
	P->ay = -I * (float)(P->dz / (4 * P->dy * P->dy * k0 * P->n_0));

	// ホスト側バッファ (E2/Eyx/b はデバイス側のみ : createDeviceStructs で確保)
	P->E1 = (floatcomplex *)malloc((P->Nx * P->Ny) * sizeof(floatcomplex));
	P->E2 = NULL;
	P->Eyx = NULL;
	P->b = NULL;
	P->Efinal = (floatcomplex *)malloc((P->Nx * P->Ny) * sizeof(floatcomplex)); // Output E field(Array)
	P->n_out = (floatcomplex *)malloc((P->Nx * P->Ny) * sizeof(floatcomplex)); // Output refractive index(Array)
	P->multiplier = (float *)malloc((P->Nx * P->Ny) * sizeof(float));

	// 初期電界 (ガウシアンビーム) と端部吸収体 (multiplier) の設定
	// - ビームウェスト : 入力 (beam =) を優先し、未指定なら領域幅の 1/4
	// - ビーム中心    : 入力 (beam = w0 x0 y0) > 給電点 (feed) > 領域中心
	// - 振幅          : 給電点の電圧 (未指定なら 1)
	double w0, xc0, yc0;
	{
		const double Lx = Xn[Nx] - Xn[0];
		const double Ly = Yn[Ny] - Yn[0];
		w0 = (BPM.w0 > 0) ? BPM.w0 : 0.25 * MIN(Lx, Ly);  // ビームウェスト
		if (BPM.center) {
			xc0 = BPM.x0;
			yc0 = BPM.y0;
		}
		else if (NFeed > 0) {
			xc0 = Xn[Feed[0].i];
			yc0 = Yn[Feed[0].j];
		}
		else {
			// 対称境界がある軸は鏡像面 (メッシュ始端) を既定の中心にする
			xc0 = BPM.symx ? Xn[0] : 0.5 * (Xn[0] + Xn[Nx]);
			yc0 = BPM.symy ? Yn[0] : 0.5 * (Yn[0] + Yn[Ny]);
		}
		const double volt = (NFeed > 0) ? Feed[0].volt : 1;
		const double xEdge = BPM.symx ? (0.90 * Lx) : (0.45 * Lx);           // 吸収体開始位置 (領域中心からの距離)
		const double yEdge = BPM.symy ? (0.90 * Ly) : (0.45 * Ly);
		const double xmid = BPM.symx ? Xn[0] : 0.5 * (Xn[0] + Xn[Nx]);
		const double ymid = BPM.symy ? Yn[0] : 0.5 * (Yn[0] + Yn[Ny]);
		const double alpha = 3.0e14;              // 吸収係数 [1/m^3] (BPM-MATLAB 既定値)
		// 入射ビームの傾き (beamtilt =) : 横方向波数の位相ランプ
		const double ktx = k0 * P->n_0 * sin(BPM.tiltx * DTOR);
		const double kty = k0 * P->n_0 * sin(BPM.tilty * DTOR);
		double power = 0;
		for (int iy = 0; iy < P->Ny; iy++) {
			for (int ix = 0; ix < P->Nx; ix++) {
				long i = ix + (long)iy * P->Nx;
				double xd = Xc[ix] - xc0;
				double yd = Yc[iy] - yc0;
				float amp = (float)(volt * exp(-(xd * xd + yd * yd) / (w0 * w0)));
				double ph = -(ktx * xd) - (kty * yd);
				floatcomplex e = {(float)(amp * cos(ph)), (float)(amp * sin(ph))};
				P->E1[i] = e;
				power += (double)amp * amp;
				double dist = MAX(0.0, MAX(fabs(Xc[ix] - xmid) - xEdge, fabs(Yc[iy] - ymid) - yEdge));
				P->multiplier[i] = (float)exp(-P->dz * dist * dist * alpha);
			}
		}
		P->precisePower = power;
	}

	// ============================================================
	// モード解析 (modes = <nModes> [excite]) : CPU 版と同一のホスト側処理
	//   (bpm/modes.cpp の虚軸伝搬法。曲げは反映せず、スカラー演算子で解析)
	// ============================================================
	int nModesFound = 0;
	std::complex<double> *modeFields = NULL;
	double *modeNeff = NULL;
	if (BPM.nmodes > 0) {
		const long NN2 = (long)P->Nx * P->Ny;
		wabpm_params Wm;
		Wm.Nx = Nx;
		Wm.Ny = Ny;
		Wm.dx = P->dx;
		Wm.dy = P->dy;
		Wm.k0 = k0;
		Wm.n0 = P->n_0;
		Wm.wideangle = 0;
		Wm.pol = 0;
		Wm.symx = BPM.symx;
		Wm.symy = BPM.symy;

		// 入力断面の比誘電率 (実部のみ : モードは無損失断面で定義する)
		std::complex<double> *n2m = new std::complex<double>[NN2];
		double nmax = P->n_0;
		for (int iy = 0; iy < P->Ny; iy++) {
			for (int ix = 0; ix < P->Nx; ix++) {
				const floatcomplex nm = n_mat[iEx[NA(ix, iy, P->iz_start)]];
				const double nr = CREALF(nm);
				n2m[ix + (long)iy * P->Nx] = std::complex<double>(nr * nr, 0.0);
				if (nr > nmax) nmax = nr;
			}
		}
		const double mu_max = k0 * k0 * (nmax * nmax - P->n_0 * P->n_0);
		Wm.dz = (mu_max > 0) ? (2.0 * k0 * P->n_0 / mu_max) : (10.0 * P->dx);

		modeFields = new std::complex<double>[(size_t)BPM.nmodes * NN2];
		modeNeff = new double[BPM.nmodes];
		nModesFound = wabpm_find_modes(&Wm, n2m, BPM.nmodes, 20000, 1e-12,
		                               modeFields, modeNeff);
		delete[] n2m;

		sprintf(str, "modes : %d / %d converged (imaginary-distance BPM, scalar)",
		        nModesFound, BPM.nmodes);
		if (io) fprintf(fp, "%s\n", str);
		fprintf(stdout, "%s\n", str);
		for (int m = 0; m < nModesFound; m++) {
			sprintf(str, "mode %d : neff = %.6f", m + 1, modeNeff[m]);
			if (io) fprintf(fp, "%s\n", str);
			fprintf(stdout, "%s\n", str);
		}

		// モード整合励振 : ピーク振幅 = 給電電圧、beamtilt の位相ランプは維持
		if (BPM.modeExcite && nModesFound > 0) {
			const double volt = (NFeed > 0) ? Feed[0].volt : 1;
			const double ktx = k0 * P->n_0 * sin(BPM.tiltx * DTOR);
			const double kty = k0 * P->n_0 * sin(BPM.tilty * DTOR);
			double mmax = 0;
			for (long i = 0; i < NN2; i++) {
				mmax = MAX(mmax, std::abs(modeFields[i]));
			}
			const double scale = (mmax > 0) ? (volt / mmax) : 1;
			double power = 0;
			for (int iy = 0; iy < P->Ny; iy++) {
				for (int ix = 0; ix < P->Nx; ix++) {
					const long i = ix + (long)iy * P->Nx;
					const std::complex<double> a = scale * modeFields[i];
					const double ph = -(ktx * (Xc[ix] - xc0)) - (kty * (Yc[iy] - yc0));
					const std::complex<double> e =
					    a * std::complex<double>(cos(ph), sin(ph));
					floatcomplex ef = {(float)e.real(), (float)e.imag()};
					P->E1[i] = ef;
					power += std::norm(e);
				}
			}
			P->precisePower = power;
			sprintf(str, "modes : mode-matched excitation with mode 1 (peak amp = %g)", volt);
			if (io) fprintf(fp, "%s\n", str);
			fprintf(stdout, "%s\n", str);
		}
	}

	P->EfieldPower = 0;
	P->precisePowerDiff = 0;

	// 出力用ホスト配列の初期化 (BPM では FDTD の反復を行わないため明示的にゼロ化する)
	memset(Eiter, 0, Iter_size);
	memset(Hiter, 0, Iter_size);
	Niter = 0;
	if (NFeed > 0) {
		memset(VFeed, 0, Feed_size);
		memset(IFeed, 0, Feed_size);
	}
	if (NPoint > 0) {
		memset(VPoint, 0, Point_size);
	}
	// CUDA 版ではホスト側の near3d 配列は memcopy3_gpu() 内で確保されるため、
	// BPM (FDTD 反復なし) では writeout 用にゼロ配列をここで確保する
	if ((NN > 0) && (NFreq2 > 0)) {
		const size_t num = (size_t)NFreq2 * NN;
		cEx_r = (float *)calloc(num, sizeof(float));
		cEx_i = (float *)calloc(num, sizeof(float));
		cEy_r = (float *)calloc(num, sizeof(float));
		cEy_i = (float *)calloc(num, sizeof(float));
		cEz_r = (float *)calloc(num, sizeof(float));
		cEz_i = (float *)calloc(num, sizeof(float));
		cHx_r = (float *)calloc(num, sizeof(float));
		cHx_i = (float *)calloc(num, sizeof(float));
		cHy_r = (float *)calloc(num, sizeof(float));
		cHy_i = (float *)calloc(num, sizeof(float));
		cHz_r = (float *)calloc(num, sizeof(float));
		cHz_i = (float *)calloc(num, sizeof(float));
	}

	sprintf(str, "lambda : %.6e [m], n_0 : %.4f", lambda, P->n_0);
	if (io) fprintf(fp, "%s\n", str);
	fprintf(stdout, "%s\n", str);
	sprintf(str, "beam : w0 = %.6e [m], center = (%.6e, %.6e) [m]", w0, xc0, yc0);
	if (io) fprintf(fp, "%s\n", str);
	fprintf(stdout, "%s\n", str);
	if (BPM.RoC > 0) {
		sprintf(str, "bend : RoC = %.6e [m], direction = %.1f [deg]", BPM.RoC, BPM.bendDir);
		if (io) fprintf(fp, "%s\n", str);
		fprintf(stdout, "%s\n", str);
	}

	// ============================================================
	// 非線形吸収 (TPA) / 入力パワー掃引 (ONN 光活性化関数) の準備
	// CPU 版 sol/solve_bpm.cpp と同一の規約 :
	//   - 物理スケーリング : 初期界を ∫∫|E|^2 dA = P_in [W] に正規化し、
	//     |E|^2 がそのまま強度 I [W/m^2] になるようにする (dA = Dx*Dy)。
	//     tpa / powersweep 未指定時は従来通り無次元の界のまま計算する (後方互換)。
	//   - 実効断面積 A_eff = (∫|E|^2 dA)^2 / ∫|E|^4 dA を初期界から実計算する
	// ============================================================
	const int    sweep    = (PowerSweep.npoints > 0);
	const int    nsweep   = sweep ? PowerSweep.npoints : 1;
	const int    physcale = (sweep || haveTpa);   // 物理スケーリングの有無
	const double E0rawpow = P->precisePower;      // 初期界の生の |E|^2 和 (面積要素なし)
	floatcomplex *E0 = NULL;                      // 初期界の保存 (掃引毎の再初期化用)
	double *sweepPin = NULL, *sweepPout = NULL;
	// TPA 係数のスライス配列 (デバイスへ転送する)
	float *betaSlice = haveTpa
	    ? (float *)malloc((size_t)P->Nx * P->Ny * sizeof(float)) : NULL;
	if (physcale) {
		E0 = (floatcomplex *)malloc((size_t)P->Nx * P->Ny * sizeof(floatcomplex));
		memcpy(E0, P->E1, (size_t)P->Nx * P->Ny * sizeof(floatcomplex));
		sweepPin  = (double *)malloc(nsweep * sizeof(double));
		sweepPout = (double *)malloc(nsweep * sizeof(double));
		for (int n = 0; n < nsweep; n++) {
			if (!sweep) {
				sweepPin[n] = 1;  // tpa のみ指定時の既定入力パワー [W]
			}
			else if (nsweep == 1) {
				sweepPin[n] = PowerSweep.pmin;
			}
			else if (PowerSweep.logscale) {
				sweepPin[n] = PowerSweep.pmin *
					pow(PowerSweep.pmax / PowerSweep.pmin, (double)n / (nsweep - 1));
			}
			else {
				sweepPin[n] = PowerSweep.pmin +
					((PowerSweep.pmax - PowerSweep.pmin) * n) / (nsweep - 1);
			}
			sweepPout[n] = 0;
		}
		// 実効断面積 (初期界から実計算。ガウシアンでは A_eff = pi*w0^2 に一致)
		double sumI2 = 0;
		for (long i = 0; i < (long)P->Nx * P->Ny; i++) {
			const floatcomplex e = E0[i];
			const double i2 = ((double)CREALF(e) * CREALF(e)) + ((double)CIMAGF(e) * CIMAGF(e));
			sumI2 += i2 * i2;
		}
		const double Aeff = (sumI2 > 0) ? (E0rawpow * E0rawpow / sumI2) * (Dx * Dy) : 0;
		sprintf(str, "ONN: A_eff = %.6e [m^2] (from input field), L = %.6e [m]",
				Aeff, Zn[Nz] - Zn[0]);
		if (io) fprintf(fp, "%s\n", str);
		fprintf(stdout, "%s\n", str);
	}
	if (haveTpa) {
		for (int n = 0; n < NTpaB; n++) {
			sprintf(str, "ONN: TPA material id = %d, beta = %g [cm/GW] = %g [m/W]",
					(int)TpaB[n].m, TpaB[n].beta, tpa_mat[TpaB[n].m]);
			if (io) fprintf(fp, "%s\n", str);
			fprintf(stdout, "%s\n", str);
		}
		if (!sweep) {
			sprintf(str, "ONN: no powersweep -> single run with P_in = %g [W]", sweepPin[0]);
			if (io) fprintf(fp, "%s\n", str);
			fprintf(stdout, "%s\n", str);
		}
	}

	// 伝搬の可視化用 : 中心行 (y = Ny/2) の強度 |E(x, z)|^2 を全ステップ記録する
	float *Ixz = (float *)malloc((size_t)P->Nx * (P->iz_end - P->iz_start) * sizeof(float));
	floatcomplex *Erow = (floatcomplex *)malloc(P->Nx * sizeof(floatcomplex));

	// |E(x,y)|^2 スナップショット (frames = interval 指定時のみ, /field/frames へ出力)
	const int  frameInterval = BPM.frames;
	const long nframes = (frameInterval > 0)
	    ? ((P->iz_end - P->iz_start - 1) / frameInterval + 1) : 0;
	float *Frames = (nframes > 0)
	    ? (float *)malloc((size_t)nframes * P->Nx * P->Ny * sizeof(float)) : NULL;
	if (nframes > 0) {
		sprintf(str, "frames : interval = %d steps, count = %ld", frameInterval, nframes);
		if (io) fprintf(fp, "%s\n", str);
		fprintf(stdout, "%s\n", str);
	}

	// HDF5ファイルの作成
	file_id = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// ============================================================
	// 入力パワー掃引ループ (powersweep 未指定時は 1 回のみ実行)
	// 各回 : 初期界を P_in [W] に正規化 -> z 伝搬 -> P_out を記録する。
	// /field 出力 (Ixz/Efinal/frames) は最終掃引点の結果になる。
	// ============================================================
	for (int isweep = 0; isweep < nsweep; isweep++) {

	// 物理スケーリング時 : 初期界を ∫∫|E|^2 dA = P_in になるよう再設定する
	if (physcale) {
		const double Pin = sweepPin[isweep];
		const float s = (float)sqrt(Pin / (E0rawpow * Dx * Dy));
		for (long i = 0; i < (long)P->Nx * P->Ny; i++) {
			P->E1[i] = E0[i] * s;
		}
		P->precisePower = E0rawpow * (double)s * s;
		P->EfieldPower = 0;
		P->precisePowerDiff = 0;
		sprintf(str, "sweep %d / %d : P_in = %g [W]", isweep + 1, nsweep, Pin);
		fprintf(stdout, "%s\n", str);
	}

	if (BPM.pol || BPM.wideangle) {
		// ============================================================
		// 拡張パス : 広角 BPM (Pade(1,1)) / 半ベクトル BPM (倍精度, bpm/wabpm.cu)
		// 屈折率項を演算子の内部に含めるため、一般化 ADI で伝搬する。
		// 曲げは等価屈折率法で n2 スライスへ反映する (近軸パスと同じ変換)。
		// ============================================================
		sprintf(str, "mode : %s, polarization = %s (CUDA)",
				(BPM.wideangle ? "wide-angle Pade(1,1)" : "paraxial"),
				(BPM.pol == 1 ? "x (semivectorial)" : (BPM.pol == 2 ? "y (semivectorial)" : "scalar")));
		if (io) fprintf(fp, "%s\n", str);
		fprintf(stdout, "%s\n", str);

		wabpm_params W;
		W.Nx = Nx;
		W.Ny = Ny;
		W.dx = P->dx;
		W.dy = P->dy;
		W.dz = P->dz;
		W.k0 = k0;
		W.n0 = P->n_0;
		W.wideangle = BPM.wideangle;
		W.pol = BPM.pol;
		W.symx = BPM.symx;
		W.symy = BPM.symy;

		wabpm_cplx *Ed = (wabpm_cplx *)malloc((size_t)Nx * Ny * sizeof(wabpm_cplx));
		wabpm_cplx *n2 = (wabpm_cplx *)malloc((size_t)Nx * Ny * sizeof(wabpm_cplx));

		// 曲げ (等価屈折率法) : n_bend = n*(1 - n^2*xb*rho_e/(2*RoC))*exp(xb/RoC)
		// (xb = x*cos(dir) + y*sin(dir), 吸収 = 虚部は変換しない)
		const int  bend = (BPM.RoC > 0);
		const double cosB = bend ? cos(BPM.bendDir * DTOR) : 1.0;
		const double sinB = bend ? sin(BPM.bendDir * DTOR) : 0.0;

		// 初期電界 (チルト位相込みの E1) を倍精度へ
		for (long i = 0; i < (long)Nx * Ny; i++) {
			Ed[i] = wabpm_cplx(CREALF(P->E1[i]), CIMAGF(P->E1[i]));
		}

		struct wabpm_gpu *G = wabpm_gpu_create(&W, Ed, P->multiplier);

		for (long iz = P->iz_start; iz < P->iz_end; iz++) {
			sprintf(str, "step iz : %ld / %ld", iz + 1, (long)Nz);
			fprintf(stdout, "%s\n", str);

			const double t0 = cputime();

			// このスライスの複素比誘電率 (n_mat は損失を +imag で保持 -> 物理符号 n = nr - i*ni)
			// テーパ/ツイスト指定時は入力断面を座標変換 (相似縮小 + 回転) して参照する
			// (CPU 版 sol/solve_bpm.cpp と同一の変換式)
			const long izr = iz - P->iz_start;
			const double xfScale = 1.0 / (1.0 - (double)P->taperPerStep * izr);
			const double xfCos = cos(-(double)P->twistPerStep * izr);
			const double xfSin = sin(-(double)P->twistPerStep * izr);
			for (int iy = 0; iy < Ny; iy++) {
				for (int ix = 0; ix < Nx; ix++) {
					long index = (long)(Nx * iy) + ix;
					double nr, nimag;
					if (xform) {
						const double x = P->dx * (ix - (Nx - 1) / 2.0);
						const double y = P->dy * (iy - (Ny - 1) / 2.0);
						const double xs = xfScale * ((xfCos * x) - (xfSin * y));
						const double ys = xfScale * ((xfSin * x) + (xfCos * y));
						double fx = (xs / P->dx) + ((Nx - 1) / 2.0);
						double fy = (ys / P->dy) + ((Ny - 1) / 2.0);
						fx = MIN(MAX(0.0, fx), (Nx - 1) * (1.0 - 1e-7));
						fy = MIN(MAX(0.0, fy), (Ny - 1) * (1.0 - 1e-7));
						const int i0 = (int)floor(fx);
						const int j0 = (int)floor(fy);
						const double tx = fx - i0;
						const double ty = fy - j0;
						const floatcomplex n00 = n_mat[iEx[NA(i0,     j0,     P->iz_start)]];
						const floatcomplex n10 = n_mat[iEx[NA(i0 + 1, j0,     P->iz_start)]];
						const floatcomplex n01 = n_mat[iEx[NA(i0,     j0 + 1, P->iz_start)]];
						const floatcomplex n11 = n_mat[iEx[NA(i0 + 1, j0 + 1, P->iz_start)]];
						nr    = ((1 - tx) * (1 - ty) * CREALF(n00)) + (tx * (1 - ty) * CREALF(n10))
						      + ((1 - tx) * ty * CREALF(n01))       + (tx * ty * CREALF(n11));
						nimag = ((1 - tx) * (1 - ty) * CIMAGF(n00)) + (tx * (1 - ty) * CIMAGF(n10))
						      + ((1 - tx) * ty * CIMAGF(n01))       + (tx * ty * CIMAGF(n11));
					}
					else {
						const floatcomplex nm = n_mat[iEx[NA(ix, iy, iz)]];
						nr    = CREALF(nm);
						nimag = CIMAGF(nm);
					}
					if (bend) {
						const double x = P->dx * (ix - (Nx - 1) / 2.0);
						const double y = P->dy * (iy - (Ny - 1) / 2.0);
						const double xb = x * cosB + y * sinB;
						nr = nr * (1.0 - nr * nr * xb * BPM.rho_e / (2.0 * BPM.RoC))
						        * exp(xb / BPM.RoC);
					}
					wabpm_cplx n(nr, -nimag);
					n2[index] = n * n;
				}
			}

			// 1 ステップ伝搬 (+ 端部吸収体, wabpm_gpu_step 内で適用)
			wabpm_gpu_step(G, n2);

			// TPA (二光子吸収) : E *= exp(-(beta/2)*I*dz) をデバイス上で適用
			if (haveTpa) {
				for (int iy = 0; iy < Ny; iy++) {
					for (int ix = 0; ix < Nx; ix++) {
						betaSlice[(long)(Nx * iy) + ix] = (float)tpa_mat[iEx[NA(ix, iy, iz)]];
					}
				}
				wabpm_gpu_tpa(G, betaSlice, Dz);
			}

			// 中心行の強度を記録
			wabpm_gpu_get_row(G, P->Ny / 2, Ed);  // Ed の先頭 Nx 要素を一時利用
			for (long ix = 0; ix < Nx; ix++) {
				Ixz[(iz - P->iz_start) * Nx + ix] = (float)thrust::norm(Ed[ix]);
			}

			// スナップショットを記録 (全電界をホストへ取得)
			if (Frames && ((iz - P->iz_start) % frameInterval == 0)) {
				wabpm_gpu_get_field(G, Ed);
				float *f = &Frames[(size_t)((iz - P->iz_start) / frameInterval) * Nx * Ny];
				for (long i = 0; i < (long)Nx * Ny; i++) {
					f[i] = (float)thrust::norm(Ed[i]);
				}
			}

			*tdft += cputime() - t0;
		}

		// 最終電界・屈折率分布・出力電力を格納
		wabpm_gpu_get_field(G, Ed);
		double power = 0;
		for (long i = 0; i < (long)Nx * Ny; i++) {
			floatcomplex e = {(float)Ed[i].real(), (float)Ed[i].imag()};
			P->Efinal[i] = e;
			power += thrust::norm(Ed[i]);
		}
		P->precisePower = power;
		for (int iy = 0; iy < Ny; iy++) {
			for (int ix = 0; ix < Nx; ix++) {
				long index = (long)(Nx * iy) + ix;
				P->n_out[index] = n_mat[iEx[NA(ix, iy, P->iz_end - 1)]];
			}
		}
		wabpm_gpu_destroy(G);
		free(Ed);
		free(n2);
	}
	else {
	// デバイス側構造体の確保と転送 (E1, multiplier, n_in を転送し E2/Eyx/b/n_out を確保)
	struct parameters *P_dev;
	struct debug D_var = {{0.0, 0.0, 0.0}, {0, 0, 0}};
	struct debug *D = &D_var;
	struct debug *D_dev;
	createDeviceStructs(P, &P_dev, D, &D_dev);

	// カーネル起動パラメータ
	int temp, nBlocks;
	CUDA_CHECK(cudaOccupancyMaxPotentialBlockSize(&nBlocks, &temp, &substep1a, 0, 0));
	dim3 blockDims(TILE_DIM, TILE_DIM, 1);

	// TPA 用のデバイスバッファ (beta スライスと電力簿記の集計)
	// applyTPA は 1 次元のグリッドストライドループのため 1D ブロックで起動する
	// (2D ブロックだと threadIdx.y 方向のスレッドが同一要素を重複更新してしまう)
	float  *d_beta = NULL;
	double *d_sums = NULL;
	const int tpaTpb = 256;
	const int tpaNb = (int)((((long)P->Nx * P->Ny) + tpaTpb - 1) / tpaTpb);
	if (haveTpa) {
		CUDA_CHECK(cudaMalloc(&d_beta, (size_t)P->Nx * P->Ny * sizeof(float)));
		CUDA_CHECK(cudaMalloc(&d_sums, 2 * sizeof(double)));
	}

	// BPM は時間反復ではなく Z 軸方向に 1 回伝搬する
	for (long iz = P->iz_start; iz < P->iz_end; iz++) {
		sprintf(str, "step iz : %ld / %ld", iz + 1, (long)Nz);
		fprintf(stdout, "%s\n", str);

		const double t0 = cputime();

		// Douglas-Gunn ADI 法による 1 ステップ伝搬
		substep1a<<<nBlocks, blockDims>>>(P_dev); // xy -> yx
		substep1b<<<nBlocks, blockDims>>>(P_dev); // yx -> yx
		substep2a<<<nBlocks, blockDims>>>(P_dev); // yx -> xy
		substep2b<<<nBlocks, blockDims>>>(P_dev); // xy -> xy
		applyMultiplier<<<nBlocks, blockDims>>>(P_dev, iz, D_dev); // xy -> xy

		CUDA_CHECK(cudaDeviceSynchronize()); // Wait until all kernels have finished

		// TPA (二光子吸収) : E2 *= exp(-(beta/2)*I*dz)
		// 次ステップの fieldCorrection が TPA 減衰を打ち消さないよう、
		// 電力簿記 precisePower も同率で減らす (CPU 版と同一の手順)
		if (haveTpa) {
			updatePrecisePower<<<1, 1>>>(P_dev);  // 係留中の precisePowerDiff を先に確定する
			for (int iy = 0; iy < P->Ny; iy++) {
				for (int ix = 0; ix < P->Nx; ix++) {
					betaSlice[(long)(P->Nx * iy) + ix] = (float)tpa_mat[iEx[NA(ix, iy, iz)]];
				}
			}
			CUDA_CHECK(cudaMemcpy(d_beta, betaSlice, (size_t)P->Nx * P->Ny * sizeof(float), cudaMemcpyHostToDevice));
			CUDA_CHECK(cudaMemset(d_sums, 0, 2 * sizeof(double)));
			applyTPA<<<tpaNb, tpaTpb>>>(P_dev, d_beta, Dz, d_sums);
			scalePrecisePowerByTPA<<<1, 1>>>(P_dev, d_sums);
			CUDA_CHECK(cudaDeviceSynchronize());
		}

		// 中心行の強度を記録 (このステップの結果はデバイス側 E2 にある)
		{
			struct parameters P_now;
			CUDA_CHECK(cudaMemcpy(&P_now, P_dev, sizeof(struct parameters), cudaMemcpyDeviceToHost));
			const long row0 = (long)(P->Ny / 2) * P->Nx;
			CUDA_CHECK(cudaMemcpy(Erow, P_now.E2 + row0, P->Nx * sizeof(floatcomplex), cudaMemcpyDeviceToHost));
			for (long ix = 0; ix < P->Nx; ix++) {
				Ixz[(iz - P->iz_start) * P->Nx + ix] =
					(CREALF(Erow[ix]) * CREALF(Erow[ix])) + (CIMAGF(Erow[ix]) * CIMAGF(Erow[ix]));
			}

			// スナップショットを記録 (デバイス側 E2 の全体をコピー)
			if (Frames && ((iz - P->iz_start) % frameInterval == 0)) {
				floatcomplex *Efull = (floatcomplex *)malloc((size_t)P->Nx * P->Ny * sizeof(floatcomplex));
				CUDA_CHECK(cudaMemcpy(Efull, P_now.E2, (size_t)P->Nx * P->Ny * sizeof(floatcomplex), cudaMemcpyDeviceToHost));
				float *f = &Frames[(size_t)((iz - P->iz_start) / frameInterval) * P->Nx * P->Ny];
				for (long i = 0; i < (long)P->Nx * P->Ny; i++) {
					f[i] = (CREALF(Efull[i]) * CREALF(Efull[i])) + (CIMAGF(Efull[i]) * CIMAGF(Efull[i]));
				}
				free(Efull);
			}
		}

		if (iz + 1 < P->iz_end) swapEPointers<<<1, 1>>>(P_dev, iz);
		updatePrecisePower<<<1, 1>>>(P_dev);
		CUDA_CHECK(cudaDeviceSynchronize());

		*tdft += cputime() - t0;
	}

	// 最終電界 (デバイス側 E2) と屈折率分布を Efinal / n_out へ回収し、デバイスメモリを解放する
	CUDA_CHECK(cudaDeviceSynchronize());
	retrieveAndFreeDeviceStructs(P, P_dev, D, D_dev);
	if (haveTpa) {
		CUDA_CHECK(cudaFree(d_beta));
		CUDA_CHECK(cudaFree(d_sums));
	}
	}

	// 物理スケーリング時 : 出力パワー P_out = ∫∫|E_end|^2 dA [W] を記録する
	if (physcale) {
		double rawpow = 0;
		for (long i = 0; i < (long)P->Nx * P->Ny; i++) {
			const floatcomplex e = P->Efinal[i];
			rawpow += ((double)CREALF(e) * CREALF(e)) + ((double)CIMAGF(e) * CIMAGF(e));
		}
		sweepPout[isweep] = rawpow * Dx * Dy;
		sprintf(str, "ONN: P_in=%g W -> P_out=%g W (T=%g)",
				sweepPin[isweep], sweepPout[isweep],
				(sweepPin[isweep] > 0) ? (sweepPout[isweep] / sweepPin[isweep]) : 0);
		if (io) {
			fprintf(fp, "%s\n", str);
			fflush(fp);
		}
		fprintf(stdout, "%s\n", str);
	}

	}  // 掃引ループ (isweep) 終了

	// 光活性化関数曲線 P_out(P_in) の CSV 出力
	if (sweep) {
		FILE *fcsv = fopen(FN_activation, "w");
		if (fcsv != NULL) {
			fprintf(fcsv, "P_in_W,P_out_W,transmission\n");
			for (int n = 0; n < nsweep; n++) {
				fprintf(fcsv, "%.6e,%.6e,%.6e\n", sweepPin[n], sweepPout[n],
						(sweepPin[n] > 0) ? (sweepPout[n] / sweepPin[n]) : 0);
			}
			fclose(fcsv);
			sprintf(str, "ONN: activation curve -> %s (%d points)", FN_activation, nsweep);
			if (io) fprintf(fp, "%s\n", str);
			fprintf(stdout, "%s\n", str);
		}
	}

	// result
	if (io) {
		sprintf(str, "    --- propagated %ld steps, output power = %.6e ---",
				(long)(P->iz_end - P->iz_start), P->precisePower);
		fprintf(fp,     "%s\n", str);
		fprintf(stdout, "%s\n", str);
		fflush(fp);
		fflush(stdout);
	}

	// time steps
	Ntime = 1;

	// 結果 (最終電界・屈折率分布) の書き込み
	{
		hid_t field_group_id = H5Gcreate(file_id, "/field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hsize_t field_dims[2] = {(hsize_t)P->Ny, (hsize_t)P->Nx};
		float *tmp = (float *)malloc(P->Nx * P->Ny * sizeof(float));
		struct {
			const char *name;
			floatcomplex *data;
			int imagpart;
		} field_data[] = {
			{"Efinal_r", P->Efinal, 0},
			{"Efinal_i", P->Efinal, 1},
			{"n_out_r",  P->n_out,  0},
			{"n_out_i",  P->n_out,  1}
		};
		for (size_t n = 0; n < sizeof(field_data) / sizeof(field_data[0]); n++) {
			for (long i = 0; i < P->Nx * P->Ny; i++) {
				tmp[i] = field_data[n].imagpart ? CIMAGF(field_data[n].data[i])
												: CREALF(field_data[n].data[i]);
			}
			dataspace_id = H5Screate_simple(2, field_dims, NULL);
			dataset_id = H5Dcreate(field_group_id, field_data[n].name, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
		free(tmp);

		// 伝搬マップ |E(x, y=Ny/2, z)|^2 の書き込み
		{
			hsize_t ixz_dims[2] = {(hsize_t)(P->iz_end - P->iz_start), (hsize_t)P->Nx};
			dataspace_id = H5Screate_simple(2, ixz_dims, NULL);
			dataset_id = H5Dcreate(field_group_id, "Ixz", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Ixz);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}

		// スナップショット |E(x,y)|^2 (nframes x Ny x Nx) の書き込み
		if (Frames) {
			hsize_t fr_dims[3] = {(hsize_t)nframes, (hsize_t)P->Ny, (hsize_t)P->Nx};
			dataspace_id = H5Screate_simple(3, fr_dims, NULL);
			dataset_id = H5Dcreate(field_group_id, "frames", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Frames);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
		H5Gclose(field_group_id);
	}

	// モード解析結果 (/modes) の書き込み : mode<i> (Ny x Nx) と neff
	if (nModesFound > 0) {
		hid_t modes_group_id = H5Gcreate(file_id, "/modes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		const long NN2 = (long)P->Nx * P->Ny;
		float *tmp = (float *)malloc((size_t)NN2 * sizeof(float));
		char dsname[32];
		for (int m = 0; m < nModesFound; m++) {
			for (long i = 0; i < NN2; i++) {
				tmp[i] = (float)modeFields[(size_t)m * NN2 + i].real();
			}
			hsize_t mdims[2] = {(hsize_t)P->Ny, (hsize_t)P->Nx};
			sprintf(dsname, "mode%d", m + 1);
			dataspace_id = H5Screate_simple(2, mdims, NULL);
			dataset_id = H5Dcreate(modes_group_id, dsname, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
		free(tmp);
		{
			hsize_t ndims[1] = {(hsize_t)nModesFound};
			dataspace_id = H5Screate_simple(1, ndims, NULL);
			dataset_id = H5Dcreate(modes_group_id, "neff", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, modeNeff);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
		H5Gclose(modes_group_id);
	}

	// メタデータの作成
	hid_t metadata_group_id = H5Gcreate(file_id, "/metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// 時間に関するメタデータの書き込み
	double time_metadata[1] = {Solver.maxiter * Dt};
	hsize_t time_count[1] = {1};
	dataspace_id = H5Screate_simple(1, time_count, NULL);
	dataset_id = H5Dcreate(metadata_group_id, "time", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_metadata);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);

	// Title
	hsize_t title_dims[1] = {256};
	dataspace_id = H5Screate_simple(1, title_dims, NULL);
	dataset_id = H5Dcreate(metadata_group_id, "Title", H5T_NATIVE_CHAR, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, Title);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);

	// 各種整数型メタデータの書き込み
	struct {
		const char *name;
		void *value;
		hid_t type;
	} metadata[] = {
		{"Nx", &Nx, H5T_NATIVE_INT},
		{"Ny", &Ny, H5T_NATIVE_INT},
		{"Nz", &Nz, H5T_NATIVE_INT},
		{"Ni", &Ni, H5T_NATIVE_INT},
		{"Nj", &Nj, H5T_NATIVE_INT},
		{"Nk", &Nk, H5T_NATIVE_INT},
		{"N0", &N0, H5T_NATIVE_INT},
		{"NN", &NN, H5T_NATIVE_INT64},
		{"NFreq1", &NFreq1, H5T_NATIVE_INT},
		{"NFreq2", &NFreq2, H5T_NATIVE_INT},
		{"NFeed", &NFeed, H5T_NATIVE_INT},
		{"NPoint", &NPoint, H5T_NATIVE_INT},
		{"Niter", &Niter, H5T_NATIVE_INT},
		{"Ntime", &Ntime, H5T_NATIVE_INT},
		{"Solver_maxiter", &Solver.maxiter, H5T_NATIVE_INT},
		{"Solver_nout", &Solver.nout, H5T_NATIVE_INT},
		{"NGline", &NGline, H5T_NATIVE_INT},
		{"IPlanewave", &IPlanewave, H5T_NATIVE_INT}
	};

	for (size_t i = 0; i < sizeof(metadata) / sizeof(metadata[0]); i++) {
		dataspace_id = H5Screate(H5S_SCALAR);
		dataset_id = H5Dcreate(metadata_group_id, metadata[i].name, metadata[i].type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, metadata[i].type, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata[i].value);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
	}

	// Dtの書き込み
	dataspace_id = H5Screate(H5S_SCALAR);
	dataset_id = H5Dcreate(metadata_group_id, "Dt", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Dt);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);

	// BPM パラメータの書き込み
	{
		const double n_0_d = P->n_0;
		const double frame_interval_d = frameInterval;
		struct {
			const char *name;
			const double *value;
		} bpm_metadata[] = {
			{"lambda",  &lambda},
			{"n_0",     &n_0_d},
			{"beam_w0", &w0},
			{"beam_x0", &xc0},
			{"beam_y0", &yc0},
			{"frame_interval", &frame_interval_d},
			{"grid_dx", &Dx},
			{"grid_dy", &Dy},
			{"grid_dz", &Dz}
		};
		for (size_t i = 0; i < sizeof(bpm_metadata) / sizeof(bpm_metadata[0]); i++) {
			dataspace_id = H5Screate(H5S_SCALAR);
			dataset_id = H5Dcreate(metadata_group_id, bpm_metadata[i].name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, bpm_metadata[i].value);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
	}

	// Planewaveの書き込み
	dataspace_id = H5Screate(H5S_SCALAR);
	dataset_id = H5Dcreate(metadata_group_id, "Planewave", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Planewave);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);

	// 配列データの書き込み
	struct {
		const char *name;
		double *data;
		size_t size;
	} arrays[] = {
		{"Xn", Xn, (size_t)(Nx + 1)},
		{"Yn", Yn, (size_t)(Ny + 1)},
		{"Zn", Zn, (size_t)(Nz + 1)},
		{"Xc", Xc, (size_t)Nx},
		{"Yc", Yc, (size_t)Ny},
		{"Zc", Zc, (size_t)Nz},
		{"Eiter", Eiter, (size_t)Niter},
		{"Hiter", Hiter, (size_t)Niter},
		{"VFeed", VFeed, (size_t)(NFeed * (Solver.maxiter + 1))},
		{"IFeed", IFeed, (size_t)(NFeed * (Solver.maxiter + 1))},
		{"VPoint", VPoint, (size_t)(NPoint * (Solver.maxiter + 1))},
		{"Freq1", Freq1, (size_t)NFreq1},
		{"Freq2", Freq2, (size_t)NFreq2},
		{"Gline", reinterpret_cast<double*>(Gline), (size_t)(NGline * 2 * 3)}
	};

	for (size_t i = 0; i < sizeof(arrays) / sizeof(arrays[0]); i++) {
		if ((arrays[i].size == 0) || (arrays[i].data == NULL)) continue;
		hsize_t array_dims[1] = {arrays[i].size};
		dataspace_id = H5Screate_simple(1, array_dims, NULL);
		dataset_id = H5Dcreate(metadata_group_id, arrays[i].name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arrays[i].data);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
	}

	// NSurfaceデータの書き込み
	dataspace_id = H5Screate(H5S_SCALAR);
	dataset_id = H5Dcreate(metadata_group_id, "NSurface", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &NSurface);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);

	// Surfaceデータの書き込み
	// surface_t構造体に対応する複合データ型を定義
	hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(surface_t));
	H5Tinsert(memtype, "nx", HOFFSET(surface_t, nx), H5T_NATIVE_DOUBLE);
	H5Tinsert(memtype, "ny", HOFFSET(surface_t, ny), H5T_NATIVE_DOUBLE);
	H5Tinsert(memtype, "nz", HOFFSET(surface_t, nz), H5T_NATIVE_DOUBLE);
	H5Tinsert(memtype, "x", HOFFSET(surface_t, x), H5T_NATIVE_DOUBLE);
	H5Tinsert(memtype, "y", HOFFSET(surface_t, y), H5T_NATIVE_DOUBLE);
	H5Tinsert(memtype, "z", HOFFSET(surface_t, z), H5T_NATIVE_DOUBLE);
	H5Tinsert(memtype, "ds", HOFFSET(surface_t, ds), H5T_NATIVE_DOUBLE);

	if ((NSurface > 0) && (Surface != NULL)) {
		hsize_t surface_dims[1] = {(hsize_t)NSurface};
		dataspace_id = H5Screate_simple(1, surface_dims, NULL);
		dataset_id = H5Dcreate(metadata_group_id, "Surface", memtype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Surface);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
	}
	H5Tclose(memtype);

	// メタデータグループのクローズ
	H5Gclose(metadata_group_id);

	status = H5Fclose(file_id);

	// free (FDTD 用デバイスメモリ)
	memfree2_gpu();
	memfree3_gpu();

	// ホスト側バッファの解放 (E2/Eyx/b はデバイス側のみで retrieveAndFreeDeviceStructs が解放済み)
	free(P->E1);
	free(P->Efinal);
	free(P->n_out);
	free(P->n_in);
	free(P->multiplier);
	free(n_mat);
	free(Ixz);
	free(Erow);
	if (Frames) free(Frames);
	free(tpa_mat);
	if (betaSlice != NULL) free(betaSlice);
	if (E0 != NULL) free(E0);
	if (sweepPin != NULL) free(sweepPin);
	if (sweepPout != NULL) free(sweepPout);
	if (modeFields) delete[] modeFields;
	if (modeNeff) delete[] modeNeff;
}
