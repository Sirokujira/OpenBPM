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
	P->xSymmetry = 0;
	P->ySymmetry = 0;

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

	P->taperPerStep = 0.0;
	P->twistPerStep = 0.0;
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
			xc0 = 0.5 * (Xn[0] + Xn[Nx]);
			yc0 = 0.5 * (Yn[0] + Yn[Ny]);
		}
		const double volt = (NFeed > 0) ? Feed[0].volt : 1;
		const double xEdge = 0.45 * Lx;           // 吸収体開始位置 (領域中心からの距離)
		const double yEdge = 0.45 * Ly;
		const double xmid = 0.5 * (Xn[0] + Xn[Nx]);
		const double ymid = 0.5 * (Yn[0] + Yn[Ny]);
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

	// TPA (tpa) / パワー掃引 (powersweep) は CUDA 版では未対応 :
	// 指定されていても線形の単発計算になるため、明示的に警告を出す
	if ((NTpaB > 0) || (PowerSweep.npoints > 0)) {
		sprintf(str, "*** warning : tpa / powersweep are not supported in CUDA version (ignored). Use CPU version (obpm) instead.");
		if (io) fprintf(fp, "%s\n", str);
		fprintf(stdout, "%s\n", str);
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
			for (int iy = 0; iy < Ny; iy++) {
				for (int ix = 0; ix < Nx; ix++) {
					long index = (long)(Nx * iy) + ix;
					const floatcomplex nm = n_mat[iEx[NA(ix, iy, iz)]];
					double nr = CREALF(nm);
					if (bend) {
						const double x = P->dx * (ix - (Nx - 1) / 2.0);
						const double y = P->dy * (iy - (Ny - 1) / 2.0);
						const double xb = x * cosB + y * sinB;
						nr = nr * (1.0 - nr * nr * xb * BPM.rho_e / (2.0 * BPM.RoC))
						        * exp(xb / BPM.RoC);
					}
					wabpm_cplx n(nr, -CIMAGF(nm));
					n2[index] = n * n;
				}
			}

			// 1 ステップ伝搬 (+ 端部吸収体, wabpm_gpu_step 内で適用)
			wabpm_gpu_step(G, n2);

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
}
