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


// ============================================================
// HDF5 の自己記述用ヘルパ (GUI・汎用ビューア向け)
//   BPM は定常解法で時間軸を持たない。行進 (marching) 軸は z であり、
//   FDTD の時間軸に相当するのは z である。この対応を読み手が推測せずに
//   済むよう、軸名・単位・座標データセットの所在を属性として書き出す。
//   参考 : 上流の BPM-Matlab も P.z / P.xzSlice / P.modeOverlaps のように
//          z を独立変数として保持する (時間軸には見せない)。
// ============================================================
static void h5AttrStr(hid_t obj, const char *name, const char *value)
{
	hid_t t = H5Tcopy(H5T_C_S1);
	H5Tset_size(t, strlen(value) + 1);
	H5Tset_strpad(t, H5T_STR_NULLTERM);
	hid_t s = H5Screate(H5S_SCALAR);
	hid_t a = H5Acreate(obj, name, t, s, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(a, t, value);
	H5Aclose(a);
	H5Sclose(s);
	H5Tclose(t);
}

static void h5AttrLong(hid_t obj, const char *name, long value)
{
	hid_t s = H5Screate(H5S_SCALAR);
	hid_t a = H5Acreate(obj, name, H5T_NATIVE_LONG, s, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(a, H5T_NATIVE_LONG, &value);
	H5Aclose(a);
	H5Sclose(s);
}

// データセットへ「何を・どの軸で・どの単位で」を付ける。
//   dims        : 軸名をカンマ区切りで並べたもの (先頭が最も外側の軸)
//   coordinates : 各軸の座標データセットのパス (空白区切り、dims と同順)
static void h5Annotate(hid_t obj, const char *long_name, const char *units,
					   const char *dims, const char *coordinates)
{
	h5AttrStr(obj, "long_name", long_name);
	h5AttrStr(obj, "units", units);
	if (dims)        h5AttrStr(obj, "dims", dims);
	if (coordinates) h5AttrStr(obj, "coordinates", coordinates);
}
// z 断面の統計量 (GUI 表示用の /trace) を集計値から算出する。
// CPU 版 sol/solve_bpm.cpp の computeTrace と同一の式。
// acc[0..5] = sum|E|^2, sum x|E|^2, sum y|E|^2, sum x^2|E|^2, sum y^2|E|^2, peak
static void traceFromAcc(const double *acc, double dA,
                         double *power, double *peak,
                         double *cx, double *cy, double *wx, double *wy)
{
	const double s = acc[0];
	*power = s * dA;
	*peak  = acc[5];
	if (s > 0) {
		const double mx = acc[1] / s, my = acc[2] / s;
		const double vx = (acc[3] / s) - (mx * mx);
		const double vy = (acc[4] / s) - (my * my);
		*cx = mx;
		*cy = my;
		*wx = 2 * sqrt(vx > 0 ? vx : 0);
		*wy = 2 * sqrt(vy > 0 ? vy : 0);
	}
	else {
		*cx = *cy = *wx = *wy = 0;
	}
}

// 各モードへのパワー占有率 eta_m = |<phi_m, E>|^2 / (||phi_m||^2 * ||E||^2)。
// CPU 版 sol/solve_bpm.cpp の computeOverlap と同一の式 (ホスト側で計算する:
// モード形状はホストにしかなく、1 ステップあたり nModes*N の演算で済むため)。
// getE(i) は格子点 i の複素電界を返す。結果は out[m*stride] へ格納する。
template <typename F>
static void computeOverlap(int nModes, long N, const std::complex<double> *modes,
                           F getE, double *out, long stride)
{
	double pe = 0;
	for (long i = 0; i < N; i++) pe += std::norm(getE(i));
	for (int m = 0; m < nModes; m++) {
		const std::complex<double> *phi = &modes[(size_t)m * N];
		std::complex<double> ip = 0;
		double pm = 0;
		for (long i = 0; i < N; i++) {
			ip += std::conj(phi[i]) * getE(i);
			pm += std::norm(phi[i]);
		}
		out[(size_t)m * stride] = ((pe > 0) && (pm > 0)) ? (std::norm(ip) / (pm * pe)) : 0.0;
	}
}

// ============================================================
// BPM 用のセル中心材料 ID (CPU 版 sol/solve_bpm.cpp と同一)
//   FDTD の iEx 等は Yee 千鳥格子 (Ex は (Xc, Yn, Zn)) でサンプリングされ、
//   BPM の界の座標 (Xc, Yc) に対して材料分布だけが y に半セルずれる。
//   setupId() と同じ描画順・既定値 (id 0)・eps で、セル中心
//   (Xc, Yc, Zc) に対して直接ラスタライズする。
// ============================================================
static id_t bpmIdPoint(double x, double y, double z)
{
	const double geps = EPS * sqrt(
		((Xn[Nx] - Xn[0]) * (Xn[Nx] - Xn[0])) +
		((Yn[Ny] - Yn[0]) * (Yn[Ny] - Yn[0])) +
		((Zn[Nz] - Zn[0]) * (Zn[Nz] - Zn[0])));
	id_t id = 0;
	for (int64_t n = 0; n < NGeometry; n++) {
		if (ingeometry(x, y, z, Geometry[n].shape, Geometry[n].g, geps)) {
			id = Geometry[n].m;
		}
	}
	return id;
}

// スライス iz のセル中心材料 ID を ids[ix + Nx*iy] へ格納する
static void bpmIdSlice(id_t *ids, int64_t iz)
{
	for (int iy = 0; iy < Ny; iy++) {
		for (int ix = 0; ix < Nx; ix++) {
			ids[ix + ((long)iy * Nx)] = bpmIdPoint(Xc[ix], Yc[iy], Zc[iz]);
		}
	}
}

// ============================================================
// 屈折率分布の直接入力 (ripfile = <path>)
//   geometry プリミティブの代わりに、任意の屈折率分布を与える。
//   上流 BPM-Matlab の「RIP を行列で与える」使い方に対応するもの。
//   - .csv        : 2D (Ny 行 x Nx 列, 実数 n)。全スライス共通 (z 不変)
//   - .h5 / .hdf5 : データセット /rip/n。2D (Ny x Nx) または 3D (Nz x Ny x Nx)
//   格納は行優先 (行 iy = 0 が y 最小)。寸法はメッシュと一致すること。
//   ripfile は屈折率のみを置き換える。TPA の β は従来どおり geometry の
//   材料 ID から決まる (両者は併用できる)。
// ============================================================
static float *ripN  = NULL;   // 実数 n。2D: Nx*Ny、3D: Nx*Ny*Nz (行優先)
static int    ripNz = 0;      // 0 = 未使用, 1 = 2D (z 不変), >1 = 3D のスライス数

// 入力 .ofd のディレクトリ基準でも解決して開く
static FILE *ripOpen(const char *path)
{
	FILE *fp = fopen(path, "rb");
	if ((fp == NULL) && (path[0] != '/') && (InputPath[0] != '\0')) {
		char buf[2048];
		strcpy(buf, InputPath);
		char *sl = strrchr(buf, '/');
		char *bs = strrchr(buf, '\\');
		char *cut = (bs > sl) ? bs : sl;
		if (cut != NULL) {
			cut[1] = '\0';
			strcat(buf, path);
			fp = fopen(buf, "rb");
		}
	}
	return fp;
}

// BPM.ripfile を読み込む (未指定なら何もしない)。失敗時はメッセージを出して終了
static void ripLoad(void)
{
	if (BPM.ripfile[0] == '\0') return;
	const char *path = BPM.ripfile;
	const char *ext = strrchr(path, '.');
	const long NN2 = (long)Nx * Ny;

	if ((ext != NULL) && ((strcmp(ext, ".h5") == 0) || (strcmp(ext, ".hdf5") == 0))) {
		// HDF5 : /rip/n (2D Ny x Nx または 3D Nz x Ny x Nx)
		hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
		if ((file < 0) && (path[0] != '/') && (InputPath[0] != '\0')) {
			char buf[2048];
			strcpy(buf, InputPath);
			char *sl = strrchr(buf, '/');
			char *bs = strrchr(buf, '\\');
			char *cut = (bs > sl) ? bs : sl;
			if (cut != NULL) {
				cut[1] = '\0';
				strcat(buf, path);
				file = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
			}
		}
		if (file < 0) {
			fprintf(stderr, "*** ripfile : cannot open \"%s\"\n", path);
			exit(1);
		}
		hid_t dset = H5Dopen(file, "/rip/n", H5P_DEFAULT);
		if (dset < 0) {
			fprintf(stderr, "*** ripfile : dataset /rip/n not found in \"%s\"\n", path);
			exit(1);
		}
		hid_t space = H5Dget_space(dset);
		const int nd = H5Sget_simple_extent_ndims(space);
		hsize_t dims[3] = {0, 0, 0};
		H5Sget_simple_extent_dims(space, dims, NULL);
		if (nd == 2) {
			if (((int)dims[0] != Ny) || ((int)dims[1] != Nx)) {
				fprintf(stderr, "*** ripfile : /rip/n is %dx%d, expected Ny x Nx = %dx%d\n",
						(int)dims[0], (int)dims[1], Ny, Nx);
				exit(1);
			}
			ripNz = 1;
		}
		else if (nd == 3) {
			if (((int)dims[0] != Nz) || ((int)dims[1] != Ny) || ((int)dims[2] != Nx)) {
				fprintf(stderr, "*** ripfile : /rip/n is %dx%dx%d, expected Nz x Ny x Nx = %dx%dx%d\n",
						(int)dims[0], (int)dims[1], (int)dims[2], Nz, Ny, Nx);
				exit(1);
			}
			ripNz = Nz;
		}
		else {
			fprintf(stderr, "*** ripfile : /rip/n must be 2D or 3D (got %dD)\n", nd);
			exit(1);
		}
		ripN = (float *)malloc((size_t)NN2 * ripNz * sizeof(float));
		if (H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ripN) < 0) {
			fprintf(stderr, "*** ripfile : read error \"%s\"\n", path);
			exit(1);
		}
		H5Sclose(space);
		H5Dclose(dset);
		H5Fclose(file);
	}
	else {
		// CSV : Ny 行 x Nx 列 (区切りはカンマ / 空白 / タブ)。要素数の厳密一致を要求
		FILE *fp = ripOpen(path);
		if (fp == NULL) {
			fprintf(stderr, "*** ripfile : cannot open \"%s\"\n", path);
			exit(1);
		}
		ripN = (float *)malloc((size_t)NN2 * sizeof(float));
		long cnt = 0;
		char line[65536];
		while (fgets(line, sizeof(line), fp) != NULL) {
			if (line[0] == '#') continue;
			char *q = line;
			char *qe;
			for (;;) {
				while ((*q == ',') || (*q == ' ') || (*q == '\t')) q++;
				const double v = strtod(q, &qe);
				if (qe == q) break;
				if (cnt >= NN2) {
					fprintf(stderr, "*** ripfile : too many values in \"%s\" (expected %ld = Nx x Ny)\n",
							path, NN2);
					exit(1);
				}
				ripN[cnt++] = (float)v;
				q = qe;
			}
		}
		fclose(fp);
		if (cnt != NN2) {
			fprintf(stderr, "*** ripfile : %ld values in \"%s\", expected %ld (Nx x Ny = %d x %d)\n",
					cnt, path, NN2, Nx, Ny);
			exit(1);
		}
		ripNz = 1;
	}
}

// セル中心の複素屈折率 : ripfile があればそれを、無ければ材料 ID から返す
// (n_mat は損失を +imag で保持する規約。rip は実数 n のみ = 無損失)
static floatcomplex bpmNPoint(const floatcomplex *n_mat, int ix, int iy, int64_t iz)
{
	if (ripN != NULL) {
		const size_t off = (ripNz > 1) ? ((size_t)iz * Nx * Ny) : 0;
		floatcomplex v = {ripN[off + ix + ((size_t)iy * Nx)], 0.0f};
		return v;
	}
	return n_mat[bpmIdPoint(Xc[ix], Yc[iy], Zc[iz])];
}

// スライス iz のセル中心複素屈折率を dst[ix + Nx*iy] へ格納する
static void bpmNSlice(floatcomplex *dst, const floatcomplex *n_mat, int64_t iz)
{
	if (ripN != NULL) {
		const size_t off = (ripNz > 1) ? ((size_t)iz * Nx * Ny) : 0;
		for (long i = 0; i < (long)Nx * Ny; i++) {
			floatcomplex v = {ripN[off + i], 0.0f};
			dst[i] = v;
		}
		return;
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int iy = 0; iy < Ny; iy++) {
		for (int ix = 0; ix < Nx; ix++) {
			dst[ix + ((long)iy * Nx)] = n_mat[bpmIdPoint(Xc[ix], Yc[iy], Zc[iz])];
		}
	}
}


// ============================================================
// PML (複素座標伸長) の係数構築 (CPU 版 sol/solve_bpm.cpp と同一)
//
// 位相規約 exp(-i k z) に対し s(x) = 1 - i*sigma(x) と伸長すると、
// 外向きの横方向波は PML 層内で |E| ~ exp(-|kx| ∫sigma dx) で減衰する。
// sigma(d) = smax*(d/W)^m (m = 3), smax は「かすめ入射 (kx = k0*n0) の
// 往復振幅反射率が R0」から smax = -(m+1)*ln(R0)/(2*k0*n0*W) と決める。
// 対称境界側は物理境界ではないので PML を置かない。
// ============================================================
static double bpmPmlSigma(double u, double lo, double hi, double Wpml, double smax, int symLow)
{
	const double m = 3.0;
	double d = 0;
	if (!symLow && (u < lo)) d = lo - u;
	else if (u > hi)         d = u - hi;
	if (d <= 0) return 0;
	if (d > Wpml) d = Wpml;
	return smax * pow(d / Wpml, m);
}

static void bpmPmlAxis(const double *node, const double *cent, int N,
                       double Wpml, double smax, int symLow,
                       std::complex<double> *gm, std::complex<double> *gp)
{
	const double lo = node[0] + Wpml;
	const double hi = node[N] - Wpml;
	for (int i = 0; i < N; i++) {
		const std::complex<double> sc(1.0, -bpmPmlSigma(cent[i],     lo, hi, Wpml, smax, symLow));
		const std::complex<double> sm(1.0, -bpmPmlSigma(node[i],     lo, hi, Wpml, smax, symLow));
		const std::complex<double> sp(1.0, -bpmPmlSigma(node[i + 1], lo, hi, Wpml, smax, symLow));
		gm[i] = 1.0 / (sc * sm);
		gp[i] = 1.0 / (sc * sp);
	}
}

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

	// BPM で無効な入力キーワードの警告 (パースはされるが伝搬には反映されない)
	{
	    struct { int cond; const char *name; } ignored[] = {
	        {IPlanewave != 0,          "planewave"},
	        {NPoint > 0,               "point (S パラメータ出力はゼロになる)"},
	        {NInductor > 0,            "load"},
	        {rFeed != 0,               "rfeed"},
	        {iABC == 1,                "abc = 1 (PML; BPM は独自の端部吸収体を使用)"},
	        {PBCx || PBCy || PBCz,     "pbc"},
	        // 波長掃引は CPU 版 (obpm) のみ実装。CUDA 版は先頭波長のみ計算する
	        {BPM.wlsweep != 0,         "wlsweep (CUDA 版は未対応 : CPU 版 obpm を使用してください)"},
	    };
	    for (size_t n = 0; n < sizeof(ignored) / sizeof(ignored[0]); n++) {
	        if (ignored[n].cond) {
	            sprintf(str, "*** warning : keyword '%s' is ignored by the BPM solver.", ignored[n].name);
	            if (io) fprintf(fp, "%s\n", str);
	            fprintf(stdout, "%s\n", str);
	        }
	    }
	    // tpa / powersweep は CUDA 版でも実装済み (無視されない) だが、
	    // GPU 実機での実行検証が未実施のため、その旨だけを伝える。
	    // 「無視される」と誤って警告すると正しい結果まで捨てられてしまう。
	    if ((NTpaB > 0) || (PowerSweep.npoints > 0)) {
	        sprintf(str, "*** note : tpa / powersweep are implemented in the CUDA solver"
	                     " but not yet validated on real GPU hardware."
	                     " Cross-check against the CPU solver (obpm) if accuracy matters.");
	        if (io) fprintf(fp, "%s\n", str);
	        fprintf(stdout, "%s\n", str);
	    }
	    // CUDA 版は波長掃引に未対応 (先頭周波数のみ使用)。
	    // wlsweep 指定時は上の一覧で専用の warning を出しているため二重に出さない。
	    if (!BPM.wlsweep && ((NFreq2 > 1) || ((NFreq2 == 0) && (NFreq1 > 1)))) {
	        sprintf(str, "*** warning : multiple frequencies are not swept by the BPM solver (only the first is used).");
	        if (io) fprintf(fp, "%s\n", str);
	        fprintf(stdout, "%s\n", str);
	    }
	}


	// 材料毎の複素屈折率テーブル : n = sqrt(epsr - i*sigma/(omega*eps0))
	// 損失を imag(n) > 0 として格納する (applyMultiplier の exp(d*imag(n)), d < 0 の規約)
	floatcomplex *n_mat = (floatcomplex *)malloc(NMaterial * sizeof(floatcomplex));

	// セル中心材料 ID のスライスバッファ (CPU 版と同一の役割)
	floatcomplex *bpmNIn = (floatcomplex *)malloc((size_t)Nx * Ny * sizeof(floatcomplex));
	floatcomplex *bpmNZ  = (floatcomplex *)malloc((size_t)Nx * Ny * sizeof(floatcomplex));
	id_t *bpmIdZ  = (id_t *)malloc((size_t)Nx * Ny * sizeof(id_t));

	// ripfile = <path> 指定時 : 屈折率分布を読み込む (CPU 版と同一)
	ripLoad();
	if (ripN != NULL) {
		sprintf(str, "ripfile : %s (%s)", BPM.ripfile,
		        (ripNz > 1) ? "3D, z-varying" : "2D, z-invariant");
		if (io) fprintf(fp, "%s\n", str);
		fprintf(stdout, "%s\n", str);
	}
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
		P->n_0 = CREALF(bpmNPoint(n_mat, Nx / 2, Ny / 2, 0));
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
		bpmNSlice(&P->n_in[(size_t)P->Nx_n * P->Ny_n * iz], n_mat, iz);
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

	// --- PML (pml = <width> [<R0>]) の係数構築 ---
	// 未指定なら全ポインタが NULL のままで、カーネルは従来の ax/ay を使う
	const int usePml = (BPM.pmlw > 0);
	P->pmlxm = P->pmlxp = P->pmlym = P->pmlyp = NULL;
	std::complex<double> *pmlgxm = NULL, *pmlgxp = NULL, *pmlgym = NULL, *pmlgyp = NULL;
	if (usePml) {
		const double m_grade = 3.0;
		const double Lx_dom = Xn[Nx] - Xn[0];
		const double Ly_dom = Yn[Ny] - Yn[0];
		if ((2 * BPM.pmlw >= Lx_dom) || (2 * BPM.pmlw >= Ly_dom)) {
			sprintf(str, "*** warning : pml width %.3e [m] is too large for the domain"
			             " (%.3e x %.3e [m]).", BPM.pmlw, Lx_dom, Ly_dom);
			if (io) fprintf(fp, "%s\n", str);
			fprintf(stdout, "%s\n", str);
		}
		const double smax = -(m_grade + 1) * log(BPM.pmlR0) / (2 * k0 * P->n_0 * BPM.pmlw);
		pmlgxm = new std::complex<double>[Nx];
		pmlgxp = new std::complex<double>[Nx];
		pmlgym = new std::complex<double>[Ny];
		pmlgyp = new std::complex<double>[Ny];
		bpmPmlAxis(Xn, Xc, Nx, BPM.pmlw, smax, BPM.symx, pmlgxm, pmlgxp);
		bpmPmlAxis(Yn, Yc, Ny, BPM.pmlw, smax, BPM.symy, pmlgym, pmlgyp);
		P->pmlxm = (floatcomplex *)malloc((size_t)Nx * sizeof(floatcomplex));
		P->pmlxp = (floatcomplex *)malloc((size_t)Nx * sizeof(floatcomplex));
		P->pmlym = (floatcomplex *)malloc((size_t)Ny * sizeof(floatcomplex));
		P->pmlyp = (floatcomplex *)malloc((size_t)Ny * sizeof(floatcomplex));
		for (int ix = 0; ix < Nx; ix++) {
			floatcomplex vm = {(float)pmlgxm[ix].real(), (float)pmlgxm[ix].imag()};
			floatcomplex vp = {(float)pmlgxp[ix].real(), (float)pmlgxp[ix].imag()};
			P->pmlxm[ix] = vm;
			P->pmlxp[ix] = vp;
		}
		for (int iy = 0; iy < Ny; iy++) {
			floatcomplex vm = {(float)pmlgym[iy].real(), (float)pmlgym[iy].imag()};
			floatcomplex vp = {(float)pmlgyp[iy].real(), (float)pmlgyp[iy].imag()};
			P->pmlym[iy] = vm;
			P->pmlyp[iy] = vp;
		}
		const int ncx = (int)(BPM.pmlw / P->dx + 0.5);
		const int ncy = (int)(BPM.pmlw / P->dy + 0.5);
		sprintf(str, "pml : width = %.3e [m] (%d x %d cells), R0 = %.1e,"
		             " sigma_max = %.3f (legacy absorber disabled)",
		        BPM.pmlw, ncx, ncy, BPM.pmlR0, smax);
		if (io) fprintf(fp, "%s\n", str);
		fprintf(stdout, "%s\n", str);
	}

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
				// PML 使用時は振幅吸収体を掛けない (境界吸収は PML が担当する)
				P->multiplier[i] = usePml ? 1.0f : (float)exp(-P->dz * dist * dist * alpha);
			}
		}
		P->precisePower = power;
	}

	// モード励振 (launch = mode <m>) : 先頭スライスの屈折率分布から導波モードを
	// 虚軸伝搬法 (bpm/modes.cpp) で求め、ガウシアンの代わりに励振する。
	// 振幅は feed 電圧 (モードは L2 正規化済みのため入力電力 = volt^2)。
	// beamtilt 指定時はガウシアンと同じ位相ランプをモード界にも掛ける。
	double mode_neff = 0;   // HDF5 /metadata/mode_neff (モード励振時のみ > 0)
	if (BPM.launchMode >= 0) {
		const long Nm = (long)Nx * Ny;
		std::complex<double> *n2m = new std::complex<double>[Nm];
		double nmax = 0;
		bpmNSlice(bpmNZ, n_mat, 0);
		for (int iy = 0; iy < Ny; iy++) {
		    for (int ix = 0; ix < Nx; ix++) {
				const floatcomplex nmv = bpmNZ[(long)Nx * iy + ix];
				std::complex<double> nn(CREALF(nmv), -CIMAGF(nmv));
				n2m[(long)Nx * iy + ix] = nn * nn;
				if (nn.real() > nmax) nmax = nn.real();
		    }
		}
		int launched = 0;
		wabpm_params Wm;
		Wm.Nx = Nx;
		Wm.Ny = Ny;
		Wm.dx = Dx;
		Wm.dy = Dy;
		Wm.k0 = k0;
		Wm.n0 = P->n_0;
		Wm.wideangle = 0;
		Wm.pol = BPM.pol;
		Wm.symx = BPM.symx;   // 対称境界はモード解析にも適用する (未設定だと未初期化)
		Wm.symy = BPM.symy;
		// 導波判定のしきい値 : 境界リングの n^2 (bpm/modes.cpp と同じ値)。
		// 参照屈折率 n0 と比較すると、GRIN のように n0 = 軸上最大値の入力で
		// 導波構造があるのに「モード無し」と誤判定する
		const double nthr2 = wabpm_guided_threshold(&Wm, n2m);
		if (((double)nmax * nmax) > nthr2) {
		    /* 重ね合わせる最大モード番号まで求める */
		    int m = 0;
		    for (int n = 0; n < BPM.launchNModes; n++) {
		        if (BPM.launchIdx[n] > m) m = BPM.launchIdx[n];
		    }
		    // a*mu_max ~ 0.5 となる虚軸ステップ幅 (tests/test_modes.cpp と同じ)。
		    // スペクトル幅は「最大屈折率 - 導波しきい値」で見積もる
		    const double mu_max = k0 * k0 * (((double)nmax * nmax) - nthr2);
		    Wm.dz = 2 * k0 * P->n_0 / mu_max;
		    std::complex<double> *modes = new std::complex<double>[(size_t)(m + 1) * Nm];
		    double *neff = new double[m + 1];
		    const int nFound = wabpm_find_modes(&Wm, n2m, m + 1, 20000, 1e-10, modes, neff);
		    if (nFound >= m + 1) {
				const double volt = (NFeed > 0) ? Feed[0].volt : 1;
				// 入射ビームの傾き (beamtilt =) : ガウシアンと同じ横方向波数の
				// 位相ランプをモード界にも適用する (位相基準はビーム中心)
				const double ktx = k0 * P->n_0 * sin(BPM.tiltx * DTOR);
				const double kty = k0 * P->n_0 * sin(BPM.tilty * DTOR);
				// モードの重ね合わせ s = sum_k coef_k * mode_{idx_k} を作り、
				// L2 ノルムで正規化してから feed 電圧を掛ける。これにより
				// モード数や係数に依らず入力電力は volt^2 に保たれる
				// (係数は分岐比のみを決める)。
				std::complex<double> *sup = new std::complex<double>[Nm];
				for (long i = 0; i < Nm; i++) sup[i] = 0.0;
				for (int n = 0; n < BPM.launchNModes; n++) {
				    const std::complex<double> *src = &modes[(size_t)BPM.launchIdx[n] * Nm];
				    const double c = BPM.launchCoef[n];
				    for (long i = 0; i < Nm; i++) sup[i] += c * src[i];
				}
				double nrm2 = 0;
				for (long i = 0; i < Nm; i++) nrm2 += std::norm(sup[i]);
				const double scale = (nrm2 > 0) ? (volt / sqrt(nrm2)) : 0;

				double power = 0;
				for (long i = 0; i < Nm; i++) {
				    const long ix = i % Nx;
				    const long iy = i / Nx;
				    const double ph = -(ktx * (Xc[ix] - xc0)) - (kty * (Yc[iy] - yc0));
				    const std::complex<double> em = scale * sup[i]
				                                  * std::complex<double>(cos(ph), sin(ph));
				    floatcomplex e = {(float)em.real(), (float)em.imag()};
				    P->E1[i] = e;
				    power += std::norm(em);
				}
				delete[] sup;
				P->precisePower = power;
				mode_neff = neff[BPM.launchIdx[0]];   /* 代表値 : 先頭に指定したモード */
				if (BPM.launchNModes == 1) {
				    sprintf(str, "launch : mode %d, neff = %.6f (found %d)",
				            BPM.launchIdx[0], neff[BPM.launchIdx[0]], nFound);
				    if (io) fprintf(fp, "%s\n", str);
				    fprintf(stdout, "%s\n", str);
				}
				else {
				    int len = sprintf(str, "launch : superposition of %d modes (found %d) :",
				                      BPM.launchNModes, nFound);
				    for (int n = 0; n < BPM.launchNModes; n++) {
				        len += sprintf(str + len, " %d(c=%g, neff=%.6f)",
				                       BPM.launchIdx[n], BPM.launchCoef[n],
				                       neff[BPM.launchIdx[n]]);
				    }
				    if (io) fprintf(fp, "%s\n", str);
				    fprintf(stdout, "%s\n", str);
				}
				if ((BPM.tiltx != 0) || (BPM.tilty != 0)) {
				    sprintf(str, "launch : beamtilt (%.3f, %.3f) [deg] applied to the mode field",
				            BPM.tiltx, BPM.tilty);
				    if (io) fprintf(fp, "%s\n", str);
				    fprintf(stdout, "%s\n", str);
				}
				launched = 1;
		    }
		    delete[] modes;
		    delete[] neff;
		}
		if (!launched) {
		    sprintf(str, "*** warning : launch = mode %d: mode not found (no guided structure or not converged). Falling back to Gaussian.", BPM.launchMode);
		    if (io) fprintf(fp, "%s\n", str);
		    fprintf(stdout, "%s\n", str);
		}
		delete[] n2m;
	}

	// ============================================================
	// モード解析 (modes = <nModes> [excite]) : CPU 版と同一のホスト側処理
	//   (bpm/modes.cpp の虚軸伝搬法。曲げは反映せず、スカラー演算子で解析)
	//   launch = mode と併用した場合、励振は launch を優先し excite は無視する
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
		bpmNSlice(bpmNZ, n_mat, P->iz_start);
		for (int iy = 0; iy < P->Ny; iy++) {
			for (int ix = 0; ix < P->Nx; ix++) {
				const floatcomplex nm = bpmNZ[ix + (long)iy * P->Nx];
				const double nr = CREALF(nm);
				n2m[ix + (long)iy * P->Nx] = std::complex<double>(nr * nr, 0.0);
				if (nr > nmax) nmax = nr;
			}
		}
		// 虚軸ステップ幅 : 最大固有値 mu_max に対し a*mu_max ~ 0.5 となるよう選ぶ。
		// スペクトル幅は「最大屈折率 - 導波しきい値 (境界リングの n)」で見積もる
		const double nthr2m = wabpm_guided_threshold(&Wm, n2m);
		const double mu_max = k0 * k0 * (nmax * nmax - nthr2m);
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
		if (BPM.modeExcite && (BPM.launchMode >= 0)) {
			sprintf(str, "*** warning : modes = ... excite is ignored because launch = mode is specified.");
			if (io) fprintf(fp, "%s\n", str);
			fprintf(stdout, "%s\n", str);
		}
		if (BPM.modeExcite && (BPM.launchMode < 0) && (nModesFound > 0)) {
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

	// 伝搬の可視化用 (追加) : 中心列 (x = Nx/2) の強度と z ごとのスカラー推移
	const long ntr = P->iz_end - P->iz_start;
	float  *Iyz     = (float *)malloc((size_t)P->Ny * ntr * sizeof(float));
	double *trZ     = (double *)malloc((size_t)ntr * sizeof(double));
	double *trPower = (double *)malloc((size_t)ntr * sizeof(double));
	double *trPeak  = (double *)malloc((size_t)ntr * sizeof(double));
	double *trCx    = (double *)malloc((size_t)ntr * sizeof(double));
	double *trCy    = (double *)malloc((size_t)ntr * sizeof(double));
	double *trWx    = (double *)malloc((size_t)ntr * sizeof(double));
	double *trWy    = (double *)malloc((size_t)ntr * sizeof(double));
	// モード結合率 (modes 指定時のみ) : eta_m(z) = |<phi_m, E>|^2 / (||phi_m||^2 ||E||^2)
	// 格納は行優先 [m*ntr + t] (HDF5 では nModesFound x ntr の 2 次元配列)
	double *trOverlap = (nModesFound > 0)
		? (double *)malloc((size_t)nModesFound * ntr * sizeof(double)) : NULL;
	floatcomplex *Ecol = (floatcomplex *)malloc(P->Ny * sizeof(floatcomplex));
	// 近軸パスのモード結合率用 : 全電界のホスト受け皿 (modes 指定時のみ確保)
	floatcomplex *Eovl = (trOverlap != NULL)
		? (floatcomplex *)malloc((size_t)P->Nx * P->Ny * sizeof(floatcomplex)) : NULL;
	// 統計量カーネル用 : セル中心座標と集計バッファをデバイスへ
	double *d_xc = NULL, *d_yc = NULL, *d_acc = NULL;
	double accHost[6];
	CUDA_CHECK(cudaMalloc(&d_xc, (size_t)P->Nx * sizeof(double)));
	CUDA_CHECK(cudaMalloc(&d_yc, (size_t)P->Ny * sizeof(double)));
	CUDA_CHECK(cudaMalloc(&d_acc, 6 * sizeof(double)));
	CUDA_CHECK(cudaMemcpy(d_xc, Xc, (size_t)P->Nx * sizeof(double), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d_yc, Yc, (size_t)P->Ny * sizeof(double), cudaMemcpyHostToDevice));

	// |E(x,y)|^2 スナップショット (frames = interval 指定時のみ, /field/frames へ出力)
	const int  frameInterval = BPM.frames;
	const long nframes = (frameInterval > 0)
	    ? ((P->iz_end - P->iz_start - 1) / frameInterval + 1) : 0;
	float *Frames = (nframes > 0)
	    ? (float *)malloc((size_t)nframes * P->Nx * P->Ny * sizeof(float)) : NULL;
	// 複素電界のスナップショット (frames = <interval> complex 指定時のみ。位相表示用)
	const int framesCplx = (nframes > 0) && BPM.framesComplex;
	float *FramesR = framesCplx ? (float *)malloc((size_t)nframes * P->Nx * P->Ny * sizeof(float)) : NULL;
	float *FramesI = framesCplx ? (float *)malloc((size_t)nframes * P->Nx * P->Ny * sizeof(float)) : NULL;
	if (nframes > 0) {
		sprintf(str, "frames : interval = %d steps, count = %ld%s", frameInterval, nframes,
		        BPM.framesComplex ? " (complex field also recorded)" : "");
		if (io) fprintf(fp, "%s\n", str);
		fprintf(stdout, "%s\n", str);
	}

	// HDF5ファイルの作成
	file_id = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// ルート属性 : このファイルが BPM (定常解法) の出力であり、系列の独立変数が
	// 時間ではなく z であることを宣言する (CPU 版 sol/solve_bpm.cpp と同一)。
	{
		hid_t root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
		h5AttrStr(root_id, "solver", "OpenBPM (beam propagation method)");
		h5AttrStr(root_id, "domain", "steady-state (single frequency, no time axis)");
		h5AttrStr(root_id, "marching_axis", "z");
		h5AttrStr(root_id, "marching_axis_values", "/trace/z");
		h5AttrLong(root_id, "marching_steps", (long)ntr);
		h5AttrLong(root_id, "time_dependent", 0);
		H5Gclose(root_id);
	}

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
		// PML (複素座標伸長) : 近軸パスと同じ係数を拡張パスにも渡す
		W.gxm = pmlgxm;
		W.gxp = pmlgxp;
		W.gym = pmlgym;
		W.gyp = pmlgyp;

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

		// テーパ/ツイスト : 座標変換の参照元となる入力断面の材料 ID を確定する
		if (xform) {
			bpmNSlice(bpmNIn, n_mat, P->iz_start);
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
			bpmNSlice(bpmNZ, n_mat, iz);            // 現スライスの屈折率 (n2 が参照)
			if (haveTpa) bpmIdSlice(bpmIdZ, iz);    // TPA の β は材料 ID から決まる
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
						const floatcomplex n00 = bpmNIn[i0     + ((long)j0       * Nx)];
						const floatcomplex n10 = bpmNIn[i0 + 1 + ((long)j0       * Nx)];
						const floatcomplex n01 = bpmNIn[i0     + ((long)(j0 + 1) * Nx)];
						const floatcomplex n11 = bpmNIn[i0 + 1 + ((long)(j0 + 1) * Nx)];
						nr    = ((1 - tx) * (1 - ty) * CREALF(n00)) + (tx * (1 - ty) * CREALF(n10))
						      + ((1 - tx) * ty * CREALF(n01))       + (tx * ty * CREALF(n11));
						nimag = ((1 - tx) * (1 - ty) * CIMAGF(n00)) + (tx * (1 - ty) * CIMAGF(n10))
						      + ((1 - tx) * ty * CIMAGF(n01))       + (tx * ty * CIMAGF(n11));
					}
					else {
						const floatcomplex nm = bpmNZ[index];
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
						betaSlice[(long)(Nx * iy) + ix] = (float)tpa_mat[bpmIdZ[(long)(Nx * iy) + ix]];
					}
				}
				wabpm_gpu_tpa(G, betaSlice, Dz);
			}

			// 中心行 / 中心列の強度と z ごとのスカラー推移を記録
			{
				const long t = iz - P->iz_start;
				wabpm_gpu_get_row(G, P->Ny / 2, Ed);  // Ed の先頭 Nx 要素を一時利用
				for (long ix = 0; ix < Nx; ix++) {
					Ixz[t * Nx + ix] = (float)thrust::norm(Ed[ix]);
				}
				wabpm_gpu_get_field(G, Ed);
				const long col0 = Nx / 2;
				for (long iy = 0; iy < Ny; iy++) {
					Iyz[t * Ny + iy] = (float)thrust::norm(Ed[col0 + iy * Nx]);
				}
				trZ[t] = Zn[0] + ((iz + 1) * Dz);
				wabpm_gpu_trace(G, d_xc, d_yc, d_acc, accHost);
				traceFromAcc(accHost, Dx * Dy,
				             &trPower[t], &trPeak[t], &trCx[t], &trCy[t], &trWx[t], &trWy[t]);
				// モード結合率 (全電界は Iyz 用に取得済みの Ed を再利用)
				if (trOverlap) {
					computeOverlap(nModesFound, (long)Nx * Ny, modeFields,
					               [&](long i) {
					                   return std::complex<double>(Ed[i].real(), Ed[i].imag());
					               },
					               &trOverlap[t], ntr);
				}
			}

			// スナップショットを記録 (全電界をホストへ取得)
			if (Frames && ((iz - P->iz_start) % frameInterval == 0)) {
				wabpm_gpu_get_field(G, Ed);
				const size_t off = (size_t)((iz - P->iz_start) / frameInterval) * Nx * Ny;
				float *f = &Frames[off];
				for (long i = 0; i < (long)Nx * Ny; i++) {
					f[i] = (float)thrust::norm(Ed[i]);
					if (framesCplx) {
						FramesR[off + i] = (float)Ed[i].real();
						FramesI[off + i] = (float)Ed[i].imag();
					}
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
		bpmNSlice(bpmNZ, n_mat, P->iz_end - 1);
		for (int iy = 0; iy < Ny; iy++) {
			for (int ix = 0; ix < Nx; ix++) {
				long index = (long)(Nx * iy) + ix;
				P->n_out[index] = bpmNZ[index];
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
			bpmIdSlice(bpmIdZ, iz);
			for (int iy = 0; iy < P->Ny; iy++) {
				for (int ix = 0; ix < P->Nx; ix++) {
					betaSlice[(long)(P->Nx * iy) + ix] = (float)tpa_mat[bpmIdZ[(long)(P->Nx * iy) + ix]];
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
			const long t = iz - P->iz_start;
			const long row0 = (long)(P->Ny / 2) * P->Nx;
			CUDA_CHECK(cudaMemcpy(Erow, P_now.E2 + row0, P->Nx * sizeof(floatcomplex), cudaMemcpyDeviceToHost));
			for (long ix = 0; ix < P->Nx; ix++) {
				Ixz[t * P->Nx + ix] =
					(CREALF(Erow[ix]) * CREALF(Erow[ix])) + (CIMAGF(Erow[ix]) * CIMAGF(Erow[ix]));
			}
			// 中心列 (x = Nx/2) : stride コピーのため 1 要素ずつ取得する
			CUDA_CHECK(cudaMemcpy2D(Ecol, sizeof(floatcomplex),
			                        P_now.E2 + (P->Nx / 2), P->Nx * sizeof(floatcomplex),
			                        sizeof(floatcomplex), P->Ny, cudaMemcpyDeviceToHost));
			for (long iy = 0; iy < P->Ny; iy++) {
				Iyz[t * P->Ny + iy] =
					(CREALF(Ecol[iy]) * CREALF(Ecol[iy])) + (CIMAGF(Ecol[iy]) * CIMAGF(Ecol[iy]));
			}
			// z ごとのスカラー推移 (デバイス側で集計)
			trZ[t] = Zn[0] + ((iz + 1) * Dz);
			CUDA_CHECK(cudaMemset(d_acc, 0, 6 * sizeof(double)));
			fieldTrace<<<tpaNb, tpaTpb>>>(P_dev, d_xc, d_yc, d_acc);
			CUDA_CHECK(cudaDeviceSynchronize());
			CUDA_CHECK(cudaMemcpy(accHost, d_acc, 6 * sizeof(double), cudaMemcpyDeviceToHost));
			traceFromAcc(accHost, Dx * Dy,
			             &trPower[t], &trPeak[t], &trCx[t], &trCy[t], &trWx[t], &trWy[t]);
			// モード結合率 : モード形状はホストにしかないため全電界を転送する。
			// modes 指定時のみのコストで、未指定時は従来通り行/列だけの転送。
			if (trOverlap) {
				CUDA_CHECK(cudaMemcpy(Eovl, P_now.E2,
				                      (size_t)P->Nx * P->Ny * sizeof(floatcomplex),
				                      cudaMemcpyDeviceToHost));
				computeOverlap(nModesFound, (long)P->Nx * P->Ny, modeFields,
				               [&](long i) {
				                   return std::complex<double>(CREALF(Eovl[i]), CIMAGF(Eovl[i]));
				               },
				               &trOverlap[t], ntr);
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


	// モード結合率の要約 : 出口断面 (z = z_end) での各モードのパワー占有率をログへ
	// (CPU 版 sol/solve_bpm.cpp と同一の出力)
	if (trOverlap && (ntr > 0)) {
		for (int m = 0; m < nModesFound; m++) {
			sprintf(str, "modes : overlap eta_%d = %.6f (z = z_end)",
			        m + 1, trOverlap[((size_t)m * ntr) + (ntr - 1)]);
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
		const char *field_label[] = {
			"final electric field (real part)",
			"final electric field (imaginary part)",
			"refractive index at output plane (real part)",
			"refractive index at output plane (imaginary part)"
		};
		for (size_t n = 0; n < sizeof(field_data) / sizeof(field_data[0]); n++) {
			for (long i = 0; i < P->Nx * P->Ny; i++) {
				tmp[i] = field_data[n].imagpart ? CIMAGF(field_data[n].data[i])
												: CREALF(field_data[n].data[i]);
			}
			dataspace_id = H5Screate_simple(2, field_dims, NULL);
			dataset_id = H5Dcreate(field_group_id, field_data[n].name, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
			h5Annotate(dataset_id, field_label[n], (n < 2) ? "arb. unit" : "1",
			           "y,x", "/metadata/Yc /metadata/Xc");
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
			h5Annotate(dataset_id, "propagation map |E(x, y=Ny/2, z)|^2", "arb. unit",
			           "z,x", "/trace/z /metadata/Xc");
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}

		// 伝搬マップ |E(x=Nx/2, y, z)|^2 の書き込み (Ixz の y 版)
		{
			hsize_t iyz_dims[2] = {(hsize_t)ntr, (hsize_t)P->Ny};
			dataspace_id = H5Screate_simple(2, iyz_dims, NULL);
			dataset_id = H5Dcreate(field_group_id, "Iyz", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Iyz);
			h5Annotate(dataset_id, "propagation map |E(x=Nx/2, y, z)|^2", "arb. unit",
			           "z,y", "/trace/z /metadata/Yc");
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}

		// スナップショット |E(x,y)|^2 (nframes x Ny x Nx) の書き込み
		if (Frames) {
			hsize_t fr_dims[3] = {(hsize_t)nframes, (hsize_t)P->Ny, (hsize_t)P->Nx};
			dataspace_id = H5Screate_simple(3, fr_dims, NULL);
			dataset_id = H5Dcreate(field_group_id, "frames", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Frames);
			h5Annotate(dataset_id, "snapshot |E(x,y)|^2 along z", "arb. unit",
			           "frame,y,x", "/field/frames_z /metadata/Yc /metadata/Xc");
			h5AttrLong(dataset_id, "frame_interval_steps", (long)frameInterval);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
			// 各フレームの z 位置 (CPU 版と同一)
			{
				double *fz = (double *)malloc((size_t)nframes * sizeof(double));
				for (long f = 0; f < nframes; f++) {
					const long t = f * frameInterval;
					fz[f] = trZ[(t < ntr) ? t : (ntr - 1)];
				}
				hsize_t fz_dims[1] = {(hsize_t)nframes};
				dataspace_id = H5Screate_simple(1, fz_dims, NULL);
				dataset_id = H5Dcreate(field_group_id, "frames_z", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fz);
				h5Annotate(dataset_id, "z position of each snapshot frame", "m",
				           "frame", NULL);
				H5Dclose(dataset_id);
				H5Sclose(dataspace_id);
				free(fz);
			}
			// 複素電界 (位相表示用)
			if (framesCplx) {
				const char *cn[2] = {"frames_r", "frames_i"};
				float *cd[2] = {FramesR, FramesI};
				for (int c = 0; c < 2; c++) {
					dataspace_id = H5Screate_simple(3, fr_dims, NULL);
					dataset_id = H5Dcreate(field_group_id, cn[c], H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
					status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cd[c]);
					h5Annotate(dataset_id,
					           (c == 0) ? "snapshot E(x,y) along z (real part)"
					                    : "snapshot E(x,y) along z (imaginary part)",
					           "arb. unit", "frame,y,x", NULL);
					H5Dclose(dataset_id);
					H5Sclose(dataspace_id);
				}
			}
		}
		H5Gclose(field_group_id);
	}

	// z ごとのスカラー推移 (/trace) の書き込み : GUI での伝搬グラフ表示用
	// (CPU 版 sol/solve_bpm.cpp と同一の内容・単位)
	{
		hid_t trace_group_id = H5Gcreate(file_id, "/trace", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		struct {
			const char *name;
			const double *data;
			const char *label;
			const char *units;
		} trace_data[] = {
			{"z",          trZ,      "cross-section position along the propagation axis", "m"},
			{"power",      trPower,  "cross-section power (W when tpa/powersweep is used)", "arb. unit"},
			{"peak",       trPeak,   "peak intensity |E|^2",                    "arb. unit"},
			{"centroid_x", trCx,     "intensity centroid x",                    "m"},
			{"centroid_y", trCy,     "intensity centroid y",                    "m"},
			{"width_x",    trWx,     "beam width 2*sigma along x",              "m"},
			{"width_y",    trWy,     "beam width 2*sigma along y",              "m"}
		};
		// このグループが「z を独立変数とする系列」であることを宣言する
		h5AttrStr(trace_group_id, "description",
		          "series along the propagation axis (BPM has no time axis)");
		h5AttrStr(trace_group_id, "axis", "z");
		h5AttrStr(trace_group_id, "axis_values", "/trace/z");
		h5AttrLong(trace_group_id, "npoints", (long)ntr);
		hsize_t tr_dims[1] = {(hsize_t)ntr};
		for (size_t n = 0; n < sizeof(trace_data) / sizeof(trace_data[0]); n++) {
			dataspace_id = H5Screate_simple(1, tr_dims, NULL);
			dataset_id = H5Dcreate(trace_group_id, trace_data[n].name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, trace_data[n].data);
			h5Annotate(dataset_id, trace_data[n].label, trace_data[n].units,
			           "z", "/trace/z");
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
		// モード結合率 (nModesFound x ntr) : overlap[m][t] は z = z[t] における
		// モード m+1 のパワー占有率。/modes/neff と同じモード順序。
		if (trOverlap) {
			hsize_t ov_dims[2] = {(hsize_t)nModesFound, (hsize_t)ntr};
			dataspace_id = H5Screate_simple(2, ov_dims, NULL);
			dataset_id = H5Dcreate(trace_group_id, "overlap", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, trOverlap);
			h5Annotate(dataset_id, "power fraction in each guided mode", "1",
			           "mode,z", "/modes/neff /trace/z");
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
		H5Gclose(trace_group_id);
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
			{"grid_dz", &Dz},
			// モード励振 (launch = mode) 時の実効屈折率。ガウシアン励振では 0
			{"mode_neff", &mode_neff}
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
	if (P->pmlxm != NULL) free(P->pmlxm);
	if (P->pmlxp != NULL) free(P->pmlxp);
	if (P->pmlym != NULL) free(P->pmlym);
	if (P->pmlyp != NULL) free(P->pmlyp);
	if (pmlgxm != NULL) { delete[] pmlgxm; delete[] pmlgxp; delete[] pmlgym; delete[] pmlgyp; }
	free(n_mat);
	free(bpmNIn);
	free(bpmNZ);
	free(bpmIdZ);
	if (ripN != NULL) { free(ripN); ripN = NULL; ripNz = 0; }
	free(Ixz);
	free(Iyz);
	free(Erow);
	free(Ecol);
	free(trZ); free(trPower); free(trPeak);
	free(trCx); free(trCy); free(trWx); free(trWy);
	if (trOverlap) free(trOverlap);
	if (Eovl) free(Eovl);
	CUDA_CHECK(cudaFree(d_xc));
	CUDA_CHECK(cudaFree(d_yc));
	CUDA_CHECK(cudaFree(d_acc));
	if (Frames) free(Frames);
	if (FramesR) free(FramesR);
	if (FramesI) free(FramesI);
	free(tpa_mat);
	if (betaSlice != NULL) free(betaSlice);
	if (E0 != NULL) free(E0);
	if (sweepPin != NULL) free(sweepPin);
	if (sweepPout != NULL) free(sweepPout);
	if (modeFields) delete[] modeFields;
	if (modeNeff) delete[] modeNeff;
}
