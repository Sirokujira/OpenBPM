/*
postbpm.c

BPM 出力 (time_series_data.h5 の /field) の可視化。
  - /field/Ixz       : 伝搬マップ |E(x, y=Ny/2, z)|^2 (Nz x Nx) → 等高線ページ
  - /field/Efinal_r/i: 最終断面電界 (Ny x Nx) → |E| 等高線ページ
  - /field/frames    : |E(x,y)|^2 スナップショット (nframes x Ny x Nx)
                       → 等間隔に最大 BPM_FRAME_PAGES 枚の等高線ページ
  - /trace/*         : z 伝搬に沿ったスカラー推移 → 折れ線グラフページ
                       (パワー・ピーク強度・ビーム幅・重心・モード結合率)
FDTD 用の post が読む obpm.out には含まれないため、ここで HDF5 を直接読む。
BPM 出力が無い場合 (純粋な FDTD 実行) は何もしない。
*/

#include "obpm.h"
#include "obpm_prototype.h"
#include "ev.h"

#include "hdf5.h"

#define BPM_H5_FILE   "time_series_data.h5"
#define BPM_CSV_FILE  "bpm_ixz.csv"
#define BPM_DB_RANGE  40.0    /* 等高線のダイナミックレンジ [dB] */
#define BPM_FRAME_PAGES 6     /* frames から描画する最大ページ数 (等間隔に間引く) */

/* 2次元データセットを読む (呼び出し側で free)。失敗時 NULL */
static double *read2d(hid_t file, const char *path, hsize_t dims[2])
{
	if (H5Lexists(file, path, H5P_DEFAULT) <= 0) return NULL;
	hid_t dset = H5Dopen(file, path, H5P_DEFAULT);
	if (dset < 0) return NULL;
	hid_t space = H5Dget_space(dset);
	if (H5Sget_simple_extent_ndims(space) != 2) {
		H5Sclose(space);
		H5Dclose(dset);
		return NULL;
	}
	H5Sget_simple_extent_dims(space, dims, NULL);
	double *buf = (double *)malloc((size_t)dims[0] * dims[1] * sizeof(double));
	if (buf != NULL) {
		if (H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		            H5P_DEFAULT, buf) < 0) {
			free(buf);
			buf = NULL;
		}
	}
	H5Sclose(space);
	H5Dclose(dset);
	return buf;
}

/* 3次元データセットを読む (呼び出し側で free)。失敗時 NULL */
static double *read3d(hid_t file, const char *path, hsize_t dims[3])
{
	if (H5Lexists(file, path, H5P_DEFAULT) <= 0) return NULL;
	hid_t dset = H5Dopen(file, path, H5P_DEFAULT);
	if (dset < 0) return NULL;
	hid_t space = H5Dget_space(dset);
	if (H5Sget_simple_extent_ndims(space) != 3) {
		H5Sclose(space);
		H5Dclose(dset);
		return NULL;
	}
	H5Sget_simple_extent_dims(space, dims, NULL);
	double *buf = (double *)malloc((size_t)dims[0] * dims[1] * dims[2] * sizeof(double));
	if (buf != NULL) {
		if (H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		            H5P_DEFAULT, buf) < 0) {
			free(buf);
			buf = NULL;
		}
	}
	H5Sclose(space);
	H5Dclose(dset);
	return buf;
}

/* スカラーのメタデータを読む (無ければ def を返す) */
static double readScalar(hid_t file, const char *path, double def)
{
	if (H5Lexists(file, path, H5P_DEFAULT) <= 0) return def;
	hid_t dset = H5Dopen(file, path, H5P_DEFAULT);
	if (dset < 0) return def;
	double v = def;
	if (H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v) < 0) {
		v = def;
	}
	H5Dclose(dset);
	return v;
}

/* 1次元データセットを読む (呼び出し側で free)。失敗時 NULL、要素数を *n へ */
static double *read1d(hid_t file, const char *path, hsize_t *n)
{
	if (H5Lexists(file, path, H5P_DEFAULT) <= 0) return NULL;
	hid_t dset = H5Dopen(file, path, H5P_DEFAULT);
	if (dset < 0) return NULL;
	hid_t space = H5Dget_space(dset);
	if (H5Sget_simple_extent_ndims(space) != 1) {
		H5Sclose(space);
		H5Dclose(dset);
		return NULL;
	}
	H5Sget_simple_extent_dims(space, n, NULL);
	double *buf = (double *)malloc((size_t)(*n) * sizeof(double));
	if (buf != NULL) {
		if (H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		            H5P_DEFAULT, buf) < 0) {
			free(buf);
			buf = NULL;
		}
	}
	H5Sclose(space);
	H5Dclose(dset);
	return buf;
}

/* /trace/* を読む。長さが expect と違う場合は破棄して NULL を返す
   (同じ横軸 /trace/z に載せられないものを描かないため) */
static double *readTrace(hid_t file, const char *path, hsize_t expect)
{
	hsize_t n = 0;
	double *v = read1d(file, path, &n);
	if ((v != NULL) && (n != expect)) {
		free(v);
		v = NULL;
	}
	return v;
}

/* z を横軸とする折れ線グラフページを1枚描く。
   nseries 本の系列を同じ縦軸 (共通の min/max) で重ねる。
   縦軸を共有するので、単位の異なる量を混ぜないこと。 */
static void tracePage(int nz, const double *z,
	int nseries, const double *const *series, const char *const *labels,
	const char *title, const char *ylabel)
{
	if ((nz < 2) || (nseries < 1)) return;

	/* 全系列の値域 (0 を含めて原点からの大小が読めるようにする) */
	double ymin = 0, ymax = 0;
	for (int m = 0; m < nseries; m++) {
		for (int i = 0; i < nz; i++) {
			const double v = series[m][i];
			if (v < ymin) ymin = v;
			if (v > ymax) ymax = v;
		}
	}
	if (ymin == ymax) {
		/* 定数系列 : 目盛りが潰れないよう幅を持たせる */
		const double d = (ymax != 0) ? (0.05 * fabs(ymax)) : 1.0;
		ymin -= d;
		ymax += d;
	}

	/* レイアウト (post/plot2dFeed.c と同じ流儀) */
	const double px0 = 0.13 * Width2d;
	const double px1 = 0.92 * Width2d;
	const double py0 = 0.13 * Height2d;
	const double py1 = 0.88 * Height2d;
	const double fs  = Fontsize2d;
	/* 系列の色 (最大 4 本。赤・青・緑・マゼンタ) */
	const unsigned char rgb[4][3] = {
		{200, 0, 0}, {0, 0, 200}, {0, 140, 0}, {180, 0, 180}
	};

	char str[BUFSIZ];
	ev2d_newPage();
	ev2d_setColor(0, 0, 0);
	ev2dlib_grid(px0, py0, px1, py1, 10, 10);

	for (int m = 0; m < nseries; m++) {
		const int c = m % 4;
		ev2d_setColor(rgb[c][0], rgb[c][1], rgb[c][2]);
		ev2dlib_func2(nz - 1, z, series[m], ymin, ymax, px0, py0, px1, py1);
		/* 凡例 (グラフ枠の上に横並び) */
		ev2d_drawString(px0 + m * 12.0 * fs, py1 + 1.9 * fs, 0.8 * fs, labels[m]);
	}

	ev2d_setColor(0, 0, 0);
	ev2d_drawString(px0, py1 + 3.2 * fs, fs, title);
	sprintf(str, "%.4e", z[0]);
	char smax[BUFSIZ];
	sprintf(smax, "%.4e", z[nz - 1]);
	ev2dlib_Xaxis(px0, px1, py0, 0.8 * fs, str, smax, "z[m]");
	sprintf(str, "%.4e", ymin);
	sprintf(smax, "%.4e", ymax);
	ev2dlib_Yaxis(py0, py1, px0, 0.8 * fs, str, smax, ylabel);
}

/* dB スケールの等高線ページを1枚描く (data は行 = 縦軸, 列 = 横軸) */
static void contourPage(const double *data, int rows, int cols,
	const char *title, const char *xlabel, const char *ylabel,
	double xmin, double xmax, double ymin, double ymax)
{
	/* ピーク基準の dB へ変換 (床 -BPM_DB_RANGE dB) */
	double peak = 0;
	for (int i = 0; i < rows * cols; i++) {
		if (data[i] > peak) peak = data[i];
	}
	if (peak <= 0) return;

	double **z = (double **)malloc((cols + 1) * sizeof(double *));
	double *zbuf = (double *)malloc((size_t)(cols + 1) * (rows + 1) * sizeof(double));
	double *xg = (double *)malloc((cols + 1) * sizeof(double));
	double *yg = (double *)malloc((rows + 1) * sizeof(double));
	if ((z == NULL) || (zbuf == NULL) || (xg == NULL) || (yg == NULL)) {
		free(z); free(zbuf); free(xg); free(yg);
		return;
	}

	/* ページ座標系: 余白付きの矩形にマップ */
	const double px0 = 0.14 * Width2d, px1 = 0.92 * Width2d;
	const double py0 = 0.16 * Height2d, py1 = 0.84 * Height2d;
	for (int i = 0; i <= cols; i++) {
		xg[i] = px0 + (px1 - px0) * i / cols;
	}
	for (int j = 0; j <= rows; j++) {
		yg[j] = py0 + (py1 - py0) * j / rows;
	}
	/* ev2dlib_contour は z[i][j] (i: 横, j: 縦) の節点値を期待する */
	for (int i = 0; i <= cols; i++) {
		z[i] = &zbuf[(size_t)i * (rows + 1)];
		for (int j = 0; j <= rows; j++) {
			const int ci = (i < cols) ? i : cols - 1;
			const int rj = (j < rows) ? j : rows - 1;
			const double v = data[(size_t)rj * cols + ci];
			double db = 10.0 * log10((v > 0 ? v : 1e-30) / peak);
			if (db < -BPM_DB_RANGE) db = -BPM_DB_RANGE;
			z[i][j] = db;
		}
	}

	ev2d_newPage();
	ev2dlib_contour(cols, rows, xg, yg, z, -BPM_DB_RANGE, 0.0, 1);
	ev2d_setColor(0, 0, 0);
	ev2d_drawRectangle(px0, py0, px1, py1);

	char str[BUFSIZ];
	const double fs = 0.025 * Width2d;
	ev2d_drawString(px0, py1 + 0.8 * fs, fs, title);
	sprintf(str, "%s : %.3e - %.3e", xlabel, xmin, xmax);
	ev2d_drawString(px0, py0 - 1.6 * fs, 0.8 * fs, str);
	sprintf(str, "%s : %.3e - %.3e", ylabel, ymin, ymax);
	ev2d_drawString(px0, py0 - 2.8 * fs, 0.8 * fs, str);
	sprintf(str, "0dB(peak) - -%.0fdB", BPM_DB_RANGE);
	ev2d_drawString(px1 - 12.0 * 0.8 * fs, py1 + 0.8 * fs, 0.8 * fs, str);

	free(z);
	free(zbuf);
	free(xg);
	free(yg);
}

int plot2dBpm(void)
{
	hid_t file = H5Fopen(BPM_H5_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file < 0) return 0;
	if (H5Lexists(file, "/field", H5P_DEFAULT) <= 0) {
		H5Fclose(file);
		return 0;
	}

	int pages = 0;

	/* ── 伝搬マップ Ixz (Nz x Nx) ── */
	hsize_t dims[2];
	double *ixz = read2d(file, "/field/Ixz", dims);
	if (ixz != NULL) {
		const int nz = (int)dims[0], nx = (int)dims[1];
		const double x0 = (nx == Nx) ? Xn[0] : 0, x1 = (nx == Nx) ? Xn[Nx] : nx;
		const double z0 = (nz <= Nz) ? Zn[0] : 0, z1 = (nz <= Nz) ? Zn[nz] : nz;
		contourPage(ixz, nz, nx, "BPM propagation map |E(x,z)|^2 [dB]",
		            "x[m]", "z[m]", x0, x1, z0, z1);
		pages++;

		/* CSV 出力 (z 行 x 列) */
		FILE *fp = fopen(BPM_CSV_FILE, "w");
		if (fp != NULL) {
			for (int iz = 0; iz < nz; iz++) {
				for (int ix = 0; ix < nx; ix++) {
					fprintf(fp, "%s%.6e", (ix ? "," : ""),
					        ixz[(size_t)iz * nx + ix]);
				}
				fprintf(fp, "\n");
			}
			fclose(fp);
		}

		/* サマリ: 入口/出口の断面パワー比 (伝搬損失の目安) */
		double pin = 0, pout = 0;
		for (int ix = 0; ix < nx; ix++) {
			pin += ixz[ix];
			pout += ixz[(size_t)(nz - 1) * nx + ix];
		}
		printf("BPM: propagation map %d x %d, P(end)/P(start) = %.4f\n",
		       nz, nx, (pin > 0) ? pout / pin : 0.0);
		free(ixz);
	}

	/* ── 最終断面 |Efinal| (Ny x Nx) ── */
	hsize_t dr[2], di[2];
	double *er = read2d(file, "/field/Efinal_r", dr);
	double *ei = read2d(file, "/field/Efinal_i", di);
	if ((er != NULL) && (ei != NULL) &&
	    (dr[0] == di[0]) && (dr[1] == di[1])) {
		const int ny = (int)dr[0], nx = (int)dr[1];
		double *mag = (double *)malloc((size_t)ny * nx * sizeof(double));
		if (mag != NULL) {
			for (size_t i = 0; i < (size_t)ny * nx; i++) {
				mag[i] = er[i] * er[i] + ei[i] * ei[i];
			}
			const double x0 = (nx == Nx) ? Xn[0] : 0;
			const double x1 = (nx == Nx) ? Xn[Nx] : nx;
			const double y0 = (ny == Ny) ? Yn[0] : 0;
			const double y1 = (ny == Ny) ? Yn[Ny] : ny;
			contourPage(mag, ny, nx, "BPM final field |E(x,y)|^2 [dB]",
			            "x[m]", "y[m]", x0, x1, y0, y1);
			pages++;
			free(mag);
		}
	}
	free(er);
	free(ei);

	/* ── 伝搬スナップショット frames (nframes x Ny x Nx) ──
	   入力に frames = <interval> を指定した場合のみ存在する。
	   全フレームを描くとページ数が過大になるため等間隔に間引き、
	   各ページのタイトルには (取得できれば) z 位置を入れる。 */
	hsize_t fd[3];
	double *frames = read3d(file, "/field/frames", fd);
	if (frames != NULL) {
		const int nf = (int)fd[0], ny = (int)fd[1], nx = (int)fd[2];
		const double interval = readScalar(file, "/metadata/frame_interval", 0);
		const double dz = readScalar(file, "/metadata/grid_dz", 0);
		const int npage = (nf < BPM_FRAME_PAGES) ? nf : BPM_FRAME_PAGES;
		const double x0 = (nx == Nx) ? Xn[0] : 0;
		const double x1 = (nx == Nx) ? Xn[Nx] : nx;
		const double y0 = (ny == Ny) ? Yn[0] : 0;
		const double y1 = (ny == Ny) ? Yn[Ny] : ny;
		for (int p = 0; p < npage; p++) {
			/* 先頭と末尾を必ず含む等間隔サンプリング */
			const int idx = (npage > 1) ? (int)((double)p * (nf - 1) / (npage - 1)) : 0;
			char title[BUFSIZ];
			if ((interval > 0) && (dz > 0)) {
				sprintf(title, "BPM snapshot |E(x,y)|^2 [dB]  z = %.6e [m] (frame %d/%d)",
				        Zn[0] + (idx * interval * dz), idx + 1, nf);
			}
			else {
				sprintf(title, "BPM snapshot |E(x,y)|^2 [dB] (frame %d/%d)", idx + 1, nf);
			}
			contourPage(&frames[(size_t)idx * ny * nx], ny, nx, title,
			            "x[m]", "y[m]", x0, x1, y0, y1);
			pages++;
		}
		printf("BPM: snapshots %d frames (%d x %d), plotted %d page(s)\n",
		       nf, ny, nx, npage);
		free(frames);
	}

	/* ── z 伝搬に沿ったスカラー推移 (/trace) ──
	   BPM は定常解法で時間軸を持たない。系列の独立変数は z である
	   (HDF5 のルート属性 marching_axis = "z" を参照)。
	   単位の異なる量を同じ縦軸に混ぜないよう、ページを分けて描く。 */
	if (H5Lexists(file, "/trace", H5P_DEFAULT) > 0) {
		hsize_t nzt = 0;
		double *z = read1d(file, "/trace/z", &nzt);
		if ((z != NULL) && (nzt >= 2)) {
			const int nz = (int)nzt;
			/* /trace/* はすべて /trace/z と同じ長さで書かれる。
			   長さが違うものは (将来の不整合に備えて) 読み捨てる。 */
			double *power = readTrace(file, "/trace/power", nzt);
			double *peak  = readTrace(file, "/trace/peak", nzt);
			double *wx    = readTrace(file, "/trace/width_x", nzt);
			double *wy    = readTrace(file, "/trace/width_y", nzt);
			double *cx    = readTrace(file, "/trace/centroid_x", nzt);
			double *cy    = readTrace(file, "/trace/centroid_y", nzt);

			if (power != NULL) {
				const double *ser[1] = {power};
				const char *lab[1] = {"power"};
				tracePage(nz, z, 1, ser, lab,
				          "BPM trace : cross-section power", "P");
				pages++;
			}
			if (peak != NULL) {
				const double *ser[1] = {peak};
				const char *lab[1] = {"peak"};
				tracePage(nz, z, 1, ser, lab,
				          "BPM trace : peak intensity |E|^2", "|E|^2");
				pages++;
			}
			if ((wx != NULL) && (wy != NULL)) {
				const double *ser[2] = {wx, wy};
				const char *lab[2] = {"width_x", "width_y"};
				tracePage(nz, z, 2, ser, lab,
				          "BPM trace : beam width 2*sigma", "w[m]");
				pages++;
			}
			if ((cx != NULL) && (cy != NULL)) {
				const double *ser[2] = {cx, cy};
				const char *lab[2] = {"centroid_x", "centroid_y"};
				tracePage(nz, z, 2, ser, lab,
				          "BPM trace : intensity centroid", "c[m]");
				pages++;
			}

			/* モード結合率 (nModes x nz) : modes 指定時のみ存在する。
			   先頭 4 モードまでを 1 ページに重ねる (単位は共通の占有率)。 */
			hsize_t ovdims[2];
			double *ov = read2d(file, "/trace/overlap", ovdims);
			if ((ov != NULL) && ((int)ovdims[1] == nz)) {
				const int nm = ((int)ovdims[0] < 4) ? (int)ovdims[0] : 4;
				const double *ser[4];
				const char *lab[4];
				char labbuf[4][32];
				for (int m = 0; m < nm; m++) {
					ser[m] = &ov[(size_t)m * nz];
					sprintf(labbuf[m], "mode %d", m + 1);
					lab[m] = labbuf[m];
				}
				tracePage(nz, z, nm, ser, lab,
				          "BPM trace : mode coupling ratio", "eta");
				pages++;
				printf("BPM: trace graphs (%d modes in overlap)\n", (int)ovdims[0]);
			}
			free(ov);

			printf("BPM: propagation traces %d points, z = %.6e .. %.6e [m]\n",
			       nz, z[0], z[nz - 1]);
			free(power); free(peak);
			free(wx); free(wy); free(cx); free(cy);
		}
		free(z);
	}

	H5Fclose(file);
	return pages;
}
