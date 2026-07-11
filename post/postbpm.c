/*
postbpm.c

BPM 出力 (time_series_data.h5 の /field) の可視化。
  - /field/Ixz       : 伝搬マップ |E(x, y=Ny/2, z)|^2 (Nz x Nx) → 等高線ページ
  - /field/Efinal_r/i: 最終断面電界 (Ny x Nx) → |E| 等高線ページ
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

	H5Fclose(file);
	return pages;
}
