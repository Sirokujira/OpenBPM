/*
input_data.c

input data
*/

#include "obpm.h"
#include "obpm_prototype.h"

#define MAXTOKEN 1000


int input_data(FILE *fp)
{
	int    ntoken, ngeom, nline;
	int    version = 0;
	int    nxr = 0, nyr = 0, nzr = 0;
	int    *dxr = NULL, *dyr = NULL, *dzr = NULL;
	double *xr = NULL, *yr = NULL, *zr = NULL;
	double *xfeed = NULL, *yfeed = NULL, *zfeed = NULL;
	double *xpoint = NULL, *ypoint = NULL, *zpoint = NULL;
	int    nload = 0;
	double *xload = NULL, *yload = NULL, *zload = NULL, *pload = NULL;
	char   *dload = NULL, *cload = NULL;
	char   prog[BUFSIZ];
	char   strline[BUFSIZ], strkey[BUFSIZ], strsave[BUFSIZ];
	char   strprop[BUFSIZ];
	char   *token[MAXTOKEN];
	const int array_inc = 10000;		// reduce malloc times
	const char sep[] = " \t";			// separator
	const char errfmt1[] = "*** too many %s data #%d\n";
	const char errfmt2[] = "*** invalid %s data\n";
	const char errfmt3[] = "*** invalid %s data #%d\n";

	// initialize

	NMaterial = 2;  // air + PEC
	Material = (material_t *)malloc(NMaterial * sizeof(material_t));
	for (int64_t m = 0; m < NMaterial; m++) {
		Material[m].type = 1;
		Material[m].epsr = 1;
		Material[m].esgm = 0;
		Material[m].amur = 1;
		Material[m].msgm = 0;
	}

	NGeometry = 0;

	NFeed = 0;
	rFeed = 0;

	IPlanewave = 0;

	NPoint = 0;

	NInductor = 0;

	iABC = 0;  // Mur-1st
	PBCx = PBCy = PBCz = 0;

	Dt = 0;
	Tw = 0;

	Solver.maxiter = 3000;
	Solver.nout = 50;
	Solver.converg = 1e-3;

	BPM.w0 = 0;       // 0 = auto
	BPM.x0 = BPM.y0 = 0;
	BPM.center = 0;
	BPM.n0 = 0;       // 0 = auto
	BPM.RoC = 0;      // 0 = straight
	BPM.bendDir = 0;
	BPM.rho_e = 0;
	BPM.pol = 0;      // scalar
	BPM.wideangle = 0;
	BPM.tiltx = BPM.tilty = 0;
	BPM.frames = 0;   // 0 = スナップショット記録なし
	BPM.framesComplex = 0;
	BPM.nmodes = 0;   // 0 = モード解析なし
	BPM.modeExcite = 0;
	BPM.taper = 1;    // 1 = テーパなし (出口/入口の横方向スケール比)
	BPM.twist = 0;    // 0 = ツイストなし [deg/m]
	BPM.symx = 0;     // 0 = 対称境界なし
	BPM.symy = 0;
	BPM.launchMode = -1;  // -1 = ガウシアン励振 (既定)
	BPM.launchNModes = 0; // 0 = ガウシアン励振
	BPM.wlsweep = 0;      // 0 = 先頭波長のみ (既定)

	// 非線形吸収 (TPA) と入力パワー掃引 (省略時は従来動作 = 線形単発)
	NTpaB = 0;
	TpaB = NULL;
	PowerSweep.npoints = 0;   // 0 = 掃引なし
	PowerSweep.pmin = PowerSweep.pmax = 0;
	PowerSweep.logscale = 0;

	NFreq1 =
	NFreq2 = 0;

	Plot3dGeom = 0;

	// read

	nline = 0;
	while (fgets(strline, sizeof(strline), fp) != NULL) {
		// skip a empty line
		if (strlen(strline) <= 1) continue;

		// skip a comment line
		if (strline[0] == '#') continue;

		// delete "\n"
		//printf("%zd\n", strlen(strline));
		if (strstr(strline, "\r\n") != NULL) {
			strline[strlen(strline) - 2] = '\0';
		}
		else if ((strstr(strline, "\r") != NULL) || (strstr(strline, "\n") != NULL)) {
			strline[strlen(strline) - 1] = '\0';
		}
		//printf("%zd\n", strlen(strline));

		// "end" -> break
		if (!strncmp(strline, "end", 3)) break;

		// save "strline"
		strcpy(strsave, strline);

		// token ("strline" is destroyed)
		ntoken = tokenize(strline, sep, token, MAXTOKEN);
		//for (int i = 0; i < ntoken; i++) printf("%d %s\n", i, token[i]);

		// check number of data and "=" (exclude header)
		if ((nline > 0) && ((ntoken < 3) || strcmp(token[1], "="))) continue;

		// keyword
		strcpy(strkey, token[0]);

		// input
		if      (nline == 0) {
			strcpy(prog, strkey);
			if (strcmp(prog, "OpenBPM") && strcmp(prog, "OpenTHFD")) {
				printf("%s\n", "*** not OpenBPM data");
				return 1;
			}
			if (ntoken < 3) {
				printf("%s\n", "*** no version data");
				return 1;
			}
			version = (10 * atoi(token[1])) + atoi(token[2]);
			nline++;
		}
		else if (!strcmp(strkey, "title")) {
			strcpy(Title, strchr(strsave, '=') + 2);
		}
		else if (!strcmp(strkey, "xmesh")) {
			if ((ntoken < 5) || (ntoken % 2 == 0)) {
				printf(errfmt2, strkey);
				return 1;
			}
			nxr = (ntoken - 3) / 2;
			xr = (double *)malloc((nxr + 1) * sizeof(double));
			dxr = (int *)malloc(nxr * sizeof(int));
			sscanf(token[2], "%lf", &xr[0]);
			for (int i = 0; i < nxr; i++) {
				sscanf(token[2 * i + 3], "%d", &dxr[i]);
				sscanf(token[2 * i + 4], "%lf", &xr[i + 1]);
			}
		}
		else if (!strcmp(strkey, "ymesh")) {
			if ((ntoken < 5) || (ntoken % 2 == 0)) {
				printf(errfmt2, strkey);
				return 1;
			}
			nyr = (ntoken - 3) / 2;
			yr = (double *)malloc((nyr + 1) * sizeof(double));
			dyr = (int *)malloc(nyr * sizeof(int));
			sscanf(token[2], "%lf", &yr[0]);
			for (int j = 0; j < nyr; j++) {
				sscanf(token[2 * j + 3], "%d", &dyr[j]);
				sscanf(token[2 * j + 4], "%lf", &yr[j + 1]);
			}
		}
		else if (!strcmp(strkey, "zmesh")) {
			if ((ntoken < 5) || (ntoken % 2 == 0)) {
				printf(errfmt2, strkey);
				return 1;
			}
			nzr = (ntoken - 3) / 2;
			zr = (double *)malloc((nzr + 1) * sizeof(double));
			dzr = (int *)malloc(nzr * sizeof(int));
			sscanf(token[2], "%lf", &zr[0]);
			for (int k = 0; k < nzr; k++) {
				sscanf(token[2 * k + 3], "%d", &dzr[k]);
				sscanf(token[2 * k + 4], "%lf", &zr[k + 1]);
			}
		}
		else if (!strcmp(strkey, "material")) {
			if (ntoken < 7) {
				printf(errfmt3, strkey, (int)NMaterial - 1);
				return 1;
			}
			if (NMaterial % array_inc == 2) {   // 2 : initial set (air + PEC)
				Material = (material_t *)realloc(Material, (NMaterial + array_inc) * sizeof(material_t));
			}

			if (NMaterial >= MAXMATERIAL) {
				printf(errfmt1, strkey, (int)NMaterial - 1);
				return 1;
			}

			int type = 1;
			double epsr = 1, esgm = 0, amur = 1, msgm = 0;
			double einf = 0, ae = 0, be = 0, ce = 0;
			if (!strcmp(prog, "OpenBPM") && (version < 22)) {
				type = 1;
				epsr = atof(token[2]);
				esgm = atof(token[3]);
				amur = atof(token[4]);
				msgm = atof(token[5]);
			}
			else if (!strcmp(token[2], "1")) {
				type = 1;
				epsr = atof(token[3]);
				esgm = atof(token[4]);
				amur = atof(token[5]);
				msgm = atof(token[6]);
			}
			else if (!strcmp(token[2], "2")) {
				type = 2;
				einf = atof(token[3]);
				ae   = atof(token[4]);
				be   = atof(token[5]);
				ce   = atof(token[6]);
				epsr = 1;
				esgm = 0;
				amur = 1;
				msgm = 0;
			}
			if ((type == 1) && ((epsr <= 0) || (esgm < 0) || (amur <= 0) || (msgm < 0))) {
				printf(errfmt3, strkey, (int)NMaterial - 1);
				return 1;
			}
			else if ((type == 2) && (einf <= 0)) {
				printf(errfmt3, strkey, (int)NMaterial - 1);
				return 1;
			}
			Material[NMaterial].type = type;
			Material[NMaterial].epsr = epsr;
			Material[NMaterial].esgm = esgm;
			Material[NMaterial].amur = amur;
			Material[NMaterial].msgm = msgm;
			Material[NMaterial].einf = einf;
			Material[NMaterial].ae   = ae;
			Material[NMaterial].be   = be;
			Material[NMaterial].ce   = ce;
			NMaterial++;
			//printf("%d\n", NMaterial);
		}
		else if (!strcmp(strkey, "geometry")) {
			if (ntoken < 4) {
				printf(errfmt3, strkey, (int)NGeometry + 1);
				return 1;
			}
			if (NGeometry % array_inc == 0) {
				Geometry = (geometry_t *)realloc(Geometry, (NGeometry + array_inc) * sizeof(geometry_t));
			}
			Geometry[NGeometry].m     = (id_t)atoi(token[2]);
			Geometry[NGeometry].shape = atoi(token[3]);
			switch (Geometry[NGeometry].shape) {
				case 1:
				case 2:
				case 11:
				case 12:
				case 13:
					ngeom = 6;
					break;
				case 31:
				case 32:
				case 33:
				case 41:
				case 42:
				case 43:
				case 51:
				case 52:
				case 53:
					ngeom = 8;
					break;
				default:
					ngeom = 0;
					break;
			}
			if (ntoken < 4 + ngeom) {
				printf(errfmt3, strkey, (int)NGeometry + 1);
				return 1;
			}
			for (int n = 0; n < ngeom; n++) {
				Geometry[NGeometry].g[n] = atof(token[4 + n]);
			}
			NGeometry++;
		}
		else if (!strcmp(strkey, "name")) {
			;
		}
		else if (!strcmp(strkey, "feed")) {
			if (ntoken > 8) {
				Feed  = (feed_t *)realloc(Feed,  (NFeed + 1) * sizeof(feed_t));
				xfeed = (double *)realloc(xfeed, (NFeed + 1) * sizeof(double));
				yfeed = (double *)realloc(yfeed, (NFeed + 1) * sizeof(double));
				zfeed = (double *)realloc(zfeed, (NFeed + 1) * sizeof(double));
				Feed[NFeed].dir   = (char)toupper((int)token[2][0]);
				xfeed[NFeed]      = atof(token[3]);
				yfeed[NFeed]      = atof(token[4]);
				zfeed[NFeed]      = atof(token[5]);
				Feed[NFeed].volt  = atof(token[6]);
				Feed[NFeed].delay = atof(token[7]);
				Feed[NFeed].z0    = atof(token[8]);
				NFeed++;
			}
			else {
				printf(errfmt3, strkey, NFeed + 1);
				return 1;
			}
		}
		else if (!strcmp(strkey, "planewave")) {
			if (ntoken < 5) {
				printf(errfmt2, strkey);
				return 1;
			}
			IPlanewave = 1;
			Planewave.theta = atof(token[2]);
			Planewave.phi   = atof(token[3]);
			Planewave.pol   = atoi(token[4]);
		}
		else if (!strcmp(strkey, "point")) {
			//if (((NPoint == 0) && (ntoken < 7)) || (ntoken < 6)) {
			if (ntoken < 6) {
				printf(errfmt3, strkey, NPoint + 1);
				return 1;
			}
			Point  = (point_t *)realloc(Point,  (NPoint + 1) * sizeof(point_t));
			xpoint = (double *) realloc(xpoint, (NPoint + 1) * sizeof(double));
			ypoint = (double *) realloc(ypoint, (NPoint + 1) * sizeof(double));
			zpoint = (double *) realloc(zpoint, (NPoint + 1) * sizeof(double));
			Point[NPoint].dir = (char)toupper((int)token[2][0]);
			xpoint[NPoint]    = atof(token[3]);
			ypoint[NPoint]    = atof(token[4]);
			zpoint[NPoint]    = atof(token[5]);
			if (NPoint == 0) {
				// propagation on port #1
				if (ntoken > 6) {
					strcpy(strprop, token[6]);
				}
				else {
					strcpy(strprop, "+X");
				}
			}
			NPoint++;
		}
		else if (!strcmp(strkey, "load")) {
			if (ntoken < 8) {
				printf(errfmt3, strkey, nload + 1);
				return 1;
			}
			dload = (char   *)realloc(dload, (nload + 1) * sizeof(char));
			xload = (double *)realloc(xload, (nload + 1) * sizeof(double));
			yload = (double *)realloc(yload, (nload + 1) * sizeof(double));
			zload = (double *)realloc(zload, (nload + 1) * sizeof(double));
			cload = (char   *)realloc(cload, (nload + 1) * sizeof(char));
			pload = (double *)realloc(pload, (nload + 1) * sizeof(double));
			dload[nload] = (char)toupper((int)token[2][0]);
			xload[nload] = atof(token[3]);
			yload[nload] = atof(token[4]);
			zload[nload] = atof(token[5]);
			cload[nload] = (char)toupper((int)token[6][0]);
			pload[nload] = atof(token[7]);
			nload++;
		}
		else if (!strcmp(strkey, "rfeed")) {
			rFeed = atof(token[2]);
		}
		else if (!strcmp(strkey, "abc")) {
			if      ((ntoken >= 3) && !strncmp(token[2], "0", 1)) {
				iABC = 0;
			}
			else if ((ntoken >= 6) && !strncmp(token[2], "1", 1)) {
				iABC = 1;
				cPML.l = atoi(token[3]);
				cPML.m = atof(token[4]);
				cPML.r0 = atof(token[5]);
			}
			else {
				printf(errfmt2, strkey);
				return 1;
			}
		}
		else if (!strcmp(strkey, "pbc")) {
			if      (ntoken >= 5) {
				PBCx = atoi(token[2]);
				PBCy = atoi(token[3]);
				PBCz = atoi(token[4]);
			}
			else {
				printf(errfmt2, strkey);
				return 1;
			}
		}
		else if (!strcmp(strkey, "frequency1")) {
			if (ntoken > 4) {
				double f0 = atof(token[2]);
				double f1 = atof(token[3]);
				int fdiv = atoi(token[4]);
				if (fdiv < 0) {
					printf(errfmt2, strkey);
					return 1;
				}
				double df = (fdiv > 0) ? (f1 - f0) / fdiv : 0;
				NFreq1 = fdiv + 1;
				Freq1 = (double *)malloc(NFreq1 * sizeof(double));
				for (int n = 0; n < NFreq1; n++) {
					Freq1[n] = f0 + (n * df);
				}
			}
		}
		else if (!strcmp(strkey, "frequency2")) {
			if (ntoken > 4) {
				double f0 = atof(token[2]);
				double f1 = atof(token[3]);
				int fdiv = atoi(token[4]);
				if (fdiv < 0) {
					printf(errfmt2, strkey);
					return 1;
				}
				double df = (fdiv > 0) ? (f1 - f0) / fdiv : 0;
				NFreq2 = fdiv + 1;
				Freq2 = (double *)malloc(NFreq2 * sizeof(double));
				for (int n = 0; n < NFreq2; n++) {
					Freq2[n] = f0 + (n * df);
				}
			}
		}
		else if (!strcmp(prog, "OpenTHFD") && !strcmp(strkey, "frequency")) {
			if (ntoken > 4) {
				double f0 = atof(token[2]);
				double f1 = atof(token[3]);
				int fdiv = atoi(token[4]);
				if (fdiv < 0) {
					printf(errfmt2, strkey);
					return 1;
				}
				double df = (fdiv > 0) ? (f1 - f0) / fdiv : 0;
				NFreq1 = NFreq2 = fdiv + 1;
				Freq1 = (double *)malloc(NFreq1 * sizeof(double));
				Freq2 = (double *)malloc(NFreq2 * sizeof(double));
				for (int n = 0; n < NFreq1; n++) {
					Freq1[n] = Freq2[n] = f0 + (n * df);
				}
			}
		}
		else if (!strcmp(strkey, "solver")) {
			if (ntoken > 4) {
				Solver.maxiter = atoi(token[2]);
				Solver.nout    = atoi(token[3]);
				Solver.converg = atof(token[4]);
				if (Solver.maxiter % Solver.nout != 0) {
					Solver.maxiter = (Solver.maxiter / Solver.nout + 1) * Solver.nout;
				}
			}
		}
		else if (!strcmp(strkey, "timestep")) {
			Dt = atof(token[2]);
		}
		else if (!strcmp(strkey, "pulsewidth")) {
			Tw = atof(token[2]);
		}
		else if (!strcmp(strkey, "plot3dgeom")) {
			Plot3dGeom = atoi(token[2]);
		}
		else if (!strcmp(strkey, "beam")) {
			// beam = <w0[m]> [<x0[m]> <y0[m]>]
			BPM.w0 = atof(token[2]);
			if (BPM.w0 <= 0) {
				printf(errfmt2, strkey);
				return 1;
			}
			if (ntoken > 4) {
				BPM.x0 = atof(token[3]);
				BPM.y0 = atof(token[4]);
				BPM.center = 1;
			}
		}
		else if (!strcmp(strkey, "refindex")) {
			// refindex = <n0> : BPM reference refractive index
			BPM.n0 = atof(token[2]);
			if (BPM.n0 <= 0) {
				printf(errfmt2, strkey);
				return 1;
			}
		}
		else if (!strcmp(strkey, "bend")) {
			// bend = <RoC[m]> [<direction[deg]> [<rho_e>]] : 曲げ半径 (等価屈折率法)
			BPM.RoC = atof(token[2]);
			if (BPM.RoC <= 0) {
				printf(errfmt2, strkey);
				return 1;
			}
			if (ntoken > 3) BPM.bendDir = atof(token[3]);
			if (ntoken > 4) BPM.rho_e = atof(token[4]);
		}
		else if (!strcmp(strkey, "polarization")) {
			// polarization = <scalar|x|y> : 半ベクトル BPM の偏波方向
			if      (!strcmp(token[2], "x")) BPM.pol = 1;
			else if (!strcmp(token[2], "y")) BPM.pol = 2;
			else if (!strcmp(token[2], "scalar") || !strcmp(token[2], "0")) BPM.pol = 0;
			else {
				printf(errfmt2, strkey);
				return 1;
			}
		}
		else if (!strcmp(strkey, "wideangle")) {
			// wideangle = <0|1> : 広角 BPM (Pade(1,1))
			BPM.wideangle = atoi(token[2]);
		}
		else if (!strcmp(strkey, "beamtilt")) {
			// beamtilt = <theta_x[deg]> [<theta_y[deg]>] : 入射ビームの傾き
			BPM.tiltx = atof(token[2]);
			if (ntoken > 3) BPM.tilty = atof(token[3]);
		}
		else if (!strcmp(strkey, "frames")) {
			// frames = <interval> [complex] : |E(x,y)|^2 スナップショットの記録間隔
			// [z ステップ] (0 = 記録なし。/field/frames へ出力)。
			// complex 指定時は複素電界も /field/frames_r, /field/frames_i へ出力する
			// (位相分布の表示用。記憶容量は 3 倍になる)
			BPM.frames = atoi(token[2]);
			BPM.framesComplex = ((ntoken > 3) && !strcmp(token[3], "complex")) ? 1 : 0;
		}
		else if (!strcmp(strkey, "taper")) {
			// taper = <ratio> : 出口での横方向スケール比 (入口 = 1)。
			// 屈折率分布を z に沿って相似縮小/拡大する (座標変換, 近軸パス由来)
			BPM.taper = atof(token[2]);
		}
		else if (!strcmp(strkey, "twist")) {
			// twist = <rate[deg/m]> : 導波路のツイスト率 (座標系を z に沿って回転)
			BPM.twist = atof(token[2]);
		}
		else if (!strcmp(strkey, "symmetry")) {
			// symmetry = <x|y|xy> [sym|anti] : 対称境界 (計算領域を 1/2〜1/4 に削減)
			// 鏡像面は指定軸のメッシュ始端 (x = xmesh の始点 / y = ymesh の始点)。
			// メッシュには半領域のみを与えること (出力も半領域になる)
			const int t = ((ntoken > 3) && !strcmp(token[3], "anti")) ? 2 : 1;
			if (strchr(token[2], 'x') != NULL) BPM.symx = t;
			if (strchr(token[2], 'y') != NULL) BPM.symy = t;
			if ((BPM.symx == 0) && (BPM.symy == 0)) {
				printf(errfmt2, strkey);
				return 1;
			}
		}
		else if (!strcmp(strkey, "modes")) {
			// modes = <nModes> [excite] : 入力断面のモード解析 (虚軸伝搬法)。
			// neff を obpm.log へ、モード形状を /modes へ出力する。
			// excite 指定で入射ビームを基本モードに置き換える (モード整合励振)
			BPM.nmodes = atoi(token[2]);
			BPM.modeExcite = ((ntoken > 3) && !strcmp(token[3], "excite")) ? 1 : 0;
		}
		else if (!strcmp(strkey, "launch")) {
			// launch = gauss | mode [<m>[:<coef>] ...] : 励振方法 (既定 gauss)
			// mode 指定時は先頭スライスの屈折率分布から導波モードを虚軸伝搬法で
			// 求めて励振する (m 省略時は基本モード 0)。複数指定すると重ね合わせ、
			// <m>:<coef> で係数を与える (省略時 1)。入力電力は feed 電圧^2 に
			// 正規化されるため、係数は分岐比のみを決める
			if      (!strcmp(token[2], "mode")) {
				BPM.launchNModes = 0;
				for (int n = 3; n < ntoken; n++) {
					if (BPM.launchNModes >= BPM_MAX_LAUNCH_MODES) {
						printf(errfmt2, strkey);
						return 1;
					}
					// "<m>" または "<m>:<coef>"
					char *sep = strchr(token[n], ':');
					double coef = 1;
					if (sep != NULL) {
						*sep = '\0';
						coef = atof(sep + 1);
					}
					const int idx = atoi(token[n]);
					if ((idx < 0) || (coef == 0)) {
						printf(errfmt2, strkey);
						return 1;
					}
					BPM.launchIdx[BPM.launchNModes] = idx;
					BPM.launchCoef[BPM.launchNModes] = coef;
					BPM.launchNModes++;
				}
				if (BPM.launchNModes == 0) {
					// 引数省略時は基本モード
					BPM.launchIdx[0] = 0;
					BPM.launchCoef[0] = 1;
					BPM.launchNModes = 1;
				}
				BPM.launchMode = BPM.launchIdx[0];
			}
			else if (!strcmp(token[2], "gauss")) {
				BPM.launchMode = -1;
				BPM.launchNModes = 0;
			}
			else {
				printf(errfmt2, strkey);
				return 1;
			}
		}
		else if (!strcmp(strkey, "wlsweep")) {
			// wlsweep = <0|1> : 波長掃引 (frequency2 の全点を順に計算する)
			// 省略時 (0) は先頭周波数のみ = 従来動作
			BPM.wlsweep = atoi(token[2]);
		}
		else if (!strcmp(strkey, "tpa")) {
			// tpa = <material id> <beta[cm/GW]> : 二光子吸収 (TPA) 係数
			// (メタマテリアル装荷 Si 導波路の実測例 : beta = 424 cm/GW,
			//  Honda, Shoji, Amemiya, Opt. Lett. 49, 5811 (2024))
			if (ntoken < 4) {
				printf(errfmt2, strkey);
				return 1;
			}
			TpaB = (tpab_t *)realloc(TpaB, (NTpaB + 1) * sizeof(tpab_t));
			TpaB[NTpaB].m    = (id_t)atoi(token[2]);
			TpaB[NTpaB].beta = atof(token[3]);
			if (TpaB[NTpaB].beta <= 0) {
				printf(errfmt2, strkey);
				return 1;
			}
			NTpaB++;
		}
		else if (!strcmp(strkey, "powersweep")) {
			// powersweep = <Pmin[W]> <Pmax[W]> <npoints> [log|lin] : 入力パワー掃引
			// (光活性化関数曲線 P_out(P_in) の抽出用。省略時は従来の単発計算)
			if (ntoken < 5) {
				printf(errfmt2, strkey);
				return 1;
			}
			PowerSweep.pmin    = atof(token[2]);
			PowerSweep.pmax    = atof(token[3]);
			PowerSweep.npoints = atoi(token[4]);
			PowerSweep.logscale = 0;  // 既定 lin
			if (ntoken > 5) {
				if      (!strcmp(token[5], "log")) PowerSweep.logscale = 1;
				else if (!strcmp(token[5], "lin")) PowerSweep.logscale = 0;
				else {
					printf(errfmt2, strkey);
					return 1;
				}
			}
			if ((PowerSweep.pmin <= 0) || (PowerSweep.pmax < PowerSweep.pmin) ||
			    (PowerSweep.npoints < 1)) {
				printf(errfmt2, strkey);
				return 1;
			}
		}
	}
/*
	// debug
	//printf("title = %s\n", Title);
	//printf("xmesh = %e", xr[0]); for (int i = 0; i < nxr; i++) printf(" %d %e", dxr[i], xr[i + 1]); printf("\n");
	//printf("ymesh = %e", yr[0]); for (int j = 0; j < nyr; j++) printf(" %d %e", dyr[j], yr[j + 1]); printf("\n");
	//printf("zmesh = %e", zr[0]); for (int k = 0; k < nzr; k++) printf(" %d %e", dzr[k], zr[k + 1]); printf("\n");
	//for (int n = 0; n < NMaterial; n++) if (Material[n].type == 1) printf("material = %d %.3e %.3e %.3e %.3e\n", Material[n].epsr, Material[n].esgm, Material[n].amur, Material[n].msgm);
	//for (int n = 0; n < NMaterial; n++) if (Material[n].type == 2) printf("material = %.3e %.3e %.3e %.3e\n", Material[n].einf, Material[n].ae, Material[n].be, Material[n].ce);
	//for (int n = 0; n < NGeometry; n++) printf("geometry = %d %.3e %.3e %.3e %.3e %.3e %.3e\n", Geometry[n].m, Geometry[n].g[0], Geometry[n].g[1], Geometry[n].g[2], Geometry[n].g[3], Geometry[n].g[4], Geometry[n].g[5]);
	//for (int n = 0; n < NFeed; n++) printf("feed = %c %e %e %e %e\n", Feed[n].dir, xfeed[n], yfeed[n], zfeed[n], Feed[n].volt);
	//for (int n = 0; n < NFreq1; n++) printf("frequency1 = %e\n", Freq1[n]);
	//for (int n = 0; n < NFreq2; n++) printf("frequency2 = %e\n", Freq2[n]);
	//printf("%d %d %e %e\n", iABC, cPML.l, cPML.m, cPML.r0);
	//printf("solver = %d %d %e\n", Solver.maxiter, Solver.nout, Solver.converg);
*/
	// error check

	if (nxr <= 0) {
		printf("%s\n", "*** no xmesh data");
		return 1;
	}
	if (nyr <= 0) {
		printf("%s\n", "*** no ymesh data");
		return 1;
	}
	if (nzr <= 0) {
		printf("%s\n", "*** no zmesh data");
		return 1;
	}
	for (int i = 0; i < nxr; i++) {
		if ((xr[i] >= xr[i + 1]) || (dxr[i] <= 0)) {
			printf("%s\n", "*** invalid xmesh data");
			return 1;
		}
	}
	for (int j = 0; j < nyr; j++) {
		if ((yr[j] >= yr[j + 1]) || (dyr[j] <= 0)) {
			printf("%s\n", "*** invalid ymesh data");
			return 1;
		}
	}
	for (int k = 0; k < nzr; k++) {
		if ((zr[k] >= zr[k + 1]) || (dzr[k] <= 0)) {
			printf("%s\n", "*** invalid zmesh data");
			return 1;
		}
	}
	if (!NFeed && !IPlanewave) {
		printf("%s\n", "*** no source");
		return 1;
	}
	if (NFeed && IPlanewave) {
		printf("%s\n", "*** feed and planewave");
		return 1;
	}
	if ((Solver.maxiter <= 0) || (Solver.nout <= 0)) {
		printf("%s\n", "*** invalid solver data");
		return 1;
	}

	for (int64_t n = 0; n < NGeometry; n++) {
		if (Geometry[n].m >= NMaterial) {
			printf("*** invalid material id of geometry data #%zd\n", n + 1);
			return 1;
		}
	}
	for (int n = 0; n < NTpaB; n++) {
		if (TpaB[n].m >= NMaterial) {
			printf("*** invalid material id of tpa data #%d\n", n + 1);
			return 1;
		}
	}
	for (int n = 0; n < NFeed; n++) {
		if ((Feed[n].dir != 'X') && (Feed[n].dir != 'Y') && (Feed[n].dir != 'Z')) {
			printf("%s #%d\n", "*** invalid feed direction", n + 1);
			return 1;
		}
	}
	for (int n = 0; n < NPoint; n++) {
		if ((Point[n].dir != 'X') && (Point[n].dir != 'Y') && (Point[n].dir != 'Z')) {
			printf("%s #%d\n", "*** invalid point direction", n + 1);
			return 1;
		}
	}
	for (int n = 0; n < nload; n++) {
		if (((dload[n] != 'X') && (dload[n] != 'Y') && (dload[n] != 'Z')) ||
		    ((cload[n] != 'R') && (cload[n] != 'C') && (cload[n] != 'L'))) {
			printf("*** invalid load parameter #%d\n", n + 1);
			return 1;
		}
	}

	// warnings
/*
	if ((NFreq1 <= 0) && (NFreq2 <= 0)) {
		printf("%s\n", "*** no frequency data");
	}
*/
	// PBC -> Mur-1st
	if ((iABC == 1) && (PBCx || PBCy || PBCz)) {
		printf("%s\n", "*** warning : PBC -> Mur-1st");
		iABC = 0;
	}

	// number of cells
	setup_cells(nxr, nyr, nzr, dxr, dyr, dzr);

	// node
	setup_node(nxr, nyr, nzr, xr, yr, zr, dxr, dyr, dzr);

	// cell center
	setup_center();

	// feed
	if (NFeed) {
		setup_feed(xfeed, yfeed, zfeed);
	}

	// plane wave
	if (IPlanewave) {
		setup_planewave();
	}

	// point
	if (NPoint) {
		setup_point(xpoint, ypoint, zpoint, strprop);
	}

	// load
	if (nload) {
		setup_load(nload, dload, xload, yload, zload, cload, pload, array_inc);
	}

	// fit geometry without thickness
	if (NGeometry) {
		fitgeometry();
	}

	// setup geometry lines (plot)
	if (NGeometry) {
		setup_geomline3d();
	}

	// free
	free(xr);
	free(yr);
	free(zr);
	free(dxr);
	free(dyr);
	free(dzr);
	if (NFeed) {
		free(xfeed);
		free(yfeed);
		free(zfeed);
	}
	if (NPoint) {
		free(xpoint);
		free(ypoint);
		free(zpoint);
	}

	return 0;
}
