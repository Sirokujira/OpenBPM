#include "obpm.h"
#include "user_define_complex.h"
#include "obpm_prototype.h"
#include "ev.h"

#include "hdf5.h"
#define FILE_NAME "time_series_data.h5"

// 光活性化関数曲線 (powersweep 指定時に出力)
#define FN_activation "activation_curve.csv"

#include "bpm/bpm_prototype.h"
#include "bpm/wabpm.h"

#include <complex>


void solve_bpm(int io, double *tdft, FILE *fp) {
    // HDF5ファイルの作成
    // 関数から?(fp の入替え?)
    hid_t file_id;
    // local
    hid_t dataset_id, dataspace_id;
    herr_t status;

    char str[BUFSIZ];

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
        };
        for (size_t n = 0; n < sizeof(ignored) / sizeof(ignored[0]); n++) {
            if (ignored[n].cond) {
                sprintf(str, "*** warning : keyword '%s' is ignored by the BPM solver.", ignored[n].name);
                if (io) fprintf(fp, "%s\n", str);
                fprintf(stdout, "%s\n", str);
            }
        }
        // 波長掃引は未対応 (先頭周波数のみ使用)
        if ((NFreq2 > 1) || ((NFreq2 == 0) && (NFreq1 > 1))) {
            sprintf(str, "*** warning : multiple frequencies are not swept by the BPM solver (only the first is used).");
            if (io) fprintf(fp, "%s\n", str);
            fprintf(stdout, "%s\n", str);
        }
    }


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
    //   1 cm/GW = 1e-2 m / 1e9 W = 1e-11 m/W (例 : 424 cm/GW = 4.24e-9 m/W)
    // 出典 : Honda, Shoji, Amemiya, Opt. Lett. 49, 5811 (2024)
    //        (メタマテリアル装荷 Si 導波路の TPA による光活性化関数)
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
    P->dx = (float)Dx; // x軸の区切りデータの1step幅
    P->dy = (float)Dy; // y軸の区切りデータの1step幅
    P->dz = (float)Dz; // z軸の区切りデータの1step幅
    P->iz_start = 0; // z軸の区切りデータ先頭(index)
    P->iz_end = Nz; // z軸の区切りデータ末尾(index)
    //enum symmetry : unsignate char
    //  enumeration
    //    NoSymmetry   (0)
    //    Symmetry     (1)
    //    AntiSymmetry (2)
    //  end
    //end
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

	// input refractive index(Array) : zスライス毎に材料 ID から更新する
    P->n_in = (floatcomplex *)malloc((P->Nx * P->Ny)*sizeof(floatcomplex));
    // 屈折率分布はスライス毎の 2D 配列として与える
    P->Nx_n = Nx;
    P->Ny_n = Ny;
    P->Nz_n = 1;
    P->dz_n = P->dz;
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

    // 電界バッファ
	P->E1 = (floatcomplex *)malloc((P->Nx * P->Ny)*sizeof(floatcomplex));
	P->E2 = (floatcomplex *)malloc((P->Nx * P->Ny)*sizeof(floatcomplex));
    P->Efinal = (floatcomplex *)malloc((P->Nx * P->Ny)*sizeof(floatcomplex)); // Output E field(Array)
    P->n_out = (floatcomplex *)malloc((P->Nx * P->Ny)*sizeof(floatcomplex)); //Output refractive index(Array)
	P->multiplier = (float *)malloc((P->Nx * P->Ny)*sizeof(float));

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

    // モード励振 (launch = mode <m>) : 先頭スライスの屈折率分布から導波モードを
    // 虚軸伝搬法 (bpm/modes.cpp) で求め、ガウシアンの代わりに励振する。
    // 振幅は feed 電圧 (モードは L2 正規化済みのため入力電力 = volt^2)。
    // beamtilt 指定時はガウシアンと同じ位相ランプをモード界にも掛ける。
    double mode_neff = 0;   // HDF5 /metadata/mode_neff (モード励振時のみ > 0)
    if (BPM.launchMode >= 0) {
        const long Nm = (long)Nx * Ny;
        std::complex<double> *n2m = new std::complex<double>[Nm];
        double nmax = 0;
        for (int iy = 0; iy < Ny; iy++) {
            for (int ix = 0; ix < Nx; ix++) {
                const floatcomplex nmv = n_mat[iEx[NA(ix, iy, 0)]];
                std::complex<double> nn(CREALF(nmv), -CIMAGF(nmv));
                n2m[(long)Nx * iy + ix] = nn * nn;
                if (nn.real() > nmax) nmax = nn.real();
            }
        }
        int launched = 0;
        if (nmax > P->n_0) {
            const int m = BPM.launchMode;
            wabpm_params Wm;
            Wm.Nx = Nx;
            Wm.Ny = Ny;
            Wm.dx = Dx;
            Wm.dy = Dy;
            Wm.k0 = k0;
            Wm.n0 = P->n_0;
            Wm.wideangle = 0;
            Wm.pol = BPM.pol;
            // a*mu_max ~ 0.5 となる虚軸ステップ幅 (tests/test_modes.cpp と同じ)
            const double mu_max = k0 * k0 * ((nmax * nmax) - ((double)P->n_0 * P->n_0));
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
                double power = 0;
                for (long i = 0; i < Nm; i++) {
                    const long ix = i % Nx;
                    const long iy = i / Nx;
                    const double ph = -(ktx * (Xc[ix] - xc0)) - (kty * (Yc[iy] - yc0));
                    const std::complex<double> em = volt * modes[(size_t)m * Nm + i]
                                                  * std::complex<double>(cos(ph), sin(ph));
                    floatcomplex e = {(float)em.real(), (float)em.imag()};
                    P->E1[i] = e;
                    power += std::norm(em);
                }
                P->precisePower = power;
                mode_neff = neff[m];
                sprintf(str, "launch : mode %d, neff = %.6f (found %d)", m, neff[m], nFound);
                if (io) fprintf(fp, "%s\n", str);
                fprintf(stdout, "%s\n", str);
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
    // 非線形吸収 (TPA) / 入力パワー掃引 (ONN 光活性化関数) の準備
    //   - 物理スケーリング : 初期界を ∫∫|E|^2 dA = P_in [W] になるよう正規化し、
    //     |E|^2 がそのまま強度 I [W/m^2] になるようにする (面積要素 dA = Dx*Dy)。
    //     tpa / powersweep 未指定時は従来通り無次元の界のまま計算する (後方互換)。
    //   - 実効断面積 A_eff = (∫|E|^2 dA)^2 / ∫|E|^4 dA を初期界から実計算する
    //     (平面波近似の解析解 T = 1/(1 + beta*(P_in/A_eff)*L) との比較用)
    // ============================================================
    const int    sweep    = (PowerSweep.npoints > 0);
    const int    nsweep   = sweep ? PowerSweep.npoints : 1;
    const int    physcale = (sweep || haveTpa);   // 物理スケーリングの有無
    const double E0rawpow = P->precisePower;      // 初期界の生の |E|^2 和 (面積要素なし)
    floatcomplex *E0 = NULL;                      // 初期界の保存 (掃引毎の再初期化用)
    double *sweepPin = NULL, *sweepPout = NULL;
    if (physcale) {
        E0 = (floatcomplex *)malloc((P->Nx * P->Ny) * sizeof(floatcomplex));
        memcpy(E0, P->E1, (P->Nx * P->Ny) * sizeof(floatcomplex));
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

    P->EfieldPower = 0;
    P->precisePowerDiff = 0;
    #ifdef _OPENMP
    bool useAllCPUs = false;
    long numThreads = useAllCPUs || omp_get_num_procs() == 1? omp_get_num_procs(): omp_get_num_procs()-1;
    #else
    long numThreads = 1;
    #endif
    P->b = (floatcomplex *)malloc(numThreads*MAX(P->Nx,P->Ny)*sizeof(floatcomplex));

    // 解放時の二重 free 防止用 (swapEPointers でポインタが入れ替わるため)
    floatcomplex *E_in = P->E1;
    floatcomplex *E_buf = P->E2;

    fprintf(stdout, "initfield\n");

    // initial field
    initfield();

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

    // 伝搬の可視化用 : 中心行 (y = Ny/2) の強度 |E(x, z)|^2 を全ステップ記録する
    float *Ixz = (float *)malloc((size_t)P->Nx * (P->iz_end - P->iz_start) * sizeof(float));

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
    // (掃引 2 回目以降のポインタ・電力簿記のリセットを兼ねる)
    if (physcale) {
        const double Pin = sweepPin[isweep];
        const float s = (float)sqrt(Pin / (E0rawpow * Dx * Dy));
        P->E1 = E_in;
        P->E2 = E_buf;
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
        // 拡張パス : 広角 BPM (Pade(1,1)) / 半ベクトル BPM (倍精度)
        // 屈折率項を演算子の内部に含めるため、位相 multiplier ではなく
        // 一般化 ADI (bpm/wabpm.cpp) で伝搬する。曲げは等価屈折率法で
        // n2 スライスへ反映する (近軸パスと同じ変換、実部のみ)。
        // ============================================================
        sprintf(str, "mode : %s, polarization = %s",
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

        std::complex<double> *Ed = new std::complex<double>[(size_t)Nx * Ny];
        std::complex<double> *n2 = new std::complex<double>[(size_t)Nx * Ny];

        // 曲げ (等価屈折率法) : n_bend = n*(1 - n^2*xb*rho_e/(2*RoC))*exp(xb/RoC)
        // (xb = x*cos(dir) + y*sin(dir), 吸収 = 虚部は変換しない)
        const int  bend = (BPM.RoC > 0);
        const double cosB = bend ? cos(BPM.bendDir * DTOR) : 1.0;
        const double sinB = bend ? sin(BPM.bendDir * DTOR) : 0.0;

        // 初期電界 (チルト位相込みの E1) を倍精度へ
        for (long i = 0; i < (long)Nx * Ny; i++) {
            Ed[i] = std::complex<double>(CREALF(P->E1[i]), CIMAGF(P->E1[i]));
        }

        for (long iz = P->iz_start; iz < P->iz_end; iz++) {
            sprintf(str, "step iz : %ld / %ld", iz + 1, (long)Nz);
            fprintf(stdout, "%s\n", str);

            const double t0 = cputime();

            // このスライスの複素比誘電率 (n_mat は損失を +imag で保持 -> 物理符号 n = nr - i*ni に戻す)
            for (int iy = 0; iy < P->Ny; iy++) {
                for (int ix = 0; ix < P->Nx; ix++) {
                    long index = (long)(P->Nx * iy) + ix;
                    int64_t nn = NA(ix, iy, iz);
                    const floatcomplex nm = n_mat[iEx[nn]];
                    double nr = CREALF(nm);
                    if (bend) {
                        const double x = P->dx * (ix - (P->Nx - 1) / 2.0);
                        const double y = P->dy * (iy - (P->Ny - 1) / 2.0);
                        const double xb = x * cosB + y * sinB;
                        nr = nr * (1.0 - nr * nr * xb * BPM.rho_e / (2.0 * BPM.RoC))
                                * exp(xb / BPM.RoC);
                    }
                    std::complex<double> n(nr, -CIMAGF(nm));
                    n2[index] = n * n;
                }
            }

            // 1 ステップ伝搬 + 端部吸収体
            wabpm_step(&W, Ed, n2);
            for (long i = 0; i < (long)Nx * Ny; i++) {
                Ed[i] *= P->multiplier[i];
            }

            // TPA (二光子吸収) : E *= exp(-(beta/2)*I*dz), I = |E|^2 [W/m^2]
            // (物理スケーリング済みの界を仮定。強度の減衰率 alpha = beta*I に対し
            //  界には alpha/2 を適用する)
            if (haveTpa) {
                for (int iy = 0; iy < P->Ny; iy++) {
                    for (int ix = 0; ix < P->Nx; ix++) {
                        const long index = (long)(P->Nx * iy) + ix;
                        const double beta = tpa_mat[iEx[NA(ix, iy, iz)]];
                        if (beta > 0) {
                            const double i2 = std::norm(Ed[index]);
                            Ed[index] *= exp(-0.5 * beta * i2 * Dz);
                        }
                    }
                }
            }

            // 中心行の強度を記録
            {
                const long row0 = (long)(P->Ny / 2) * P->Nx;
                for (long ix = 0; ix < P->Nx; ix++) {
                    Ixz[(iz - P->iz_start) * P->Nx + ix] = (float)std::norm(Ed[row0 + ix]);
                }
            }

            // スナップショットを記録
            if (Frames && ((iz - P->iz_start) % frameInterval == 0)) {
                float *f = &Frames[(size_t)((iz - P->iz_start) / frameInterval) * P->Nx * P->Ny];
                for (long i = 0; i < (long)Nx * Ny; i++) {
                    f[i] = (float)std::norm(Ed[i]);
                }
            }

            *tdft += cputime() - t0;
        }

        // 最終電界・屈折率分布・出力電力を格納
        double power = 0;
        for (long i = 0; i < (long)Nx * Ny; i++) {
            floatcomplex e = {(float)Ed[i].real(), (float)Ed[i].imag()};
            P->Efinal[i] = e;
            power += std::norm(Ed[i]);
        }
        P->precisePower = power;
        for (int iy = 0; iy < P->Ny; iy++) {
            for (int ix = 0; ix < P->Nx; ix++) {
                long index = (long)(P->Nx * iy) + ix;
                P->n_out[index] = n_mat[iEx[NA(ix, iy, P->iz_end - 1)]];
            }
        }
        delete[] Ed;
        delete[] n2;
    }
    else {
    // BPM は時間反復ではなく Z 軸方向に 1 回伝搬する
    for(long iz = P->iz_start; iz < P->iz_end; iz++) {
        sprintf(str, "step iz : %ld / %ld", iz + 1, (long)Nz);
        fprintf(stdout, "%s\n", str);

        const double t0 = cputime();

        // このスライスの複素屈折率分布を材料 ID から設定 (虚部 = 導電率による吸収)
        for (int iy = 0; iy < P->Ny; iy++)
        {
            for (int ix = 0; ix < P->Nx; ix++)
            {
                long index = (long)(P->Nx * iy) + ix;
                int64_t nn = NA(ix, iy, iz);
                P->n_in[index] = n_mat[iEx[nn]];
            }
        }

        // Douglas-Gunn ADI 法による 1 ステップ伝搬
        #ifdef _OPENMP
        #pragma omp parallel num_threads(numThreads)
        #endif
        {
            substep1a(P);
            substep1b(P);
            substep2a(P);
            substep2b(P);
            applyMultiplier(P, iz, NULL);
        }

        // TPA (二光子吸収) : E *= exp(-(beta/2)*I*dz), I = |E|^2 [W/m^2]
        // (物理スケーリング済みの界を仮定。強度の減衰率 alpha = beta*I に対し
        //  界には alpha/2 を適用する)
        // 次ステップの fieldCorrection (= sqrt(precisePower/EfieldPower)) が
        // TPA 減衰を打ち消さないよう、電力簿記 precisePower も同率で減らす。
        if (haveTpa) {
            updatePrecisePower(P);  // 係留中の precisePowerDiff を先に確定する
            double pb = 0, pa = 0;
            for (int iy = 0; iy < P->Ny; iy++) {
                for (int ix = 0; ix < P->Nx; ix++) {
                    const long index = (long)(P->Nx * iy) + ix;
                    const floatcomplex e = P->E2[index];
                    const double i2 = ((double)CREALF(e) * CREALF(e)) + ((double)CIMAGF(e) * CIMAGF(e));
                    const double beta = tpa_mat[iEx[NA(ix, iy, iz)]];
                    pb += i2;
                    if (beta > 0) {
                        const float g = (float)exp(-0.5 * beta * i2 * Dz);
                        P->E2[index] = e * g;
                        pa += i2 * (double)g * g;
                    }
                    else {
                        pa += i2;
                    }
                }
            }
            if (pb > 0) P->precisePower *= pa / pb;
        }

        // 中心行の強度を記録 (このステップの結果は E2 にある)
        {
            const long row0 = (long)(P->Ny / 2) * P->Nx;
            for (long ix = 0; ix < P->Nx; ix++) {
                const floatcomplex e = P->E2[row0 + ix];
                Ixz[(iz - P->iz_start) * P->Nx + ix] =
                    (CREALF(e) * CREALF(e)) + (CIMAGF(e) * CIMAGF(e));
            }
        }

        // スナップショットを記録
        if (Frames && ((iz - P->iz_start) % frameInterval == 0)) {
            float *f = &Frames[(size_t)((iz - P->iz_start) / frameInterval) * P->Nx * P->Ny];
            for (long i = 0; i < (long)P->Nx * P->Ny; i++) {
                const floatcomplex e = P->E2[i];
                f[i] = (CREALF(e) * CREALF(e)) + (CIMAGF(e) * CIMAGF(e));
            }
        }

        // E1/E2 の入替え : swapEPointers はステップ数が奇数のとき新規バッファを
        // malloc するため、パワー掃引の複数回実行ではバッファが増えないよう
        // 単純スワップに置き換える (最終結果はループ後に Efinal へコピーされる。
        // 演算内容は同一で、単発実行の結果も変わらない)
        if (iz + 1 < P->iz_end) {
            P->EfieldPower = 0;
            floatcomplex *etmp = P->E1;
            P->E1 = P->E2;
            P->E2 = etmp;
        }
        updatePrecisePower(P);

        *tdft += cputime() - t0;
    }

    // 最終電界を Efinal へ格納
    if (P->E2 != P->Efinal) {
        memcpy(P->Efinal, P->E2, P->Nx * P->Ny * sizeof(floatcomplex));
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
    // (小パワーで線形透過、大パワーで TPA 飽和 -> ReLU 相当の応答曲線)
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

    // メタデータの作成
    hid_t metadata_group_id = H5Gcreate(file_id, "/metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // 時間に関するメタデータの書き込み(収束条件で終了時の対応)
    //int maxIter = MIN(Solver.maxiter, Ntime);
    //int maxNOut = MIN(Solver.nout, Niter);
    //double time_metadata[1] = {maxIter * Dt};
    double time_metadata[1] = {Solver.maxiter * Dt};
    //dataspace_id = H5Screate_simple(1, count, NULL);
    hsize_t time_count[1] = {1};
    dataspace_id = H5Screate_simple(1, time_count, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "time", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_metadata);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // グリッドに関するメタデータの書き込み
    //double grid_metadata[3] = {Dx, Dy, Dz};
    //hsize_t grid_count[1] = {3};
    //dataspace_id = H5Screate_simple(1, grid_count, NULL);
    //dataset_id = H5Dcreate(metadata_group_id, "grid", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid_metadata);
    //H5Dclose(dataset_id);
    //H5Sclose(dataspace_id);

    // その他のメタデータの書き込み
/*
    // title, dt, source, fPlanewave, z0, Ni, Nj, Nk, N0, NN
    const char *title = Title;
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "title", H5T_C_S1, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_C_S1, H5S_ALL, H5S_ALL, H5P_DEFAULT, title);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // 平面波が伝わるときのインピーダンスの値
    //double metadata_values[8] = {Dt, Planewave.z0, Ni, Nj, Nk, N0, NN};
    double metadata_values[8] = {Dt, 0.0, Ni, Nj, Nk, N0, NN};
    hsize_t metadata_count[1] = {8};
    dataspace_id = H5Screate_simple(1, metadata_count, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "metadata_values", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata_values);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // 配列に関するメタデータの書き込み (Xn, Yn, Zn, Freq1, Freq2)
    hsize_t array_count[1];
    double *arrays[] = {Xn, Yn, Zn, Freq1, Freq2};
    const char *array_names[] = {"Xn", "Yn", "Zn", "Freq1", "Freq2"};
    size_t array_sizes[] = {Nx + 1, Ny + 1, Nz + 1, NFreq1, NFreq2};

    for (int i = 0; i < 5; i++) {
        array_count[0] = array_sizes[i];
        dataspace_id = H5Screate_simple(1, array_count, NULL);
        dataset_id = H5Dcreate(metadata_group_id, array_names[i], H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arrays[i]);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
*/
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

    for (int i = 0; i < sizeof(metadata) / sizeof(metadata[0]); i++) {
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
        //{"Gline", Gline, NGline * 2 * 3}
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

    // free
    memfree2();

    // swapEPointers でポインタが入れ替わるため、重複を除外して解放する
    {
        floatcomplex *bufs[] = {E_in, E_buf, P->E1, P->E2, P->Efinal};
        const int nbuf = sizeof(bufs) / sizeof(bufs[0]);
        for (int i = 0; i < nbuf; i++) {
            int duplicated = 0;
            for (int j = 0; j < i; j++) {
                if (bufs[j] == bufs[i]) {
                    duplicated = 1;
                    break;
                }
            }
            if (!duplicated && (bufs[i] != NULL)) free(bufs[i]);
        }
    }
    free(P->n_out);
    free(P->n_in);
    free(P->multiplier);
    free(P->b);
    free(n_mat);
    free(tpa_mat);
    if (E0 != NULL) free(E0);
    if (sweepPin != NULL) free(sweepPin);
    if (sweepPout != NULL) free(sweepPout);
    free(Ixz);
    if (Frames) free(Frames);

    return;
}

