#include "obpm.h"
#include "user_define_complex.h"
#include "obpm_prototype.h"
#include "ev.h"

#include "hdf5.h"
#define FILE_NAME "time_series_data.h5"

// 光活性化関数曲線 (powersweep 指定時に出力)
#define FN_activation "activation_curve.csv"

#include "bpm/bpm_prototype.h"
#include "bpm/ElectricFieldProfile.hpp"
#include "bpm/RefractiveIndexProfile.hpp"
#include "bpm/wabpm.h"

#include <complex>



// z 断面の統計量 (GUI 表示用の /trace) を計算する。
// getI(ix, iy) は |E|^2 を返す呼び出し可能オブジェクト。
//   power  : 断面の総パワー (面積要素 dA を掛ける。物理スケーリング時は [W])
//   peak   : |E|^2 の最大値
//   cx, cy : 強度重心 [m]
//   wx, wy : 強度の 2 次モーメント幅 2*sigma [m]
template <typename F>
static void computeTrace(int Nx, int Ny, const double *xc, const double *yc,
                         double dA, F getI,
                         double *power, double *peak,
                         double *cx, double *cy, double *wx, double *wy)
{
    double s = 0, sx = 0, sy = 0, sxx = 0, syy = 0, pk = 0;
    for (int iy = 0; iy < Ny; iy++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double iv = getI(ix, iy);   // I は虚数単位マクロのため別名にする
            s   += iv;
            sx  += xc[ix] * iv;
            sy  += yc[iy] * iv;
            sxx += xc[ix] * xc[ix] * iv;
            syy += yc[iy] * yc[iy] * iv;
            if (iv > pk) pk = iv;
        }
    }
    *power = s * dA;
    *peak  = pk;
    if (s > 0) {
        const double mx = sx / s, my = sy / s;
        const double vx = (sxx / s) - (mx * mx);
        const double vy = (syy / s) - (my * my);
        *cx = mx;
        *cy = my;
        *wx = 2 * sqrt(vx > 0 ? vx : 0);
        *wy = 2 * sqrt(vy > 0 ? vy : 0);
    }
    else {
        *cx = *cy = *wx = *wy = 0;
    }
}

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

	// input refractive index(Array) : zスライス毎に材料 ID から更新する
    P->n_in = (floatcomplex *)malloc((P->Nx * P->Ny)*sizeof(floatcomplex));
    // 屈折率分布はスライス毎の 2D 配列として与える
    P->Nx_n = Nx;
    P->Ny_n = Ny;
    P->Nz_n = 1;
    P->dz_n = P->dz;
    // テーパ (taper = 出口/入口の横方向スケール比) とツイスト (twist = [deg/m]) :
    // 屈折率分布を座標変換 (相似縮小 + 回転) で z に沿って変化させる。
    //   カーネル規約 : スケール係数 = 1 - taperPerStep*iz、回転角 = twistPerStep*iz [rad]
    // 変換時の屈折率分布は**入力断面 (iz_start) の 2D 分布**を参照する
    // (geometry 自体の z 変化とは併用しない)。
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
    // モード解析 (modes = <nModes> [excite])
    //   入力断面 (iz_start) の屈折率分布から虚軸伝搬法 (bpm/modes.cpp) で
    //   導波モードを求め、neff をログへ、モード形状を HDF5 (/modes) へ出力する。
    //   excite 指定時は入射ビームを基本モードで置き換える (モード整合励振)。
    //   曲げ (bend) は反映しない (直線導波路の入力断面のモード)。
    //   半ベクトルの deflation は近似のためスカラー演算子で解析する。
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
        // 虚軸ステップ幅 : 最大固有値 mu_max に対し a*mu_max ~ 0.5 となるよう選ぶ
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

        // モード整合励振 : 入射ビームを基本モードで置き換える。
        // 振幅はガウシアン励振と同じ規約 (ピーク振幅 = 給電電圧)、
        // beamtilt の位相ランプは維持する。
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

    // 伝搬の可視化用 (追加) : 中心列 (x = Nx/2) の強度 |E(y, z)|^2 と、
    // z ごとのスカラー推移 (GUI でのグラフ表示用, /trace へ出力)
    const long ntr = P->iz_end - P->iz_start;
    float  *Iyz = (float *)malloc((size_t)P->Ny * ntr * sizeof(float));
    double *trZ      = (double *)malloc((size_t)ntr * sizeof(double));
    double *trPower  = (double *)malloc((size_t)ntr * sizeof(double));
    double *trPeak   = (double *)malloc((size_t)ntr * sizeof(double));
    double *trCx     = (double *)malloc((size_t)ntr * sizeof(double));
    double *trCy     = (double *)malloc((size_t)ntr * sizeof(double));
    double *trWx     = (double *)malloc((size_t)ntr * sizeof(double));
    double *trWy     = (double *)malloc((size_t)ntr * sizeof(double));

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
        W.symx = BPM.symx;
        W.symy = BPM.symy;

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
            // テーパ/ツイスト指定時は入力断面を座標変換 (相似縮小 + 回転) して参照する
            // (近軸カーネル bpm/FDBPMpropagator.c と同一の変換式)
            const long izr = iz - P->iz_start;
            const double xfScale = 1.0 / (1.0 - (double)P->taperPerStep * izr);
            const double xfCos = cos(-(double)P->twistPerStep * izr);
            const double xfSin = sin(-(double)P->twistPerStep * izr);
            for (int iy = 0; iy < P->Ny; iy++) {
                for (int ix = 0; ix < P->Nx; ix++) {
                    long index = (long)(P->Nx * iy) + ix;
                    double nr, nimag;
                    if (xform) {
                        // 出力格子点に対応する入力断面上の位置を双一次補間で取得
                        const double x = P->dx * (ix - (P->Nx - 1) / 2.0);
                        const double y = P->dy * (iy - (P->Ny - 1) / 2.0);
                        const double xs = xfScale * ((xfCos * x) - (xfSin * y));
                        const double ys = xfScale * ((xfSin * x) + (xfCos * y));
                        double fx = (xs / P->dx) + ((P->Nx - 1) / 2.0);
                        double fy = (ys / P->dy) + ((P->Ny - 1) / 2.0);
                        fx = MIN(MAX(0.0, fx), (P->Nx - 1) * (1.0 - 1e-7));
                        fy = MIN(MAX(0.0, fy), (P->Ny - 1) * (1.0 - 1e-7));
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
                        const double x = P->dx * (ix - (P->Nx - 1) / 2.0);
                        const double y = P->dy * (iy - (P->Ny - 1) / 2.0);
                        const double xb = x * cosB + y * sinB;
                        nr = nr * (1.0 - nr * nr * xb * BPM.rho_e / (2.0 * BPM.RoC))
                                * exp(xb / BPM.RoC);
                    }
                    std::complex<double> n(nr, -nimag);
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

            // 中心行 / 中心列の強度と z ごとのスカラー推移を記録
            {
                const long t = iz - P->iz_start;
                const long row0 = (long)(P->Ny / 2) * P->Nx;
                for (long ix = 0; ix < P->Nx; ix++) {
                    Ixz[t * P->Nx + ix] = (float)std::norm(Ed[row0 + ix]);
                }
                const long col0 = P->Nx / 2;
                for (long iy = 0; iy < P->Ny; iy++) {
                    Iyz[t * P->Ny + iy] = (float)std::norm(Ed[col0 + iy * P->Nx]);
                }
                trZ[t] = Zn[0] + ((iz + 1) * Dz);
                computeTrace(P->Nx, P->Ny, Xc, Yc, Dx * Dy,
                             [&](int ix, int iy) { return std::norm(Ed[ix + (long)iy * P->Nx]); },
                             &trPower[t], &trPeak[t], &trCx[t], &trCy[t], &trWx[t], &trWy[t]);
            }

            // スナップショットを記録
            if (Frames && ((iz - P->iz_start) % frameInterval == 0)) {
                const size_t off = (size_t)((iz - P->iz_start) / frameInterval) * P->Nx * P->Ny;
                float *f = &Frames[off];
                for (long i = 0; i < (long)Nx * Ny; i++) {
                    f[i] = (float)std::norm(Ed[i]);
                }
                if (framesCplx) {
                    for (long i = 0; i < (long)Nx * Ny; i++) {
                        FramesR[off + i] = (float)Ed[i].real();
                        FramesI[off + i] = (float)Ed[i].imag();
                    }
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
        // テーパ/ツイスト指定時は入力断面を渡し、z 変化はカーネル側の座標変換に任せる
        const int64_t iz_n = xform ? P->iz_start : iz;
        for (int iy = 0; iy < P->Ny; iy++)
        {
            for (int ix = 0; ix < P->Nx; ix++)
            {
                long index = (long)(P->Nx * iy) + ix;
                int64_t nn = NA(ix, iy, iz_n);
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

        // 中心行 / 中心列の強度と z ごとのスカラー推移を記録 (結果は E2 にある)
        {
            const long t = iz - P->iz_start;
            const long row0 = (long)(P->Ny / 2) * P->Nx;
            for (long ix = 0; ix < P->Nx; ix++) {
                const floatcomplex e = P->E2[row0 + ix];
                Ixz[t * P->Nx + ix] = (CREALF(e) * CREALF(e)) + (CIMAGF(e) * CIMAGF(e));
            }
            const long col0 = P->Nx / 2;
            for (long iy = 0; iy < P->Ny; iy++) {
                const floatcomplex e = P->E2[col0 + iy * P->Nx];
                Iyz[t * P->Ny + iy] = (CREALF(e) * CREALF(e)) + (CIMAGF(e) * CIMAGF(e));
            }
            trZ[t] = Zn[0] + ((iz + 1) * Dz);
            computeTrace(P->Nx, P->Ny, Xc, Yc, Dx * Dy,
                         [&](int ix, int iy) {
                             const floatcomplex e = P->E2[ix + (long)iy * P->Nx];
                             return ((double)CREALF(e) * CREALF(e)) + ((double)CIMAGF(e) * CIMAGF(e));
                         },
                         &trPower[t], &trPeak[t], &trCx[t], &trCy[t], &trWx[t], &trWy[t]);
        }

        // スナップショットを記録
        if (Frames && ((iz - P->iz_start) % frameInterval == 0)) {
            const size_t off = (size_t)((iz - P->iz_start) / frameInterval) * P->Nx * P->Ny;
            float *f = &Frames[off];
            for (long i = 0; i < (long)P->Nx * P->Ny; i++) {
                const floatcomplex e = P->E2[i];
                f[i] = (CREALF(e) * CREALF(e)) + (CIMAGF(e) * CIMAGF(e));
                if (framesCplx) {
                    FramesR[off + i] = CREALF(e);
                    FramesI[off + i] = CIMAGF(e);
                }
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

        // 伝搬マップ |E(x=Nx/2, y, z)|^2 の書き込み (Ixz の y 版)
        {
            hsize_t iyz_dims[2] = {(hsize_t)ntr, (hsize_t)P->Ny};
            dataspace_id = H5Screate_simple(2, iyz_dims, NULL);
            dataset_id = H5Dcreate(field_group_id, "Iyz", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Iyz);
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
            // 複素電界 (位相表示用)
            if (framesCplx) {
                const char *cn[2] = {"frames_r", "frames_i"};
                float *cd[2] = {FramesR, FramesI};
                for (int c = 0; c < 2; c++) {
                    dataspace_id = H5Screate_simple(3, fr_dims, NULL);
                    dataset_id = H5Dcreate(field_group_id, cn[c], H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cd[c]);
                    H5Dclose(dataset_id);
                    H5Sclose(dataspace_id);
                }
            }
        }
        H5Gclose(field_group_id);
    }

    // z ごとのスカラー推移 (/trace) の書き込み : GUI での伝搬グラフ表示用
    // すべて長さ ntr (= 伝搬ステップ数) の 1 次元配列で、trace[t] は
    // z = trace/z[t] における値 (t 番目のステップ実行後の断面)
    {
        hid_t trace_group_id = H5Gcreate(file_id, "/trace", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        struct {
            const char *name;
            const double *data;
        } trace_data[] = {
            {"z",          trZ},       // 断面位置 [m]
            {"power",      trPower},   // 断面の総パワー (物理スケーリング時は [W])
            {"peak",       trPeak},    // |E|^2 の最大値
            {"centroid_x", trCx},      // 強度重心 x [m]
            {"centroid_y", trCy},      // 強度重心 y [m]
            {"width_x",    trWx},      // 強度の 2 次モーメント幅 2*sigma x [m]
            {"width_y",    trWy}       // 同 y [m]
        };
        hsize_t tr_dims[1] = {(hsize_t)ntr};
        for (size_t n = 0; n < sizeof(trace_data) / sizeof(trace_data[0]); n++) {
            dataspace_id = H5Screate_simple(1, tr_dims, NULL);
            dataset_id = H5Dcreate(trace_group_id, trace_data[n].name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, trace_data[n].data);
            H5Dclose(dataset_id);
            H5Sclose(dataspace_id);
        }
        H5Gclose(trace_group_id);
    }

    // モード解析結果 (/modes) の書き込み : mode<i> (nModesFound x Ny x Nx) と neff
    if (nModesFound > 0) {
        hid_t modes_group_id = H5Gcreate(file_id, "/modes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        const long NN2 = (long)P->Nx * P->Ny;
        float *tmp = (float *)malloc((size_t)NN2 * sizeof(float));
        char dsname[32];
        for (int m = 0; m < nModesFound; m++) {
            // モードは実数 (無損失断面) だが規約に合わせ実部を格納する
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
    free(Iyz);
    free(trZ); free(trPower); free(trPeak);
    free(trCx); free(trCy); free(trWx); free(trWy);
    if (Frames) free(Frames);
    if (FramesR) free(FramesR);
    if (FramesI) free(FramesI);
    if (modeFields) delete[] modeFields;
    if (modeNeff) delete[] modeNeff;

    return;
}

