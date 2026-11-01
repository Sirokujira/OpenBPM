#!/usr/bin/env python3
"""
check_fullvector_he11.py

完全ベクトル BPM の差分式 (bpm/wabpm.cpp の wabpm_fv_*) が、格子を細かくすると
ステップインデックスファイバの**ベクトル厳密解 HE11** に収束することを確かめる
収束検証スクリプト (CI では実行しない。要 numpy / scipy)。

BPM の虚軸伝搬では細かい格子が遅すぎるため、同じ差分式 (Stern の対角ブロック +
4 隅の交差項) を疎行列として組み、固有値問題を shift-invert で直接解く。
円周の階段近似の雑音を除くため、n^2 はセル内のサブピクセル平均を使う。

既定のケース (tests/test_fullvec.cpp と同じ):
  n1 = 1.45 / n2 = 1.0 (空気), a = 0.5 um, lambda = 1.55 um (V = 2.13)
  HE11 = 1.176447 (ベクトル厳密解), LP01 = 1.225891 (スカラー, 弱導波近似)

実測 (sub = 8):
  dx [um]   完全ベクトル (HE11 との差)   スカラー (LP01 との差)
  0.05      1.178368 (+1.92e-3)         1.225903 (+1.2e-5)
  0.025     1.177564 (+1.12e-3)         1.225921 (+3.0e-5)
  0.0125    1.177081 (+6.3e-4)          1.225871 (-2.0e-5)
  0.00833   1.176906 (+4.6e-4)          1.225874 (-1.7e-5)
  -> 完全ベクトルは HE11 へ 1 次収束 (線形外挿 1.176556, 差 +1.1e-4)。
     スカラーは LP01 に一致し、HE11 とは 0.049 離れたまま (偏波補正を表せない)。

使い方:
  python3 tools/check_fullvector_he11.py [dx_um ...]    (既定 0.05 0.025 0.0125)
"""

import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.optimize import brentq
from scipy.special import jv, jvp, kv, kvp

LAM = 1.55e-6
K0 = 2 * np.pi / LAM
N1, N2, A = 1.45, 1.0, 0.5e-6
HALF = 2e-6          # 計算領域の半幅 (Dirichlet 境界)


def exact_roots():
    """HE11 (ベクトル) と LP01 (スカラー) の厳密 neff"""
    def uw(ne):
        return A * K0 * np.sqrt(N1**2 - ne**2), A * K0 * np.sqrt(ne**2 - N2**2)

    def he(ne):
        U, W = uw(ne)
        jt = jvp(1, U) / (U * jv(1, U))
        kt = kvp(1, W) / (W * kv(1, W))
        r = (N2 / N1)**2
        return (jt + kt) * (jt + r * kt) - (1 / U**2 + 1 / W**2) * (1 / U**2 + r / W**2)

    def lp(ne):
        U, W = uw(ne)
        return U * jv(1, U) / jv(0, U) - W * kv(1, W) / kv(0, W)

    out = []
    for f in (he, lp):
        xs = np.linspace(N1 - 1e-9, N2 + 1e-9, 20001)
        fs = [f(x) for x in xs]
        for i in range(len(xs) - 1):
            if np.sign(fs[i]) != np.sign(fs[i + 1]) and abs(fs[i]) < 1e3 and abs(fs[i + 1]) < 1e3:
                out.append(brentq(f, xs[i + 1], xs[i]))
                break
    return out


def index_map(dx, sub):
    n = int(round(2 * HALF / dx))
    c = (np.arange(n) - (n - 1) / 2) * dx
    X, Y = np.meshgrid(c, c)                       # [iy, ix]
    o = (np.arange(sub) + 0.5) / sub - 0.5
    inside = np.zeros_like(X)
    for u in o:
        for v in o:
            inside += ((X + u * dx)**2 + (Y + v * dx)**2) < A * A
    inside /= sub * sub
    return n, N2**2 + (N1**2 - N2**2) * inside


def second_diff(n, n2, dx, along_x, stern):
    """1 方向の 2 階差分 (stern = True なら d/du[(1/n^2) d(n^2 E)/du])"""
    rows, cols, vals = [], [], []
    for iy in range(n):
        for ix in range(n):
            i = ix + iy * n
            nc = n2[iy, ix]
            for s in (-1, 1):
                jx, jy = (ix + s, iy) if along_x else (ix, iy + s)
                inside = (0 <= jx < n) and (0 <= jy < n)
                nn = n2[jy, jx] if inside else nc
                w, wb = (2 * nn / (nc + nn), -2 * nc / (nc + nn)) if stern else (1.0, -1.0)
                rows.append(i); cols.append(i); vals.append(wb / dx**2)
                if inside:
                    rows.append(i); cols.append(jx + jy * n); vals.append(w / dx**2)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n * n, n * n))


def cross_terms(n, n2, dx):
    """Pxy, Pyx (bpm/wabpm.cpp の apply_cross と同じ 4 隅の中心差分)"""
    rows, cols, vxy, vyx = [], [], [], []
    for iy in range(n):
        for ix in range(n):
            for sx in (-1, 1):
                jx = ix + sx
                if not 0 <= jx < n:
                    continue
                for sy in (-1, 1):
                    jy = iy + sy
                    if not 0 <= jy < n:
                        continue
                    nk = n2[jy, jx]
                    rows.append(ix + iy * n); cols.append(jx + jy * n)
                    vxy.append(sx * sy * (nk / n2[iy, jx] - 1) / (4 * dx * dx))
                    vyx.append(sx * sy * (nk / n2[jy, ix] - 1) / (4 * dx * dx))
    shape = (n * n, n * n)
    return sp.csr_matrix((vxy, (rows, cols)), shape=shape), sp.csr_matrix((vyx, (rows, cols)), shape=shape)


def solve(dx, sub, full):
    n, n2 = index_map(dx, sub)
    V = sp.diags(K0**2 * (n2.ravel() - N2**2))
    if full:
        pxx = second_diff(n, n2, dx, True, True) + second_diff(n, n2, dx, False, False) + V
        pyy = second_diff(n, n2, dx, True, False) + second_diff(n, n2, dx, False, True) + V
        pxy, pyx = cross_terms(n, n2, dx)
        P = sp.bmat([[pxx, pxy], [pyx, pyy]]).tocsc()
        guess = 1.18
    else:
        P = (second_diff(n, n2, dx, True, False) + second_diff(n, n2, dx, False, False) + V).tocsc()
        guess = 1.226
    w = sla.eigs(P, k=2, sigma=K0**2 * (guess**2 - N2**2), which='LM', return_eigenvectors=False)
    return float(np.max(np.sqrt(N2**2 + np.real(w) / K0**2)))


def main():
    dxs = [float(v) for v in sys.argv[1:]] or [0.05, 0.025, 0.0125]
    he11, lp01 = exact_roots()
    print(f"exact : HE11 = {he11:.6f}, LP01 = {lp01:.6f} (polarization correction {he11 - lp01:+.5f})")
    for d in dxs:
        fv = solve(d * 1e-6, 8, True)
        sc = solve(d * 1e-6, 8, False)
        print(f"dx = {d:.5f} um : full-vector {fv:.6f} ({fv - he11:+.2e} vs HE11), "
              f"scalar {sc:.6f} ({sc - lp01:+.2e} vs LP01)", flush=True)


if __name__ == "__main__":
    main()
