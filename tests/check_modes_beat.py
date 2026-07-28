#!/usr/bin/env python3
"""
多モード重ね合わせ励振 (launch = mode <m1> <m2>) の定量検証。

2 モードの重ね合わせは実効屈折率差でビート (モードビート) を起こし、
ビーム重心が周期 L = lambda / |neff_0 - neff_1| で横方向に振動する。
z=0 で 2 モードが同位相 (重心が振幅の極値) なので

    c(z) = c(0) * cos(2*pi*z / L)   ->   L = 2*pi*z_end / arccos(c(z_end)/c(0))

として伝搬マップ /field/Ixz の重心軌跡から L を実測し、ログに出力された
neff から求まる理論値と比較する。

あわせて電力保存 (係数に依らず入力電力に正規化されること) も検査する。

使い方:
  tests/check_modes_beat.py --h5 time_series_data.h5 --log obpm.log
                            [--dx-total <領域幅[m]>] [--z-total <伝搬長[m]>]
                            [--rtol 0.10]

標準ライブラリ + h5py/numpy を使用 (伝搬マップの読み出しに必要)。
"""
import argparse
import math
import os
import re
import sys

try:
    import h5py
    import numpy as np
except ImportError:
    sys.exit("h5py and numpy are required (pip install h5py numpy)")


def parse_log(path):
    """obpm.log から重ね合わせた各モードの neff を取り出す。"""
    neffs = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if "superposition of" in line:
                neffs = [float(v) for v in re.findall(r"neff=([0-9.]+)", line)]
    return neffs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5", default="time_series_data.h5")
    ap.add_argument("--log", default="obpm.log")
    ap.add_argument("--rtol", type=float, default=0.10)
    args = ap.parse_args()

    for p in (args.h5, args.log):
        if not os.path.exists(p):
            raise SystemExit(f"not found: {p}")

    neffs = parse_log(args.log)
    if len(neffs) < 2:
        raise SystemExit(f"{args.log}: 重ね合わせ (2 モード以上) のログが見つかりません")

    with h5py.File(args.h5, "r") as f:
        Ixz = f["/field/Ixz"][:]                   # (Nz, Nx)
        md = f["/metadata"]
        lam = float(md["lambda"][()])
        xc = md["Xc"][:] if "Xc" in md else None
        zc = md["Zc"][:] if "Zc" in md else None

    Nz, Nx = Ixz.shape
    x = xc if xc is not None else np.arange(Nx, dtype=float)
    # 各 z における重心
    p = Ixz.sum(axis=1)
    if np.any(p <= 0):
        raise SystemExit("伝搬マップに強度ゼロの断面があります")
    cen = (Ixz * x).sum(axis=1) / p

    z_end = float(zc[-1] - zc[0]) if zc is not None else float(Nz - 1)
    c0, c1 = cen[0], cen[-1]
    if abs(c0) < 1e-12:
        raise SystemExit("初期重心が 0 です (重ね合わせが非対称になっていません)")

    ratio = c1 / c0
    if not (-1.0 < ratio <= 1.0):
        raise SystemExit(f"重心比が範囲外です: {ratio}")

    # 位相 phi = arccos(ratio)。ratio ~ 1 は「重心が動いていない」= ビートが
    # 起きていないことを意味するので、無限大にせず明示的に失敗させる。
    phi = math.acos(ratio)
    if phi < 1e-6:
        print(f"FAIL: ビートが検出できません (重心が動いていない: "
              f"{c0 * 1e6:+.3f} -> {c1 * 1e6:+.3f} um)")
        print("=== FAILED (1 failure(s)) ===")
        return 1

    # 実測ビート長 (1 周期未満の区間から位相を逆算する)
    L_meas = 2 * math.pi * z_end / phi
    L_theory = lam / abs(neffs[0] - neffs[1])
    err = abs(L_meas - L_theory) / L_theory

    print(f"neff        : {', '.join(f'{v:.6f}' for v in neffs)}")
    print(f"lambda      : {lam:.6e} m")
    print(f"centroid    : {c0 * 1e6:+.3f} um (z=0) -> {c1 * 1e6:+.3f} um (z={z_end * 1e6:.0f} um)")
    print(f"beat length : 実測 {L_meas * 1e6:.0f} um / 理論 {L_theory * 1e6:.0f} um  "
          f"(誤差 {err:.2%}, tol {args.rtol:.0%})")

    failures = 0
    if err > args.rtol:
        print("FAIL: ビート長が理論値と一致しません")
        failures += 1
    # 電力保存 (重ね合わせは L2 正規化 -> 伝搬による損失のみ)
    ptot = Ixz.sum(axis=1)
    drift = abs(ptot[-1] - ptot[0]) / ptot[0]
    print(f"power drift : {drift:.3%} (伝搬マップ断面の総和)")
    if drift > 0.05:
        print("FAIL: 電力が保存していません")
        failures += 1

    print("=== %s (%d failure(s)) ===" % ("PASSED" if failures == 0 else "FAILED", failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
