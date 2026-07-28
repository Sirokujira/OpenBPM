#!/usr/bin/env python3
"""
波長掃引 (wlsweep = 1) の検証。

掃引した各波長は「その波長を単独指定した実行」と厳密に等価でなければならない
(波長ごとに複素屈折率・参照屈折率・ADI 係数・励振界を組み直しているため、
どれか一つでも取りこぼすと結果がずれる)。本テストは:

  1. 掃引を実行し spectrum.csv の点数と波長 (= c/frequency) の整合を確認
  2. 掃引の **最終波長** を単独指定した入力を作って実行し、
     /field/Efinal が掃引結果と完全一致することを確認
     (/field は最終波長の結果、という仕様の裏取りも兼ねる)

使い方:
  tests/check_wlsweep.py --bin <obpm> --ofd <sweep.ofd> [--workdir DIR]
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

try:
    import h5py
    import numpy as np
except ImportError:
    sys.exit("h5py and numpy are required (pip install h5py numpy)")

C0 = 2.99792458e8


def run_solver(binary, ofd, workdir):
    r = subprocess.run([binary, "-n", "2", "-no-fdtd-out", ofd],
                       cwd=workdir, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:])
        raise SystemExit(f"solver failed (exit {r.returncode}): {ofd}")


def read_field(path):
    with h5py.File(path, "r") as f:
        E = f["/field/Efinal_r"][:] + 1j * f["/field/Efinal_i"][:]
        lam = float(f["/metadata/lambda"][()])
    return E, lam


def read_spectrum(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        head = f.readline()
        if "lambda_m" not in head:
            raise SystemExit(f"{path}: unexpected header: {head.strip()}")
        for line in f:
            if line.strip():
                rows.append([float(v) for v in line.split(",")])
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--ofd", required=True)
    ap.add_argument("--workdir", default=None)
    args = ap.parse_args()

    binary = os.path.abspath(args.bin)
    ofd = os.path.abspath(args.ofd)
    tmp = args.workdir or tempfile.mkdtemp(prefix="wlsweep_")
    os.makedirs(tmp, exist_ok=True)

    failures = 0

    # --- 1) 掃引を実行 ---
    run_solver(binary, ofd, tmp)
    spec_path = os.path.join(tmp, "spectrum.csv")
    if not os.path.exists(spec_path):
        raise SystemExit("spectrum.csv が生成されていません (wlsweep = 1 の入力ですか)")
    rows = read_spectrum(spec_path)
    print(f"spectrum.csv : {len(rows)} points")
    if len(rows) < 2:
        print("FAIL: 掃引点が 2 点未満です")
        failures += 1

    # lambda と frequency の整合 (lambda = c / f)
    for lam, freq, t in rows:
        # CSV は %.6e で丸めて書かれるため許容は 1e-6 (桁違いの取り違えを検出する目的)
        if abs(lam * freq - C0) / C0 > 1e-6:
            print(f"FAIL: lambda*frequency != c ({lam:.6e} * {freq:.6e})")
            failures += 1
        if not (0.0 <= t <= 1.0 + 1e-9):
            print(f"FAIL: transmission が範囲外です: {t}")
            failures += 1

    sweep_h5 = os.path.join(tmp, "sweep.h5")
    shutil.copy(os.path.join(tmp, "time_series_data.h5"), sweep_h5)
    E_sweep, lam_sweep = read_field(sweep_h5)

    # --- 2) 最終波長を単独指定した実行と比較 ---
    last_lam, last_freq = rows[-1][0], rows[-1][1]
    if abs(lam_sweep - last_lam) / last_lam > 1e-6:   # CSV の丸め分を許容
        print(f"FAIL: /metadata/lambda ({lam_sweep:.6e}) が spectrum 末尾 "
              f"({last_lam:.6e}) と一致しません")
        failures += 1

    src = open(ofd, encoding="utf-8").read().splitlines()
    out = []
    for line in src:
        if line.startswith("frequency2"):
            out.append(f"frequency2 = {last_freq:.10e} {last_freq:.10e} 0")
        elif line.startswith("wlsweep"):
            continue          # 単独実行 (掃引なし)
        else:
            out.append(line)
    single_ofd = os.path.join(tmp, "single.ofd")
    open(single_ofd, "w", encoding="utf-8").write("\n".join(out) + "\n")

    run_solver(binary, single_ofd, tmp)
    E_single, lam_single = read_field(os.path.join(tmp, "time_series_data.h5"))

    if E_sweep.shape != E_single.shape:
        print(f"FAIL: 形状が違います {E_sweep.shape} vs {E_single.shape}")
        failures += 1
    else:
        dmax = float(np.abs(E_sweep - E_single).max())
        scale = float(np.abs(E_single).max())
        rel = dmax / scale if scale > 0 else dmax
        print(f"最終波長 {last_lam * 1e6:.4f} um : 掃引 vs 単独実行 "
              f"最大差 {dmax:.3e} (相対 {rel:.2e})")
        if rel > 1e-9:
            print("FAIL: 掃引結果が単独実行と一致しません "
                  "(波長ごとの再設定に取りこぼしがあります)")
            failures += 1

    print("=== %s (%d failure(s)) ===" % ("PASSED" if failures == 0 else "FAILED", failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
