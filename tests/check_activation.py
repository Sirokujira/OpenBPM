#!/usr/bin/env python3
"""
ONN 光活性化関数 (TPA + powersweep) の定量検証。

activation_curve.csv と obpm.log を読み、透過率 T(P_in) を平面波近似の
解析解と比較する:

    dI/dz = -beta * I^2  ->  I(z) = I0 / (1 + beta*I0*z)
    T(P_in) = 1 / (1 + beta * (P_in / A_eff) * L)

beta / A_eff / L はソルバーが obpm.log に出力した値を使う (テスト側で
入力を二重定義しない)。ガウシアン断面では径方向で飽和度が異なるため
平面波近似からは深飽和側でずれる (既知: 最大 +4% 程度) ので、既定許容は
相対 7%。あわせて単調非増加と飽和の有無も検査する。

標準ライブラリのみ使用 (CI ランナーに numpy 等を要求しない)。

使い方:
  tests/check_activation.py [--dir DIR] [--rtol 0.07]
"""
import argparse
import math
import os
import re
import sys

# Windows のコンソール/パイプは cp1252 等になり日本語出力で UnicodeEncodeError に
# なるため UTF-8 に強制する (CI の ctest 経由で実際に発生した)
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def parse_log(path):
    """obpm.log から beta [m/W]、A_eff [m^2]、L [m] を取り出す。"""
    beta = aeff = length = None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = re.search(r"A_eff\s*=\s*([0-9.eE+-]+).*?L\s*=\s*([0-9.eE+-]+)", line)
            if m:
                aeff, length = float(m.group(1)), float(m.group(2))
            m = re.search(r"beta\s*=\s*[0-9.eE+-]+\s*\[cm/GW\]\s*=\s*([0-9.eE+-]+)\s*\[m/W\]", line)
            if m:
                beta = float(m.group(1))
    return beta, aeff, length


def parse_curve(path):
    """activation_curve.csv を [(P_in, P_out, T), ...] として読む。"""
    rows = []
    with open(path, encoding="utf-8") as f:
        header = f.readline()
        if "P_in" not in header:
            raise SystemExit(f"{path}: unexpected header: {header.strip()}")
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) < 3:
                raise SystemExit(f"{path}: malformed row: {line}")
            rows.append((float(parts[0]), float(parts[1]), float(parts[2])))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".", help="obpm.log / activation_curve.csv のあるディレクトリ")
    ap.add_argument("--rtol", type=float, default=0.07, help="解析解との相対許容誤差")
    ap.add_argument("--min-points", type=int, default=4, help="必要な掃引点数")
    args = ap.parse_args()

    log_path = os.path.join(args.dir, "obpm.log")
    csv_path = os.path.join(args.dir, "activation_curve.csv")
    for p in (log_path, csv_path):
        if not os.path.exists(p):
            raise SystemExit(f"not found: {p}")

    beta, aeff, length = parse_log(log_path)
    if beta is None or aeff is None or length is None:
        raise SystemExit(f"{log_path}: beta / A_eff / L を読み取れませんでした")
    rows = parse_curve(csv_path)

    print(f"beta = {beta:.6g} [m/W], A_eff = {aeff:.6g} [m^2], L = {length:.6g} [m]")
    print(f"{'P_in[W]':>12} {'T(BPM)':>10} {'T(理論)':>10} {'誤差':>9}")

    failures = 0
    prev_t = None
    max_err = 0.0
    for p_in, p_out, t in rows:
        t_ana = 1.0 / (1.0 + beta * (p_in / aeff) * length)
        err = (t - t_ana) / t_ana
        max_err = max(max_err, abs(err))
        flag = " " if abs(err) <= args.rtol else "*"
        print(f"{p_in:12.4e} {t:10.5f} {t_ana:10.5f} {err:+8.2%}{flag}")
        if abs(err) > args.rtol:
            failures += 1
        # 単調非増加 (飽和吸収なので入力を上げると透過率は下がる)
        if prev_t is not None and t > prev_t * (1 + 1e-6):
            print(f"  FAIL: transmission is not monotonically non-increasing "
                  f"({prev_t:.6f} -> {t:.6f})")
            failures += 1
        # 出力パワーは入力を超えない (受動素子)
        if p_out > p_in * (1 + 1e-6):
            print(f"  FAIL: P_out > P_in ({p_out:.6g} > {p_in:.6g})")
            failures += 1
        prev_t = t

    if len(rows) < args.min_points:
        print(f"FAIL: too few sweep points: {len(rows)} < {args.min_points}")
        failures += 1
    # 掃引範囲で実際に飽和していること (テストとして意味のある範囲か)
    if rows and not (rows[0][2] > rows[-1][2] * 1.5):
        print(f"FAIL: no saturation over the sweep "
              f"(T: {rows[0][2]:.4f} -> {rows[-1][2]:.4f})")
        failures += 1

    print(f"max |error| = {max_err:.2%} (tol = {args.rtol:.0%})")
    print("=== %s (%d failure(s)) ===" % ("PASSED" if failures == 0 else "FAILED", failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
