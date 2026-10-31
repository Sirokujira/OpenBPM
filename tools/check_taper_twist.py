#!/usr/bin/env python3
"""
check_taper_twist.py

テーパ / ツイスト (座標変換) の幾何判定。出口の屈折率分布 /field/n_out_r から
コア (n > n_clad) の 2 次モーメントを取り、主軸の角度と広がりを期待値と比べる。

  主軸の角度   : ツイストの合計回転角 (反時計回りが正, 0..180 deg で比較)
  長辺の 2σ    : 一様な矩形 (辺 L) の 2σ = L/sqrt(3)。テーパ後の長辺 L で比べる

判定幅の根拠 (data/sample/taper_twist.ofd, dx = 1/3 um):
  - 角度 : 幅 1um (3 セル) の回転した矩形のラスタ化で ±1 deg 程度の測定バイアスが
           出る (ステップ数を変えても 90.6〜91.3 deg を揺れる) ため ±2 deg
  - 2σ   : 長辺 4um の理論 2.309um に対し実測 2.33〜2.35um のため ±3%
           (テーパが効かなければ 4.62um、ツイストが逆向きなら角度が 180 - θ になる)

使い方:
  python3 tools/check_taper_twist.py time_series_data.h5 --angle 90 --long 4e-6
  (n_clad, 領域は HDF5 の /metadata から取る)
"""

import argparse
import math
import sys

import h5py
import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("h5", help="time_series_data.h5")
    ap.add_argument("--angle", type=float, required=True, help="期待する主軸角 [deg]")
    ap.add_argument("--long", type=float, required=True, help="テーパ後の長辺 [m]")
    ap.add_argument("--nclad", type=float, default=1.4440, help="クラッド屈折率")
    ap.add_argument("--angle-tol", type=float, default=2.0, help="角度の許容 [deg]")
    ap.add_argument("--long-rtol", type=float, default=0.03, help="長辺 2σ の相対許容")
    a = ap.parse_args()

    f = h5py.File(a.h5, "r")
    n = f["/field/n_out_r"][:]
    xc = f["/metadata/Xc"][:] if "/metadata/Xc" in f else None
    yc = f["/metadata/Yc"][:] if "/metadata/Yc" in f else None
    ny, nx = n.shape
    if xc is None or len(xc) != nx or yc is None or len(yc) != ny:
        print("FAIL: /metadata/Xc, Yc が無いか寸法が合わない")
        return 1
    X, Y = np.meshgrid(xc, yc)
    w = np.clip(n - a.nclad, 0, None)          # コアだけに重み
    s = w.sum()
    if s <= 0:
        print("FAIL: n_out にコアが見つからない")
        return 1
    mx, my = (w * X).sum() / s, (w * Y).sum() / s
    cxx = (w * (X - mx) ** 2).sum() / s
    cyy = (w * (Y - my) ** 2).sum() / s
    cxy = (w * (X - mx) * (Y - my)).sum() / s
    ang = (0.5 * math.degrees(math.atan2(2 * cxy, cxx - cyy))) % 180.0
    ev = np.sort(np.linalg.eigvalsh([[cxx, cxy], [cxy, cyy]]))
    two_sigma_long = 2 * math.sqrt(ev[1])

    exp_ang = a.angle % 180.0
    dang = abs(ang - exp_ang)
    dang = min(dang, 180.0 - dang)
    exp_long = a.long / math.sqrt(3.0)
    rlong = two_sigma_long / exp_long - 1

    ok = (dang <= a.angle_tol) and (abs(rlong) <= a.long_rtol)
    print(f"{'OK  ' if ok else 'FAIL'} 主軸 {ang:.2f} deg (期待 {exp_ang:.2f}, 差 {dang:.2f}, 許容 {a.angle_tol}) / "
          f"長辺 2σ {two_sigma_long * 1e6:.3f} um (期待 {exp_long * 1e6:.3f}, {rlong * 100:+.1f}%, 許容 {a.long_rtol * 100:.0f}%)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
