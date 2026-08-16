#!/usr/bin/env python3
"""data/sample/grin_rip.csv (放物型 GRIN の屈折率分布) を生成するスクリプト。

n(r) = n0 * sqrt(1 - 2*Delta*(r/a)^2)   (r < a)
     = n0 * sqrt(1 - 2*Delta)           (r >= a, クラッド)

この分布では近軸光線がちょうど調和振動し、自己集束ピッチは
  Lambda = 2*pi*a / sqrt(2*Delta)
となる (data/sample/grin.ofd の検証ポイント)。

メッシュは grin.ofd と一致させること (Ny 行 x Nx 列、行 0 が y 最小、
セル中心 Xc/Yc でサンプリング)。
"""
import numpy as np

n0, delta, a = 1.47, 0.02, 10e-6
Nx = Ny = 90
L = 30e-6                       # 領域幅 [-15um, +15um]
dx = L / Nx
xc = -L / 2 + (np.arange(Nx) + 0.5) * dx
yc = xc.copy()

r2 = xc[None, :] ** 2 + yc[:, None] ** 2
arg = np.maximum(1.0 - 2.0 * delta * r2 / (a * a), 1.0 - 2.0 * delta)
n = n0 * np.sqrt(arg)

np.savetxt("data/sample/grin_rip.csv", n, fmt="%.6f", delimiter=",")
print(f"wrote data/sample/grin_rip.csv ({Ny} x {Nx})")
print(f"pitch Lambda = 2*pi*a/sqrt(2*Delta) = {2*np.pi*a/np.sqrt(2*delta)*1e6:.2f} um")
