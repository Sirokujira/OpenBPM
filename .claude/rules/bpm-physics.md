---
paths:
  - "bpm/**"
  - "sol/solve_bpm.cpp"
  - "cuda/solve_bpm.cu"
  - "include/bpm/**"
---

# BPM 数値計算の規約

## 符号・単位の約束

- 複素屈折率: 物理符号は `n = nr - i*ni` (ni > 0 が損失)。
  材料テーブル `n_mat` は損失を **+imag** で保持しているため、
  伝搬計算へ渡す前に虚部の符号を反転する (`-CIMAGF(nm)`)。
- 比誘電率→屈折率: `n = sqrt(epsr)`、導電率は `ni = sigma/(2*omega*eps0*nr)` 相当で虚部へ。
- TPA: 強度減衰率 alpha = beta*I に対し **界には alpha/2** を適用
  (`E *= exp(-(beta/2)*I*dz)`)。単位変換 1 cm/GW = 1e-11 m/W。
- 物理スケーリング: `tpa`/`powersweep` 指定時のみ `∫∫|E|^2 dA = P_in [W]` に正規化。
  未指定時は無次元の界 (後方互換を壊さない)。
- 座標: 格子中心原点 `x = dx*(ix - (Nx-1)/2)`。配列は行優先 `[iy*Nx + ix]`。
- 曲げ (等価屈折率法): `n_bend = n*(1 - n^2*xb*rho_e/(2*RoC))*exp(xb/RoC)`、
  `xb = x*cos(dir) + y*sin(dir)`。**実部のみ**変換し吸収 (虚部) は変えない。

## 伝搬エンジンの使い分け

- 近軸スカラー: `bpm/FDBPMpropagator.c(.cu)` — 位相 multiplier + ADI。
- 広角 Pade(1,1) / 半ベクトル: `bpm/wabpm.cpp(.cu)` — 屈折率項を演算子に含む一般化 ADI
  (`(1 + d+P)E' = (1 + d-P)E`)。倍精度。
- モード解析: `bpm/modes.cpp` — 虚軸伝搬 + Gram-Schmidt。neff は Rayleigh 商。
- 半ベクトル (Stern) 差分は非対称。均一媒質で標準ラプラシアンに一致することが検算点。

## 変更時の必須検証

数値カーネルに触れたら、該当する解析解テストを必ず実行して確認する:

- 回折: `w(z) = w0*sqrt(1+(z/zR)^2)`
- 吸収: `P(z)/P(0) = exp(-2*k0*n''*z)`
- 曲げ偏向 (一様媒質): `<x>(z) = z^2/(2*RoC)`
- モード: LP 分散方程式 `u J_{l+1}(u)/J_l(u) = w K_{l+1}(w)/K_l(w)`
- TPA: `T = 1/(1 + beta*(P_in/A_eff)*L)`

許容誤差を緩めてテストを通すことは禁止。誤差が悪化した場合は原因を特定する。
