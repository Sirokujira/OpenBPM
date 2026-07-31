---
description: 解析解との定量比較でソルバーの物理を検証する (回帰確認の決定版)
---

.claude/rules/physics-validation.md の解析解セットを使って、現在のビルドの
物理検証を実行してください:

1. ctest (test_wabpm / test_modes / test_allset) を実行
2. 追加で以下のサンプル実行結果を解析解と突き合わせる (h5py 使用):
   - freespace.ofd : ガウシアン回折幅 w(z) vs w0*sqrt(1+(z/zR)^2)
   - lossy.ofd : 電力減衰 vs exp(-2*k0*n''*z)
   - fiber_mode.ofd : 電力保存 (~1.0) とプロファイル定常性
   - tilt_wideangle.ofd : 横変位 (近軸 z*sinθ / 広角 z*tanθ への接近)
3. 各項目を「実測 / 理論 / 誤差%」の表で報告
4. 誤差が既知の水準 (回折 ~1% / 減衰 ~0.1% / モード電力 ~1e-4) を
   明確に超えるものは、メッシュ半分で再実行して離散化誤差か実装退行かを
   切り分ける

$ARGUMENTS で特定の検証項目のみ指定可能。
