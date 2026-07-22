---
description: サンプル入力を実行して理論値・解析解と照合する (e2e 検証)
allowed-tools: Bash(cmake:*), Bash(ctest:*), Bash(./bin/obpm:*), Bash(bin/obpm:*), Bash(sh tools/check_activation.sh:*), Bash(python3:*)
---

OpenBPM の e2e 検証を実行してください。ビルド済みの `bin/obpm` を使い、
一時ディレクトリ (スクラッチパッド) でサンプルを実行して理論値と照合します。

対象 (引数で個別指定がなければ以下の代表セットを実行):

1. **ONN 活性化曲線**: `data/sample/onn_activation.ofd` を実行し、
   `sh tools/check_activation.sh activation_curve.csv obpm.log` で
   解析解 `T = 1/(1 + beta*(P_in/A_eff)*L)` との一致 (±8%) を判定する
2. **自由空間回折**: `data/freespace.ofd` を実行し、最終ビーム幅を
   `w(z) = w0*sqrt(1+(z/zR)^2)` と比較する (HDF5 を Python/h5py で読む)
3. **吸収媒質**: `data/lossy.ofd` を実行し、出力電力を `exp(-2*k0*n''*z)` と比較する
4. **曲げファイバ**: `data/fiber_bend.ofd` を実行し、モード重心の外側シフトと
   電力減衰 (README 記載の期待値: 重心 +5.8um, P/P0 = 0.988) を確認する

各 `.ofd` ファイル冒頭のコメントに「期待される動作」と理論値が記載されているので、
それを判定基準として使うこと。

報告: 各項目の 実測値 / 理論値 / 相対誤差 を表にまとめ、逸脱があれば原因を調査する。

$ARGUMENTS
