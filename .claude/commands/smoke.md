---
description: data/ の全サンプル .ofd を実行し、正常終了と HDF5 出力を検証する (CI と同じ手順)
---

data/ 直下の全サンプル .ofd に対して CI と同一のスモークを実行してください:

1. bin/obpm が無ければ先にビルド (/build-test 相当)
2. 一時ディレクトリで各サンプルを `bin/obpm -n 2 -no-fdtd-out <sample>` で実行
3. 各ケースで obpm.log の "normal end" と time_series_data.h5 の生成を確認
4. 出力電力 (ログの "output power") を一覧表示し、直近の既知値から
   大きく変わったサンプルがあれば指摘する
5. FAIL したサンプルはログ末尾を添えて報告

$ARGUMENTS にサンプル名 (拡張子なし) が指定された場合はそれのみ実行する。
