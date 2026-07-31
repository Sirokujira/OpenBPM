---
description: CPU版をビルドし、ctest (解析解検証) と fiber スモークまで一気に実行する
---

OpenBPM の標準ビルド+テストを実行してください:

1. `cmake -S . -B build -DWITH_CUDA=OFF -DWITH_TESTS=ON` で構成
   (依存欠如で失敗したら `libhdf5-dev` `libeigen3-dev` の導入を試みる)
2. `cmake --build build -j` でビルド (error があれば修正方針を報告)
3. `ctest --test-dir build --output-on-failure` を実行 (3/3 PASS が条件)
4. 一時ディレクトリで `bin/obpm -n 2 data/sample/fiber.ofd` を実行し、
   obpm.log の "normal end" と time_series_data.h5 の生成を確認
5. 結果を要約報告 (ビルド警告のうち新規のものがあれば列挙)

$ARGUMENTS が "cuda" を含む場合は `-DWITH_CUDA=ON` で `obpm_cuda` の
コンパイル/リンク確認も追加で行う (GPU が無い環境では実行はしない)。
