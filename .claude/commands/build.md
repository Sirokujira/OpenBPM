---
description: ビルドと単体テストを実行して結果を要約する
allowed-tools: Bash(cmake:*), Bash(ctest:*), Bash(git status:*)
---

OpenBPM をビルドして単体テストを実行し、結果を要約してください。

手順:

1. `cmake -S . -B build -DWITH_TESTS=ON` で configure する
   (既に build/ があれば再利用してよい)
2. `cmake --build build -j` でビルドする
3. `ctest --test-dir build --output-on-failure` で単体テストを実行する
4. CUDA 関連ファイル (`cuda/`, `bpm/*.cu`) に未コミットの変更がある場合は
   `cmake -S . -B build-cuda -DWITH_CUDA=ON -DWITH_TESTS=ON` +
   `cmake --build build-cuda -j` でコンパイル検証も行う

報告条件:
- エラー・警告 (新規のもの)・テスト失敗があれば、原因箇所を特定して報告する
- すべて成功なら「ビルド成功・テスト N 件通過」と簡潔に報告する
- 失敗時に許容誤差を緩めて通そうとしない (.claude/rules/testing.md 参照)

$ARGUMENTS
