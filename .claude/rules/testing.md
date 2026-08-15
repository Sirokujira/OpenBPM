---
paths:
  - "tests/**"
  - ".github/workflows/**"
  - "tools/check_activation.sh"
---

# テスト・CI の規約

## 単体テスト (tests/)

- 外部テストフレームワークは使わない。自己完結ハーネス
  (`check_close` / `check_true` + 終了コード) の既存パターンを踏襲する。
- テストは必ず**解析解**または**厳密な数学的性質** (エネルギー保存・直交性・
  スカラー帰着) と比較する。実装同士の比較だけのテストは回帰検知にしかならない。
- 依存は最小に保つ: `test_wabpm` / `test_modes` は OpenMP のみ、
  `test_allset` のみ Eigen 可。HDF5 に依存するテストは作らない (e2e は CI スモークで行う)。
- 許容誤差は物理的根拠を持って設定する (離散化誤差の見積りをコメントに書く)。
  マージンを不必要に広げない。
- 追加したテストは `CMakeLists.txt` の `WITH_TESTS` ブロックに登録し、
  `ctest` で通ることを確認してからコミットする。

## CI (.github/workflows/ci.yml)

- ジョブ構成: Linux (gcc) / macOS (AppleClang + libomp) / Windows (MSVC + vcpkg)。
- Linux/macOS は単体テスト (ctest) + サンプルスモーク + ONN 解析解判定
  (`tools/check_activation.sh`, 許容 ±8%) を実行する。
- スモークの判定基準を弱める変更 (許容誤差の拡大・判定の削除) は、
  物理的な理由と根拠を PR に明記しない限り行わない。
- YAML 変更後は `python3 -c "import yaml; yaml.safe_load(open('.github/workflows/ci.yml'))"`
  で構文検証してからコミットする。

## 再現性

BPM は OpenMP で並列実行する。**同一スレッド数なら結果はビット単位で再現する**
(集約はスレッド番号順、集約に寄与するループは `schedule(static)`)。
スレッド数が変わると電力集約の部分和の個数が変わるため最終桁が動く。
したがって:

- md5 での回帰比較は**同じスレッド数**同士で行う
- CI の判定は `output power` の 7 桁など、スレッド数に依らない量で行う
- `wlsweep_check` のような「実行間の完全一致」を要求するテストは、
  並列関連の変更後に必ず実行する (再現性の退行を最初に検知する)
