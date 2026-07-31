# OpenBPM エージェント向けガイド (Codex / その他 AI エージェント用)

このファイルは OpenAI Codex など AGENTS.md を読むエージェント向け。
Claude Code 向けの CLAUDE.md と内容を同期させている (人手で編集する場合は両方更新)。

OpenFDTD のコード構造を土台にしたビーム伝搬法 (BPM) 光導波路ソルバー (C/C++)。
伝搬カーネルは BPM-MATLAB の FDBPM (Douglas-Gunn ADI) 移植 + 独自拡張
(広角 Pade(1,1) / 半ベクトル / モードソルバ / TPA 非線形 / 波長掃引)。

## ビルド・テスト (まずこれが通ることを確認)

```sh
# 依存: cmake, gcc/g++, libhdf5-dev, libeigen3-dev (+ OpenMP, テストに h5py/numpy)
cmake -S . -B build -DWITH_CUDA=OFF -DWITH_MPI=OFF -DWITH_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure   # 解析解との比較検証 9 本
```

回帰 (必ず確認):
- `bin/obpm -n 2 data/sample/fiber.ofd` → obpm.log に `normal end`、
  `output power = 3.122518e+02` が**不変**であること
- `obpm_post` 後の `bpm_ixz.csv` が変更前と md5 一致であること

## 物理検証ファースト (最重要ルール)

これは数値物理ソルバーである。「ビルドが通る・落ちない」は合格条件ではない。
**数値が物理として正しいことを示して初めて完了**とする。詳細な検証マトリクスと
使える解析解 (ガウシアン回折 w(z)、吸収 exp(-2k0n''z)、スラブ分散方程式、
LP モード、GI 自己集束ピッチ、チルト変位) は `.claude/rules/physics-validation.md`
を参照。理論との差が出たらメッシュ/ステップを半分にして収束方向を確認し、
離散化誤差と実装バグを区別する。

`bpm/FDBPMpropagator.c/.cu` は BPM-MATLAB (github.com/ankrh/BPM-Matlab) の移植。
挙動を変える修正は上流と照合してから行う。

## 入力キーワード追加時の配線チェックリスト

.ofd キーワードの追加/変更は以下を**全部**揃えて 1 コミットにする
(詳細: `.claude/rules/keyword-wiring.md`):

1. `include/obpm.h` — BPM 構造体へフィールド追加 (単位と既定値をコメント)
2. `sol/input_data.c` — 既定値の初期化 + パース。**省略時は従来動作とバイト一致**
3. `sol/solve_bpm.cpp` — CPU 実装 (近軸パスと wabpm パスの**両方**)
4. `cuda/solve_bpm.cu` — CUDA 実装。未対応なら**実行時 warning 必須** (サイレント無視は禁止)
5. `data/` — 動作サンプル .ofd (期待される物理挙動をコメントで明記)
6. `ReadMe.md` — キーワード表へ行追加 (CUDA 対応状況も明記)
7. CI — 必要ならスモーク追加

## アーキテクチャ地図

| パス | 役割 |
|---|---|
| `src/sol_Main.cpp` | `obpm` エントリ。入力読込 → `solve_bpm()` → 出力 |
| `sol/solve_bpm.cpp` | BPM 本体 (CPU)。屈折率構築→励振→伝搬 (近軸/広角)→HDF5 出力 |
| `cuda/solve_bpm.cu` | 同 CUDA 版。**sol 側と鏡写しに保つこと** |
| `bpm/FDBPMpropagator.c/.cu` | スカラー近軸カーネル (BPM-MATLAB 移植) |
| `bpm/wabpm.cpp/.cu` | 広角/半ベクトルの一般化 ADI (倍精度) |
| `bpm/modes.cpp` | モードソルバ (虚軸伝搬法)。`launch = mode` から使用 |
| `sol/input_data.c` | .ofd パーサ。BPM 拡張キーワードもここ |
| `post/postbpm.c` | obpm_post の BPM (/field) 可視化 |
| `tests/` | ctest 9 本 (単体: wabpm/modes/allset/fdbpm、結合: ONN/モードビート/波長掃引) |

## 移植性の絶対規則 (Windows CI で実際に踏んだもの)

- **C99 VLA 禁止** (MSVC)。`malloc` + 明示インデックスで書く。
- 複素数は `CREALF` / `CIMAGF` マクロ経由 (MSVC では C ソースが
  `std::complex<float>` としてビルドされる)。`complex × double` の直接乗算禁止。
- libm は `MATH_LIB` 変数経由。C は C99、C++ は C++17。ソースは UTF-8。

## ハマりどころ

- **main が並行して進む**: 作業開始時と push 前に `git fetch origin main` で分岐確認。
- **obpm_post は obpm.out を必須入力にする**: BPM の既定出力から消さない。
- **`I` マクロ**: `bpm_prototype.h` が虚数単位 `I` を定義。変数名に使わない。
- **Dt/Tw は setup() が自動計算**: 「ユーザー指定か」の判定に `Dt != 0` は不可。
- **メッシュ配列**: `Xn/Yn/Zn` は節点 (N+1)、`Xc/Yc/Zc` はセル中心。
  セル幅は `(Xn[Nx]-Xn[0])/Nx` (括弧の位置に注意)。
- **並列 HDF5 は mpi.h を要求**: `WITH_MPI=ON` 時は CMake が `HDF5::HDF5` に
  `MPI::MPI_C` を伝播させる。壊すと obpm/obpm_post のビルドが落ちる。
- **CI 全ジョブが数秒で failure ならコードではない**: private リポジトリの
  Actions 無料枠枯渇 (`runner_id: 0`・ログ 404 が症状)。`Ping` ワークフローの
  手動実行で切り分け可能。

## コード規約

- 既存ファイルのインデント流儀に合わせる (sol/ · cuda/ はタブ、bpm/ は 4 空白)。
- コメントは日本語で可 (既存に合わせる)。
- コミットメッセージは日本語、末尾に検証内容 (何をどう確認したか) を書く。
