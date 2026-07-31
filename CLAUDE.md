# OpenBPM プロジェクトメモリ

OpenFDTD のコード構造 (入力・メッシュ・材料・後処理) を元にしたビーム伝搬法 (BPM)
ソルバー。伝搬カーネルは BPM-MATLAB の FDBPM (Douglas-Gunn ADI) を移植したもの。
OpenFDTD-X (GUI) から QProcess で起動される処理カーネルでもある。
ドキュメント・コメント・コミットメッセージは**日本語**で書く。

> 他エージェント (Codex 等) 向けに同内容を集約した `AGENTS.md` がある。
> 本ファイルまたは `.claude/rules/*.md` の規約を変更したら `AGENTS.md` も更新すること。

## ビルド / テスト / 実行

```sh
# ビルド (実行ファイルは bin/ に出力: obpm, obpm_post)
cmake -S . -B build -DWITH_TESTS=ON
cmake --build build -j

# 単体テスト (解析解との比較検証)
ctest --test-dir build --output-on-failure

# サンプル実行 (結果は time_series_data.h5)
./bin/obpm data/fiber.ofd

# 可視化 (PNG / GIF)
python3 tools/plot_ixz.py time_series_data.h5

# ONN 活性化曲線の解析解検証 (CI と同じ判定)
sh tools/check_activation.sh activation_curve.csv obpm.log
```

- 依存: CMake >= 3.18, C/C++ コンパイラ, Eigen3, HDF5, OpenMP
- オプション: `-DWITH_CUDA=ON` (CUDA 版 obpm_cuda), `-DWITH_MPI=ON` (FDTD のみ)
- 回帰の基準値: `data/sample/fiber.ofd` は output power = 3.122518e+02 が不変で
  あること (機能追加で既存入力の結果を変えない)

## ディレクトリ構成

| パス | 内容 |
|---|---|
| `sol/` | CPU ソルバー本体。BPM の入口は `sol/solve_bpm.cpp` |
| `bpm/` | BPM カーネル。近軸 `FDBPMpropagator.c`、広角/半ベクトル `wabpm.cpp`、モードソルバ `modes.cpp` |
| `cuda/` | CUDA 版。BPM は `cuda/solve_bpm.cu` + `bpm/*.cu` |
| `mpi/`, `cuda_mpi/` | MPI 版 (FDTD のみ、BPM 未対応) |
| `include/` | ヘッダ。グローバル状態は `obpm.h` (EXTERN パターン)、BPM API は `include/bpm/` |
| `post/` | 後処理 (obpm_post)。BPM 出力は `postbpm.c` |
| `tests/` | 単体テスト (フレームワーク非依存の自己完結ハーネス) |
| `tools/` | Python 可視化・CI 検証スクリプト |
| `data/` | OpenFDTD 形式 (.ofd) のサンプル入力。理論値との比較ポイントをコメントに記載 |
| `docs/implementation-checklist.md` | 実装漏れ監査と対応状況の台帳 |

## 重要な規約 (詳細は .claude/rules/ 参照)

- **検証第一**: 数値カーネルの変更は必ず解析解 (回折・減衰・分散方程式など) と
  比較して検証する。既存テストの許容誤差を緩めて通すのは禁止。
- **後方互換**: 入力キーワードの省略時は従来動作と一致させる (fiber 回帰を壊さない)。
- **複素屈折率の符号**: `n_mat` は損失を +imag で保持。物理符号は `n = nr - i*ni`
  (伝搬計算前に符号を反転する)。混同すると増幅になる。
- **CPU 両経路 + CUDA パリティ**: 物理の追加は CPU の近軸/拡張の**両経路**に入れる
  (片方だけは不可)。CUDA 未対応の場合は実行時警告 + ReadMe への明記が必須
  (サイレント無視は禁止)。
- **移植性**: C99 VLA 禁止・複素数は `CREALF`/`CIMAGF` マクロ経由など、
  MSVC/macOS で実際に踏んだ規則を `.claude/rules/portability.md` に集約。
- **入力キーワード追加**: `include/obpm.h` の BPM 構造体 → `sol/input_data.c` の
  既定値と解析 → ReadMe.md のキーワード表、の 3 点セット。新機能には
  `data/sample/` の検証ケースと CI スモーク (3 OS) を付ける。
- **HDF5 出力**: `/field/*` (Ny×Nx 行優先) と `/metadata/*` (スカラー)。
  新規データセット追加時は `tools/plot_ixz.py` と `post/postbpm.c` の対応も検討。

## 作業の進め方

- 実装漏れ対応は `docs/implementation-checklist.md` を起点にし、対応後は
  同ファイルの状態 (✅/現状/検証内容) を更新する。
- GPU 実機はこの開発環境にないため、CUDA 変更は「CUDA 12.0 でのコンパイル検証」まで。
  その旨をコミット/チェックリストに明記する。
- コミットメッセージは日本語で「何を・なぜ」を要約し、検証結果 (テスト通過・
  解析解との誤差) を本文に含める。
