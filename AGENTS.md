# OpenBPM — エージェント向けガイド (AGENTS.md)

OpenFDTD のコード構造 (入力・メッシュ・材料・後処理) を元にしたビーム伝搬法 (BPM)
ソルバー。伝搬カーネルは BPM-MATLAB の FDBPM (Douglas-Gunn ADI) を移植したもの。
OpenFDTD-X (GUI) から QProcess で起動される処理カーネルでもある。
ドキュメント・コメント・コミットメッセージ・レビューコメントは**日本語**で書く。

> Claude Code 用の `CLAUDE.md` と `.claude/rules/*.md` に同じ規約がある。
> 本ファイルはそれらを 1 つに集約した自己完結版。規約を変更する場合は両方を更新すること。

## 環境セットアップ

```sh
apt-get update && apt-get install -y cmake g++ libeigen3-dev libhdf5-dev
```

## ビルド / テスト / 実行

```sh
# ビルド (実行ファイルは bin/ に出力: obpm, obpm_post)
cmake -S . -B build -DWITH_TESTS=ON
cmake --build build -j

# 単体テスト (解析解との比較検証) — 変更後は必ず実行する
ctest --test-dir build --output-on-failure

# サンプル実行 (結果は time_series_data.h5)
./bin/obpm data/fiber.ofd

# ONN 活性化曲線の解析解検証 (CI と同じ判定)
sh tools/check_activation.sh activation_curve.csv obpm.log
```

- 回帰の基準値: `data/sample/fiber.ofd` は output power = 3.122518e+02 が不変で
  あること (機能追加で既存入力の結果を変えない)
- CUDA 版は `-DWITH_CUDA=ON` (この環境に GPU がない場合はコンパイル検証まで。
  その旨をコミット・PR に明記する)

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

## 数値計算の規約 (レビュー時の重点確認項目)

- **検証第一**: 数値カーネルの変更は必ず解析解と比較して検証する。
  回折 `w(z)=w0√(1+(z/zR)²)`、吸収 `exp(-2k0n''z)`、曲げ偏向 `<x>=z²/(2RoC)`、
  LP モード分散方程式、TPA `T=1/(1+β(P/A_eff)L)`。
  **既存テストの許容誤差を緩めて通すのは禁止。**
- **複素屈折率の符号**: 材料テーブル `n_mat` は損失を **+imag** で保持。物理符号は
  `n = nr - i*ni` のため伝搬計算前に虚部を反転する (`-CIMAGF(nm)`)。混同すると増幅になる。
- **TPA**: 強度減衰率 α = βI に対し界には α/2 を適用 (`E *= exp(-(β/2)·I·dz)`)。
  1 cm/GW = 1e-11 m/W。物理スケーリング (∫∫|E|²dA = P_in) は tpa/powersweep 指定時のみ。
- **曲げ (等価屈折率法)**: `n_bend = n·(1 - n²·xb·ρe/(2RoC))·exp(xb/RoC)`。
  実部のみ変換し吸収 (虚部) は変えない。
- **座標**: 格子中心原点 `x = dx*(ix-(Nx-1)/2)`、配列は行優先 `[iy*Nx+ix]`。
- **後方互換**: 入力キーワード省略時は従来動作と一致 (fiber 回帰を壊さない)。

## 機能追加の規則

- **CPU 両経路 + CUDA パリティ**: 物理の追加は CPU の近軸 (`FDBPMpropagator.c`) と
  拡張 (`wabpm.cpp`) の**両経路**に入れる。CUDA 未対応の場合は実行時警告 +
  ReadMe への明記が必須 (サイレント無視は禁止)。
- **入力キーワード追加**: `include/obpm.h` の BPM 構造体 → `sol/input_data.c` の
  既定値と解析 → ReadMe.md のキーワード表、の 3 点セット。
  新機能には `data/sample/` の検証ケースと CI スモークを付ける。
- **HDF5 出力**: `/field/*` (Ny×Nx 行優先) と `/metadata/*` (スカラー)。
  新規データセット追加時は `tools/plot_ixz.py` と `post/postbpm.c` の対応も検討。

## 移植性の絶対規則 (Windows/macOS CI で実際に踏んだもの)

- **C99 VLA 禁止** (MSVC 非対応)。`malloc` + 明示インデックスで書く。
- 複素数は `CREALF` / `CIMAGF` マクロ経由でアクセスする (MSVC では C ソースが
  C++ / `std::complex<float>` としてビルドされる)。`complex × double` の直接乗算を
  書かない (`complex × float` は可)。
- 暗黙の関数宣言・ナローイング変換を残さない (macOS AppleClang はエラー扱い)。
- libm は `MATH_LIB` 変数経由。3 OS の CI (Linux/macOS/Windows) が通ることを確認する。

## テストの規約

- 外部テストフレームワークは使わない。`tests/` の自己完結ハーネス
  (`check_close`/`check_true` + 終了コード) の既存パターンを踏襲する。
- テストは解析解または厳密な数学的性質 (エネルギー保存・直交性) と比較する。
- 許容誤差は物理的根拠 (離散化誤差の見積り) をコメントに書いて設定する。
- CI の判定基準を弱める変更 (許容誤差の拡大・判定の削除) は物理的根拠を PR に
  明記しない限り行わない。

## 作業の進め方

- 実装漏れ対応は `docs/implementation-checklist.md` を起点にし、対応後は状態を更新する。
- コミットメッセージは日本語で「何を・なぜ」を要約し、検証結果 (テスト通過・
  解析解との誤差) を本文に含める。
