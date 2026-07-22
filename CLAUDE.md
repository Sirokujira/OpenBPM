# OpenBPM 開発ガイド (Claude Code 用)

OpenFDTD のコード構造を土台にしたビーム伝搬法 (BPM) ソルバー。
伝搬カーネルは BPM-MATLAB の FDBPM (Douglas-Gunn ADI) 移植 + 独自拡張
(広角 Pade(1,1) / 半ベクトル / モードソルバ / TPA 非線形)。

## ビルド・テスト (まずこれが通ることを確認)

```sh
# 依存: cmake, gcc/g++, libhdf5-dev, libeigen3-dev (+ OpenMP)
cmake -S . -B build -DWITH_CUDA=OFF -DWITH_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure   # 解析解との比較検証 (3本)
```

```sh
# スモーク (CI と同じ): サンプル実行 + 正常終了 + HDF5 出力確認
cd /tmp && bin/obpm -n 2 -no-fdtd-out data/sample/fiber.ofd
grep "normal end" obpm.log && test -s time_series_data.h5
```

- CUDA 版: `-DWITH_CUDA=ON` → `obpm_cuda` (nvcc 必要。GPU 実機がない環境では
  コンパイル/リンク確認まで)
- MPI 版: `-DWITH_MPI=ON` (要 libopenmpi-dev + **libhdf5-openmpi-dev**)。
  `obpm_mpi` は **FDTD ソルバーであり BPM ではない**

## アーキテクチャ地図

| パス | 役割 |
|---|---|
| `sol/solve_bpm.cpp` | BPM 本体 (CPU)。入力→屈折率→励振→伝搬→HDF5 出力 |
| `cuda/solve_bpm.cu` | 同 CUDA 版。**sol 側と鏡写しに保つこと** (下記ルール参照) |
| `bpm/FDBPMpropagator.c/.cu` | スカラー近軸カーネル (BPM-MATLAB 移植、改変は上流照合必須) |
| `bpm/wabpm.cpp/.cu` | 広角/半ベクトルの一般化 ADI (倍精度) |
| `bpm/modes.cpp` | モードソルバ (虚軸伝搬法)。`launch = mode` から使用 |
| `sol/input_data.c` | .ofd パーサ。BPM 拡張キーワードもここ |
| `post/postbpm.c` | obpm_post の BPM (/field) 可視化 |
| `tests/` | ctest (wabpm/modes/allset)。`tools/plot_ixz.py` は伝搬マップ描画 |

出力: `time_series_data.h5` (`/field/Efinal_*`, `/field/Ixz`, `/field/frames`,
`/metadata/*`) + `obpm.out` (FDTD 形式。obpm_post の入力。`-no-fdtd-out` で省略可)

## 重要ルール (詳細は @.claude/rules/ を参照)

- @.claude/rules/physics-validation.md
- @.claude/rules/keyword-wiring.md

## Gotchas (このリポジトリ固有のハマりどころ)

- **main が並行して進む**: 別セッションの PR が main へ随時マージされる。
  作業開始時と push 前に必ず `git fetch origin main` して分岐を確認する
  (`/sync-check` コマンド参照)。ブランチが古い main 基準のときは
  最新 main に rebase/マージしてから作業する。
- **obpm_post は obpm.out を必須入力にする**: BPM の既定出力から obpm.out を
  消してはいけない (97MB でも既定 ON。省略はユーザーの `-no-fdtd-out` 判断)。
- **`I` マクロ**: `bpm_prototype.h` が複素虚数単位 `I` を定義している。
  変数名に `I` を使わない。虚数単位をキャストで自作しない (過去に実数 1.0 に
  なる事故があった)。
- **Dt/Tw は setup() が自動計算する**: 「ユーザーが指定したか」の判定に
  `Dt != 0` 等は使えない。
- **CI は obpm.log の "normal end" と HDF5 の存在を検証する**: ソルバーの
  終了メッセージやファイル名を変えるときは CI (.github/workflows/ci.yml,
  linux/macos/windows の 3 ジョブ) も併せて更新する。
- **メッシュ配列**: `Xn/Yn/Zn` は節点 (要素数 N+1)、`Xc/Yc/Zc` はセル中心。
  セル幅は `(Xn[Nx]-Xn[0])/Nx` (演算子優先順位のバグが過去に複数あった)。

## コード規約

- 既存ファイルのインデント流儀に合わせる (sol/ · cuda/ はタブ、bpm/ は 4 空白)。
- コメントは日本語で可 (既存に合わせる)。ソースは UTF-8 (MSVC は /utf-8)。
- C は C99、C++ は C++17。MSVC も CI 対象 (VLA・GNU 拡張を新規に持ち込まない)。
- コミットメッセージは日本語、末尾に検証内容 (何をどう確認したか) を書く。
