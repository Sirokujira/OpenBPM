# OpenBPM 開発ガイド (Claude Code 用)

OpenFDTD のコード構造を土台にしたビーム伝搬法 (BPM) 光導波路ソルバー (C/C++)。
伝搬カーネルは BPM-MATLAB の FDBPM (Douglas-Gunn ADI) 移植 + 独自拡張
(広角 Pade(1,1) / 半ベクトル / モードソルバ / TPA 非線形 / 波長掃引)。
OpenFDTD-X (GUI) から QProcess で起動される処理カーネルでもある。

## ビルド・テスト (まずこれが通ることを確認)

```sh
# 依存: cmake, gcc/g++, libhdf5-dev, libeigen3-dev (+ OpenMP)
cmake -S . -B build -DWITH_CUDA=OFF -DWITH_MPI=OFF -DWITH_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure   # 解析解との比較検証 (9本、要 h5py/numpy)
```

```sh
# 回帰 (fiber): output power = 3.122518e+02 が不変、
# obpm_post 後の bpm_ixz.csv が変更前と md5 一致であること
mkdir -p /tmp/smoke && cp data/sample/fiber.ofd /tmp/smoke/ && cd /tmp/smoke
$OLDPWD/bin/obpm -n 2 fiber.ofd && $OLDPWD/bin/obpm_post -n 2 fiber.ofd
grep "normal end" obpm.log

# ONN 活性化 (TPA + powersweep): activation_curve.csv が単調非増加で
# 解析解 T = 1/(1 + β(P/A_eff)L) と ±7% (A_eff・L はログ出力を使う)
cp data/sample/onn_activation.ofd /tmp/smoke/ && $OLDPWD/bin/obpm -n 2 onn_activation.ofd
```

- CUDA 版: `-DWITH_CUDA=ON` → `obpm_cuda` (nvcc 必要。GPU 実機がない環境では
  コンパイル/リンク確認まで)
- MPI 版: `-DWITH_MPI=ON` (要 libopenmpi-dev + **libhdf5-openmpi-dev**)。
  `obpm_mpi` は **FDTD ソルバーであり BPM ではない**

## アーキテクチャ地図

| パス | 役割 |
|---|---|
| `src/sol_Main.cpp` | `obpm` エントリ。入力読込 → `solve_bpm()` → `outputChars()`/`writeout()` |
| `sol/solve_bpm.cpp` | BPM 本体 (CPU)。屈折率構築→励振→伝搬 (近軸/広角)→HDF5 出力 |
| `cuda/solve_bpm.cu` | 同 CUDA 版。**sol 側と鏡写しに保つこと** (下記ルール参照) |
| `bpm/FDBPMpropagator.c/.cu` | スカラー近軸カーネル (BPM-MATLAB 移植、改変は上流照合必須) |
| `bpm/wabpm.cpp/.cu` | 広角/半ベクトルの一般化 ADI (倍精度) |
| `bpm/modes.cpp` | モードソルバ (虚軸伝搬法)。`launch = mode` から使用 |
| `sol/input_data.c` | .ofd パーサ。BPM 拡張キーワードもここ |
| `post/postbpm.c` | obpm_post の BPM (/field) 可視化 |
| `tests/` | ctest (wabpm/modes/allset)。`tools/plot_ixz.py` は伝搬マップ描画 |

出力: `time_series_data.h5` (`/field/Efinal_*`, `/field/Ixz`, `/field/frames`,
`/metadata/*`) + `obpm.out` (FDTD 形式。obpm_post の入力。`-no-fdtd-out` で省略可)
+ `activation_curve.csv` (`powersweep` 指定時)

## 重要ルール (詳細は @.claude/rules/ を参照)

- @.claude/rules/physics-validation.md
- @.claude/rules/keyword-wiring.md

## 移植性の絶対規則 (Windows CI で実際に踏んだもの)

- **C99 VLA 禁止** (MSVC)。`malloc` + 明示インデックスで書く。
- 複素数は `CREALF` / `CIMAGF` マクロ経由でアクセスする
  (MSVC では C ソースが C++ / `std::complex<float>` としてビルドされる)。
  `complex × double` の直接乗算を書かない (`complex × float` は可)。
- libm は `MATH_LIB` 変数経由。MSVC フラグは既存 CMake ブロックに従う。
- C は C99、C++ は C++17。ソースは UTF-8 (MSVC は /utf-8)。

## 機能追加の規則

- 入力キーは `sol/input_data.c` に追加し、**省略時は従来動作とバイト一致**。
- 物理スケーリング規約: `tpa`/`powersweep` 使用時のみ場を
  ∫∫|E|²dA = P_in [W] に正規化 (|E|² = 強度 I)。未使用時は従来の
  無次元場のまま (fiber 回帰を壊さない)。
- CPU の近軸/広角の**両経路**に同じ物理を入れる (片方だけの実装は不可)。
- CUDA 版 (`obpm_cuda`) の対応状況を ReadMe.md に明記し、未対応キーワードは
  **実行時 warning を出す** (サイレント無視は禁止)。
- 新機能には data/ の検証ケース + CI スモークを付ける。

## Gotchas (このリポジトリ固有のハマりどころ)

- **main が並行して進む**: 別セッションの PR が main へ随時マージされる。
  作業開始時と push 前に必ず `git fetch origin main` して分岐を確認する
  (`/sync-check` コマンド、SessionStart フックが自動表示)。
- **obpm_post は obpm.out を必須入力にする**: BPM の既定出力から obpm.out を
  消してはいけない (97MB でも既定 ON。省略はユーザーの `-no-fdtd-out` 判断)。
- **`I` マクロ**: `bpm_prototype.h` が複素虚数単位 `I` を定義している。
  変数名に `I` を使わない。虚数単位をキャストで自作しない (過去に実数 1.0 に
  なる事故があった)。
- **Dt/Tw は setup() が自動計算する**: 「ユーザーが指定したか」の判定に
  `Dt != 0` 等は使えない。
- **CI は obpm.log の "normal end" と HDF5 の存在を検証する**: 終了メッセージや
  ファイル名を変えるときは CI (linux/macOS/Windows の 3 ジョブ) も併せて更新する。
- **メッシュ配列**: `Xn/Yn/Zn` は節点 (要素数 N+1)、`Xc/Yc/Zc` はセル中心。
  セル幅は `(Xn[Nx]-Xn[0])/Nx` (演算子優先順位のバグが過去に複数あった)。

## コード規約

- 既存ファイルのインデント流儀に合わせる (sol/ · cuda/ はタブ、bpm/ は 4 空白)。
- コメントは日本語で可 (既存に合わせる)。
- コミットメッセージは日本語、末尾に検証内容 (何をどう確認したか) を書く。

## CI

`.github/workflows/ci.yml`: Linux / macOS / Windows (MSVC + Ninja +
vcpkg `hdf5[core,zlib]:x64-windows-static-md`) + CUDA ビルド検証。
Linux では ctest と data/ 全サンプルのスモークも実行。タグ `v*` で Release 添付。
