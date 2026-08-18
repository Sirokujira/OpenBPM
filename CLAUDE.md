# OpenBPM 開発ガイド (Claude Code 用)

OpenFDTD のコード構造 (入力・メッシュ・材料・後処理) を土台にしたビーム伝搬法 (BPM)
光導波路ソルバー (C/C++)。伝搬カーネルは BPM-MATLAB の FDBPM (Douglas-Gunn ADI)
移植 + 独自拡張 (広角 Pade(1,1) / 半ベクトル / モードソルバ / TPA 非線形 /
波長掃引 / テーパ・ツイスト / 対称境界)。
OpenFDTD-X (GUI) から QProcess で起動される処理カーネルでもある。
ドキュメント・コメント・コミットメッセージは**日本語**で書く。

> 他エージェント (Codex 等) 向けに同内容を集約した `AGENTS.md` がある。
> 本ファイルまたは `.claude/rules/*.md` の規約を変更したら `AGENTS.md` も更新すること。

## ビルド・テスト (まずこれが通ることを確認)

```sh
# 依存: cmake >= 3.18, gcc/g++, libhdf5-dev, libeigen3-dev (+ OpenMP)
# macOS: brew install libomp hdf5 eigen (AppleClang は OpenMP を同梱しない)
# 実行ファイルは bin/ に出力 (obpm, obpm_post)
cmake -S . -B build -DWITH_CUDA=OFF -DWITH_MPI=OFF -DWITH_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure   # 解析解との比較検証 (10本、要 h5py/numpy)
```

```sh
# 回帰 (fiber): output power = 3.122506e+02 が不変、
# obpm_post 後の bpm_ixz.csv が変更前と md5 一致であること
# (md5 比較は「同じスレッド数」同士で行う。BPM は並列実行するため、
#  電力集約の部分和の個数がスレッド数で変わり最終桁が動く。
#  output power の 7 桁はスレッド数に依らず不変)
mkdir -p /tmp/smoke && cp data/sample/fiber.ofd /tmp/smoke/ && cd /tmp/smoke
$OLDPWD/bin/obpm -n 2 fiber.ofd && $OLDPWD/bin/obpm_post -n 2 fiber.ofd
grep "normal end" obpm.log

# ONN 活性化 (TPA + powersweep): activation_curve.csv が単調非増加で
# 解析解 T = 1/(1 + β(P/A_eff)L) と ±7% (A_eff・L はログ出力を使う)
cp data/sample/onn_activation.ofd /tmp/smoke/ && $OLDPWD/bin/obpm -n 2 onn_activation.ofd
sh $OLDPWD/tools/check_activation.sh activation_curve.csv obpm.log

# 可視化 (PNG / GIF): 伝搬マップ・/trace グラフ・モード形状・結合率
python3 tools/plot_ixz.py time_series_data.h5
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
| `bpm/modes.cpp` | モードソルバ (虚軸伝搬法)。`launch = mode` と `modes = <n>` から使用 |
| `sol/input_data.c` | .ofd パーサ。BPM 拡張キーワードもここ |
| `include/` | ヘッダ。グローバル状態は `obpm.h` (EXTERN パターン)、BPM API は `include/bpm/` |
| `post/postbpm.c` | obpm_post の BPM 可視化 (/field 等高線 + /trace 折れ線グラフ) |
| `mpi/`, `cuda_mpi/` | MPI 版 (FDTD のみ、BPM 未対応) |
| `tests/` | ctest 10 本 (単体: wabpm/modes/allset/fdbpm/pml、結合: ONN 活性化 / モードビート / 波長掃引) |
| `tools/` | Python 可視化 (`plot_ixz.py`)・CI 検証スクリプト (`check_activation.sh`) |
| `data/` | OpenFDTD 形式 (.ofd) のサンプル入力。理論値との比較ポイントをコメントに記載 |
| `docs/implementation-checklist.md` | 実装漏れ監査と対応状況の台帳 |

出力: `time_series_data.h5` (`/field/Efinal_*`, `/field/Ixz`, `/field/Iyz`,
`/field/frames`, `/field/frames_z`, `/trace/*`, `/modes/*`, `/metadata/*`)
+ `obpm.out` (FDTD 形式。obpm_post の入力。`-no-fdtd-out` で省略可)
+ `activation_curve.csv` (`powersweep` 指定時) / `spectrum.csv` (`wlsweep` 指定時)

## 重要ルール (詳細は @.claude/rules/ を参照)

- @.claude/rules/physics-validation.md
- @.claude/rules/keyword-wiring.md
- @.claude/rules/bpm-physics.md
- @.claude/rules/cuda.md
- @.claude/rules/testing.md
- @.claude/rules/portability.md

## 移植性の絶対規則 (Windows/macOS CI で実際に踏んだもの)

- **C99 VLA 禁止** (MSVC)。`malloc` + 明示インデックスで書く。
- 複素数は `CREALF` / `CIMAGF` マクロ経由でアクセスする
  (MSVC では C ソースが C++ / `std::complex<float>` としてビルドされる)。
  `complex × double` の直接乗算を書かない (`complex × float` は可)。
- libm は `MATH_LIB` 変数経由。MSVC フラグは既存 CMake ブロックに従う。
- C は C99、C++ は C++17。ソースは UTF-8 (MSVC は /utf-8)。
- 暗黙の関数宣言を残さない (macOS AppleClang はエラー扱い)。

## 重要な規約

- **検証第一**: 数値カーネルの変更は必ず解析解 (回折・減衰・分散方程式など) と
  比較して検証する。既存テストの許容誤差を緩めて通すのは禁止。
- **後方互換**: 入力キーは `sol/input_data.c` に追加し、**省略時は従来動作と
  バイト一致** (fiber 回帰を壊さない)。
- **複素屈折率の符号**: `n_mat` は損失を +imag で保持。物理符号は `n = nr - i*ni`
  (伝搬計算前に符号を反転する)。混同すると増幅になる。
- **物理スケーリング規約**: `tpa`/`powersweep` 使用時のみ場を
  ∫∫|E|²dA = P_in [W] に正規化 (|E|² = 強度 I)。未使用時は従来の
  無次元場のまま。
- **CPU 両経路 + CUDA パリティ**: 物理の追加は CPU の近軸/拡張の**両経路**に入れる
  (片方だけは不可)。CUDA 版 (`obpm_cuda`) の対応状況を ReadMe.md に明記し、
  未対応キーワードは**実行時 warning を出す** (サイレント無視は禁止)。
- **入力キーワード追加**: `include/obpm.h` の BPM 構造体 → `sol/input_data.c` の
  既定値と解析 → ReadMe.md のキーワード表、の 3 点セット。新機能には
  `data/sample/` の検証ケースと CI スモーク (3 OS) を付ける。
- **HDF5 出力**: `/field/*` (Ny×Nx 行優先) と `/metadata/*` (スカラー)。
  新規データセット追加時は `tools/plot_ixz.py` と `post/postbpm.c` の対応も検討。

## Gotchas (このリポジトリ固有のハマりどころ)

- **main が並行して進む**: 別セッションの PR が main へ随時マージされる。
  作業開始時と push 前に必ず `git fetch origin main` して分岐を確認する
  (`/sync-check` コマンド、SessionStart フックが自動表示)。
- **obpm_post は obpm.out を必須入力にする**: BPM の既定出力から obpm.out を
  消してはいけない (97MB でも既定 ON。省略はユーザーの `-no-fdtd-out` 判断)。
- **`I` マクロ**: `bpm_prototype.h` が複素虚数単位 `I` を定義している。
  変数名に `I` を使わない。虚数単位をキャストで自作しない (過去に実数 1.0 に
  なる事故があった)。
- **励振キーワードは 2 系統ある**: `launch = mode <m> [coef...]` (重ね合わせ可) と
  `modes = <n> [excite]` (解析 + 基本モード励振)。併用時は **launch を優先**し
  `excite` を無視する (警告を出す)。
- **BPM に時間軸は無い**: `time_series_data.h5` は OpenFDTD 由来の名前で、
  `Ntime = 1` 固定・`/metadata/time` も 1 点。系列の独立変数は **z** で、
  ルート属性 `marching_axis = "z"` / `marching_axis_values = "/trace/z"` が
  それを宣言している (上流 BPM-Matlab の `P.z` と同じ規約)。z を時間軸に
  見せかける実装はしない。属性を消すと GUI が横軸を決められなくなる
  (CI の「HDF5 axis declaration」スモークが検知する)。
- **Dt/Tw は setup() が自動計算する**: 「ユーザーが指定したか」の判定に
  `Dt != 0` 等は使えない。
- **CI は obpm.log の "normal end" と HDF5 の存在を検証する**: 終了メッセージや
  ファイル名を変えるときは CI (linux/macOS/Windows の 3 OS ジョブ) も併せて更新する。
- **CI 全ジョブが数秒で failure ならコードではない**: private リポジトリの
  Actions 無料枠 (月 2,000 分、macOS ×10 換算) の枯渇。症状は `runner_id: 0`・
  ステップ実行なし・ログ 404。`Ping` ワークフロー (echo のみ) の手動実行で
  1 分未満で切り分けられる。ワークフロー YAML の構文エラーなら conclusion が
  `startup_failure` になるので区別できる。
- **並列 HDF5 のヘッダは mpi.h を要求する**: `WITH_MPI=ON` では MPI を使わない
  `obpm`/`obpm_post` もビルドが壊れるため、CMakeLists が `HDF5::HDF5` に
  `MPI::MPI_C` を伝播させている。CI の build-mpi は全ターゲットをビルドして
  これを検知する (`--target obpm_mpi` だけでは素通りする)。
- **メッシュ配列**: `Xn/Yn/Zn` は節点 (要素数 N+1)、`Xc/Yc/Zc` はセル中心。
  セル幅は `(Xn[Nx]-Xn[0])/Nx` (演算子優先順位のバグが過去に複数あった)。
- **OpenMP は C と CXX の両方をリンクする**: `OpenMP::OpenMP_C` だけだと
  `-fopenmp` が `$<COMPILE_LANGUAGE:C>` で守られ、C++ ソース
  (`sol/solve_bpm.cpp` 等) に渡らない。さらに `bpm/FDBPMpropagator.c` の
  `#pragma omp for` は `solve_bpm.cpp` 側の parallel 領域に束縛されるため、
  C++ 側が OpenMP 無しだと orphaned となり **BPM 全体が直列実行**になる
  (バナーは `CPU+OpenMP`、`thread=N` と表示するので気づきにくい)。
- **集約に寄与するループは `schedule(static)`**: `schedule(dynamic)` だと
  各スレッドが担当する行が実行ごとに変わり、部分和の項の集合が変化するため
  同一入力の再実行で結果が一致しなくなる (`wlsweep_check` が落ちる)。
  集約自体もスレッド番号順に行って順序を固定している。

## コード規約

- 既存ファイルのインデント流儀に合わせる (sol/ · cuda/ はタブ、bpm/ は 4 空白)。
- コメントは日本語で可 (既存に合わせる)。
- コミットメッセージは日本語、末尾に検証内容 (何をどう確認したか) を書く。

## 作業の進め方

- 実装漏れ対応は `docs/implementation-checklist.md` を起点にし、対応後は
  同ファイルの状態 (✅/現状/検証内容) を更新する。
- GPU 実機はこの開発環境にないため、CUDA 変更は「CUDA 12.0 でのコンパイル検証」まで。
  その旨をコミット/チェックリストに明記する。
- コミットメッセージは日本語で「何を・なぜ」を要約し、検証結果 (テスト通過・
  解析解との誤差) を本文に含める。

## CI

`.github/workflows/ci.yml`: Linux / macOS / Windows (MSVC + Ninja +
vcpkg `hdf5[core,zlib]:x64-windows-static-md`) の 3 OS + CUDA / MPI ビルド検証の
計 5 ジョブ。3 OS とも ctest を実行、Linux は data/ 全サンプルのスモークも実行。
タグ `v*` で Release 添付。

課金対策 (private リポジトリの無料枠節約) のため:
- push トリガーは main とタグのみ。作業ブランチは PR イベントでカバー
  (二重実行の防止)。PR の無いブランチは workflow_dispatch で手動実行
- build-macos (課金 ×10) は PR では skip (main push / タグ / 手動のみ)
- `ping.yml` (Ping) はランナー割り当て診断用の最小ワークフロー (echo のみ)
