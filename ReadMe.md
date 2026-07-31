# OpenBPM

OpenFDTD のコード構造（入力データ・メッシュ・材料・後処理）を元にした
ビーム伝搬法 (BPM) ソルバーです。伝搬カーネルは BPM-MATLAB の
FDBPM (Douglas-Gunn ADI 法) を移植しています。

## 処理部の構成

| ディレクトリ | 役割 |
|---|---|
| `src/sol_Main.cpp` | ソルバー `obpm` のエントリ。入力読込 → `solve_bpm()` → `outputChars()`/`writeout()` |
| `sol/solve_bpm.cpp` | BPM 処理本体: 複素屈折率分布の構築、ガウシアンビーム励振、Douglas-Gunn ADI 伝搬、端部吸収体、HDF5 出力 (`/field`, `/metadata`) |
| `bpm/FDBPMpropagator.c(u)` | 近軸スカラー BPM の伝搬カーネル (BPM-MATLAB 移植、CPU/CUDA) |
| `bpm/wabpm.cpp` / `bpm/wabpm.cu` | 広角 (Pade(1,1)) / 半ベクトル BPM の一般化 ADI エンジン (倍精度、CPU/CUDA) |
| `sol/` (その他) | OpenFDTD 由来の入力パーサ (`input_data.c`)・メッシュ・材料・出力群 |
| `post/` | ポストプロセッサ `obpm_post` (ev2/ev3、収束・波形プロット)。※下記の既知事項参照 |
| `include/bpm/` | BPM 共有ヘッダ (`bpm_prototype.h` のプラットフォーム別 complex マクロ等) |

> ℹ **既知事項**:
> - `obpm_post` は OpenFDTD 由来の FDTD 量の描画 (obpm.out ベース) に加え、
>   BPM が `time_series_data.h5` の `/field` に書く伝搬マップ (`Ixz`)・
>   最終電界 (`Efinal_*`)・伝搬スナップショット (`frames`、等間隔に最大 6 枚) の
>   可視化にも対応しています (`post/postbpm.c`)。
>   HDF5 を Python (h5py) で直接読む・GUI (OpenFDTD-X) の H5 ビューアで
>   参照することも可能です。
> - obpm.out (FDTD 形式、BPM ではほぼゼロの大容量バイナリ) が不要な場合は
>   ソルバーに `-no-fdtd-out` を付けるとスキップできます
>   (その場合 `obpm_post` の FDTD 量描画は不可、BPM の H5 出力は影響なし)。
> - BPM ソルバーで無効な入力キーワード (`planewave`/`point`/`load`/`rfeed`/
>   `abc = 1`/`pbc`) は実行時に warning を表示して無視します。
>   波長掃引は `wlsweep = 1` で有効 (省略時は `frequency2` の先頭周波数のみ使用)。
>   分散性材料 (type 2) は einf 近似です。

## ビルド

必要パッケージ: CMake (>=3.18), C/C++ コンパイラ, Eigen3, HDF5, OpenMP

```sh
cmake -S . -B build
cmake --build build -j
```

実行ファイルは `bin/` に出力されます (`obpm`, `obpm_post`)。

- CUDA 版: `-DWITH_CUDA=ON` (要 CUDA Toolkit)
  - `obpm_cuda` は CPU 版と同じ BPM 機能 (beam/refindex/導電率吸収/メッシュ警告/bend/
    beamtilt、広角 `wideangle`・半ベクトル `polarization` を含む) に対応しています。
    スカラー近軸は `bpm/FDBPMpropagator.cu`、広角/半ベクトルは `bpm/wabpm.cu` を使用します。
  - **`tpa` / `powersweep` (ONN 光活性化関数) と `wlsweep` (波長掃引) は CPU 版のみ**です。CUDA 版で
    指定した場合は実行時に warning を表示して無視します (CPU 版 `obpm` を使用してください)。
  - ビルドは CUDA 12.0 で確認済みです (GPU 実機での実行検証は未実施)。
- MPI 版: `-DWITH_MPI=ON` (要 MPI + 並列 HDF5。例: `libopenmpi-dev libhdf5-openmpi-dev`)
  - `obpm_mpi` は OpenFDTD 由来の **FDTD 時間領域ソルバー (BPM ではない)**。
    `mpirun -np <N> ./bin/obpm_mpi <input>.ofd` で実行します。

## 可視化

出力 HDF5 を PNG / GIF に変換する後処理スクリプトを用意しています
(要 h5py, numpy, matplotlib):

```sh
python3 tools/plot_ixz.py time_series_data.h5 [--db] [--prefix out] [--fps 15]
```

- `*_ixz.png` : 伝搬マップ `/field/Ixz` のヒートマップ
- `*_final.png` : 最終電界 `|E(x,y)|^2` と屈折率分布
- `*_prop.gif` : `|E(x,y)|^2` の伝搬アニメーション (入力に `frames = <interval>` を
  指定して `/field/frames` を記録した場合のみ)

## モード解析

虚軸伝搬法 (imaginary-distance BPM) + Gram-Schmidt 直交化によるモードソルバを
`bpm/modes.cpp` (`wabpm_find_modes`) に実装しています。実効屈折率の降順に導波モードを
求め、neff は Rayleigh 商から算出します。`include/bpm/allset.hpp` の `findModes()`
(BPM-MATLAB 互換 API) はこのソルバを使用します。

## テスト

解析解と比較する単体テストを `tests/` に用意しています。

```sh
pip install h5py numpy   # 結合テスト (mode_beat / wlsweep) が HDF5 を読むため
cmake -S . -B build -DWITH_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

| テスト | 検証項目 |
|---|---|
| `test_wabpm` | 自由空間ガウシアン回折 `w(z)=w0*sqrt(1+(z/zR)^2)` (近軸/広角)、エネルギー保存、吸収減衰 `exp(-2*k0*n''*z)`、曲げ偏向 `<x>=z^2/(2*RoC)`、半ベクトル差分の一様媒質整合 |
| `test_modes` | ステップインデックスファイバ LP01/LP11 の実効屈折率 (分散方程式の厳密解と比較)、モード直交性 |
| `test_allset` | `findModes` / `modeSuperposition` / `offsetField` / `tiltField` (P_Struct API) |
| `test_fdbpm` | スカラー近軸カーネル (既定経路, `bpm/FDBPMpropagator.c`) の自由空間回折 `w(z)`、エネルギー保存、一様吸収 `exp(-2*k0*n''*z)`、ビーム中心保持 |
| `mode_beat_solver` / `_check` | 多モード重ね合わせ励振の結合テスト。`data/fiber_mode_superposition.ofd` を実行し、伝搬マップの重心軌跡から求めたモードビート長を理論値 `lambda/\|neff0-neff1\|` と比較 (実測 5710um / 理論 5762um, 誤差 0.9%)。電力保存も検査 (`tests/check_modes_beat.py`) |
| `wlsweep_check` | 波長掃引の結合テスト。`data/spectrum_sweep.ofd` の掃引結果が、最終波長を単独指定した実行と**完全一致**することを確認 (波長ごとの再設定に取りこぼしが無いことの検証)。`spectrum.csv` の点数・`lambda*f = c` も検査 (`tests/check_wlsweep.py`) |
| `onn_activation_solver` / `_check` | ONN 光活性化関数の結合テスト。`data/sample/onn_activation.ofd` を実行し、`activation_curve.csv` を解析解 `T = 1/(1 + beta*(P_in/A_eff)*L)` と比較 (相対 7% 以内)。単調非増加・飽和・`P_out <= P_in` も検査 (`tests/check_activation.py`、CI の 3 OS でも同一スクリプトを使用) |

## 実行

入力データは OpenFDTD 形式です。BPM を適用したサンプル (波長 1.55um):

| ファイル | 内容 |
|---|---|
| `data/fiber.ofd` | 光ファイバ (SMF-28 相当, コア半径 4.1um, V=2.26 シングルモード)。LP01 モードの閉じ込め (Marcuse 近似 w=4.7um と比較可能) |
| `data/fiber_offset.ofd` | 光ファイバへの軸ずれ入射 (2um オフセット)。コア中心への引き込みと結合損失 |
| `data/spectrum_sweep.ofd` | 波長掃引 (`wlsweep = 1`、1.50-1.60um の 5 点)。`spectrum.csv` に透過率スペクトルを出力 |
| `data/fiber_mode.ofd` | SMF の基本モード励振 (`launch = mode 0`)。虚軸伝搬法で求めた LP01 を初期界にし、ガウシアン励振より放射損の少ない定常伝搬 |
| `data/fiber_mode_superposition.ofd` | MMF の LP01+LP11 重ね合わせ励振 (`launch = mode 0 1`)。モードビート (実測周期 5710um / 理論 λ/Δneff = 5762um) |
| `data/fiber_mmf.ofd` | マルチモードファイバ (コア径 50um, NA=0.2, V=20.3)。オフセット励振によるコア内反射と閉じ込め |
| `data/fiber_splice.ofd` | コア径不整合の融着接続 (4.1um -> 2.0um)。接続点でのモード変換と放射 (理論 T=0.83) |
| `data/fiber_gap.ofd` | コネクタ間隙 (40um)。間隙での回折とファイバ再入射時の LP01 再整定 |
| `data/fiber_coupler.ofd` | 2 コア方向性結合器 (コア間ギャップ 2um)。エバネッセント結合による電力移行 (62% @ 500um) |
| `data/fiber_attenuator.ofd` | 減衰ファイバ (コアに sigma=13.3 S/m)。コア閉じ込め率を考慮した実効減衰 (0.578 @ 200um, 理論 ~0.57) |
| `data/fiber_bend.ofd` | 曲げファイバ (RoC=2mm, 等価屈折率法)。モードの外側シフトと放射漏れ |
| `data/fiber_gi_mmf.ofd` | GI 型マルチモードファイバ (放物分布の 8 層近似)。自己集束のピッチが理論値と一致 (実測 1150um / 理論 1145um) |
| `data/fiber_taper.ofd` | テーパファイバ (コア半径 4.1 -> 2.0um を階段近似)。断熱モード変換 (P/P0=0.998, 急峻接続との対比) |
| `data/tilt_wideangle.ofd` | 広角 BPM (Pade(1,1))。20 度チルトビームの横変位が厳密値 z*tan(20)=21.8um へ収束 (近軸は z*sin(20)=20.5um で頭打ち) |
| `data/slab_polarization.ofd` | 半ベクトル BPM。高コントラストスラブ (n=2.0/1.0) の界面電界ジャンプ TE 0.875/解析0.876, TM 3.14/解析3.23 |
| `data/waveguide.ofd` | ステップインデックス導波路 (コア 4um 角 n=1.450 / クラッド n=1.444)。コアへのビーム閉じ込め |
| `data/freespace.ofd` | 均一媒質中のガウシアンビーム回折。解析解 w(z)=w0*sqrt(1+(z/zR)^2) と比較可能 |
| `data/lossy.ofd` | 導電率 sigma=100 S/m の吸収媒質。解析解 exp(-2*k0*n''*z) と比較可能 |

```sh
./bin/obpm data/waveguide.ofd
```

- 周波数 (`frequency2`) から波長を決定し、材料 (`material`/`geometry`) の
  比誘電率と導電率から複素屈折率分布 (実部=屈折率, 虚部=吸収) を構築して
  z 軸方向にビームを伝搬します。
- 結果は `time_series_data.h5` (HDF5) に出力されます
  (`/field/Efinal_r`, `/field/Efinal_i` : 最終電界、`/field/n_out_r` : 屈折率分布、
  `/field/Ixz` : 伝搬マップ |E(x, y=Ny/2, z)|^2 (Nz x Nx)、
  `/metadata/lambda`, `/metadata/n_0`, `/metadata/beam_w0` : BPM パラメータ、
  `/metadata/mode_neff` : モード励振 (`launch = mode`) 時の実効屈折率
  (ガウシアン励振では 0))。

### BPM 用キーワード (OpenFDTD 形式の拡張)

| キーワード | 形式 | 意味 |
|---|---|---|
| `beam` | `beam = <w0[m]> [<x0[m]> <y0[m]>]` | 入射ガウシアンビームのウェスト (1/e^2 強度半径) と中心。中心省略時は `feed` 位置 (なければ領域中心) |
| `refindex` | `refindex = <n0>` | BPM の参照屈折率。省略時は領域中心の材料から自動取得 |
| `bend` | `bend = <RoC[m]> [<dir[deg]> [<rho_e>]]` | 曲げ半径 (等価屈折率法 n_eq = n*exp(x/RoC))。dir は曲げ方向 (0 = +x)、rho_e は弾性光学係数。省略時は直線 |
| `polarization` | `polarization = <scalar\|x\|y>` | 半ベクトル BPM の偏波方向。x/y 指定時は偏波方向の界面で D 法線成分連続 (Stern 差分) を反映。省略時はスカラー |
| `wideangle` | `wideangle = <0\|1>` | 広角 BPM (Pade(1,1))。屈折率項を演算子内部に含めた一般化 ADI で伝搬。省略時は近軸 |
| `beamtilt` | `beamtilt = <tx[deg]> [<ty[deg]>]` | 入射ビームの傾き (横方向波数の位相ランプ)。`launch = mode` と併用した場合はモード界にも同じ位相ランプが適用される |
| `frames` | `frames = <interval>` | `|E(x,y)|^2` スナップショットの記録間隔 (z ステップ単位)。`/field/frames` (nframes x Ny x Nx) へ出力。省略時 (0) は記録なし |
| `launch` | `launch = <gauss\|mode [<m>[:<coef>] ...]>` | 励振方法。`mode` は先頭スライスの屈折率分布から導波モードを虚軸伝搬法 (`bpm/modes.cpp`) で求めて励振する (引数省略時は基本モード 0)。複数指定すると**重ね合わせ**、`<m>:<coef>` で係数を与える (省略時 1、最大 8 モード)。入力電力は係数に依らず `feed` 電圧^2 に正規化されるため係数は分岐比のみを決める。モードが見つからない場合はガウシアンへフォールバック。省略時は `gauss` |
| `wlsweep` | `wlsweep = <0\|1>` | 波長掃引。`frequency2` の全点を順に計算し、各波長の透過率を `spectrum.csv` (lambda_m, frequency_Hz, transmission) へ出力する。`/field` と HDF5 メタデータは**最終波長**の結果。省略時 (0) は先頭波長のみ = 従来動作。**CPU 版のみ** |
| `tpa` | `tpa = <material id> <beta[cm/GW]>` | 材料に二光子吸収 (TPA) 係数 beta を付与 (複数行可)。省略時は線形計算 |
| `powersweep` | `powersweep = <Pmin[W]> <Pmax[W]> <npoints> [log\|lin]` | 入力パワー掃引 (既定 lin)。`activation_curve.csv` に P_out(P_in) を出力。省略時は従来の単発計算 |

- `polarization` / `wideangle` 指定時は倍精度の一般化伝搬エンジンを使用します
  (CPU 版 `bpm/wabpm.cpp` / CUDA 版 `bpm/wabpm.cu`、同一アルゴリズム)。
  `bend` との併用も可能です (等価屈折率法を屈折率分布に反映)。

- `feed` の電圧はビーム振幅として使用されます。
- メッシュは等間隔を推奨します (不均一の場合は平均セル幅で計算し警告を表示)。

### ONN 光活性化関数 (非線形吸収 TPA + パワー掃引)

光ニューラルネットワーク (ONN) の光活性化関数を BPM で軽量に設計反復するための
機能です。メタマテリアル装荷 Si 導波路の二光子吸収 (TPA) による P_out(P_in) の
飽和特性 (ReLU 相当) を対象としています。

出典: K. Honda, Y. Shoji, T. Amemiya, Opt. Lett. 49, 5811 (2024)
(TPA 係数 beta = 424 cm/GW、P_out(P_in) の飽和特性で ReLU 相当の活性化関数)

- `tpa = <material id> <beta[cm/GW]>` を指定すると、その材料のセルで各伝搬
  ステップに `E *= exp(-(beta/2)*I*dz)` の非線形減衰を適用します
  (強度の減衰率 alpha = beta*I に対し界は alpha/2。
  単位変換は 1 cm/GW = 1e-11 m/W、例: 424 cm/GW = 4.24e-9 m/W)。
- `tpa` / `powersweep` 指定時は初期界を `∫∫|E|^2 dA = P_in [W]` に正規化し、
  `|E|^2` がそのまま強度 I [W/m^2] になります (物理スケーリング)。
  未指定時は従来通り無次元の界で計算します (後方互換)。
- `powersweep = <Pmin[W]> <Pmax[W]> <npoints> [log|lin]` で入力パワーを掃引し、
  各点の出力パワー `P_out = ∫∫|E_end|^2 dA` を **activation_curve.csv**
  (列: `P_in_W, P_out_W, transmission`) に出力します。obpm.log には
  `ONN: P_in=... W -> P_out=... W (T=...)` の行が残ります。
  `powersweep` 省略で `tpa` のみの場合は P_in = 1 W の単発計算になります。
- 実効断面積 `A_eff = (∫|E|^2 dA)^2 / ∫|E|^4 dA` を初期界から実計算して
  obpm.log に出力します。一様断面の直線導波路では平面波近似の解析解
  `T(P_in) = 1 / (1 + beta*(P_in/A_eff)*L)` と比較できます。
- サンプル: `data/sample/onn_activation.ofd` (fiber.ofd ベースの直線導波路 +
  `tpa = 3 424` + log 掃引 8 点)。全掃引点で上記解析解と ±7% 以内で一致し
  (最大 +4.1% @ 最深飽和点)、小パワーで線形透過 (T = 0.98)、大パワーで飽和
  (T = 0.30) の単調非増加曲線になります。深い飽和 (beta*P*L/A_eff > 3 程度)
  では TPA がモード中心を選択的に焼いて界形状が平坦化するため、固定 A_eff の
  解析解との差が拡大します (平面波近似の限界)。
- `/field` 出力 (Ixz/Efinal/frames) は最終掃引点 (最大パワー) の結果です。
- `wideangle` / `polarization` の一般化伝搬エンジンとも併用できます。
  CUDA 版 (`obpm_cuda`) は現状 `tpa` / `powersweep` に未対応です (無視されます)。

## CI / Release

GitHub Actions で以下のジョブを実行します。トリガーは **main への push・タグ push・
PR (base = main)・手動実行 (workflow_dispatch)** です。

| ジョブ | 実行条件 | 内容 |
|---|---|---|
| `build-cpu` (Linux/gcc) | すべて | CPU ビルド + `ctest` + `fiber.ofd` スモーク + `data/` 全サンプルスモーク + ONN 活性化関数の定量検証 |
| `build-macos` (AppleClang + Homebrew libomp/Eigen) | **PR 以外** (main push / タグ / 手動) | CPU ビルド + `ctest` + スモーク + ONN 定量検証 |
| `build-windows` (MSVC + Ninja + vcpkg) | すべて | CPU ビルド + `ctest` + スモーク + ONN 定量検証 |
| `build-cuda` | すべて | CUDA Toolkit を導入し `obpm_cuda` のコンパイル・リンクを検証 (ランナーに GPU は無いため実行はしない) |
| `build-mpi` | すべて | MPI + 並列 HDF5 を導入し全ターゲット (`obpm` / `obpm_post` / `obpm_mpi`) のビルドを検証 |

- ONN 活性化関数の検証は各 OS とも `tests/check_activation.py` (解析解との
  相対 7% 比較) を使用します
- ビルド成果物は artifact (`obpm-linux-x64` / `obpm-macos-arm64` /
  `obpm-windows-x64`) に保存
- `v*` タグを push すると GitHub Release にバイナリが自動添付されます

> ℹ **Actions 課金についての注記** (private リポジトリの無料枠は月 2,000 分、
> macOS ×10 / Windows ×2 換算):
> - push トリガーを main とタグに限定しています。作業ブランチは PR イベント側で
>   カバーされるため、PR を持つブランチでの二重実行 (push + pull_request) は
>   発生しません。PR の無いブランチは workflow_dispatch で手動実行できます。
> - 同一ブランチへの連続 push は `concurrency` により古い run を自動キャンセルします。
> - 課金レートが最も高い macOS ジョブは PR では実行しません (上表)。
> - **全ジョブが数秒で failure になる場合** (`runner_id: 0`、ステップ実行なし、
>   ログ 404) はコードではなく無料枠の枯渇・課金設定が原因です。最小ワークフロー
>   `Ping` (`.github/workflows/ping.yml`、echo のみ) を手動実行すると 1 分未満で
>   切り分けられます。

## 姉妹リポジトリ

| リポジトリ | 手法 |
|---|---|
| [OpenFDTD](https://github.com/Sirokujira/OpenFDTD) | 電磁 FDTD |
| [OpenRCWA](https://github.com/Sirokujira/OpenRCWA) | 周期構造 RCWA |
| [OpenFDTD-X](https://github.com/Sirokujira/OpenFDTD-X) | Qt6 GUI フロントエンド |

# Reference
OpenFDTD
http://www.e-em.co.jp/OpenFDTD/

BPM-MATLAB
https://gitlab.gbar.dtu.dk/biophotonics/Programs/BPM-Matlab

BPM
https://github.com/warthan07/Nemaktis

# 屈折率資料
https://github.com/polyanskiy/refractiveindex.info-database
