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

> ⚠ **既知事項**: `obpm_post` は OpenFDTD 由来のポスト処理
> (obpm.out ベースの FDTD 量の描画) のままで、BPM が
> `time_series_data.h5` の `/field` に書く伝搬マップ (`Ixz`) や
> 最終電界 (`Efinal_*`) の可視化には未対応です。当面は HDF5 を
> Python (h5py) などで直接読んでください。GUI (OpenFDTD-X) の
> H5 ビューアからも参照できます。

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
  - ビルドは CUDA 12.0 で確認済みです (GPU 実機での実行検証は未実施)。
- MPI 版: `-DWITH_MPI=ON` (要 MPI)
  - **MPI 版は OpenFDTD 由来の FDTD (時間領域) ソルバのみで、BPM は未対応です**。
    BPM の並列化は OpenMP (CPU 版 `obpm`) または CUDA 版 `obpm_cuda` を使用してください
    (BPM の ADI は z 方向に逐次のため、MPI 領域分割は転置通信が必要で費用対効果が低い)。

## 可視化

出力 HDF5 を PNG / GIF に変換する後処理スクリプトを用意しています
(要 h5py, numpy, matplotlib):

```sh
python3 tools/plot_ixz.py time_series_data.h5 [--db] [--prefix out] [--fps 15]
```

- `*_ixz.png` : 伝搬マップ `/field/Ixz` のヒートマップ
- `*_final.png` : 最終電界 `|E(x,y)|^2` と屈折率分布
- `*_modes.png` : 導波モード形状と neff (入力に `modes = <nModes>` を指定して
  `/modes` を出力した場合のみ)
- `*_prop.gif` : `|E(x,y)|^2` の伝搬アニメーション (入力に `frames = <interval>` を
  指定して `/field/frames` を記録した場合のみ)

## モード解析

虚軸伝搬法 (imaginary-distance BPM) + Gram-Schmidt 直交化によるモードソルバを
`bpm/modes.cpp` (`wabpm_find_modes`) に実装しています。実効屈折率の降順に導波モード
(n_clad < neff を満たすもののみ) を求め、neff は Rayleigh 商から算出します。

- 入力キーワード **`modes = <nModes> [excite]`** で `obpm` / `obpm_cuda` から実行できます
  (CUDA 版もホスト側で同一の解析を実行)。結果は obpm.log (`mode <i> : neff = ...`) と
  HDF5 (`/modes/mode<i>` : Ny×Nx 実数場、`/modes/neff`) に出力されます。
  解析は入力断面 (z 始端) の屈折率分布に対して行い、曲げ (bend) は反映しません。
  半ベクトルの直交化は近似となるためスカラー演算子で解析します。
- `excite` を付けると入射ビームがガウシアンではなく数値的に厳密な基本モードに
  置き換わります (モード整合励振)。サンプル `data/sample/fiber_modes.ofd` では
  LP01 の neff = 1.447167 が分散方程式の厳密解と一致し、電力保存が
  T = 0.99999 (ガウシアン励振では 0.99988) に向上します。
- `include/bpm/allset.hpp` の `findModes()` (BPM-MATLAB 互換 API) も同じソルバを
  使用します。

## テスト

解析解と比較する単体テストを `tests/` に用意しています。

```sh
cmake -S . -B build -DWITH_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

| テスト | 検証項目 |
|---|---|
| `test_wabpm` | 自由空間ガウシアン回折 `w(z)=w0*sqrt(1+(z/zR)^2)` (近軸/広角)、エネルギー保存、吸収減衰 `exp(-2*k0*n''*z)`、曲げ偏向 `<x>=z^2/(2*RoC)`、半ベクトル差分の一様媒質整合 |
| `test_modes` | ステップインデックスファイバ LP01/LP11 の実効屈折率 (分散方程式の厳密解と比較)、モード直交性 |
| `test_allset` | `findModes` / `modeSuperposition` / `offsetField` / `tiltField` (P_Struct API) |

## 実行

入力データは OpenFDTD 形式です。BPM を適用したサンプル (波長 1.55um):

| ファイル | 内容 |
|---|---|
| `data/fiber.ofd` | 光ファイバ (SMF-28 相当, コア半径 4.1um, V=2.26 シングルモード)。LP01 モードの閉じ込め (Marcuse 近似 w=4.7um と比較可能) |
| `data/fiber_offset.ofd` | 光ファイバへの軸ずれ入射 (2um オフセット)。コア中心への引き込みと結合損失 |
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
  `/metadata/lambda`, `/metadata/n_0`, `/metadata/beam_w0` : BPM パラメータ)。

### BPM 用キーワード (OpenFDTD 形式の拡張)

| キーワード | 形式 | 意味 |
|---|---|---|
| `beam` | `beam = <w0[m]> [<x0[m]> <y0[m]>]` | 入射ガウシアンビームのウェスト (1/e^2 強度半径) と中心。中心省略時は `feed` 位置 (なければ領域中心) |
| `refindex` | `refindex = <n0>` | BPM の参照屈折率。省略時は領域中心の材料から自動取得 |
| `bend` | `bend = <RoC[m]> [<dir[deg]> [<rho_e>]]` | 曲げ半径 (等価屈折率法 n_eq = n*exp(x/RoC))。dir は曲げ方向 (0 = +x)、rho_e は弾性光学係数。省略時は直線 |
| `polarization` | `polarization = <scalar\|x\|y>` | 半ベクトル BPM の偏波方向。x/y 指定時は偏波方向の界面で D 法線成分連続 (Stern 差分) を反映。省略時はスカラー |
| `wideangle` | `wideangle = <0\|1>` | 広角 BPM (Pade(1,1))。屈折率項を演算子内部に含めた一般化 ADI で伝搬。省略時は近軸 |
| `beamtilt` | `beamtilt = <tx[deg]> [<ty[deg]>]` | 入射ビームの傾き (横方向波数の位相ランプ) |
| `frames` | `frames = <interval>` | `|E(x,y)|^2` スナップショットの記録間隔 (z ステップ単位)。`/field/frames` (nframes x Ny x Nx) へ出力。省略時 (0) は記録なし |
| `taper` | `taper = <ratio>` | 出口での横方向スケール比 (入口 = 1)。屈折率分布を z に沿って相似縮小/拡大する (座標変換)。省略時 1 |
| `twist` | `twist = <rate[deg/m]>` | 導波路のツイスト率。屈折率分布を z に沿って回転させる。省略時 0 |
| `modes` | `modes = <nModes> [excite]` | 入力断面のモード解析 (虚軸伝搬法)。導波条件 (n_clad < neff) を満たすモードを neff 降順に求め、`/modes/mode<i>` と `/modes/neff` に出力。`excite` 指定で入射ビームを基本モードに置き換え (モード整合励振)。省略時は解析なし |
| `tpa` | `tpa = <material id> <beta[cm/GW]>` | 材料に二光子吸収 (TPA) 係数 beta を付与 (複数行可)。省略時は線形計算 |
| `powersweep` | `powersweep = <Pmin[W]> <Pmax[W]> <npoints> [log\|lin]` | 入力パワー掃引 (既定 lin)。`activation_curve.csv` に P_out(P_in) を出力。省略時は従来の単発計算 |

- `polarization` / `wideangle` 指定時は倍精度の一般化伝搬エンジンを使用します
  (CPU 版 `bpm/wabpm.cpp` / CUDA 版 `bpm/wabpm.cu`、同一アルゴリズム)。
  `bend` との併用も可能です (等価屈折率法を屈折率分布に反映)。

- `feed` の電圧はビーム振幅として使用されます。
- メッシュは等間隔を推奨します (不均一の場合は平均セル幅で計算し警告を表示)。

### テーパ / ツイスト (座標変換)

`taper` / `twist` を指定すると、屈折率分布を伝搬に沿って相似縮小・回転させます
(BPM-MATLAB の座標変換方式)。`geometry` を z 方向に階段近似する必要がありません
(対比: `data/fiber_taper.ofd` は階段近似)。

- **指定時は入力断面 (z 始端) の 2D 分布のみを参照します**。`geometry` 自体の
  z 変化とは併用できません (入力断面が変換されて全 z に適用されます)。
- スケール係数は `1 - (1-ratio)*(iz/Nz)`、回転角は `twist * z` です。
  `bend` / `wideangle` / `polarization` とは併用可能で、CPU の近軸・拡張の
  両経路と CUDA 版に実装されています。
- サンプル: `data/sample/taper_twist.ofd` (矩形コア 8x2um を 0.5 倍に縮小しつつ
  200um で 90 度回転)。出力屈折率分布 `/field/n_out_r` が実際に縮小・回転したもの
  になり (主軸 89.1 度 / 理論 90 度)、変換は断熱的で電力が保存します (P/P0 = 0.9985)。
- 検証: 一様媒質では変換を掛けても結果不変 (拡張パスは倍精度で厳密一致)、
  円形コアの taper では出力コア半径が `a*ratio` に一致 (2.08um / 理論 2.05um)。

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
  CUDA 版 (`obpm_cuda`) は現状 `tpa` / `powersweep` に未対応です
  (指定時は警告を表示して線形の単発計算になります。CPU 版 `obpm` を使用してください)。

## CI / Release

- push / PR ごとに Linux (gcc) と macOS (AppleClang + Homebrew libomp/Eigen)
  で CPU ビルド + `data/sample/fiber.ofd` のスモーク実行 (`normal end` 判定)
- 解析解と比較する単体テスト (`ctest` : test_wabpm / test_modes / test_allset) を
  Linux / macOS で実行
- `data/sample/onn_activation.ofd` の ONN 活性化関数スモーク
  (`activation_curve.csv` の単調飽和 + 平面波近似の解析解
  `T = 1/(1 + beta*(P_in/A_eff)*L)` との ±8% 一致判定, `tools/check_activation.sh`)
  を 3 OS で実行 (解析解判定は Linux / macOS)
- ビルド成果物は artifact (`obpm-linux-x64` / `obpm-macos-arm64`) に保存
- `v*` タグを push すると GitHub Release にバイナリが自動添付されます

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
