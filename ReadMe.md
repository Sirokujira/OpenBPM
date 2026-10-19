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
>   BPM が `time_series_data.h5` に書く伝搬マップ (`/field/Ixz`)・
>   最終電界 (`Efinal_*`)・伝搬スナップショット (`frames`、等間隔に最大 6 枚)・
>   **z 伝搬に沿ったスカラー推移 (`/trace`) の折れ線グラフ**の可視化にも
>   対応しています (`post/postbpm.c`)。`/trace` はパワー・ピーク強度・
>   ビーム幅・重心・モード結合率を 1 量ずつページに分けて描きます
>   (単位の異なる量を同じ縦軸に混ぜないため)。
>   モード形状 (`/modes`) の表示や GIF アニメーションは
>   `tools/plot_ixz.py` を使うか、HDF5 を Python (h5py) で直接読んでください。
>   GUI (OpenFDTD-X) の H5 ビューアからも参照できます。
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
# Ubuntu / Debian
sudo apt install cmake g++ libhdf5-dev libeigen3-dev

# macOS (AppleClang は OpenMP を同梱しないため libomp が必須)
brew install libomp hdf5 eigen
```

```sh
cmake -S . -B build
cmake --build build -j
```

実行ファイルは `bin/` に出力されます (`obpm`, `obpm_post`)。

> **ビルドタイプ**: `CMAKE_BUILD_TYPE` を指定しない場合は **Release**
> (最適化あり) を既定にします。以前は未指定だと GCC/Clang が `-O0` で
> コンパイルするため、同じソースでも 2.4〜2.7 倍遅いバイナリになっていました
> (実測: fiber 0.76 → 0.28 秒、slab_polarization 51.5 → 21.1 秒。
> 数値は `-O0` と Release でビット一致するので純粋な速度差)。
> デバッグビルドは `-DCMAKE_BUILD_TYPE=Debug` を明示してください。

> **並列実行について**: BPM の伝搬は OpenMP で並列化されています
> (4 コアでの実測: 近軸 1.9 倍、広角 2.4 倍、半ベクトル 2.3 倍)。
> 同一スレッド数なら結果はビット単位で再現しますが、スレッド数が変わると
> 電力集約の部分和の個数が変わるため最終桁が動きます
> (`output power` の 7 桁はスレッド数に依らず不変)。

> macOS では CMake が Homebrew の `libomp` を自動検出します
> (`brew --prefix libomp`)。Homebrew を使っていない場合や別の場所に
> インストールしている場合は、OpenMP のパスを明示指定してください:
>
> ```sh
> OMP=/path/to/libomp
> cmake -S . -B build \
>   -DOpenMP_C_FLAGS="-Xclang -fopenmp -I${OMP}/include" \
>   -DOpenMP_C_LIB_NAMES=omp \
>   -DOpenMP_CXX_FLAGS="-Xclang -fopenmp -I${OMP}/include" \
>   -DOpenMP_CXX_LIB_NAMES=omp \
>   -DOpenMP_omp_LIBRARY="${OMP}/lib/libomp.dylib"
> ```

- CUDA 版: `-DWITH_CUDA=ON` (要 CUDA Toolkit)
  - `obpm_cuda` は CPU 版と同じ BPM 機能 (beam/refindex/導電率吸収/メッシュ警告/bend/
    beamtilt、広角 `wideangle`・半ベクトル `polarization` を含む) に対応しています。
    スカラー近軸は `bpm/FDBPMpropagator.cu`、広角/半ベクトルは `bpm/wabpm.cu` を使用します。
  - `tpa` / `powersweep` (ONN 光活性化関数) も CUDA 版に実装済みです
    (両伝搬パスの TPA カーネル・物理スケーリング・A_eff 算出・掃引ループ・
    `activation_curve.csv` 出力)。ただし **GPU 実機での実行検証は未実施**のため、
    実行時に「実機未検証」の note を表示します。精度が要る用途では CPU 版
    (`obpm`) と突き合わせてください。
  - `pml` (PML 吸収境界) も CUDA 版に実装済みです (近軸・拡張の両カーネルで
    面ごとの伸長係数を使用)。こちらも **GPU 実機での実行検証は未実施**です。
  - **`wlsweep` (波長掃引) は CPU 版のみ**です。CUDA 版で指定した場合は実行時に
    warning を表示して無視します (先頭波長のみ計算されます)。
  - ビルドは CUDA 12.0 で確認済みです (GPU 実機での実行検証は未実施)。
- MPI 版: `-DWITH_MPI=ON` (要 MPI + 並列 HDF5。例: `libopenmpi-dev libhdf5-openmpi-dev`)
  - `obpm_mpi` は OpenFDTD 由来の **FDTD 時間領域ソルバー (BPM ではない)**。
    `mpirun -np <N> ./bin/obpm_mpi <input>.ofd` で実行します。
  - **BPM は MPI 未対応**です。BPM の並列化は OpenMP (CPU 版 `obpm`) または
    CUDA 版 `obpm_cuda` を使用してください (BPM の ADI は z 方向に逐次のため、
    MPI 領域分割は転置通信が必要で費用対効果が低い)。

## 可視化

出力 HDF5 を PNG / GIF に変換する後処理スクリプトを用意しています
(要 h5py, numpy, matplotlib):

```sh
python3 tools/plot_ixz.py time_series_data.h5 [--db] [--prefix out] [--fps 15]
```

- `*_ixz.png` : 伝搬マップ `/field/Ixz` のヒートマップ
- `*_final.png` : 最終電界 `|E(x,y)|^2` と屈折率分布
- `*_trace.png` : z ごとのスカラー推移 (`/trace` : 断面パワー・ピーク強度・重心・ビーム幅)
- `*_overlap.png` : 各導波モードへのパワー占有率の z 推移 (`/trace/overlap`。
  入力に `modes = <nModes>` を指定した場合のみ)
- `*_modes.png` : 導波モード形状と neff (入力に `modes = <nModes>` を指定して
  `/modes` を出力した場合のみ)
- `*_prop.gif` : `|E(x,y)|^2` の伝搬アニメーション (入力に `frames = <interval>` を
  指定して `/field/frames` を記録した場合のみ)
- `*_phase.gif` : 位相 `arg(E)` の伝搬アニメーション (`frames = <interval> complex` で
  `/field/frames_r`, `/field/frames_i` を記録した場合のみ)

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
  T = 0.99999 (ガウシアン励振では 0.99987) に向上します。
- `modes` 指定時は伝搬中の**モード結合率** η_m(z) も `/trace/overlap` に記録されます。
  「基本モードにどれだけパワーが残っているか」を z の関数として GUI 表示できます
  (`tools/plot_ixz.py` の `*_overlap.png`)。`fiber_modes.ofd` では
  `excite` 指定時 η_1 ≈ 0.9993、ガウシアン励振では η_1 ≈ 0.9944 となります。
- `include/bpm/allset.hpp` の `findModes()` (BPM-MATLAB 互換 API) も同じソルバを
  使用します。

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

入力データは OpenFDTD 形式です。1 行目のプログラム名は
`OpenBPM` / `OpenFDTD` / `OpenTHFD` のいずれでも受け付けます
(OpenFDTD-X (GUI) が書き出す `.ofd` は `OpenFDTD` ヘッダです)。
それ以外の名前だった場合は、実際に読んだ名前を添えて

```
*** not OpenBPM data : 1st line = "..." (expected OpenBPM / OpenFDTD / OpenTHFD)
```

と表示して終了します。

BPM を適用したサンプル (波長 1.55um):

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
  `/field/Iyz` : 伝搬マップ |E(x=Nx/2, y, z)|^2 (Nz x Ny)、
  `/metadata/lambda`, `/metadata/n_0`, `/metadata/beam_w0` : BPM パラメータ、
  `/metadata/mode_neff` : モード励振 (`launch = mode`) 時の実効屈折率
  (ガウシアン励振では 0))。

#### 系列データの独立変数は「時間」ではなく z

**BPM は定常 (単一周波数) 解法なので時間軸を持ちません。** 出力ファイル名
`time_series_data.h5` は OpenFDTD 由来の名残で、`/metadata/Ntime` は 1 固定、
`/metadata/time` も 1 点しかありません。FDTD の時間軸に相当する「行進 (marching)
軸」は **z (伝搬距離)** です。

この対応は上流の [BPM-Matlab](https://github.com/ankrh/BPM-Matlab) と同じ規約です
(BPM-Matlab も `P.z` / `P.xzSlice` / `P.modeOverlaps` のように z を独立変数として
保持し、`imagesc(P.z, P.x, abs(P.xzSlice).^2)` のように z 軸で描画します)。
RSoft BeamPROP もモニタ値を Z 軸に対して表示します。z を時間軸に見せかける
実装は行いません。

読み手が推測せずに済むよう、HDF5 側に軸の宣言を書き出しています:

| 属性の場所 | 属性 | 値 |
|---|---|---|
| ルート `/` | `solver` | `OpenBPM (beam propagation method)` |
| ルート `/` | `domain` | `steady-state (single frequency, no time axis)` |
| ルート `/` | `marching_axis` | `z` |
| ルート `/` | `marching_axis_values` | `/trace/z` |
| ルート `/` | `marching_steps` | 伝搬ステップ数 |
| ルート `/` | `time_dependent` | `0` |
| `/trace` グループ | `axis` / `axis_values` / `npoints` | `z` / `/trace/z` / ステップ数 |
| 各データセット | `long_name` / `units` / `dims` / `coordinates` | 例: `z,x` と `/trace/z /metadata/Xc` |

GUI や汎用ビューアは `marching_axis` と `marching_axis_values` を見れば、
`Ntime` を参照せずに横軸を決められます。

```python
import h5py
f = h5py.File("time_series_data.h5", "r")
print(f.attrs["marching_axis"])          # b'z'
zpath = f.attrs["marching_axis_values"]  # b'/trace/z'
z = f[zpath.decode()][:]                 # 横軸 [m]
p = f["/trace/power"][:]                 # 縦軸
print(f["/trace/power"].attrs["units"])  # b'arb. unit'
```

#### 伝搬に沿った系列データ (`/trace`) — GUI 表示用

z ステップごとの断面統計量を 1 次元配列 (長さ = 伝搬ステップ数) で出力します。
画面側でそのままグラフ表示できます (`tools/plot_ixz.py` の `*_trace.png` も参照)。

| データセット | 内容 |
|---|---|
| `/trace/z` | 断面位置 z [m] (t 番目のステップ実行後) |
| `/trace/power` | 断面の総パワー ∫∫\|E\|^2 dA (`tpa`/`powersweep` 使用時は [W]) |
| `/trace/peak` | \|E\|^2 の最大値 |
| `/trace/centroid_x`, `/trace/centroid_y` | 強度重心 [m] |
| `/trace/width_x`, `/trace/width_y` | 強度の 2 次モーメント幅 2σ [m] |
| `/trace/overlap` | 各導波モードへのパワー占有率 η_m(z) (nModes x Nz)。`modes = <n>` 指定時のみ |

`/trace/overlap` は η_m = \|<φ_m, E>\|^2 / (\|\|φ_m\|\|^2 \|\|E\|\|^2) で、
モード順序は `/modes/neff` と同じ (実効屈折率の降順)。直線導波路では各 η_m は
z によらずほぼ一定になり、テーパ/ツイスト/曲げではモード間の変換が推移として現れます。
Σ_m η_m < 1 のぶんが放射モード (非導波成分) です。

座標軸は `/metadata/Xc`, `Yc`, `Zc` (セル中心) と `Xn`, `Yn`, `Zn` (節点) を使用できます。
`frames` を記録した場合は各フレームの z 位置が `/field/frames_z` (nframes) に
出力されるため、`frame_interval` と `grid_dz` から計算し直す必要はありません。
`/trace` と `/field/Iyz` は常に出力され、CPU 版・CUDA 版で同一内容です
(CUDA 版はデバイス上で集計するため転送コストは断面 1 行/列ぶんに収まります)。

### BPM 用キーワード (OpenFDTD 形式の拡張)

| キーワード | 形式 | 意味 |
|---|---|---|
| `beam` | `beam = <w0[m]> [<x0[m]> <y0[m]>]` | 入射ガウシアンビームのウェスト (1/e^2 強度半径) と中心。中心省略時は `feed` 位置 (なければ領域中心) |
| `refindex` | `refindex = <n0>` | BPM の参照屈折率。省略時は領域中心の材料から自動取得 |
| `bend` | `bend = <RoC[m]> [<dir[deg]> [<rho_e>]]` | 曲げ半径 (等価屈折率法 n_eq = n*exp(x/RoC))。dir は曲げ方向 (0 = +x)、rho_e は弾性光学係数。省略時は直線 |
| `polarization` | `polarization = <scalar\|x\|y>` | 半ベクトル BPM の偏波方向。x/y 指定時は偏波方向の界面で D 法線成分連続 (Stern 差分) を反映。省略時はスカラー |
| `wideangle` | `wideangle = <0\|1>` | 広角 BPM (Pade(1,1))。屈折率項を演算子内部に含めた一般化 ADI で伝搬。省略時は近軸 |
| `beamtilt` | `beamtilt = <tx[deg]> [<ty[deg]>]` | 入射ビームの傾き (横方向波数の位相ランプ)。`launch = mode` / `modes = ... excite` と併用した場合はモード界にも同じ位相ランプが適用される |
| `frames` | `frames = <interval> [complex]` | `|E(x,y)|^2` スナップショットの記録間隔 (z ステップ単位)。`/field/frames` (nframes x Ny x Nx) へ出力。`complex` 指定時は複素電界も `/field/frames_r`, `/field/frames_i` へ出力 (位相表示用、容量 3 倍)。省略時 (0) は記録なし |
| `taper` | `taper = <ratio>` | 出口での横方向スケール比 (入口 = 1)。屈折率分布を z に沿って相似縮小/拡大する (座標変換)。省略時 1 |
| `twist` | `twist = <rate[deg/m]>` | 導波路のツイスト率。屈折率分布を z に沿って回転させる。省略時 0 |
| `symmetry` | `symmetry = <x\|y\|xy> [sym\|anti]` | 対称境界。指定軸のメッシュ始端を鏡像面とし計算領域を 1/2〜1/4 に削減 (既定 sym)。**メッシュ・出力とも半領域**になる。省略時は対称境界なし |
| `ripfile` | `ripfile = <path>` | 屈折率分布の直接入力 (geometry の代わりに使用)。`.csv` は 2D (Ny 行 x Nx 列の実数 n、z 不変)、`.h5`/`.hdf5` はデータセット `/rip/n` (2D: Ny x Nx / 3D: Nz x Ny x Nx)。寸法はメッシュと厳密一致。相対パスは入力 .ofd のディレクトリ基準でも解決。TPA の β は従来どおり geometry の材料 ID から決まる。省略時は geometry から構築 (従来動作) |
| `pml` | `pml = <width[m]> [<R0>]` | 横方向境界の PML (複素座標伸長)。`width` = 層厚 [m] (4 辺の内側に確保)、`R0` = 目標往復振幅反射率 (既定 `1e-20`、かすめ入射 kx = k0·n0 基準で `sigma_max = -(m+1)ln(R0)/(2 k0 n0 W)`, m = 3)。指定時は従来の端部吸収体を無効化する。省略時は従来の吸収体 (後方互換) |
| `modes` | `modes = <nModes> [excite]` | 入力断面のモード解析 (虚軸伝搬法)。導波条件 (n_clad < neff) を満たすモードを neff 降順に求め、`/modes/mode<i>` と `/modes/neff` に出力。`excite` 指定で入射ビームを基本モードに置き換え (モード整合励振)。省略時は解析なし |
| `launch` | `launch = <gauss\|mode [<m>[:<coef>] ...]>` | 励振方法。`mode` は先頭スライスの屈折率分布から導波モードを虚軸伝搬法 (`bpm/modes.cpp`) で求めて励振する (引数省略時は基本モード 0)。複数指定すると**重ね合わせ**、`<m>:<coef>` で係数を与える (省略時 1、最大 8 モード)。入力電力は係数に依らず `feed` 電圧^2 に正規化されるため係数は分岐比のみを決める。モードが見つからない場合はガウシアンへフォールバック。省略時は `gauss` |
| `wlsweep` | `wlsweep = <0\|1>` | 波長掃引。`frequency2` の全点を順に計算し、各波長の透過率を `spectrum.csv` (lambda_m, frequency_Hz, transmission) へ出力する。`/field` と HDF5 メタデータは**最終波長**の結果。省略時 (0) は先頭波長のみ = 従来動作。**CPU 版のみ** |
| `tpa` | `tpa = <material id> <beta[cm/GW]>` | 材料に二光子吸収 (TPA) 係数 beta を付与 (複数行可)。省略時は線形計算 |
| `powersweep` | `powersweep = <Pmin[W]> <Pmax[W]> <npoints> [log\|lin]` | 入力パワー掃引 (既定 lin)。`activation_curve.csv` に P_out(P_in) を出力。省略時は従来の単発計算 |

- `polarization` / `wideangle` 指定時は倍精度の一般化伝搬エンジンを使用します
  (CPU 版 `bpm/wabpm.cpp` / CUDA 版 `bpm/wabpm.cu`、同一アルゴリズム)。
  `bend` との併用も可能です (等価屈折率法を屈折率分布に反映)。

- `feed` の電圧はビーム振幅として使用されます。
- メッシュは等間隔を推奨します (不均一の場合は平均セル幅で計算し警告を表示)。

### 対称境界 (計算領域の削減)

`symmetry` を指定すると、指定軸のメッシュ始端を鏡像面として計算領域を 1/2〜1/4 に
削減できます (物理は変わらず計算量だけが減ります)。

- **メッシュには半領域のみを与え、出力も半領域になります** (例: `xmesh = 0 45 15e-6`)。
  出力パワーも全領域の 1/2 (両軸なら 1/4) になります。
- 鏡像面はセル中心の半セル手前 (= メッシュ始端) です。セル中心が始端から dx/2, 3dx/2, ...
  に並ぶため、半領域の格子は全領域の該当半分と厳密に一致します。
- `sym` (既定) は対称、`anti` は反対称 (鏡像面で符号反転) です。
- CPU の近軸・拡張 (広角/半ベクトル) の両経路と CUDA 版に実装されています。
- 検証: 一様媒質 (厳密に対称な構造) で半領域/1/4 領域の結果が全領域の該当部分と一致
  (近軸 2e-6 = float32 の丸め、拡張パスは倍精度で**厳密に 0**)。
  サンプル `data/sample/fiber_symmetry.ofd` (ファイバを x 半領域で計算、
  全領域の右半分と 3.5e-6 で一致し出力パワーは正確に 1/2)。
- **離散化後の**構造が鏡像面について対称である必要があります。BPM はセル中心
  (Xc, Yc) で材料をラスタライズするため、対称なメッシュと形状なら x/y どちらの
  対称面も使えます (円形コアで `symmetry = y` も全領域のちょうど 1/2 の
  パワーになることを確認済み)。

### 屈折率分布の直接入力 (ripfile)

geometry プリミティブの代わりに、任意の屈折率分布をファイルで与えられます
(上流 BPM-Matlab の「RIP を行列で与える」使い方に対応)。

```
ripfile = grin_rip.csv
```

- `.csv` : Ny 行 × Nx 列の実数 n (区切りはカンマ/空白/タブ、`#` 行はコメント)。
  行 0 が y 最小。全スライス共通 (z 不変)
- `.h5` / `.hdf5` : データセット `/rip/n`。2D (Ny × Nx) または
  3D (Nz × Ny × Nx) の実数 n。3D は z 方向に変化する構造を直接与えられる
- 寸法はメッシュ (`xmesh`/`ymesh`/`zmesh`) と**厳密一致**が必要
  (不一致は期待寸法つきのエラーで停止)
- 相対パスは実行時のカレントディレクトリに加え、**入力 .ofd の
  ディレクトリ基準**でも解決されます (GUI からの起動を想定)
- `ripfile` が置き換えるのは屈折率のみです。`tpa` の β は従来どおり
  geometry の材料 ID から決まるため、両者は併用できます
- 検証サンプル `data/sample/grin.ofd`: 放物型 GRIN
  (`tools/make_grin_rip.py` で生成した `grin_rip.csv`) にオフセット入射した
  ビーム重心の自己集束ピッチが、理論値 Λ = 2πa/√(2Δ) = 314.16 µm に対して
  コサインフィットで 315.6 µm (**誤差 0.46%**)。3D HDF5 経路は、fiber.ofd の
  geometry と同じ n 値を 3D `/rip/n` で与えた場合に全出力が**ビット一致**する
  ことで検証済み

### PML 吸収境界 (複素座標伸長)

横方向境界に PML (Perfectly Matched Layer) を置き、領域外へ出ていく光の
境界反射を抑えます。

```
pml = 4e-6          # 層厚 4um (この例では 16 セル)
pml = 4e-6 1e-8     # 目標往復振幅反射率 R0 を指定 (既定 1e-20)
```

- 横方向の微分を `∂/∂x -> (1/s(x)) ∂/∂x` に置き換えます。位相規約が
  `exp(-i k z)` なので `s(x) = 1 - i·sigma(x)` とすると、外向き波は層内で
  `|E| ~ exp(-|kx|·∫sigma dx)` で減衰します (屈折率を変えないため、
  垂直入射でもインピーダンス不整合による反射が生じません)
- `sigma(d) = sigma_max·(d/W)^3` (d = 層内での深さ)。`sigma_max` は
  「かすめ入射 (kx = k0·n0) での往復振幅反射率が R0」から決めます:
  `sigma_max = -(m+1)·ln(R0) / (2·k0·n0·W)`, m = 3
- 入射角 θ (kx = k0·n0·sinθ) での往復**電力**反射率は `R0^(2 sinθ)` になります
- 指定時は従来の端部吸収体 (BPM-MATLAB 由来の振幅 multiplier) を無効化します。
  **省略時は従来どおり**で、出力は変更前とビット一致します
- CPU の近軸・拡張 (広角/半ベクトル) の両経路と CUDA 版に実装されています
  (CUDA はコンパイル検証のみ、GPU 実機未検証)
- 対称境界 (`symmetry`) を使う軸の鏡像面側には PML を置きません

検証サンプル `data/sample/pml_tilt.ofd` (15° 傾けたガウシアンを側壁へ入射):

| 境界処理 | 出口電力 P/P0 | 壁なし参照との差 (= 反射) |
|---|---|---|
| 従来の端部吸収体 | 0.982 | 入射電力の **73%** が反射 |
| `pml = 4e-6` (近軸) | 7.67e-05 | **9.5e-08** |
| `pml = 4e-6` (広角パス) | — | **1.6e-09** |

R0 を変えたときの反射率は理論式 `R0^(2 sinθ)` に沿って下がり
(R0 = 1e-4 で実測 9.7e-3 / 理論 8.5e-3)、R0 = 1e-20 付近で離散化由来の
反射床 (~1e-7) に飽和します。

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
- CUDA 版 (`obpm_cuda`) も `tpa` / `powersweep` に対応しています
  (近軸パスは TPA カーネル + 電力簿記の補正、拡張パスは `bpm/wabpm.cu` の TPA カーネル。
  掃引ループ・物理スケーリング・A_eff・CSV 出力はホスト側で CPU 版と同一処理)。
  **ビルドは CUDA 12.0 で確認済みですが、GPU 実機での実行検証は未実施です**
  (数値検証済みなのは CPU 版。実機がある場合は CPU 版との一致確認を推奨)。

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
- 3 OS とも以下の解析解スモークを実行します:
  - モード解析 (`fiber_modes.ofd`): LP01 の neff が分散方程式の厳密解
    1.447167 と ±5e-4、モード整合励振の結合率 η_1 > 0.99
  - テーパ/ツイスト (`taper_twist.ofd`): 合計回転角のログと断熱的な電力保存
  - 対称境界 (`fiber_symmetry.ofd`): 半領域計算の出力パワー
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
