# OpenBPM

OpenFDTD のコード構造（入力データ・メッシュ・材料・後処理）を元にした
ビーム伝搬法 (BPM) ソルバーです。伝搬カーネルは BPM-MATLAB の
FDBPM (Douglas-Gunn ADI 法) を移植しています。

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

- `polarization` / `wideangle` 指定時は倍精度の一般化伝搬エンジンを使用します
  (CPU 版 `bpm/wabpm.cpp` / CUDA 版 `bpm/wabpm.cu`、同一アルゴリズム)。
  `bend` との併用も可能です (等価屈折率法を屈折率分布に反映)。

- `feed` の電圧はビーム振幅として使用されます。
- メッシュは等間隔を推奨します (不均一の場合は平均セル幅で計算し警告を表示)。

# Reference
OpenFDTD
http://www.e-em.co.jp/OpenFDTD/

BPM-MATLAB
https://gitlab.gbar.dtu.dk/biophotonics/Programs/BPM-Matlab

BPM
https://github.com/warthan07/Nemaktis

# 屈折率資料
https://github.com/polyanskiy/refractiveindex.info-database
