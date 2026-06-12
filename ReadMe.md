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
  - `obpm_cuda` は CPU 版と同じ BPM 機能 (beam/refindex/導電率吸収/メッシュ警告) に対応しています。
    屈折率分布は 3D 配列として一括でデバイスへ転送されます。
  - ビルドは CUDA 12.0 で確認済みです (GPU 実機での実行検証は未実施)。
- MPI 版: `-DWITH_MPI=ON` (要 MPI)

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

- `polarization` / `wideangle` 指定時は倍精度の一般化伝搬エンジン (`bpm/wabpm.cpp`) を使用します
  (CPU 版のみ。CUDA 版は未対応でエラー終了します。`bend` との併用は不可)。

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
