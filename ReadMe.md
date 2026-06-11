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
- MPI 版: `-DWITH_MPI=ON` (要 MPI)

## 実行

入力データは OpenFDTD 形式です。サンプル: `data/waveguide.ofd`
(ステップインデックス導波路, 波長 1.55um)

```sh
./bin/obpm data/waveguide.ofd
```

- 周波数 (`frequency2`) から波長を決定し、材料 (`material`/`geometry`) の
  比誘電率から屈折率分布を構築して z 軸方向にビームを伝搬します。
- 結果は `time_series_data.h5` (HDF5) に出力されます
  (`/field/Efinal_r`, `/field/Efinal_i` : 最終電界、`/field/n_out_r` : 屈折率分布)。

# Reference
OpenFDTD
http://www.e-em.co.jp/OpenFDTD/

BPM-MATLAB
https://gitlab.gbar.dtu.dk/biophotonics/Programs/BPM-Matlab

BPM
https://github.com/warthan07/Nemaktis

# 屈折率資料
https://github.com/polyanskiy/refractiveindex.info-database
