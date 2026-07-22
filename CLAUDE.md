# OpenBPM

ビーム伝搬法 (BPM) 光導波路ソルバー (C/C++)。OpenFDTD-X (GUI) から
QProcess で起動される処理カーネル。近軸 ADI と広角 (wabpm) の 2 経路、
CPU / MPI / CUDA 実装を持つ。

## ビルド / テスト

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=OFF
cmake --build build -j"$(nproc)"

# 回帰 (fiber): output power = 3.122518e+02 が不変、
# obpm_post 後の bpm_ixz.csv が変更前と md5 一致であること
mkdir -p /tmp/smoke && cp data/sample/fiber.ofd /tmp/smoke/ && cd /tmp/smoke
$OLDPWD/bin/obpm fiber.ofd && $OLDPWD/bin/obpm_post fiber.ofd

# ONN 活性化 (TPA + powersweep): activation_curve.csv が単調非増加で
# 解析解 T=1/(1+β(P/A_eff)L) と ±7% (A_eff・L はログ出力を使う)
cp data/sample/onn_activation.ofd /tmp/smoke/ && cd /tmp/smoke
$OLDPWD/bin/obpm onn_activation.ofd
```

## 移植性の絶対規則 (Windows CI で実際に踏んだもの)

- **C99 VLA 禁止** (MSVC)。`malloc` + 明示インデックスで書く。
- 複素数は `CREALF` / `CIMAGF` マクロ経由でアクセスする
  (MSVC では C ソースが C++ / `std::complex<float>` としてビルドされる)。
  `complex × double` の直接乗算を書かない (`complex × float` は可)。
- libm は `MATH_LIB` 変数経由。MSVC フラグは既存 CMake ブロックに従う。

## 機能追加の規則

- 入力キーは `sol/input_data.c` に追加し、省略時は従来動作とバイト一致。
- 物理スケーリング規約: `tpa`/`powersweep` 使用時のみ場を
  ∫∫|E|²dA = P_in [W] に正規化 (|E|² = 強度 I)。未使用時は従来の
  無次元場のまま (fiber 回帰を壊さない)。
- CPU の近軸/広角の**両経路**に同じ物理を入れる (片方だけの実装は不可)。
- CUDA 版 (`obpm_cuda`) の対応状況を ReadMe.md に明記する。
  現状: `tpa` / `powersweep` は CPU のみ、CUDA 未対応 (無視される)。
- 新機能には data/sample/ の検証ケース + CI スモーク (3 OS) を付ける。

## CI

`.github/workflows/ci.yml`: Linux / macOS / Windows (MSVC + Ninja +
vcpkg `hdf5[core,zlib]:x64-windows-static-md`)。タグ `v*` で Release 添付。
