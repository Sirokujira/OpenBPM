---
paths:
  - "sol/**"
  - "bpm/**"
  - "post/**"
  - "src/**"
  - "include/**"
  - "CMakeLists.txt"
---

# 移植性の絶対規則 (Windows/macOS CI で実際に踏んだもの)

- **C99 VLA 禁止** (MSVC 非対応)。`malloc` + 明示インデックスで書く。
- 複素数は `CREALF` / `CIMAGF` マクロ経由でアクセスする
  (MSVC では C ソースが C++ / `std::complex<float>` としてビルドされる)。
  `complex × double` の直接乗算を書かない (`complex × float` は可)。
- libm のリンクは `MATH_LIB` 変数経由。MSVC 固有フラグは既存の CMake ブロックに従う。
- 暗黙の関数宣言を残さない (macOS AppleClang はエラー扱い)。
  ナローイング変換 (`{}` 初期化子内の int→size_t 等) にも注意。
- Eigen のバージョン制約・OpenMP の検出方法はプラットフォームで異なるため、
  `CMakeLists.txt` の既存の分岐を壊さない。
- **macOS の AppleClang は OpenMP を同梱しない**。素の `find_package(OpenMP)` は
  必ず失敗するため、`CMakeLists.txt` が `brew --prefix libomp` から
  `OpenMP_C_FLAGS` 等を自動設定している (`if(APPLE ...)` ブロック)。
  このブロックを消すと ReadMe が案内する `cmake -S . -B build` が macOS で
  通らなくなる。CI の build-macos は `-DOpenMP_*` を渡さずに configure して
  この自動検出を検証している。
- 変更後は 3 OS の CI (Linux / macOS / Windows) が通ることを確認する。
