---
paths:
  - "cuda/**"
  - "cuda_mpi/**"
  - "bpm/*.cu"
  - "src/cuda_Main.cu"
---

# CUDA 開発の規約

- **CI の開発環境に GPU 実機はない**。CUDA 変更の検証は
  `cmake -S . -B build-cuda -DWITH_CUDA=ON` でのコンパイル (CUDA 12.0) まで。
  実行検証が未実施であることをコミットメッセージとチェックリストに明記する。
  (2026-08-21 に Windows + CUDA 13.1 + RTX 3060 で `obpm_cuda` の実行を確認。
  nvcc 13 は sm_60/70 を受けないので CMake は nvcc の版で既定 arch を切り替える。
  HDM (既定) では device メモリを host から読めない — OpenFDTD で HDF5
  スナップショットが実際に落ちた。`-cpu` / UM で出ないバグがある。)
- **CPU/CUDA パリティ**: CPU 版 (`sol/solve_bpm.cpp`) に機能を追加したら、
  CUDA 版 (`cuda/solve_bpm.cu`) にも同じ変更を反映する。すぐに実装できない場合は
  未対応の**実行時警告**を必ず入れる (サイレント無視は禁止)。
  現状の既知の未対応: `wlsweep` (警告実装済み)。
- **逆に「未対応」の警告を残したままにしない**: 実装したのに「無視される」と
  警告し続けると、正しい結果まで捨てられてしまう (過去に `tpa` / `powersweep` で
  発生)。実装済みだが実機未検証の場合は「無視される」ではなく
  「実機未検証」である旨の note にする。
- グローバル変数 (`NTpaB`, `PowerSweep`, `BPM` 等) は `obpm.h` の EXTERN パターン。
  実体は `#define MAIN` する `src/cuda_Main.cu` 側で確保される。
- 拡張パス (広角/半ベクトル) はデバイス常駐コンテキスト `wabpm_gpu_*` を使う。
  ホストへの全電界コピー (`wabpm_gpu_get_field`) は高コストなのでフレーム記録など
  必要時のみに限定する。
- CUDA エラーは `CUDA_CHECK()` マクロで検査する。
