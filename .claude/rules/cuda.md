---
paths:
  - "cuda/**"
  - "cuda_mpi/**"
  - "bpm/*.cu"
  - "src/cuda_Main.cu"
---

# CUDA 開発の規約

- **この開発環境に GPU 実機はない**。CUDA 変更の検証は
  `cmake -S . -B build-cuda -DWITH_CUDA=ON` でのコンパイル (CUDA 12.0) まで。
  実行検証が未実施であることをコミットメッセージとチェックリストに明記する。
- **CPU/CUDA パリティ**: CPU 版 (`sol/solve_bpm.cpp`) に機能を追加したら、
  CUDA 版 (`cuda/solve_bpm.cu`) にも同じ変更を反映する。すぐに実装できない場合は
  未対応の**実行時警告**を必ず入れる (サイレント無視は禁止)。
  現状の既知の未対応: `tpa` / `powersweep` (警告実装済み)。
- グローバル変数 (`NTpaB`, `PowerSweep`, `BPM` 等) は `obpm.h` の EXTERN パターン。
  実体は `#define MAIN` する `src/cuda_Main.cu` 側で確保される。
- 拡張パス (広角/半ベクトル) はデバイス常駐コンテキスト `wabpm_gpu_*` を使う。
  ホストへの全電界コピー (`wabpm_gpu_get_field`) は高コストなのでフレーム記録など
  必要時のみに限定する。
- CUDA エラーは `CUDA_CHECK()` マクロで検査する。
