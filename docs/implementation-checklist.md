# 実装漏れ 対応チェックリスト

OpenBPM のコードベース調査により判明した「実装漏れ・未実装・暫定実装（プレースホルダ）・
未対応の組み合わせ」を一覧化したもの。優先度順にチェックボックスで管理する。

凡例: 優先度 **高** = 機能の正しさ・公開機能に直結 / **中** = 機能拡張・整合性 / **低** = 整備・検証

---

## 高優先度

- [x] **モードソルバ `findModes()` がプレースホルダ** ✅ 対応済み
  - 場所: `include/bpm/allset.hpp` / 新規 `bpm/modes.cpp`
  - 対応内容: 虚軸伝搬法 (imaginary-distance BPM) + Gram-Schmidt 直交化 (deflation) の
    べき乗法によるモードソルバ `wabpm_find_modes()` を `bpm/modes.cpp` に実装。
    実効屈折率は Rayleigh 商 `neff = sqrt(n0² + <E,PE>/<E,E>/k0²)` から算出。
    `allset.hpp` の `findModes()` はこれを呼び出す実装に置き換えた
    (使用時は `bpm/modes.cpp` と `bpm/wabpm.cpp` をリンク)。
  - 検証: ステップインデックスファイバの LP01/LP11 実効屈折率が分散方程式
    (Bessel) の厳密解と誤差 < 5e-5 で一致、モード間直交性 < 1e-6
    (`tests/test_modes.cpp`, `tests/test_allset.cpp`)。
  - 付随修正: `allset.hpp` に既存のコンパイルエラーとバグを発見し修正
    - 電界 `E.field` / `Mode_Struct.field` が実数型 `MatrixXf` で複素位相を保持できず、
      `tiltField()` が `complex<double> * float` の型不整合でコンパイル不能 → `MatrixXcf` 化
    - `offsetField()` がゼロ初期化した自分自身から補間 (常に 0) かつ物理座標を
      行列インデックスとして使用 → 元の界のコピーからインデックス空間で双一次補間する実装に修正

- [x] **`bend` と `wideangle` / `polarization` の併用が未対応** ✅ 対応済み
  - 場所: `sol/solve_bpm.cpp` / `cuda/solve_bpm.cu`
  - 対応内容: 拡張パス (広角 Pade(1,1) / 半ベクトル) の n2 スライス構築時に
    等価屈折率法 `n_bend = n·(1 - n²·xb·ρe/(2RoC))·exp(xb/RoC)`
    (xb = x·cosθ + y·sinθ、実部のみ、吸収 = 虚部は不変) を適用。近軸パスと同一の変換式。
    「bend is ignored」の警告は削除。CPU / CUDA 両対応。
  - 検証:
    - 単体: 一様媒質 + 曲げでビーム重心が円弧 (レイ) `<x>(z) = z²/(2RoC)` に一致
      (誤差 2%, `tests/test_wabpm.cpp` の `bend_deflection`)
    - E2E: `fiber_bend.ofd` + `wideangle = 1` で外側 (+x) への モードシフトを確認
      (コア近傍重心 +2.0um、摂動論の平衡変位 δ = k0²n̄²w⁴/(4RoC) ≈ 2.1um と整合)。
      直線 + 広角では重心 0.00um (擬似シフトなし)
    - 参考: 近軸パス (`fiber_bend.ofd`) の重心 +4um は変位振動 (周期 ≈ 407um ≈ 伝搬距離)
      の途中をサンプルしたもので、放射減衰の扱いの差によりエンジン間で位相が異なる

## 中優先度

- [x] **動画出力 `finalizeVideo()` がプレースホルダ** ✅ 対応済み (後処理集約方針で実装)
  - 場所: `bpm/model.cpp` / `sol/solve_bpm.cpp` / `cuda/solve_bpm.cu` / `tools/plot_ixz.py`
  - 対応内容: 動画出力はソルバ側では行わず後処理に集約する方針を採用し、実データで実装した。
    - ソルバ側: 入力キーワード `frames = <interval>` を追加し、`|E(x,y)|²` スナップショットを
      `/field/frames` (nframes × Ny × Nx) に記録 (CPU 近軸/拡張・CUDA 近軸/拡張の全 4 パス)。
      メタデータ `/metadata/frame_interval`, `grid_dx/dy/dz` も追加
    - 後処理側: `tools/plot_ixz.py` が伝搬マップ PNG (`/field/Ixz`)・最終電界/屈折率 PNG・
      伝搬アニメーション GIF (`/field/frames`) を生成
    - `finalizeVideo()` の誤解を招く出力 (`"Finalizing video..."`) は削除し、方針を記述した
      no-op に変更
  - 検証: `fiber_offset.ofd` + `frames = 20` で 10 フレーム (90×90) の記録・
    メタデータ・パワー保存を確認し、PNG / GIF の生成まで動作確認済み

- [x] **CUDA DFT (近傍界 3D) の CPU フォールバック経路が TODO** ✅ 確認済み (コード変更不要)
  - 場所: `cuda/dftNear3d.cu`
  - 調査結果: CPU フォールバック関数は E/H 全 6 成分とも実装済みで、ループ範囲も
    GPU カーネルと一致 (E: `imin`〜 + 成分別の上限、H: 成分軸以外 `imin-1`〜) していた。
    `// CPU TODO` は古い残骸のため、実状を表すコメントに置き換えた。

## 低優先度（検証・整備）

- [ ] **CUDA 版の GPU 実機実行検証が未実施** (ビルド検証は完了)
  - 場所: ReadMe.md（「GPU 実機での実行検証は未実施」）
  - 現状: 本対応で追加した曲げ × 広角/半ベクトル・frames 記録を含む `obpm_cuda` の
    **CUDA 12.0 でのフルビルドは検証済み**。実機 (GPU) での実行と CPU 版との数値一致
    確認のみ GPU 環境待ち。
  - 対応案: GPU 環境で `data/*.ofd` サンプルを実行し、CPU 版 `obpm` の出力と数値一致を確認。

- [x] **CMake のオブジェクトファイル絞り込み TODO** ✅ 精査済み (現状維持と判断)
  - 場所: `CMakeLists.txt`
  - 調査結果: CUDA ターゲットは既に明示リスト (`SOURCES2`) で必要ソースに絞られており、
    TODO は古い記述だった。CPU 版 `obpm` の `file(GLOB SOURCES "sol/*.c")` は全ソースが
    現状リンクされて動作しており、絞り込みはビルド時間の削減にしかならず退行リスクが
    上回るため意図的に維持する。コメントを実状を表す記述に更新した。

---

## テスト

`-DWITH_TESTS=ON` で以下の 9 テストがビルド/登録される (`ctest` で実行。
結合テストは h5py/numpy が必要)。最新の一覧と検証内容は ReadMe.md の
「テスト」節が正となる。

| テスト | 対象 | 検証内容 |
|---|---|---|
| `test_wabpm` | `bpm/wabpm.cpp` | 回折 (近軸/広角)・エネルギー保存・吸収減衰・曲げ偏向・半ベクトル整合 |
| `test_modes` | `bpm/modes.cpp` | LP01/LP11 実効屈折率 (分散方程式の厳密解と比較)・直交性 |
| `test_allset` | `include/bpm/allset.hpp` | findModes / modeSuperposition / offsetField / tiltField |
| `test_fdbpm` | `bpm/FDBPMpropagator.c` | スカラー近軸カーネルの回折・エネルギー保存・吸収・ビーム中心保持 |
| `onn_activation_solver` / `_check` | TPA + powersweep | `activation_curve.csv` を解析解 `T = 1/(1+β(P/A_eff)L)` と比較 (±7%) |
| `mode_beat_solver` / `_check` | `launch = mode <m1> <m2>` | 伝搬マップ重心のモードビート長を理論値 `λ/Δneff` と比較 (誤差 0.9%) |
| `wlsweep_check` | `wlsweep = 1` | 掃引結果が最終波長の単独実行と完全一致することを確認 |

## 補足: 確認したが「対応不要」と判断した項目

- `sol/`, `mpi/`, `cuda/`, `post/` に多数ある `fprintf(stderr, ...)` は正常なエラーハンドリング
  （ファイルオープン失敗・HDF5 書き込み失敗等）であり、実装漏れではない。
- `setNCladding()` / `setShapes()` / `setDisplayScaling()`（`bpm/model.cpp`）は
  旧 API を意図的に無効化する `throw` 実装であり、仕様どおり。
