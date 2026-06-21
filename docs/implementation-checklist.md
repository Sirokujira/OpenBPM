# 実装漏れ 対応チェックリスト

OpenBPM のコードベース調査により判明した「実装漏れ・未実装・暫定実装（プレースホルダ）・
未対応の組み合わせ」を一覧化したもの。優先度順にチェックボックスで管理する。

凡例: 優先度 **高** = 機能の正しさ・公開機能に直結 / **中** = 機能拡張・整合性 / **低** = 整備・検証

---

## 高優先度

- [ ] **モードソルバ `findModes()` がプレースホルダ**
  - 場所: `include/bpm/allset.hpp:166`
  - 現状: `MatrixXf::Random(...)` で乱数フィールドを返し、`neff = 1.5 + 0.01*i` の固定ダミー値を設定するだけ。
    実際の固有モード計算（虚軸伝搬法・固有値解法など）が未実装。
  - 影響: `modeSuperposition()`(`allset.hpp:146`) / `getLabeledModeIndex()`(`allset.hpp:185`) も
    実モードに依存するため、モード関連機能は全て使用不能。
  - 対応案: 虚軸伝搬法 (imaginary-distance BPM) もしくは有限差分固有値解法で実モードを算出。
    BPM-MATLAB の `findModes` を参考に neff / フィールド / 損失を正しく返す。
  - 注意: 現在 `findModes` 等は `solve_bpm` の実行経路には組み込まれておらず（呼び出しは
    `allset.hpp` のコメントアウトされた例のみ）、`obpm` 実行ファイルにもリンクされていない。
    公開機能として復活させる場合はビルド統合も併せて必要。

- [ ] **`bend` と `wideangle` / `polarization` の併用が未対応**
  - 場所: `sol/solve_bpm.cpp:245` 付近 / `cuda/solve_bpm.cu:277` / ReadMe.md:71
  - 現状: 拡張パス（広角 Pade(1,1) / 半ベクトル）では曲げ (RoC) を無視し警告を出すのみ
    （`*** warning : bend is ignored in wideangle/polarization mode.`）。
  - 対応案: 一般化 ADI エンジン `bpm/wabpm.{cpp,cu}` の複素比誘電率 `n2` に等価屈折率法
    （`n_eq = n*exp(x/RoC)`）を反映させ、曲げと広角/半ベクトルを共存させる。

## 中優先度

- [ ] **動画出力 `finalizeVideo()` がプレースホルダ**
  - 場所: `bpm/model.cpp:81`
  - 現状: `"Finalizing video..."` を出力するだけで、動画ハンドルのクローズ・書き出し処理が無い。
  - 影響: `bpm/model.cpp` / `include/bpm/model.hpp` の `BPMModel` 群は現状 `obpm` の
    ビルド対象に明示列挙されておらず未リンク（`CMakeLists.txt:408`）。可視化機能として
    使うなら実装とビルド統合の双方が必要。
  - 対応案: 伝搬マップ `Ixz` を用いたフレーム書き出し（PNG 連番 / 動画）を実装、もしくは
    可視化を後処理ツール `obpm_post` 側へ集約する方針を決める。

- [ ] **CUDA DFT (近傍界 3D) の CPU フォールバック経路が TODO**
  - 場所: `cuda/dftNear3d.cu:325`（`// CPU TODO`）
  - 現状: GPU 経路は実装済みだが、CPU フォールバック分岐が暫定。OpenFDTD 由来の後処理部。
  - 対応案: CPU 版 `dft_near3d*_cpu` 呼び出しの整合性確認、または使用しないなら分岐を整理。

## 低優先度（検証・整備）

- [ ] **CUDA 版の GPU 実機実行検証が未実施**
  - 場所: ReadMe.md:22（「GPU 実機での実行検証は未実施」）
  - 現状: ビルドは CUDA 12.0 で確認済みだが、`obpm_cuda` の実機実行結果が CPU 版と一致するか未検証。
  - 対応案: GPU 環境で `data/*.ofd` サンプルを実行し、CPU 版 `obpm` の出力（`time_series_data.h5`）と
    数値一致を確認（特に `wideangle`/`polarization` パス `bpm/wabpm.cu`）。

- [ ] **CMake のオブジェクトファイル絞り込み TODO**
  - 場所: `CMakeLists.txt:121`（`#TODO: Makefile の記述を元に必要なオブジェクトファイルのみに絞る?`）
  - 現状: 不要オブジェクトを含めてビルドしている可能性。
  - 対応案: 各ターゲット（`obpm` / `obpm_cuda` / `obpm_post`）に必要なソースのみへ整理。

---

## 補足: 確認したが「対応不要」と判断した項目

- `sol/`, `mpi/`, `cuda/`, `post/` に多数ある `fprintf(stderr, ...)` は正常なエラーハンドリング
  （ファイルオープン失敗・HDF5 書き込み失敗等）であり、実装漏れではない。
- `setNCladding()` / `setShapes()` / `setDisplayScaling()`（`bpm/model.cpp`）は
  旧 API を意図的に無効化する `throw` 実装であり、仕様どおり。
