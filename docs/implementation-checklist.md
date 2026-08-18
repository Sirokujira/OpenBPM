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

# 第 2 回監査 (main = c9a7d94 時点)

TPA 非線形・パワー掃引 (PR #8)、obpm_post の BPM 可視化 (PR #6)、
CI (Linux/macOS/Windows, PR #3/#5/#7) が追加された現行 main に対する再監査。

## 高優先度

- [x] **CUDA 版が `tpa` / `powersweep` 未対応** ✅ 対応済み (⚠ 実機検証のみ未実施)
  - 場所: `cuda/solve_bpm.cu` / `bpm/FDBPMpropagator.cu` / `bpm/wabpm.cu`
  - 対応内容:
    1. 実行時警告 (サイレント無視の解消) — 先行対応済み。本対応で警告は不要になり削除
    2. TPA カーネルと掃引ループを実装し、CPU 版と同等の機能に到達:
       - 近軸: `applyTPA` (要素毎の `E *= exp(-β/2·|E|²·Δz)` + 電力和の二重精度
         atomicAdd による集計) と `scalePrecisePowerByTPA` (電力簿記 `precisePower` を
         同率で補正し、次ステップの `fieldCorrection` が TPA 減衰を打ち消すのを防ぐ)
       - 拡張 (広角/半ベクトル): `bpm/wabpm.cu` に `k_tpa` / `wabpm_gpu_tpa` を追加
       - ホスト側: 掃引ループ・物理スケーリング・A_eff 算出・`activation_curve.csv` 出力を
         CPU 版から移植 (掃引毎に `createDeviceStructs`/`retrieveAndFreeDeviceStructs` を
         対で呼び、初期界を再スケールする)
  - 検証状況 (**この開発環境に GPU がないため実行検証は未実施**):
    - CUDA 12.0 でのコンパイル検証 ✅
    - 移植したホスト側ロジック (掃引点生成 / A_eff / P_out 記録 / CSV 出力) が
      CPU 版と**文字列レベルで一致**することを機械的に照合 ✅
    - TPA 式の精度を CPU と揃えた (倍精度で `exp` を評価してから float 化) ✅
    - `applyTPA` は 1 次元グリッドストライドのため 1D ブロックで起動する
      (2D ブロックだと `threadIdx.y` 方向のスレッドが同一要素を重複更新する) ✅
    - ⚠ **未実施**: GPU 実機での実行と CPU 版との数値一致確認。実機入手時に
      `data/sample/onn_activation.ofd` を両版で実行し `activation_curve.csv` を
      比較すること (低優先度項目「GPU 実機実行検証」と併せて対応)

- [x] **CI が単体テストを実行していない** ✅ 対応済み
  - 場所: `.github/workflows/ci.yml`
  - 対応内容: Linux / macOS ジョブの configure に `-DWITH_TESTS=ON` を追加し、
    ビルド直後に `ctest --output-on-failure` (test_wabpm / test_modes / test_allset の
    解析解検証) を実行するステップを追加した。

## 中優先度

- [x] **MPI 版は BPM 未対応 (FDTD ソルバのみ)** ✅ 「FDTD のみ」と明記で解決
  - 場所: `src/mpi_Main.c` (`solve()` のみ呼び出し、`solve_bpm` への分岐なし)、`mpi/`
  - 判断: BPM の ADI は z 方向に逐次で、MPI 領域分割には転置通信が必要となり
    費用対効果が低い。実装せず ReadMe に「MPI 版は FDTD のみ、BPM の並列化は
    OpenMP (CPU 版) / CUDA 版を使用」と理由付きで明記した (第 6 回後の対応)。
    将来需要が生じたら再検討する。

- [x] **TPA の解析解との数値回帰テストがない** ✅ 対応済み
  - 場所: `.github/workflows/ci.yml` / `tools/check_activation.sh` (新規)
  - 対応内容: `activation_curve.csv` の全掃引点を平面波近似の解析解
    `T = 1/(1 + β·(P_in/A_eff)·L)` と比較する検証スクリプトを追加し (A_eff / L は
    obpm.log から自動取得、許容 ±8%)、Linux / macOS の ONN スモークを単調飽和判定から
    このスクリプト呼び出しに置き換えた。ローカル実行で全 8 点一致 (最大誤差 4.1%) を確認。

## 低優先度

- [x] **obpm_post が `/field/frames` 未対応** ✅ 対応済み
  - 場所: `post/postbpm.c`
  - 対応内容: 3 次元データセット読み込み (`read3d`) とスカラーメタデータ読み込み
    (`readScalar`) を追加し、`/field/frames` の各フレームを 1 ページずつ等高線出力する。
    `/metadata/frame_interval` と `grid_dz` から各ページのタイトルに z 位置を表示する。
  - 検証: `fiber.ofd` + `frames = 40` で 5 フレーム (90x90) のページ生成を確認。
    `frames` 無指定時はデータセットが無いためページは増えず、従来出力と同一。

- [x] **半ベクトル (pol=x/y) 指定時のモードソルバの直交化が近似**
      ✅ 実測により「随伴 deflation の実装は不要」と判断 (別のバグを発見・修正)
  - 場所: `bpm/modes.cpp`
  - 調査 (高コントラスト矩形コア n=2.0/1.0、スカラーと半ベクトルで比較):
    - deflation 自体は機能している : 既出モードとの重なりは約 1e-16 (機械精度)
    - 固有ベクトル残差 `||P φ − μ φ|| / ||P φ||` は 4e-2〜1.9e-1。ただし
      **対称演算子であるスカラー (標準 deflation が厳密に正しい) でも同程度**
      (6e-2〜1.9e-1) であり、半ベクトル側がむしろ小さい。
      → 残差の支配要因は随伴性ではなく**べき乗法の収束率**である
    - 収束判定を残差ベース (1e-4) に変えると 60000 反復でも収束せず、
      精度向上には固有値解法自体の変更 (シフト反転・Arnoldi 等) が必要と判明
    - 主要な出力である neff の精度は良好 (LP01 で厳密解と 2e-7、LP11 で 3e-5)
  - 結論: 随伴系 deflation を実装しても測定可能な改善は得られないため実装しない。
    将来精度を上げる場合は deflation の内積ではなく固有値解法を見直すこと。
  - **副産物 (重要なバグ修正)**: 調査中に、導波フィルタ導入時に入れた収束判定の
    不具合を発見・修正した。非導波状態の Rayleigh 商を neff = 0 にクランプしていたため、
    2 回連続で 0 になった時点で「収束」と誤判定し、**モード 2 以降の探索が即座に
    打ち切られていた** (高コントラスト構造で 4 モード中 1 モードしか返らない)。
    判定をクランプしない生の `v = n0^2 + mu/k0^2` で行うよう修正し、
    再発防止の回帰テスト (`test_modes.cpp` の `multimode`) を追加した。

---

# 第 3 回監査 (機能 × 実装パスの網羅精査)

第 2 回監査の対応 (CI テスト組込・TPA 解析解検証・CUDA 警告) 後の現行コードに対し、
カーネルが持つ機能と入力キーワード・実装パスの対応関係を精査した追加の洗い出し。

## 中優先度

- [x] **モードソルバがソルバ入力から利用できない (ライブラリのみ)** ✅ 対応済み
  - 場所: `bpm/modes.cpp` / `sol/input_data.c` / `sol/solve_bpm.cpp` / `cuda/solve_bpm.cu`
  - 対応内容: 入力キーワード `modes = <nModes> [excite]` を追加 (3 点セット +
    CUDA パリティ)。入力断面のモードを虚軸伝搬法で解析し、neff を obpm.log、
    モード形状を HDF5 `/modes/mode<i>` + `/modes/neff` に出力。`excite` 指定で
    入射ビームを基本モードに置き換える (モード整合励振、beamtilt 位相ランプは維持)。
    CUDA 版もホスト側で同一処理 (bpm/modes.cpp + wabpm.cpp を obpm_cuda にリンク)。
    あわせて `wabpm_find_modes` に導波条件フィルタ (n0 < neff、非導波状態の早期打切)
    を追加し、シングルモードファイバで `modes = 3` 指定時の無意味な解と
    無駄な反復 (25 秒 → 2.6 秒) を解消。
  - 検証: `data/sample/fiber_modes.ofd` (新規) で LP01 neff = 1.4471672 が
    分散方程式の厳密解 1.447167 と一致 (誤差 2e-7)。モード整合励振で電力保存
    T = 0.99999 (ガウシアン 0.99988)。既存回帰 (fiber output power = 3.122518e+02、
    ctest 3 スイート) は不変。CI スモーク (Linux/macOS, neff ±5e-4 判定) を追加。
    CUDA 12.0 コンパイル検証済み (実機実行は GPU 待ち)。

- [x] **テーパ / ツイスト機能が未公開 (カーネルは対応済み)** ✅ 対応済み
  - 場所: `include/obpm.h` / `sol/input_data.c` / `sol/solve_bpm.cpp` / `cuda/solve_bpm.cu`
  - 対応内容: 入力キーワード `taper = <ratio>` (出口/入口の横方向スケール比) と
    `twist = <rate[deg/m]>` を追加 (3 点セット)。近軸パスはカーネルの座標変換
    (`taperPerStep`/`twistPerStep`) を有効化し、変換時は入力断面を渡すよう
    n_in の充填を `iz_start` 固定に切り替え。**拡張パス (広角/半ベクトル) にも
    同一の変換 (相似縮小 + 回転 + 双一次補間) を n2 構築時に実装**し、当初「近軸限定」と
    見込んでいた制約を解消した。CUDA 版も近軸 (カーネル既存 branch が先頭スライスを
    参照) / 拡張 (CPU と同一式) の両方に対応。
  - 検証 (いずれも解析的に厳密な性質と比較):
    - 一様媒質では変換を掛けても結果不変 : 拡張パスは倍精度で**厳密に 0**、
      近軸パスは float32 双一次補間の丸めで 3.7e-6 (60 ステップ蓄積、float32 eps 相当)
    - 幾何の直接検証 (出力 RIP `/field/n_out_r`) : 矩形コアの 90 度ツイストで
      主軸 89.6 度 (理論 90)、円形コアの `taper=0.5` で等価コア半径 2.08um (理論 2.05)
    - 円形コア (回転対称) のツイスト不変性 : 階段近似した円の離散化で有限差が残るが、
      格子細分で単調減少 (dx=1/3um: 7.2e-2 → 1/6um: 4.0e-2 → 1/9um: 3.4e-2)。
      全パワー差 3e-6、neff 差 9e-6 で、形状のみの離散化誤差と確認
    - 断熱性 : 4.1→2.0um テーパ (600um) で P/P0 = 0.9998、V 低下に伴うモード拡大も確認
    - 既存回帰不変 (fiber output power = 3.122518e+02、ctest 3 スイート)、
      CUDA 12.0 コンパイル検証済み
  - 付随: 検証サンプル `data/sample/taper_twist.ofd` と CI スモーク (3 OS) を追加。

## 低優先度

- [x] **対称境界 (xSymmetry / ySymmetry) が未公開** ✅ 対応済み
  - 場所: `include/obpm.h` / `sol/input_data.c` / `sol/solve_bpm.cpp` / `cuda/solve_bpm.cu` /
    `bpm/wabpm.cpp` / `bpm/wabpm.cu` / `include/bpm/wabpm.h`
  - 対応内容: 入力キーワード `symmetry = <x|y|xy> [sym|anti]` を追加 (3 点セット)。
    指定軸のメッシュ始端を鏡像面とし、計算領域を 1/2〜1/4 に削減する。
    - 近軸パス: カーネル既存の対称境界を有効化 (`P->ySymmetry` = x 鏡像、
      `P->xSymmetry` = y 鏡像。カーネルの命名は「対称面の軸」基準)
    - **拡張パス (広角/半ベクトル) には対称境界が無かったため新規実装**:
      `wabpm_params` に `symx`/`symy` を追加し、陽的スイープ・Thomas 法・
      `wabpm_apply_P` の index 0 で鏡像セル (E[-1] = ±E[0]) を対角へ畳み込む
      (CPU/CUDA 両方)。これで規約「CPU 両経路 + CUDA パリティ」を満たす
    - ソルバ設定: ビーム中心の既定値を鏡像面に、端部吸収体を鏡像面側で無効化
  - 検証 (半領域/1/4 領域 vs 全領域の該当部分):
    - 一様媒質 (厳密に対称な構造) : 近軸は x/y/xy いずれも 2e-6〜3e-6
      (float32 の丸め相当)、**拡張パスは倍精度で厳密に 0**
    - ファイバ (x 対称) : 3.5e-6 で一致し、出力パワーは正確に 1/2 (0.499999)
    - 既存回帰不変 (fiber = 3.122518e+02、ctest 3 スイート、ONN、taper/twist)
    - CUDA 12.0 コンパイル検証済み (実機実行は GPU 待ち)
  - **調査で判明した既存の性質 (重要)**: OpenFDTD 由来の形状ラスタライズは
    y 方向に半セルの非対称性を生じる (円形コアでも上下反転と 6.4e-3 の差、
    矩形コアでも同様)。このため y 鏡像面を使うと「より対称な別構造」を
    解くことになり全領域計算と 6〜8% ずれる。対称境界の実装バグではなく
    入力構造側の性質のため、ReadMe とサンプルに注意として明記した。
  - 付随: 検証サンプル `data/sample/fiber_symmetry.ofd` と CI スモーク (3 OS) を追加。

---

# 第 4 回監査 (main = 5d63eda 時点)

前回監査以降の main の変更は PR #9 (Claude Code 設定の追加) のみで、
ソルバ・テスト・CI のコード変更はなし。マーカー走査 (TODO/FIXME/未実装/
placeholder 等) も 0 件で、**新規の未実装項目はなし**。

- 対応: PR #9 の設定 (CLAUDE.md / .claude/settings.json) と本ブランチの設定が
  重複していたため統合した。CLAUDE.md は両者の内容を統合 (移植性規則・fiber 回帰
  基準値・バイト一致の後方互換規約を取り込み)、移植性規則は
  `.claude/rules/portability.md` として独立させた。settings.json は許可リストの
  和集合 + force-push 拒否 + `$schema` を採用。

## 未対応項目の一覧 (第 4 回時点の残件サマリ)

| # | 項目 | 優先度 | 規模 | 備考 |
|---|---|---|---|---|
| 1 | CUDA 版の TPA/掃引フル実装 (警告は対応済み) | 中 | 大 | GPU 実機がなく検証はコンパイルまで。要方針判断 |
| 2 | ~~モードソルバの入力キーワード接続~~ | - | - | ✅ 対応済み (第 3 回監査の項を参照) |
| 3 | テーパ/ツイストの公開 (カーネル対応済み) | 中 | 中 | 近軸パス限定。解析解 (断熱定理) で検証可能 |
| 4 | MPI 版の BPM 対応 or「FDTD のみ」と明記 | 中 | 大/小 | 明記のみなら ReadMe 1 行 |
| 5 | GPU 実機での実行検証 | 低 | - | 環境要因 (ビルド検証は完了) |
| 6 | obpm_post の `/field/frames` 対応 | 低 | 小 | 現状 tools/plot_ixz.py で充足 |
| 7 | 半ベクトルモードの随伴直交化 | 低 | 中 | スカラー使用では問題なし |
| 8 | 対称境界 (x/ySymmetry) の公開 | 低 | 小 | 計算量最適化 |
| 9 | Windows CI のテストカバレッジ拡充 | 低 | 小 | ctest + ONN 解析解判定を追加 |

---

# 第 5 回監査 (modes キーワード対応後)

`modes = <nModes> [excite]` 対応 (第 3 回監査 #2 の解消) 後の再監査。
main に新規変更なし、マーカー走査 0 件。E2E 検証 (ONN/回折/吸収/曲げ/モード) は
全項目が理論値と一致 (最大誤差: ONN 深飽和点 4.1%、他は 0.5% 以下)。

## 新規: 低優先度

- [x] **モード形状 (`/modes`) の可視化が未対応** ✅ 対応済み (第 6 回後)
  - 場所: `tools/plot_ixz.py` / `post/postbpm.c` (いずれも `/modes` を読まない)
  - 現状: `modes` キーワードで HDF5 に出力される `/modes/mode<i>` と `/modes/neff` を
    可視化する経路がない (HDF5 規約「新規データセット追加時は plot_ixz.py と
    postbpm.c の対応も検討」への追随漏れ)。数値 (neff) は obpm.log で確認可能。
  - 対応案: `tools/plot_ixz.py` に `/modes` があればモード形状のページを追加する
    (neff をタイトルに表示)。obpm_post 側は需要が生じた時点で対応。

## 未対応項目の一覧 (第 5 回時点の残件サマリ = 9 件)

| # | 項目 | 優先度 | 規模 |
|---|---|---|---|
| 1 | ~~CUDA 版の TPA/掃引フル実装~~ | - | - | ✅ 対応済み (実機検証のみ #4 に集約) |
| 2 | ~~テーパ/ツイストの公開~~ | - | - | ✅ 対応済み |
| 3 | ~~MPI 版の BPM 対応 or 明記~~ | - | - | ✅ 「FDTD のみ」と明記で解決 |
| 4 | GPU 実機での実行検証 (TPA/掃引・テーパ/ツイスト含む) | 低 | - (環境要因) |
| 5 | obpm_post の `/field/frames` 対応 | 低 | 小 |
| 6 | ~~モード形状 (`/modes`) の可視化~~ | - | - | ✅ 対応済み |
| 7 | 半ベクトルモードの随伴直交化 | 低 | 中 |
| 8 | 対称境界 (x/ySymmetry) の公開 | 低 | 小 |
| 9 | ~~Windows CI のテストカバレッジ拡充~~ | - | - | ✅ 対応済み |

---

# 第 6 回監査 (AGENTS.md 追加後)

前回 (第 5 回) 以降の変更はドキュメントのみ (AGENTS.md 新規 + CLAUDE.md 相互参照、
コード変更なし)。main にも新規コミットなし。マーカー走査はコード 0 件。
**新規の未実装項目はなし。残件は第 5 回のサマリ (9 件) から変化なし。**

## 第 6 回後のクリーンアップ対応 (3 件)

- [x] **モード形状 (`/modes`) の可視化** ✅ 対応済み (第 5 回 #6)
  - `tools/plot_ixz.py` に `/modes` 描画を追加 (`*_modes.png`)。符号付き実数場を
    発散配色 (0 中心、最大振幅点が正になるよう正規化) でグリッド表示し、
    neff をタイトルに表示。
  - 検証: SMF (LP01 1 モード) と V≈5 数モードファイバ (LP01/LP11a/LP11b/LP21 の
    4 モード) で描画を確認。モード形状・縮退ペアの直交方位とも物理的に正しい。

- [x] **MPI 版の BPM 方針** ✅ 「FDTD のみ」と明記で解決 (第 5 回 #3)
  - 方針: BPM の ADI は z 逐次で MPI 領域分割は転置通信が必要になり費用対効果が
    低いため、実装せず ReadMe に「MPI 版は FDTD のみ、BPM は OpenMP/CUDA を使用」と
    明記した。将来需要が生じたら再検討。

- [x] **Windows CI のテストカバレッジ拡充** ✅ 対応済み (第 5 回 #9)
  - Windows ジョブに `-DWITH_TESTS=ON` + `ctest` を追加し、Git Bash ステップで
    ONN 解析解判定 (`tools/check_activation.sh`) とモード解析判定 (neff ±5e-4) を
    Linux/macOS と同一のスクリプト・基準で実行するようにした。
  - 注意: Windows ランナーでの実挙動は CI 実行で確認する (ローカルは Linux のため
    YAML 構文検証と判定ロジックの Linux 動作確認まで)。

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

---

# 第 7 回: GUI 表示用の時系列 (z 伝搬) データ出力

OpenFDTD-X (GUI) 側で伝搬の様子をグラフ表示できるよう、HDF5 に z 伝搬に沿った
データを追加した (キーワード不要・常に出力)。

- [x] **z 伝搬に沿った時系列データが不足** ✅ 対応済み
  - 場所: `sol/solve_bpm.cpp` / `cuda/solve_bpm.cu` / `bpm/FDBPMpropagator.cu` /
    `bpm/wabpm.cu` / `tools/plot_ixz.py`
  - 追加内容:
    - `/field/Iyz` (Nz x Ny) : 中心列の伝搬マップ (既存 `/field/Ixz` の y 版)
    - `/trace/*` (長さ = 伝搬ステップ数の 1 次元配列) :
      `z` (断面位置 [m])、`power` (断面総パワー)、`peak` (|E|^2 最大値)、
      `centroid_x`/`centroid_y` (強度重心 [m])、`width_x`/`width_y` (2σ 幅 [m])
    - 座標軸は既存の `/metadata/Xc,Yc,Zc`(セル中心) と `Xn,Yn,Zn`(節点) を利用可能
    - `tools/plot_ixz.py` に `*_trace.png` (4 段グラフ) を追加
  - 実装: CPU は `computeTrace` テンプレートで近軸/拡張の両パス共通。
    CUDA は `fieldTrace` (近軸) と `k_trace` (拡張, `wabpm_gpu_trace`) の
    デバイス集計カーネルで求め、ホストへは 6 要素のみ転送する
    (毎ステップの全電界転送を避けるため)。集計値から統計量を出す式
    (`traceFromAcc`) は CPU 版 `computeTrace` と同一。
  - 検証: 自由空間ガウシアン回折 (`freespace.ofd`) で
    `width_x` が解析解 `w(z) = w0*sqrt(1+(z/zR)^2)` と相対誤差 0.5% 以内で一致、
    `power` はほぼ保存、`centroid` は約 0 (対称ビーム)、`peak` は回折で単調減少。
    対称ビームで `Ixz` と `Iyz` が一致 (差 1.1e-6)。拡張パス (広角) でも同一の
    データセットが出力されることを確認。既存回帰はすべて不変
    (fiber 3.122518e+02、ctest 3 スイート、ONN、symmetry、taper/twist)。
    CUDA 12.0 コンパイル検証済み (実機実行は GPU 待ち)。

- [x] **位相表示用の複素電界スナップショットが無い** ✅ 対応済み (第 7 回の続き)
  - 場所: `include/obpm.h` / `sol/input_data.c` / `sol/solve_bpm.cpp` /
    `cuda/solve_bpm.cu` / `tools/plot_ixz.py`
  - 現状 (対応前): `/field/frames` は強度のみで、GUI で位相分布を表示できなかった
    (複素電界は最終断面 `Efinal_r/i` のみ)。
  - 対応内容: `frames = <interval> [complex]` の `complex` オプションを追加し、
    複素電界を `/field/frames_r`, `/field/frames_i` (nframes x Ny x Nx) に出力する。
    既定 (complex なし) は従来どおり強度のみで、容量も変わらない。
    `tools/plot_ixz.py` に位相アニメーション `*_phase.gif` を追加
    (振幅で不透明度を落とし、強度の低い領域の無意味な位相を目立たせない)。
  - 検証: `beamtilt = 3.0` の傾斜ビームで、第 0 フレームの x 方向位相勾配が
    解析値 `k0*n0*sin(3deg)` と相対誤差 0.1% で一致。
    `/field/frames` と `|frames_r + i*frames_i|^2` の整合性 1.8e-7。
    `complex` 無指定時は `frames_r` が出力されないこと (後方互換) を確認。
    既存回帰はすべて不変。CUDA 12.0 コンパイル検証済み。

- [x] **モード結合率の z 推移が無い** ✅ 対応済み (第 7 回の続き)
  - 場所: `sol/solve_bpm.cpp` / `cuda/solve_bpm.cu` / `tools/plot_ixz.py` /
    `ReadMe.md` / `data/sample/fiber_modes.ofd` / `.github/workflows/ci.yml`
  - 現状 (対応前): `/modes` にモード形状と neff は出力していたが、伝搬中に
    「どのモードにどれだけパワーが残っているか」は追えなかった。GUI で
    テーパ/ツイスト/曲げによるモード変換を可視化するには z 推移が必要。
  - 対応内容: `modes = <n>` 指定時のみ `/trace/overlap` (nModes x Nz) を出力する。
    η_m(z) = |<φ_m, E>|^2 / (||φ_m||^2 ||E||^2) で、モード順序は `/modes/neff` と同じ。
    Σ_m η_m < 1 のぶんが放射モード。出口断面の値は
    `modes : overlap eta_<m> = ...` として obpm.log にも出力する
    (HDF5 を開かずに CI/GUI から確認できるようにするため)。
    `tools/plot_ixz.py` に `*_overlap.png` (各モードの推移 + 合計) を追加。
  - 実装: CPU は `computeOverlap` テンプレートで近軸/拡張の両パス共通。
    CUDA も同じ式をホスト側で計算する (モード形状はホストにしかなく、
    1 ステップあたり nModes*N の演算で済むため)。拡張パスは `Iyz` 用に
    取得済みの全電界を再利用し、近軸パスは `modes` 指定時のみ全電界を転送する
    (未指定時は従来どおり行/列だけの転送でコスト増なし)。
  - 検証: `fiber_modes.ofd` でモード整合励振 η_1 = 0.9999 に対し
    ガウシアン励振では η_1 = 0.9938 と、電力保存 (T = 0.99999 vs 0.99988) と
    整合する差が出ることを確認。定量検証として `beamtilt = 1 0` を加えた場合の
    η_1 = 0.9457 が、ガウシアン近似の解析値
    η = exp(-(k0*n0*sin(θ)*w0/2)^2) = 0.9440 と相対誤差 0.2% で一致。
    多モード構造 (n_core = 1.56) では中心ガウシアン励振に対し
    対称モード (LP01 = 0.77, LP02 系 = 3e-4) のみが励振され、反対称の LP11 系が
    ~0 になるという対称性の要請を満たすこと、直線導波路では各 η_m が z に
    よらずほぼ一定 (固有状態) であることを確認。
    既存回帰はすべて不変 (fiber 3.122518e+02、ctest 3 スイート)。
    CUDA 12.0 コンパイル検証済み (実機実行は GPU 待ち)。

---

# main とのマージ (並行ブランチ `claude/openfdtd-implementation-fg2n84` の取り込み)

別セッションが同じ領域を並行実装したため、main とのマージで 11 ファイルが衝突した。
機能を落とさない方針で以下のとおり解決した。

| ファイル | 解決方針 |
|---|---|
| `include/obpm.h` / `sol/input_data.c` | 両側の BPM 構造体フィールドと既定値をすべて残す |
| `sol/solve_bpm.cpp` | main の波長掃引 (`wlsweep`) ループ構造を採用し、本ブランチの `/trace`・`/field/Iyz`・複素 frames・モード解析の各バッファを掃引ループの**外**へ移設。モード解析は波長依存のため掃引ループ内で毎回解き直す (前回分を `delete[]` してから再確保) |
| `cuda/solve_bpm.cu` | `launch = mode` (main) と `modes = <n>` (本ブランチ) の両ブロックを併存させる |
| `post/postbpm.c` | main のページ間引き (最大 6 枚) + 本ブランチの z 位置つきタイトルを合成 |
| `CLAUDE.md` / `AGENTS.md` | main の構成 (アーキテクチャ地図・Gotchas・CI 節) を土台に、本ブランチの規約・rules 参照・作業の進め方を統合 |
| `ReadMe.md` | 両側のキーワード行・メタデータ・CI 説明をすべて残す |
| `.github/workflows/ci.yml` | ONN 検証は main の `tests/check_activation.py` に一本化し、本ブランチの modes/taper/symmetry スモークを 3 OS に追加 |
| `CMakeLists.txt` / `.claude/settings.json` | 双方の内容を統合 |

- **励振キーワードの競合**: `launch = mode <m> [coef...]` (重ね合わせ可) と
  `modes = <n> [excite]` (解析 + 基本モード励振) は独立に実装されていた。
  両方を残し、**併用時は launch を優先**して `excite` を警告つきで無視する
  (CPU / CUDA とも同一の判定)。解析と `/modes`・`/trace/overlap` の出力は
  `launch` 併用時も行う。
- 検証: ctest 9 本 (両ブランチのスイート) 全通過、`data/` 全 18 サンプルが
  正常終了 + HDF5 出力あり、fiber 回帰 3.122518e+02 が不変。
  両系統の併用 (`launch` + `modes ... excite` → 警告 + launch 優先) と
  `wlsweep` + `modes` (波長ごとに neff が 1.447033 → 1.447279 と正しく変化)
  を実行確認。CUDA 12.0 コンパイル検証済み。

- [x] **obpm_post が `/trace` (z 伝搬の推移) 未対応** ✅ 対応済み
  - 場所: `post/postbpm.c` / `.github/workflows/ci.yml` / `ReadMe.md`
  - 現状 (対応前): `obpm_post` は `/field` の等高線 (Ixz・Efinal・frames) のみ
    描画し、`/trace` のスカラー推移は描けなかった。GUI (OpenFDTD-X) は
    `obpm_post` の ev2 出力を表示するため、伝搬に沿ったグラフを見るには
    `tools/plot_ixz.py` (Python + matplotlib) が必要だった。
  - 対応内容: `/trace` の折れ線グラフページを追加した。
    - `read1d` / `readTrace` : 1 次元データセットの読み出し
      (`/trace/z` と長さが違うものは横軸に載せられないため読み捨てる)
    - `tracePage` : z を横軸とする折れ線グラフを 1 枚描く
      (`ev2dlib_grid` / `ev2dlib_func2` / `ev2dlib_Xaxis` / `ev2dlib_Yaxis`。
       レイアウトは `post/plot2dFeed.c` の流儀に合わせた)
    - ページ構成 : パワー / ピーク強度 / ビーム幅 (x,y) / 重心 (x,y) /
      モード結合率 (先頭 4 モード)。**単位の異なる量を同じ縦軸に混ぜない**
      ため 1 量ずつページを分ける。縦軸は 0 を含む範囲にして原点からの
      大小が読めるようにした。
  - 検証: `fiber_modes.ofd` (Nz=200, modes=3 excite) で ev.ev2 が
    2 ページ → 7 ページに増加。描画された軸範囲を HDF5 と突き合わせ、
    power 3.4637e-11 / peak 1.0234 / width max 4.7997e-06 (= max(wx,wy)) /
    centroid max 1.6733e-07 / eta max 9.9998e-01 がすべて一致することを確認。
    ビーム幅 4.80um は SMF-28 の LP01 (Marcuse 近似 4.7um) と整合。
    後方互換 : `/trace` を削除した HDF5 でも落ちずに従来どおり描画すること、
    `data/` 全 18 サンプルで `obpm_post` が正常終了することを確認。
    fiber 回帰 (output power 3.122518e+02、bpm_ixz.csv の md5
    71567bae452c8b134f0a6c201bd8a69c) は不変。ctest 9 本通過。
  - CI: 3 OS の post スモークに ev2 のページ数判定 (>= 6) と
    `BPM: propagation traces` のログ判定を追加した。

---

# 第 8 回監査: CUDA の「未対応」警告が実装状況と食い違っていた

- [x] **CUDA 版で `tpa` / `powersweep` が実装済みなのに「無視される」と警告していた** ✅ 対応済み
  - 場所: `cuda/solve_bpm.cu` / `ReadMe.md` / `.claude/rules/cuda.md` / `AGENTS.md`
  - 現状 (対応前): 第 4 回で CUDA 版に TPA / powersweep を実装した際、
    実装前に入れていた「無視される」警告ブロックを消し忘れていた。
    `obpm_cuda` で `tpa` を指定すると
    `*** warning : keyword 'tpa (CUDA 版は未対応 : CPU 版 obpm を使用してください)'
    is ignored by the BPM solver.` と表示されるが、**実際には適用されている**。
    ReadMe と `.claude/rules/cuda.md` も「CPU 版のみ」のままだった。
  - 影響: サイレント無視の逆で、**正しい結果を「無視された」と誤認させる**。
    利用者は有効な計算結果を捨てるか、不要に CPU 版へ切り替えることになる。
  - 実装状況の確認 (削除前に CUDA 側の完全性を検証):
    - TPA カーネル : 近軸パス (`applyTPA` + `scalePrecisePowerByTPA`) と
      拡張パス (`wabpm_gpu_tpa`) の**両方**に存在
    - 物理スケーリング : `physcale` / `E0rawpow` / 初期界の再正規化あり
    - 実効断面積 `A_eff` の実計算とログ出力あり
    - 入力パワー掃引ループ (`for isweep`) と `sweepPin`/`sweepPout` の簿記あり
    - `activation_curve.csv` の出力あり
    → CPU 版と同じ機能一式が揃っており、「未対応」は誤り。
  - 対応内容: `ignored[]` から `tpa` / `powersweep` を削除し、代わりに
    **「実装済みだが GPU 実機未検証」の note** を出すようにした
    (`*** note : tpa / powersweep are implemented in the CUDA solver but not yet
    validated on real GPU hardware. Cross-check against the CPU solver (obpm)
    if accuracy matters.`)。`wlsweep` は実際に未対応なので警告のまま残す。
    ReadMe / `.claude/rules/cuda.md` / `AGENTS.md` の記述も実態に合わせ、
    「実装したら『無視される』警告を必ず消す」ことを規約に追加した。
  - 検証: CUDA 12.0 でのコンパイル検証済み。CUDA 側で実際に TPA が適用される
    コードパス (両伝搬パス・掃引ループ・CSV 出力) をソース上で確認。
    CPU 版の挙動・回帰には影響しない (CPU 側のコード変更なし)。

---

# 第 9 回: 性能 — BPM が並列実行されていなかった

- [x] **C++ ソースに `-fopenmp` が渡らず BPM が完全に直列実行** ✅ 対応済み
  - 場所: `CMakeLists.txt` / `bpm/FDBPMpropagator.c`
  - 現状 (対応前): `CMakeLists.txt` が `OpenMP::OpenMP_C` のみをリンクしており、
    その `-fopenmp` は `$<COMPILE_LANGUAGE:C>` で守られるため C++ ソース
    (`sol/solve_bpm.cpp` / `bpm/wabpm.cpp` / `bpm/modes.cpp`) に渡っていなかった。

    ```
    C_FLAGS   = -std=gnu99 -fopenmp     ← C には付く
    CXX_FLAGS = -std=gnu++17            ← C++ には付かない
    ```

    さらに `bpm/FDBPMpropagator.c` の `#pragma omp for` は `solve_bpm.cpp` 側の
    parallel 領域に束縛されるため、C++ 側が OpenMP 無しだと **orphaned** となり
    直列化する。結果として **BPM は全経路が単一スレッド**で動いていた。
    バナーは `CPU+OpenMP`、`-n 4` も受け付けて `thread=4` と表示するため
    気づきにくい。
  - プロファイル (400x400x100, 4 コア機): 伝搬が全体の 90%
    (ADI 6.07s / 屈折率スライス再構築 0.52s / 記録 0.17s、計 6.76s)。
  - 対応内容:
    1. `CMakeLists.txt` の `obpm` / `obpm_cuda` に `OpenMP::OpenMP_CXX` を追加
    2. **決定性の確保** — 有効化しただけでは実行ごとに結果が変わり
       `wlsweep_check` が落ちた。原因は 2 つあり両方を修正:
       - 集約が `#pragma omp atomic` で、スレッドの加算順が実行ごとに変わる
         → スレッド番号順の加算に変更
       - 集約に寄与するループが `schedule(dynamic)` で、各スレッドが担当する
         行が実行ごとに変わり部分和の項の集合が変化する
         → `schedule(static)` に変更 (行あたりの計算量は一定なので不利なし)
  - 効果 (4 コア、`taskset` で実 CPU 数を制御して実測):

    | ケース | 1 CPU | 4 CPU | 向上 |
    |---|---|---|---|
    | 近軸 400x400x100 | 6.47s | 3.43s | 1.89x |
    | 広角 400x400x100 | 10.66s | 4.38s | 2.44x |
    | 半ベクトル slab | 3.26s | 1.43s | 2.28x |

  - 検証:
    - 再現性: 同一入力を 3 回実行し `spectrum.csv` が完全一致、
      HDF5 も**全データセットがビット一致** (md5 差はコンテナ内部バイトのみ)
    - `ctest` 9 本通過 (`wlsweep_check` = 掃引と単独実行の完全一致を含む)
    - fiber 回帰 `output power = 3.122518e+02` がスレッド数 1/1/3 で不変
    - `data/` 全 18 サンプルが正常終了
    - CUDA 12.0 コンパイル検証済み
  - 既知の制約: スレッド数が変わると電力集約の部分和の個数が変わるため
    `bpm_ixz.csv` の最終桁が動く (1 スレッド `71567bae...` /
    3 スレッド `22d46894...`)。`output power` の 7 桁は不変。
    md5 での回帰比較は同一スレッド数同士で行うことを CLAUDE.md /
    `.claude/rules/testing.md` に明記した。
    完全なスレッド数非依存にするには行ごとの部分和配列を用意して
    固定順で総和する必要があり、別途検討とする。

- [x] **`/trace` の power/peak が数値ノイズに拡大表示される** ✅ 対応済み
  - 場所: `tools/plot_ixz.py`
  - 現状 (対応前): 縦軸が `1e-15 + 3.463e-11` のようなオフセット表示になり、
    実際は相対 2e-4 の平坦な線が激しいギザギザに見えていた。
  - 対応内容: power / peak を入口値で正規化し、0 を含む縦軸で描く
    (伝搬損失 P(z)/P(0) として読める形にする)。`obpm_post` 側は元から
    縦軸に 0 を含めるため変更不要。

---

# 第 10 回: 物理 — 屈折率分布が y に半セルずれていた (Yee 格子 ID の流用)

- [x] **BPM の屈折率分布が界の座標に対して y に半セルずれていた** ✅ 対応済み
  - 場所: `sol/solve_bpm.cpp` / `cuda/solve_bpm.cu`
  - 現状 (対応前): BPM は材料 ID を FDTD の `iEx` から読んでいたが、`iEx` は
    Yee 千鳥格子の **(Xc, Yn, Zn)** (y と z は節点) でサンプリングされる
    (`sol/setupId.c:137`)。一方 BPM の界の座標は **(Xc, Yc)** (セル中心)。
    つまり屈折率分布だけが y に半セルずれていた。
  - 実測 (fiber_modes, 90x90 格子, 円形コア r=4.1um):
    - メッシュ Xc/Yc は完全対称 (1e-21) なのにコア形状は上下反転と 48 セル不一致
    - コア重心 y = +0.167um = ちょうど dy/2
    - 実際のマスクは (Xc, Yn) サンプリングの予測と**完全一致** (476 セル)。
      正しい (Xc, Yc) なら 468 セル・非対称 0
  - 影響: 円形ファイバで centroid_y ≠ 0 / width_x ≠ width_y、`symmetry = y` が
    使えない (第 7 回で「y 対称面は要注意」と文書化していた制約の正体)
  - 対応内容: セル中心 (Xc, Yc, Zc) で直接ラスタライズするヘルパ
    `bpmIdPoint` / `bpmIdSlice` を追加 (setupId() と同じ描画順 = 後の geometry が
    上書き・同じ既定値 id 0・同じ eps)。BPM 内の全 iEx 参照 (n_0 自動検出・
    励振モード解析・近軸 n_in・拡張 n2・テーパ/ツイスト補間・TPA・n_out) を
    置き換えた。CPU / CUDA 両方に同一実装。z も節点 Zn ではなくセル中心 Zc で
    サンプリングする (z 不変構造では同一、階段構造では境界がセル中点になる)。
  - 検証:
    - 対称性の回復: 上下反転不一致 48 → **0**、コア重心 (0, 0)、
      centroid_y ~1e-13 m、width_x/width_y 終端 4.7739/4.7729um (差 0.02%)。
      マスクは理想の (Xc, Yc) ラスタライズと**厳密一致**
    - `symmetry = y` が全領域のちょうど 1/2 (1.561253e+02 = 3.122506e+02 / 2) に
      なることを確認 → 第 7 回の y 対称面の制約を解消 (ReadMe とサンプルの
      注意書きを更新)
    - 解析解: LP01 neff = 1.447135 (厳密解 1.447167 と 3.2e-5、CI 許容 ±5e-4 内)、
      チルト結合率 0.9454 vs 解析値 0.9440 (0.15%、修正前 0.18% から改善)
    - ctest 9 本通過 (ONN 解析解 / モードビート理論値 / 波長掃引一致を含む)
    - `data/` 全 18 サンプル正常終了。CUDA 12.0 コンパイル検証済み
  - **基準値の更新** (屈折率分布が変わるため参照値が動く):
    - fiber: output power 3.122518e+02 → **3.122506e+02** (スレッド数 1/2/4 で不変)
    - taper_twist: 2.258562e+02 → 2.258039e+02 (CI 窓 [220,230] 内)
    - fiber_symmetry: 1.561259e+02 → 1.561253e+02 (CI 窓 [155.0,157.5] 内)
    - fiber_modes: ガウシアン T 0.99988 → 0.99987、eta_1 0.9938 → 0.9944、
      チルト eta_1 0.9457 → 0.9454
    - CI の grep・CLAUDE.md・AGENTS.md・サンプルコメントを一括更新した

---

# 第 11 回: 機能追加 — 屈折率分布の直接入力 (ripfile)

- [x] **任意の屈折率分布を与えられない (geometry プリミティブのみ)** ✅ 対応済み
  - 場所: `include/obpm.h` / `sol/input_data.c` / `sol/solve_bpm.cpp` /
    `cuda/solve_bpm.cu` / `data/sample/grin.ofd` / `tools/make_grin_rip.py`
  - 現状 (対応前): 上流 BPM-Matlab は RIP (屈折率分布) を任意の行列で与えられる
    のに対し、移植では OpenFDTD の geometry プリミティブ経由に限定されていた。
    GRIN のような滑らかな分布は階段近似 (fiber_gi_mmf.ofd の 8 層) しか
    できなかった。上流比の実装漏れ。
  - 対応内容: `ripfile = <path>` キーワードを追加。
    - `.csv` : 2D (Ny 行 x Nx 列, 実数 n)。z 不変
    - `.h5`/`.hdf5` : データセット `/rip/n`。2D または 3D (Nz x Ny x Nx)
    - 寸法はメッシュと厳密一致 (不一致は期待寸法つきエラーで停止)
    - 相対パスは入力 .ofd のディレクトリ基準でも解決 (`InputPath` グローバルを
      追加し main で設定。GUI からの起動を想定)
    - 実装は屈折率スライス構築を `bpmNSlice` (rip があればそれ、無ければ
      材料 ID から) に一本化。全 n 参照サイト (n_0 検出・励振/モード解析・
      近軸 n_in・拡張 n2・テーパ/ツイスト補間・n_out) が共通経路を通る。
      TPA の β は従来どおり geometry の材料 ID (rip とは独立、併用可)。
      CPU / CUDA 両対応 (CUDA は 12.0 コンパイル検証のみ)
  - 検証:
    - **後方互換**: ripfile 未指定で fiber の出力が main と**ビット一致**
      (output power 3.122506e+02、bpm_ixz.csv md5 同一スレッド数で一致)
    - **解析解 (2D CSV)**: 放物型 GRIN `n(r) = 1.47*sqrt(1-2*0.02*(r/10um)^2)`
      にオフセット入射 (3um) したビーム重心の自己集束ピッチが、理論値
      Λ = 2πa/√(2Δ) = 314.16um に対しコサインフィットで 315.60um
      (**誤差 0.46%**)。放物分布では近軸ビームが厳密に調和振動するため、
      これは実装の直接検証になる
    - **3D HDF5**: fiber.ofd の geometry と同じ n 値を 3D `/rip/n` で与えた
      実行が、geometry 版と全出力 (Efinal/Ixz/power/n_out) **ビット一致**
    - エラー処理: ファイル無し・要素数不一致で期待寸法つきメッセージを出して停止
    - ctest 9 本通過。CI に GRIN ピッチ判定 (±2%) のスモークを追加 (Linux/macOS)

# 第 12 回: 物理バグ — モードソルバの導波判定が参照屈折率 n0 と比較していた

- [x] **GRIN (n0 = 軸上最大屈折率) でモードが 1 本も見つからない** ✅ 対応済み
  - 場所: `bpm/modes.cpp` / `include/bpm/wabpm.h` / `sol/solve_bpm.cpp` /
    `cuda/solve_bpm.cu`
  - 現状 (対応前): 導波判定のしきい値が `v0 = n0^2` (参照屈折率)。n0 は
    位相基準にすぎず、GRIN のように n0 = 軸上最大値と置く入力では
    neff < n0 が常に成り立つため、導波構造があっても
    「modes : 0 / 1 converged」となった (grin.ofd + modes = 1 で確認)。
    solve_bpm 側の励振ゲート `if (nmax > n0)` も同じ誤りで、
    launch = mode 自体がスキップされた。
  - 対応内容: しきい値を「計算領域の境界リング上の max Re(n^2)」に変更
    (`wabpm_guided_threshold` を bpm/modes.cpp に新設し公開)。物理的には
    「クラッド (境界) より速い波は放射モード」という導波条件そのもの。
    対称境界 (symx/symy) の鏡像エッジ (ix=0 / iy=0) はリングから除外。
    虚軸ステップ幅 dz の見積り (mu_max) も同じしきい値を使う。
- [x] **launch = mode ブロックの Wm.symx / Wm.symy が未初期化 (UB)** ✅ 対応済み
  - 場所: `sol/solve_bpm.cpp` / `cuda/solve_bpm.cu` (launch = mode ブロック)
  - 現状 (対応前): モード解析 (`modes =`) 側は初期化していたが、
    launch = mode 側の `wabpm_params Wm` はスタック値のまま
    `wabpm_find_modes` に渡していた (対称境界と組み合わせると不定動作)。
  - 対応内容: 両ファイルの launch ブロックで `Wm.symx = BPM.symx` /
    `Wm.symy = BPM.symy` を明示。CPU / CUDA 鏡写し
    (CUDA は 12.0 コンパイル検証のみ、実機未検証)。
  - 検証:
    - **解析解 (GRIN LG00)**: 放物分布の基本モード neff = n0 − √(2Δ)/(a k0)
      = 1.465066 に対し、grin.ofd + modes = 1 で **neff = 1.465062
      (誤差 4e-6)**。対応前は 0 本 (検出不能) だった
    - **回帰**: fiber 3.122506e+02 不変 / fiber_modes.ofd neff 1.447135・
      eta_1 0.999265・output power 3.100507e+02 不変 / ctest 9 本通過
      (mode_beat・wlsweep_check 含む)
    - fiber_mode.ofd (launch = mode) の neff は 1.447167 → 1.447135 に変化
      (dz の見積り変更による収束点の差)。厳密解 1.447167 との誤差は
      3.2e-5 で許容 (±5e-4) 内。modes = 経路の値と一致するようになった
    - CI: GRIN スモークに neff 判定 (解析解 1.465066 ± 5e-4) を追加。
      判定退行 (n0 比較に戻ると neff 行が消える) を検知する
