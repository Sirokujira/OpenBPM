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

- [ ] **CUDA 版が `tpa` / `powersweep` 未対応** (⚠ サイレント無視は解消済み)
  - 場所: `cuda/solve_bpm.cu`
  - 対応状況:
    1. ✅【小】実行時警告を追加済み (「tpa / powersweep are not supported in CUDA version」
       を表示して線形の単発計算にフォールバック)。ReadMe も更新。CUDA 12.0 でコンパイル検証済み
    2. 【大・未対応】TPA 減衰カーネル (E *= exp(-β/2·|E|²·Δz), 要素毎) を `bpm/wabpm.cu` と
       近軸 CUDA パスに実装し、掃引ループを移植する (GPU 実機がないため実行検証は不可、
       コンパイル検証のみ可能。着手前に要方針判断)

- [x] **CI が単体テストを実行していない** ✅ 対応済み
  - 場所: `.github/workflows/ci.yml`
  - 対応内容: Linux / macOS ジョブの configure に `-DWITH_TESTS=ON` を追加し、
    ビルド直後に `ctest --output-on-failure` (test_wabpm / test_modes / test_allset の
    解析解検証) を実行するステップを追加した。

## 中優先度

- [ ] **MPI 版は BPM 未対応 (FDTD ソルバのみ)**
  - 場所: `src/mpi_Main.c` (`solve()` のみ呼び出し、`solve_bpm` への分岐なし)、`mpi/`
  - 現状: `-DWITH_MPI=ON` でビルドされる MPI 版は OpenFDTD 由来の時間領域ソルバで、
    BPM 伝搬は実行できない。
  - 対応案: ADI の領域分割 (行/列方向の転置通信) が必要で規模が大きい。
    需要の有無を判断のうえ、対応しない場合は ReadMe に「MPI 版は FDTD のみ」と明記する。

- [x] **TPA の解析解との数値回帰テストがない** ✅ 対応済み
  - 場所: `.github/workflows/ci.yml` / `tools/check_activation.sh` (新規)
  - 対応内容: `activation_curve.csv` の全掃引点を平面波近似の解析解
    `T = 1/(1 + β·(P_in/A_eff)·L)` と比較する検証スクリプトを追加し (A_eff / L は
    obpm.log から自動取得、許容 ±8%)、Linux / macOS の ONN スモークを単調飽和判定から
    このスクリプト呼び出しに置き換えた。ローカル実行で全 8 点一致 (最大誤差 4.1%) を確認。

## 低優先度

- [ ] **obpm_post が `/field/frames` 未対応**
  - 場所: `post/postbpm.c` (Ixz と Efinal のみ描画)
  - 現状: `frames = <interval>` で記録したスナップショットは `tools/plot_ixz.py` (GIF)
    でのみ可視化でき、obpm_post のページ出力には含まれない。
  - 対応案: 需要があれば postbpm.c にフレームのページ出力を追加 (現状は Python 側で充足)。

- [ ] **半ベクトル (pol=x/y) 指定時のモードソルバの直交化が近似**
  - 場所: `bpm/modes.cpp` (冒頭コメントに記載)
  - 現状: Stern 差分は非対称のため、半ベクトル演算子での Gram-Schmidt deflation は
    厳密には随伴モードを使うべき。スカラー使用では問題なし。
  - 対応案: 半ベクトルモード解析の需要が生じた時点で随伴系の deflation を実装。

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

- [ ] **テーパ / ツイスト機能が未公開 (カーネルは対応済み)**
  - 場所: `bpm/FDBPMpropagator.c` (`taperPerStep` / `twistPerStep`、BPM-MATLAB 由来) /
    `sol/solve_bpm.cpp:146-147` (常に 0 を設定)
  - 現状: 近軸カーネルは導波路のテーパ (スケーリング) とツイスト (回転) の座標変換に
    対応しているが、入力キーワードがなく常に無効。現状テーパは `data/fiber_taper.ofd` の
    ように geometry の階段近似で代用している。
  - 対応案: 入力キーワード (例: `taper = <rate>` / `twist = <rad/m>`) を追加して
    カーネル機能を公開する (2D RIP 前提のため近軸パス限定。拡張パスとの併用は不可と明記)。

## 低優先度

- [ ] **対称境界 (xSymmetry / ySymmetry) が未公開 (カーネルは対応済み)**
  - 場所: `bpm/FDBPMpropagator.c` (対称境界対応) / `sol/solve_bpm.cpp:123-124` (常に 0)
  - 現状: 対称構造で計算領域を 1/2〜1/4 にできるカーネル機能が使われていない。
  - 対応案: 入力キーワードで公開する (最適化のため機能追加の優先度は低)。

- [ ] **Windows CI のテストカバレッジが Linux/macOS より狭い**
  - 場所: `.github/workflows/ci.yml` (build-windows ジョブ)
  - 現状: Windows は単体テスト (ctest) を実行せず、ONN スモークも単調性判定のみ
    (解析解 ±8% 判定は Linux/macOS のみ)。
  - 対応案: Windows ジョブにも `-DWITH_TESTS=ON` + ctest を追加し、
    ONN 判定を PowerShell 版解析解比較 (または Git Bash で `check_activation.sh`) に揃える。

## 補足: 精査したが「対応済み / 対応不要」と判断した項目

- 分散性材料 (`material = 2 <einf> ...`) は BPM では `einf` 近似で処理される
  (単一波長計算のため妥当、`sol/solve_bpm.cpp:75-77` にコメントあり)。
- `beamtilt` は拡張パス (wabpm) でも初期界の位相ランプとして反映済み。
- `powersweep` の掃引ループは近軸/拡張の両パスを内包しており併用可能。
- `bend` は近軸 (CPU/CUDA)・拡張 (CPU/CUDA) の全パスで等価屈折率法に対応済み。

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

- [ ] **モード形状 (`/modes`) の可視化が未対応**
  - 場所: `tools/plot_ixz.py` / `post/postbpm.c` (いずれも `/modes` を読まない)
  - 現状: `modes` キーワードで HDF5 に出力される `/modes/mode<i>` と `/modes/neff` を
    可視化する経路がない (HDF5 規約「新規データセット追加時は plot_ixz.py と
    postbpm.c の対応も検討」への追随漏れ)。数値 (neff) は obpm.log で確認可能。
  - 対応案: `tools/plot_ixz.py` に `/modes` があればモード形状のページを追加する
    (neff をタイトルに表示)。obpm_post 側は需要が生じた時点で対応。

## 未対応項目の一覧 (第 5 回時点の残件サマリ = 9 件)

| # | 項目 | 優先度 | 規模 |
|---|---|---|---|
| 1 | CUDA 版の TPA/掃引フル実装 (警告は対応済み) | 中 | 大 (要方針判断) |
| 2 | テーパ/ツイストの公開 (カーネル対応済み) | 中 | 中 |
| 3 | MPI 版の BPM 対応 or「FDTD のみ」と明記 | 中 | 大 / 小 |
| 4 | GPU 実機での実行検証 | 低 | - (環境要因) |
| 5 | obpm_post の `/field/frames` 対応 | 低 | 小 |
| 6 | モード形状 (`/modes`) の可視化 (新規) | 低 | 小 |
| 7 | 半ベクトルモードの随伴直交化 | 低 | 中 |
| 8 | 対称境界 (x/ySymmetry) の公開 | 低 | 小 |
| 9 | Windows CI のテストカバレッジ拡充 | 低 | 小 |

---

# 第 6 回監査 (AGENTS.md 追加後)

前回 (第 5 回) 以降の変更はドキュメントのみ (AGENTS.md 新規 + CLAUDE.md 相互参照、
コード変更なし)。main にも新規コミットなし。マーカー走査はコード 0 件。
**新規の未実装項目はなし。残件は第 5 回のサマリ (9 件) から変化なし。**

---

## テスト

`-DWITH_TESTS=ON` で以下の単体テストがビルドされる (`ctest` で実行):

| テスト | 対象 | 検証内容 |
|---|---|---|
| `test_wabpm` | `bpm/wabpm.cpp` | 回折 (近軸/広角)・エネルギー保存・吸収減衰・曲げ偏向・半ベクトル整合 |
| `test_modes` | `bpm/modes.cpp` | LP01/LP11 実効屈折率 (分散方程式の厳密解と比較)・直交性 |
| `test_allset` | `include/bpm/allset.hpp` | findModes / modeSuperposition / offsetField / tiltField |

## 補足: 確認したが「対応不要」と判断した項目

- `sol/`, `mpi/`, `cuda/`, `post/` に多数ある `fprintf(stderr, ...)` は正常なエラーハンドリング
  （ファイルオープン失敗・HDF5 書き込み失敗等）であり、実装漏れではない。
- `setNCladding()` / `setShapes()` / `setDisplayScaling()`（`bpm/model.cpp`）は
  旧 API を意図的に無効化する `throw` 実装であり、仕様どおり。
