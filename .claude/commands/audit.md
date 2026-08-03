---
description: 実装漏れ監査 — 未実装・プレースホルダ・パス間の機能差を洗い出して台帳を更新する
---

OpenBPM の実装漏れ監査を実行してください。結果は `docs/implementation-checklist.md`
に「第 N 回監査」として追記し、日本語で報告します。

## 監査手順

1. **台帳の確認**: `docs/implementation-checklist.md` で既知の未対応項目と
   前回監査時点のコミットを把握する
2. **差分の把握**: 前回監査以降に main に入った変更 (`git log`) を確認し、
   新機能を監査対象に加える
3. **マーカー走査**: `TODO / FIXME / 未実装 / 未対応 / placeholder / 暫定 /
   not implemented / stub` を全ソースから grep する
   (エラーハンドリングの fprintf(stderr) は対象外)
4. **機能 × パスのマトリクス精査**: 各入力キーワード (beam / refindex / bend /
   polarization / wideangle / beamtilt / frames / tpa / powersweep / ...) が
   「CPU 近軸 / CPU 拡張 / CUDA 近軸 / CUDA 拡張 / MPI」の各パスで
   対応済みか・未対応警告があるか・サイレント無視かを確認する
5. **カーネル機能の公開状況**: カーネルが持つが入力キーワードに接続されていない
   機能 (例: taper/twist, 対称境界) がないか確認する
6. **周辺整合**: テスト・CI・ReadMe・可視化 (tools/, post/) が新機能に
   追随しているか確認する

## 報告形式

- 新規項目は 優先度 (高/中/低) と規模 (小/中/大) を付けて台帳に追記する
- 各項目に「場所 / 現状 / 対応案」を記載する
- 「精査したが対応不要」と判断した項目も根拠付きで記録する (再監査の重複を防ぐ)
- 対応済み項目には検証内容 (テスト・解析解との誤差) を残す

$ARGUMENTS
