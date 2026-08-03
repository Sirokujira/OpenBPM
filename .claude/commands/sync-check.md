---
description: main との分岐状況を確認する (このリポジトリは並行セッションで main が進む)
---

このリポジトリは別セッションの PR が main へ随時マージされるため、
作業前・push 前に必ず分岐状況を確認してください:

1. `git fetch origin main <現在のブランチ>` で最新化
2. HEAD / origin/main / origin/<ブランチ> の各先端と ahead/behind を表示
3. main が先行している場合:
   - `git log --oneline HEAD..origin/main` で新規コミットを列挙
   - `git diff --name-only HEAD origin/main` と自ブランチの変更ファイルの
     重なり (衝突候補) を表示
   - 取り込み方針 (merge か rebase か) を提案する。自ブランチが
     main 由来の履歴しか持たない場合は `git checkout -B <branch> origin/main`
     での再スタートも選択肢に含める
4. 作業ツリーに未コミット変更がある場合は先に退避/コミットを促す

結論として「このまま作業してよいか / 先に統合が必要か」を明確に述べること。
