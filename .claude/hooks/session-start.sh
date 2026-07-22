#!/usr/bin/env bash
# SessionStart hook: main との分岐状況をセッション冒頭に知らせる。
# (このリポジトリは別セッションの PR で main が随時進むため、
#  古い main 基準で作業を始める事故を防ぐ)
# 失敗してもセッションを止めない (すべて best-effort)。
set +e

cd "$(git rev-parse --show-toplevel 2>/dev/null)" || exit 0

git fetch --quiet origin main 2>/dev/null

branch=$(git branch --show-current 2>/dev/null)
head=$(git rev-parse --short HEAD 2>/dev/null)
main=$(git rev-parse --short origin/main 2>/dev/null)

echo "[session-start] branch: ${branch:-?} @ ${head:-?}, origin/main @ ${main:-?}"

if [ -n "$main" ]; then
    behind=$(git rev-list --count HEAD..origin/main 2>/dev/null)
    ahead=$(git rev-list --count origin/main..HEAD 2>/dev/null)
    if [ "${behind:-0}" -gt 0 ]; then
        echo "[session-start] *** main が ${behind} コミット先行しています。作業前に取り込みを検討してください (/sync-check 参照)。"
        git log --oneline HEAD..origin/main 2>/dev/null | head -5
    fi
    if [ "${ahead:-0}" -gt 0 ]; then
        echo "[session-start] このブランチは main より ${ahead} コミット進んでいます (未マージ)。"
    fi
fi

dirty=$(git status --short 2>/dev/null | wc -l)
if [ "${dirty:-0}" -gt 0 ]; then
    echo "[session-start] 注意: 作業ツリーに ${dirty} 件の未コミット変更があります。"
fi

exit 0
