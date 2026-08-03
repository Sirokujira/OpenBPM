#!/bin/sh
# ONN 活性化曲線 (activation_curve.csv) の検証スクリプト (CI スモーク用)
#
# 検証内容:
#   1. 掃引点数 (8 点以上)
#   2. 透過率 T の単調非増加と飽和 (初点 T > 終点 T * 1.5)
#   3. 平面波近似の解析解 T = 1/(1 + beta*(P_in/A_eff)*L) との数値一致 (既定 ±8%)
#      (A_eff と L は obpm.log の "ONN: A_eff = ..." 行から取得)
#
# 使い方: check_activation.sh <activation_curve.csv> <obpm.log> [beta_m_per_W] [rtol]
#   beta の既定は data/sample/onn_activation.ofd の 424 cm/GW = 4.24e-9 m/W
set -e

CSV="${1:?usage: check_activation.sh <csv> <log> [beta] [rtol]}"
LOG="${2:?usage: check_activation.sh <csv> <log> [beta] [rtol]}"
BETA="${3:-4.24e-9}"
RTOL="${4:-0.08}"

AEFF=$(sed -n 's/.*ONN: A_eff = \([0-9.eE+-]*\) \[m\^2\].*/\1/p' "$LOG" | head -1)
LLEN=$(sed -n 's/.*L = \([0-9.eE+-]*\) \[m\].*/\1/p' "$LOG" | head -1)
if [ -z "$AEFF" ] || [ -z "$LLEN" ]; then
    echo "error: A_eff / L not found in $LOG" >&2
    exit 1
fi
echo "A_eff = $AEFF m^2, L = $LLEN m, beta = $BETA m/W, rtol = $RTOL"

awk -F, -v aeff="$AEFF" -v llen="$LLEN" -v beta="$BETA" -v rtol="$RTOL" '
NR > 1 {
    pin = $1; t = $3
    # 単調非増加
    if (prev != "" && t > prev * 1.000001) { print "FAIL: not monotone at line " NR; exit 1 }
    if (first == "") first = t
    prev = t
    # 解析解との比較
    tana = 1.0 / (1.0 + beta * (pin / aeff) * llen)
    err = (t - tana) / tana; if (err < 0) err = -err
    printf "P_in = %10.4g W : T = %.4f  analytic = %.4f  relerr = %.3f\n", pin, t, tana, err
    if (err > rtol) { printf "FAIL: relerr %.3f > %.3f at P_in = %g\n", err, rtol, pin; exit 1 }
}
END {
    if (NR < 9) { print "FAIL: too few sweep points: " NR - 1; exit 1 }
    if (!(first > prev * 1.5)) { print "FAIL: no saturation: first = " first ", last = " prev; exit 1 }
    print "activation curve OK: T " first " -> " prev " (all points within analytic tolerance)"
}' "$CSV"
