#!/bin/bash
# ================================================================
# run_4D_rules.sh
# 4D Apple-Peel ペーリング R1/R2/R3 一括実行スクリプト
#
# 使い方:
#   bash run_4D_rules.sh              # 全6胞体
#   bash run_4D_rules.sh 5 8          # 引数で胞体番号を指定
#   bash run_4D_rules.sh 5 8 > out.log 2>&1 &  # バックグラウンド実行
# ================================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WL_SCRIPT="$SCRIPT_DIR/run4DPeelingAllRules.m"
LOG="$SCRIPT_DIR/run_4D_rules_$(date +%Y%m%d_%H%M%S).log"

# 胞体番号が引数で与えられた場合は nList を上書き
if [ $# -gt 0 ]; then
  NLIST="{$(IFS=,; echo "$*")}"
  echo "nList overridden: $NLIST"
  OVERRIDE="nList = $NLIST;"
else
  OVERRIDE=""
fi

echo "========================================"
echo " 4D Apple-Peel Rule Comparison"
echo " Start: $(date)"
echo " Log  : $LOG"
echo "========================================"

if [ -n "$OVERRIDE" ]; then
  # nList を上書きする場合はラッパーファイルを一時生成
  # scriptDir を明示的に渡してパス解決を保証する
  WRAPPER=$(mktemp /tmp/run4D_XXXXXX.m)
  cat > "$WRAPPER" <<EOF
scriptDir = "$SCRIPT_DIR/";
$OVERRIDE
Get["$WL_SCRIPT"];
EOF
  wolframscript -file "$WRAPPER" 2>&1 | tee "$LOG"
  rm -f "$WRAPPER"
else
  wolframscript -file "$WL_SCRIPT" 2>&1 | tee "$LOG"
fi

echo ""
echo "Finished: $(date)"
echo "Log saved to: $LOG"
