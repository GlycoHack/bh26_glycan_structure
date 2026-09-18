#!/bin/bash
# WURCS -> SMILES 一括変換 (並列版)
#   入力: input.txt   "WURCS=..."<TAB>"G12345XX" 形式 (2列目のIDは省略可)
#   出力: wurcs-smi-id.txt   ID<TAB>WURCS<TAB>SMILES   (入力と同じ行順)
#         変換に失敗した行は SMILES 列が空になる
#
# 使い方:
#   ./convert_wurcs.sh [入力ファイル] [出力ファイル]
#   JOBS=6 ./convert_wurcs.sh           # 並列数を変える (既定 4)
#   ONLY_SUCCESS=1 ./convert_wurcs.sh   # 変換できた行のみ出力する
#
# 注意: molwurcs.jar の SMILES 出力は canonical ではなく、
#       同じ入力でも実行ごとに原子の並び順が変わる (分子としては同一)。
#       これは並列化とは無関係な jar 側の仕様。
set -uo pipefail

INPUT_FILE="${1:-input.txt}"
OUTPUT_FILE="${2:-wurcs-smi-id.txt}"
JAR_PATH="./molwurcs.jar"
JOBS="${JOBS:-4}"                   # 並列JVM数。1 JVM が約3コア使うので 4〜6 が最速
CHUNK_LINES="${CHUNK_LINES:-8000}"  # 1 JVM あたりの処理行数
ONLY_SUCCESS="${ONLY_SUCCESS:-0}"
ERR_LOG="convert_wurcs.err"

if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: $INPUT_FILE not found" >&2
    exit 1
fi
if [ ! -f "$JAR_PATH" ]; then
    echo "ERROR: $JAR_PATH not found" >&2
    exit 1
fi

WORK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/convert_wurcs.XXXXXX") || exit 1
trap 'rm -rf "$WORK_DIR"' EXIT

# 1) ID と WURCS を組で取り出す。ダブルクォートと前後空白を除去し、
#    空行と # 行は捨てる。ID も同じ行で捨てるので対応がずれない。
awk -F'\t' 'BEGIN{OFS="\t"} {
    w=$1; id=(NF>=2 ? $2 : "")
    gsub(/"/, "", w);  gsub(/"/, "", id)
    gsub(/^[ \t]+|[ \t]+$/, "", w);  gsub(/^[ \t]+|[ \t]+$/, "", id)
    if (w == "" || w ~ /^#/) next
    print id, w
}' "$INPUT_FILE" > "$WORK_DIR/pairs.txt"

TOTAL=$(wc -l < "$WORK_DIR/pairs.txt" | tr -d ' ')
if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: no WURCS found in $INPUT_FILE" >&2
    exit 1
fi
cut -f2 "$WORK_DIR/pairs.txt" > "$WORK_DIR/wurcs.txt"

# 2) チャンクに分割し、JOBS 個の JVM で並列変換する。
#    -n (--output-no-result) により失敗行も空行として出力されるので、
#    各チャンクの出力行数は入力行数と一致し、連結すれば元の行順が保たれる。
split -a 4 -l "$CHUNK_LINES" "$WORK_DIR/wurcs.txt" "$WORK_DIR/chunk."
NCHUNK=$(ls "$WORK_DIR"/chunk.* | wc -l | tr -d ' ')
echo "Starting conversion: $TOTAL WURCS / $NCHUNK chunks / $JOBS parallel JVMs"

: > "$ERR_LOG"
cat > "$WORK_DIR/worker.sh" <<'WORKER'
#!/bin/bash
chunk="$1"; jar="$2"; errlog="$3"
if java -jar "$jar" -i wurcs -o smi -n < "$chunk" 2>>"$errlog" > "$chunk.smi"; then
    echo "  done $(basename "$chunk") ($(wc -l < "$chunk" | tr -d ' ') lines)"
else
    echo "  FAILED $(basename "$chunk")" >&2
    exit 1
fi
WORKER
chmod +x "$WORK_DIR/worker.sh"

if ! ls "$WORK_DIR"/chunk.* | xargs -P "$JOBS" -I{} "$WORK_DIR/worker.sh" {} "$JAR_PATH" "$ERR_LOG"; then
    echo "ERROR: one or more chunks failed; see $ERR_LOG. $OUTPUT_FILE was not written." >&2
    exit 1
fi

# 3) チャンク出力を名前順 (= 元の行順) に連結し、行数が揃っているか確認してから貼り合わせる
cat $(ls "$WORK_DIR"/chunk.*.smi | sort) > "$WORK_DIR/smi.txt"
SMI_LINES=$(wc -l < "$WORK_DIR/smi.txt" | tr -d ' ')
if [ "$SMI_LINES" -ne "$TOTAL" ]; then
    echo "ERROR: line mismatch (wurcs=$TOTAL, smi=$SMI_LINES). $OUTPUT_FILE was not written." >&2
    exit 1
fi

# ID <TAB> WURCS <TAB> SMILES
if [ "$ONLY_SUCCESS" = "1" ]; then
    paste -d'\t' "$WORK_DIR/pairs.txt" "$WORK_DIR/smi.txt" | awk -F'\t' '$3 != ""' > "$OUTPUT_FILE"
else
    paste -d'\t' "$WORK_DIR/pairs.txt" "$WORK_DIR/smi.txt" > "$OUTPUT_FILE"
fi

OK=$(awk -F'\t' '$3 != ""' "$OUTPUT_FILE" | wc -l | tr -d ' ')
echo "Finished! Output: $OUTPUT_FILE  (ID<TAB>WURCS<TAB>SMILES)"
echo "  converted: $OK / $TOTAL  (failed: $((TOTAL - OK)), log: $ERR_LOG)"
