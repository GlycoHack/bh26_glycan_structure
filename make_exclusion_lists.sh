#!/bin/bash
# 除外候補リストの生成
#   入力: <prefix>skeletoncode-features.txt と <prefix>rescode-unique.txt
#         (analyze_rescode.sh の出力)
#   出力: <prefix>no-anomeric-residues.txt      rescode / glycans / smiles
#         <prefix>bicyclic-residues.txt         rescode / glycans / type / smiles
#         <prefix>no-anomeric-residues-rescode.smi    CDK Depict 用 (ラベル=ResCode)
#         <prefix>no-anomeric-residues-skeleton.smi   CDK Depict 用 (ラベル=SkeletonCode)
#         <prefix>bicyclic-residues-rescode.smi
#         <prefix>bicyclic-residues-skeleton.smi
#
# 使い方:
#   ./make_exclusion_lists.sh storm-
#
# 注意: SMILES は molwurcs.jar が毎回生成するため、実行ごとに原子の並び順が変わる
#       (分子としては同一)。ResCode 列と件数は決定的に再現される。
set -uo pipefail

PRE="${1:-}"
JAR_PATH="./molwurcs.jar"
SKF="${PRE}skeletoncode-features.txt"
RCU="${PRE}rescode-unique.txt"

[ -f "$SKF" ]      || { echo "ERROR: $SKF not found (run analyze_rescode.sh first)" >&2; exit 1; }
[ -f "$RCU" ]      || { echo "ERROR: $RCU not found (run analyze_rescode.sh first)" >&2; exit 1; }
[ -f "$JAR_PATH" ] || { echo "ERROR: $JAR_PATH not found" >&2; exit 1; }

WORK=$(mktemp -d "${TMPDIR:-/tmp}/exclusion.XXXXXX") || exit 1
trap 'rm -rf "$WORK"' EXIT

# ---- 1) アノマー位になりうる炭素を持たない残基 -------------------------------
# アノマー位を供給できるのは a (アノマー) / o (アルデヒド) / O (ケトン) のみ。
# いずれも持たない = グリコシド結合を形成できない。
awk -F'\t' 'NR>1 && $7==0 && $12==0 {print $1}' "$SKF" | sort -u > "$WORK/noano.sk"

awk -F'\t' 'NR==FNR{s[$1]=1; next}
{ r=$1; p=index(r,"_"); k=(p?substr(r,1,p-1):r); q=index(k,"-"); if(q) k=substr(k,1,q-1)
  if (k in s) print $2"\t"$1 }' "$WORK/noano.sk" "$RCU" \
  | sort -t$'\t' -k1,1nr -k2,2 > "$WORK/noano.res"

# ---- 2) 2環状残基 ------------------------------------------------------------
# (A) 骨格の環closure が2つ以上   (B) 2結合点をもつ MAP が架橋環を作る
awk -F'\t' 'BEGIN{OFS="\t"}
{ r=$1; n=split(r,f,"_"); c=0; br=0
  for(i=2;i<=n;i++){
    if (f[i] !~ /\*/ && f[i] ~ /^[0-9?]+-[0-9?]+$/) c++
    else { q=index(f[i],"*"); if(q){ m=substr(f[i],q); if(gsub(/\*/,"*",m)>=2) br=1 } } }
  if      (c>=2) print $2,$1,"backbone_bicyclic"
  else if (br)   print $2,$1,"bridged_ring" }' "$RCU" \
  | sort -t$'\t' -k3,3 -k1,1nr -k2,2 > "$WORK/bi.res"

# ---- 3) 各 ResCode を単残基 WURCS にして SMILES 化 ---------------------------
smiles_for() {   # $1: ResCode を第2列に持つファイル
    awk -F'\t' '{print "WURCS=2.0/1,1,0/["$2"]/1/"}' "$1" \
      | java -jar "$JAR_PATH" -i wurcs -o smi -n 2>/dev/null
}
smiles_for "$WORK/noano.res" > "$WORK/noano.smi"
smiles_for "$WORK/bi.res"    > "$WORK/bi.smi"

for pair in "noano.res:noano.smi" "bi.res:bi.smi"; do
    a="$WORK/${pair%%:*}"; b="$WORK/${pair##*:}"
    [ "$(wc -l < "$a")" = "$(wc -l < "$b")" ] || {
        echo "ERROR: line mismatch between $a and $b" >&2; exit 1; }
done

{ printf 'rescode\tglycans\tsmiles\n'
  paste -d'\t' <(cut -f2 "$WORK/noano.res") <(cut -f1 "$WORK/noano.res") "$WORK/noano.smi"
} > "${PRE}no-anomeric-residues.txt"

{ printf 'rescode\tglycans\ttype\tsmiles\n'
  paste -d'\t' <(cut -f2 "$WORK/bi.res") <(cut -f1 "$WORK/bi.res") <(cut -f3 "$WORK/bi.res") "$WORK/bi.smi"
} > "${PRE}bicyclic-residues.txt"

# ---- 4) CDK Depict 用 "<SMILES> <label>" ------------------------------------
# https://www.simolecule.com/cdkdepict/depict.html
awk -F'\t' 'NR>1 && $3!="" { print $3" "$1 }' "${PRE}no-anomeric-residues.txt" > "${PRE}no-anomeric-residues-rescode.smi"
awk -F'\t' 'NR>1 && $3!="" { sk=$1; sub(/[-_].*$/,"",sk); print $3" "sk }' "${PRE}no-anomeric-residues.txt" > "${PRE}no-anomeric-residues-skeleton.smi"
awk -F'\t' 'NR>1 && $4!="" { print $4" "$1 }' "${PRE}bicyclic-residues.txt" > "${PRE}bicyclic-residues-rescode.smi"
awk -F'\t' 'NR>1 && $4!="" { sk=$1; sub(/[-_].*$/,"",sk); print $4" "sk }' "${PRE}bicyclic-residues.txt" > "${PRE}bicyclic-residues-skeleton.smi"

printf "Done.\n"
printf "  %-46s %4s residues, %6s occurrences\n" \
  "${PRE}no-anomeric-residues.txt" \
  "$(($(wc -l < "${PRE}no-anomeric-residues.txt")-1))" \
  "$(awk -F'\t' 'NR>1{g+=$2} END{print g+0}' "${PRE}no-anomeric-residues.txt")"
printf "  %-46s %4s residues, %6s occurrences\n" \
  "${PRE}bicyclic-residues.txt" \
  "$(($(wc -l < "${PRE}bicyclic-residues.txt")-1))" \
  "$(awk -F'\t' 'NR>1{g+=$2} END{print g+0}' "${PRE}bicyclic-residues.txt")"
