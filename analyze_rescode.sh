#!/bin/bash
# ResCode / SkeletonCode / MAP code の抽出と MAP 特徴分類
#   入力: ID<TAB>WURCS<TAB>SMILES (SMILES が空でない行のみを対象にする)
#   出力: <prefix>rescode-by-id.txt, <prefix>rescode-unique.txt,
#         <prefix>skeletoncode-unique.txt, <prefix>mapcode-unique.txt,
#         <prefix>mapcode-features.txt, <prefix>mapcode-class.txt
#
# 使い方:
#   ./analyze_rescode.sh wurcs-smi-id-ok.txt ""          # 接頭辞なし
#   ./analyze_rescode.sh wurcs-smi-id-storm-ok.txt storm- # storm-rescode-by-id.txt など
set -uo pipefail

IN="${1:?usage: analyze_rescode.sh <ID/WURCS/SMILES file> [prefix]}"
PRE="${2:-}"
JAR_PATH="./molwurcs.jar"

TMP_SC=$(mktemp "${TMPDIR:-/tmp}/_sc.XXXXXX")
TMP_MP=$(mktemp "${TMPDIR:-/tmp}/_mp.XXXXXX")
TMP_FEAT=$(mktemp "${TMPDIR:-/tmp}/_feat.XXXXXX")
trap 'rm -f "$TMP_SC" "$TMP_MP" "$TMP_FEAT"' EXIT

[ -f "$IN" ]       || { echo "ERROR: $IN not found" >&2; exit 1; }
[ -f "$JAR_PATH" ] || { echo "ERROR: $JAR_PATH not found" >&2; exit 1; }
[ -f MapFeatures.java ] || { echo "ERROR: MapFeatures.java not found" >&2; exit 1; }

# MapFeatures をビルド (クラスが入力より古ければ作り直す)
if [ ! -f MapFeatures.class ] || [ MapFeatures.java -nt MapFeatures.class ]; then
    javac -cp "$JAR_PATH" MapFeatures.java || exit 1
fi

# 1) ResCode: WURCS の [...] を取り出す。ResCode 内に / と * が出るため / では分割できない。
#    括弧は入れ子にならないので非貪欲マッチで安全。
awk -F'\t' 'BEGIN{OFS="\t"} $3 != "" {
  s = $2
  while (match(s, /\[[^]]*\]/)) {
    print $1, substr(s, RSTART+1, RLENGTH-2)
    s = substr(s, RSTART+RLENGTH)
  }
}' "$IN" > "${PRE}rescode-by-id.txt"

cut -f2 "${PRE}rescode-by-id.txt" | sort | uniq -c | sort -rn \
  | awk '{c=$1; $1=""; sub(/^ /,""); print $0 "\t" c}' > "${PRE}rescode-unique.txt"

# 2) SkeletonCode と MAP code
#    WURCSImporter.extractMS()/extractMOD() に準拠:
#      ResCode を _ で分割 -> 第1フィールドの最初の - より前が SkeletonCode
#      第2以降のフィールドは MOD。* を含むものだけが MAP を持ち、最初の * 以降が MAP code
awk -F'\t' 'BEGIN{OFS="\t"}
{
  res=$1; g=$2
  nf = split(res, f, "_")
  p = index(f[1], "-")
  sc = (p ? substr(f[1], 1, p-1) : f[1])
  skel[sc] += g; skel_n[sc]++
  for (i = 2; i <= nf; i++) {
    q = index(f[i], "*")
    if (q == 0) continue
    map = substr(f[i], q)
    mp[map] += g; mp_n[map]++
  }
}
END{
  for (k in skel) print k, skel[k], skel_n[k] > SC
  for (k in mp)   print k, mp[k],   mp_n[k]   > MP
}' SC="$TMP_SC" MP="$TMP_MP" "${PRE}rescode-unique.txt"
sort -t$'\t' -k2,2nr -k1,1 "$TMP_SC" > "${PRE}skeletoncode-unique.txt"
sort -t$'\t' -k2,2nr -k1,1 "$TMP_MP" > "${PRE}mapcode-unique.txt"

# 3) MAP 特徴量 (フレームワークのパーサでグラフ解析)
cut -f1 "${PRE}mapcode-unique.txt" | java -cp "$JAR_PATH:." MapFeatures 2>/dev/null > "$TMP_FEAT"
NBAD=$(tail -n +2 "$TMP_FEAT" | awk -F'\t' '$3!="OK"' | wc -l | tr -d ' ')
[ "$NBAD" -gt 0 ] && echo "WARNING: $NBAD MAP codes failed to parse" >&2
{ printf 'map\tglycans\trescodes\tnAtom\tnC\tnO\tnN\tnS\tnP\tnHal\tnAromAtom\tnRing\tnAromRing\tringSizes\tnCarbonyl\tnEsterC\tmaxChainC\tattach\n'
  paste -d'\t' <(tail -n +2 "$TMP_FEAT" | cut -f1) \
               <(cut -f2,3 "${PRE}mapcode-unique.txt") \
               <(tail -n +2 "$TMP_FEAT" | cut -f4-18)
} > "${PRE}mapcode-features.txt"
rm -f "$TMP_FEAT"

# 4) 特徴分類 (多ラベルフラグ + primary_class)
awk -F'\t' 'BEGIN{OFS="\t"}
NR==1 { print "map","glycans","primary_class","aromatic","lipid","multi_acyl","alicyclic","macrocycle","sulfate","phosphate","halogen","acyl","maxRing","maxChainC","attach"; next }
{
  map=$1; g=$2; nS=$8; nP=$9; nHal=$10; nRing=$12; nArR=$13; rs=$14; nCO=$15; chain=$17; at=$18
  maxr=0; if(rs!="-"){n=split(rs,a,","); for(i=1;i<=n;i++) if(a[i]+0>maxr) maxr=a[i]+0}
  arom=(nArR>=1); lip=(chain>=8); multi=(nCO>=2 && chain>=4)
  alic=(nRing>=1 && nArR==0); macro=(maxr>=10)
  sul=(nS>=1); pho=(nP>=1); hal=(nHal>=1); acyl=(nCO>=1)
  if      (macro) c="macrocycle"
  else if (lip)   c="lipid"
  else if (arom)  c="aromatic"
  else if (alic)  c="alicyclic"
  else if (pho)   c="phosphate"
  else if (sul)   c="sulfate"
  else if (multi) c="multi_acyl"
  else if (acyl)  c="acyl_short"
  else if (hal)   c="halogen"
  else            c="simple"
  print map,g,c,arom?1:0,lip?1:0,multi?1:0,alic?1:0,macro?1:0,sul?1:0,pho?1:0,hal?1:0,acyl?1:0,maxr,chain,at
}' "${PRE}mapcode-features.txt" > "${PRE}mapcode-class.txt"

echo "Done. Input: $IN"
printf "  %-34s %s\n" \
  "${PRE}rescode-by-id.txt"        "$(wc -l < "${PRE}rescode-by-id.txt" | tr -d ' ') rows" \
  "${PRE}rescode-unique.txt"       "$(wc -l < "${PRE}rescode-unique.txt" | tr -d ' ') ResCodes" \
  "${PRE}skeletoncode-unique.txt"  "$(wc -l < "${PRE}skeletoncode-unique.txt" | tr -d ' ') SkeletonCodes" \
  "${PRE}mapcode-unique.txt"       "$(wc -l < "${PRE}mapcode-unique.txt" | tr -d ' ') MAP codes"
