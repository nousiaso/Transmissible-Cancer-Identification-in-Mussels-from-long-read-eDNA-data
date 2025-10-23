#!/usr/bin/env bash
# Marker-only competitive mapping of BTN2 markers (plus optional extra panel)
# Exposes multi-mappers and applies strict uniqueness gates.
set -euo pipefail
IFS=$'\n\t'
cd "$(dirname "$0")/.."

# --------------------- config ---------------------
source config.env
[[ -f config.local ]] && source config.local

MARKER_FA="refs/marker_competitive.fasta"
MMI="refs/marker_competitive.mmi"

ALL_BAM="align_markers/markers.sorted.bam"        # secondaries kept at map-time
PASS_BAM="align_markers/markers.sorted.pass.bam"  # primaries + gates
RAW_STATS="counts_markers/raw.idxstats.tsv"
PASS_STATS="counts_markers/passing.idxstats.tsv"

LOG="logs/02_markers.log"
mkdir -p refs align_markers counts_markers logs
# --------------------------------------------------

need(){ command -v "$1" >/dev/null || { echo "[ERROR] Missing: $1"; exit 2; }; }
need minimap2; need samtools; need awk
[[ -s "$READS" ]] || { echo "[ERROR] Reads missing: $READS" | tee -a "$LOG"; exit 2; }
[[ -s "$BTN_PANEL" ]] || { echo "[ERROR] BTN panel missing: $BTN_PANEL" | tee -a "$LOG"; exit 2; }

echo "[INFO] building marker-only FASTA" | tee "$LOG"
cat "$BTN_PANEL" > "$MARKER_FA"
# Optionally augment with broader set (e.g., extra mtCOI)
[[ -s "$EXTRA_MARKERS" ]] && cat "$EXTRA_MARKERS" >> "$MARKER_FA"
grep -q '^>' "$MARKER_FA" || { echo "[ERROR] marker FASTA empty"; exit 2; }

echo "[INFO] indexing markers (k=$INDEX_K w=$INDEX_W)" | tee -a "$LOG"
minimap2 -k"$INDEX_K" -w"$INDEX_W" -d "$MMI" "$MARKER_FA" >/dev/null

echo "[INFO] mapping with secondaries kept (-N 10)" | tee -a "$LOG"
minimap2 -t "$THREADS" -x map-ont -a -N 10 "$MMI" "$READS" \
  | samtools sort -@ "$SORT_THREADS" -o "$ALL_BAM"
samtools index -@ "$SORT_THREADS" "$ALL_BAM"
samtools idxstats "$ALL_BAM" > "$RAW_STATS"

echo "[INFO] filtering: primaries ; MAPQ>=$MAPQ_MIN ; len>=... ; ID>=$ID_MIN" | tee -a "$LOG"
samtools view -h -@ "$SORT_THREADS" -F 0x904 "$ALL_BAM" \
| awk -v Q="$MAPQ_MIN" -v Ln="$MINLEN_NUC" -v Lm="$MINLEN_MITO" -v ID="$ID_MIN" '
  BEGIN{OFS="\t"}
  /^@/ {print; next}
  {
    if ($5 < Q) next
    r=$3; c=$6; Lref=0; t=c
    while (match(t,/[0-9]+[MIDNSHP=X]/)){
      seg=substr(t,RSTART,RLENGTH); op=substr(seg,length(seg),1)
      n=substr(seg,1,length(seg)-1)+0
      if(op=="M"||op=="="||op=="X") Lref+=n
      t=substr(t,RSTART+RLENGTH)
    }
    need = (r ~ /^(EF1A_|H4_)/) ? Ln : Lm
    if (Lref < need) next
    nm=0; for(i=12;i<=NF;i++) if($i ~ /^NM:i:/){split($i,a,":"); nm=a[3]+0; break}
    idn = (Lref>0)?(1.0 - nm/Lref):0.0
    if (ID>0 && idn < ID) next
    print
  }' \
| samtools view -@ "$SORT_THREADS" -Sb -o "$PASS_BAM"
samtools index -@ "$SORT_THREADS" "$PASS_BAM"
samtools idxstats "$PASS_BAM" > "$PASS_STATS"

# ---------------- audit: BTN2-labelled reads & their alternative placements ---
echo "[INFO] audit BTN2-labelled reads (any hint of BTN2 in ALL hits)" | tee -a "$LOG"
BTN2_NAMES="logs/btn2_hits.read_names.txt"
ALL_FOR_BTN2="logs/btn2_hits.all_alignments.tsv"

# 1) gather read names that align to BTN2-labelled refs at any score (ALL_BAM)
samtools view "$ALL_BAM" \
| awk '$3 ~ /^(EF1A_|H4_|mtCR_|mtCOI_)/ {print $1}' \
| sort -u > "$BTN2_NAMES"

# 2) dump *all* alignments for those reads to reveal multi-mapping
echo -e "qname\trname\tflag\tmapq\taligned_ref\tNM\tidentity" > "$ALL_FOR_BTN2"
if [[ -s "$BTN2_NAMES" ]]; then
  samtools view "$ALL_BAM" \
  | awk 'BEGIN{while((getline r < "logs/btn2_hits.read_names.txt")>0) keep[r]=1}
         {if(keep[$1]) print}' \
  | awk 'BEGIN{OFS="\t"}
        {
          q=$1; r=$3; flag=$2+0; mapq=$5+0; cig=$6; t=cig; L=0
          while (match(t,/[0-9]+[MIDNSHP=X]/)){
            seg=substr(t,RSTART,RLENGTH); op=substr(seg,length(seg),1)
            n=substr(seg,1,length(seg)-1)+0
            if(op=="M"||op=="="||op=="X") L+=n
            t=substr(t,RSTART+RLENGTH)
          }
          nm=0; for(i=12;i<=NF;i++) if($i ~ /^NM:i:/){split($i,a,":"); nm=a[3]+0; break}
          idn=(L>0)?(1.0-nm/L):0.0
          print q,r,flag,mapq,L,nm,idn
        }' >> "$ALL_FOR_BTN2"
else
  echo "[INFO] no BTN2-labelled reads at any score (on marker panel)" | tee -a "$LOG"
fi

# ---------------- small per-alignment audit for PASS bam ----------------------
echo -e "rname\tqname\tmapq\taligned_ref\tNM\tidentity" > logs/markers.pass.audit.tsv
samtools view "$PASS_BAM" \
| awk 'BEGIN{OFS="\t"}
  {
    r=$3; q=$1; mapq=$5; t=$6; L=0
    while (match(t,/[0-9]+[MIDNSHP=X]/)){
      seg=substr(t,RSTART,RLENGTH); op=substr(seg,length(seg),1)
      n=substr(seg,1,length(seg)-1)+0
      if(op=="M"||op=="="||op=="X") L+=n
      t=substr(t,RSTART+RLENGTH)
    }
    nm=0; for(i=12;i<=NF;i++) if($i ~ /^NM:i:/){split($i,a,":"); nm=a[3]+0; break}
    idn=(L>0)?(1.0-nm/L):0.0
    print r,q,mapq,L,nm,idn
  }' >> logs/markers.pass.audit.tsv

echo "[OK] marker-only mapping complete. PASS stats: $PASS_STATS" | tee -a "$LOG"
