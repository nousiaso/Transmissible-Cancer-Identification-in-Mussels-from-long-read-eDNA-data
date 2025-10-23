#!/usr/bin/env bash
# Composite mapping: BTN2 panel + non-cancer mitogenomes
set -euo pipefail
IFS=$'\n\t'
cd "$(dirname "$0")/.."

# --------------------- config ---------------------
source config.env
[[ -f config.local ]] && source config.local

REF="refs/mussel_btn_composite.fasta"
MMI="refs/mussel_btn_composite.mmi"

RAW_BAM="align/btn_composite.allprim.bam"   # primaries only at map-time
RAW_STATS="counts/composite_raw.idxstats.tsv"

PASS_BAM="align/btn_composite.passing.bam"  # post-filtered by MAPQ + len (and optional ID)
PASS_STATS="counts/composite_passing.idxstats.tsv"

LOG="logs/01_composite.log"
mkdir -p align counts logs refs
# --------------------------------------------------

need(){ command -v "$1" >/dev/null || { echo "[ERROR] Missing: $1"; exit 2; }; }
need minimap2; need samtools; need awk

[[ -s "$READS" ]] || { echo "[ERROR] Reads missing: $READS" | tee -a "$LOG"; exit 2; }
[[ -s "$BTN_PANEL" ]] || { echo "[ERROR] BTN panel missing: $BTN_PANEL" | tee -a "$LOG"; exit 2; }
[[ -s "$NONCANCER_MITO" ]] || { echo "[ERROR] Non-cancer mito missing: $NONCANCER_MITO" | tee -a "$LOG"; exit 2; }

echo "[INFO] $(date) building composite reference" | tee "$LOG"
cat "$BTN_PANEL" "$NONCANCER_MITO" > "$REF" #create the composite fasta by combining the two
grep -q '^>' "$REF" || { echo "[ERROR] Composite FASTA empty"; exit 2; }

if [[ ! -s "$MMI" || "$MMI" -ot "$REF" ]]; then
  echo "[INFO] indexing -> $MMI" | tee -a "$LOG"
  minimap2 -d "$MMI" "$REF" >/dev/null
fi

echo "[INFO] mapping ONT (primaries only)..." | tee -a "$LOG"
# --secondary=no keeps only primaries; we filter further downstream
minimap2 -t "$THREADS" -x map-ont -a --secondary=no "$MMI" "$READS" \
  | samtools sort -@ "$SORT_THREADS" -o "$RAW_BAM"
samtools index -@ "$SORT_THREADS" "$RAW_BAM"
samtools idxstats "$RAW_BAM" > "$RAW_STATS"

# ----------- strict post-map filtering (MAPQ + len [+ identity]) -------------
# -F 0x904 : drop unmapped(0x4), secondary(0x100), supplementary(0x800)
TMP_SAM="$(mktemp -p . btn2.tmp.XXXXXX.sam)"
trap 'rm -f "$TMP_SAM"' EXIT

echo "[INFO] filtering: MAPQ>=$MAPQ_MIN ; len >=${MINLEN_MITO}bp (mt) / >=${MINLEN_NUC}bp (EF1A/H4) ; ID>=${ID_MIN}" | tee -a "$LOG"
samtools view -h -F 0x904 -q "$MAPQ_MIN" "$RAW_BAM" \
| awk -v Ln="$MINLEN_NUC" -v Lm="$MINLEN_MITO" -v ID="$ID_MIN" 'BEGIN{OFS="\t"}
  /^@/ {print; next}
  {
    r=$3; c=$6; Lref=0; t=c
    # reference-aligned length = sum of M/= /X
    while (match(t,/[0-9]+[MIDNSHP=X]/)){
      seg=substr(t,RSTART,RLENGTH); op=substr(seg,length(seg),1)
      n=substr(seg,1,length(seg)-1)+0
      if(op=="M"||op=="="||op=="X") Lref+=n
      t=substr(t,RSTART+RLENGTH)
    }
    need = (r ~ /^(EF1A_|H4_)/) ? Ln : Lm
    if (Lref < need) next
    # identity from NM (optional; set ID=0 to disable)
    nm=0; for(i=12;i<=NF;i++) if($i ~ /^NM:i:/){split($i,a,":"); nm=a[3]+0; break}
    idn = (Lref>0)?(1.0 - nm/Lref):0.0
    if (ID>0 && idn < ID) next
    print
  }' > "$TMP_SAM"

samtools view -bS "$TMP_SAM" | samtools sort -@ "$SORT_THREADS" -o "$PASS_BAM"
samtools index -@ "$SORT_THREADS" "$PASS_BAM"
samtools idxstats "$PASS_BAM" > "$PASS_STATS"
[[ -s "$PASS_STATS" ]] || { echo "[ERROR] No PASS stats"; exit 2; }

echo "[INFO] BTN2-only summary (PASS):" | tee -a "$LOG"
grep -E '^(EF1A_|H4_|mtCR_|mtCOI_)' "$PASS_STATS" | sort -k3,3nr | tee -a "$LOG"

echo "[OK] done. RAW=$RAW_STATS ; PASS=$PASS_STATS" | tee -a "$LOG"
