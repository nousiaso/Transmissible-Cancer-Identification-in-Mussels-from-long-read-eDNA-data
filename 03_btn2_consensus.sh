#!/usr/bin/env bash
# Optional per-target consensus polishing when coverage exists.
# Two rounds of racon against single-target references.
set -euo pipefail
IFS=$'\n\t'
cd "$(dirname "$0")/.."

# --------------------- config ---------------------
source config.env
[[ -f config.local ]] && source config.local

REF="refs/mussel_btn_composite.fasta"     # from step 01
PASS_BAM="align/btn_composite.passing.bam"
PASS_STATS="counts/composite_passing.idxstats.tsv"

OUTDIR="consensus"
LOG="logs/03_consensus.log"
mkdir -p "$OUTDIR" logs
# --------------------------------------------------

need(){ command -v "$1" >/dev/null || { echo "[ERROR] Missing: $1"; exit 2; }; }
need minimap2; need samtools; need awk; need racon

[[ -s "$REF" && -s "$PASS_BAM" ]] || { echo "[ERROR] Need composite ref & PASS bam (run step 01)"; exit 2; }
[[ -s "$PASS_STATS" ]] || { echo "[ERROR] Missing $PASS_STATS"; exit 2; }

# Helper: extract a single FASTA record by exact header
extract_fa() {
  local name="$1"; local ofa="$2"
  awk -v key=">$name" 'BEGIN{p=0}
    $0==key {p=1; print; next}
    /^>/ && p==1 {p=0}
    p==1 {print}
  ' "$REF" > "$ofa"
  [[ -s "$ofa" ]] || { echo "[ERROR] Could not extract $name"; return 1; }
}

# Choose targets:
CORE=( EF1A_BTN2_H EF1A_BTN2_Hprime mtCR_BTN2_D_2copy mtCOI_BTN2_B mtCOI_BTN2_Q )
H4_TOP=($(grep -E '^H4_BTN2_candidate_' "$PASS_STATS" | sort -k3,3nr | head -6 | awk '{print $1}'))
TARGETS=()

# Add core targets if present in PASS stats (i.e., have any coverage)
for t in "${CORE[@]}"; do
  if grep -q "^${t}\b" "$PASS_STATS"; then TARGETS+=("$t"); fi
done
TARGETS+=("${H4_TOP[@]}")

echo "[INFO] Targets selected: ${#TARGETS[@]}" | tee "$LOG"
printf ' - %s\n' "${TARGETS[@]}" | tee -a "$LOG"

for t in "${TARGETS[@]}"; do
  TBAM="consensus/${t}.bam"
  TFASTQ="consensus/${t}.reads.fastq"
  TREF="consensus/${t}.ref.fa"
  POL1="consensus/${t}.polish1.fa"
  POL2="consensus/${t}.polish2.fa"

  # Extract read subset
  samtools view -b "$PASS_BAM" "$t" > "$TBAM" || { echo "[WARN] No reads for $t"; continue; }
  if [[ $(samtools view -c "$TBAM") -lt 3 ]]; then
    echo "[WARN] <$t> has <3 reads; skipping consensus." | tee -a "$LOG"
    rm -f "$TBAM"
    continue
  fi
  samtools fastq "$TBAM" > "$TFASTQ"
  extract_fa "$t" "$TREF"

  # Map → Racon ×2
  minimap2 -x map-ont -a "$TREF" "$TFASTQ" | samtools sort -o "${TBAM%.bam}.map1.bam"
  samtools index "${TBAM%.bam}.map1.bam"
  racon -m 8 -x -6 -g -8 -w 500 "$TFASTQ" "${TBAM%.bam}.map1.bam" "$TREF" > "$POL1"

  minimap2 -x map-ont -a "$POL1" "$TFASTQ" | samtools sort -o "${TBAM%.bam}.map2.bam"
  samtools index "${TBAM%.bam}.map2.bam"
  racon -m 8 -x -6 -g -8 -w 500 "$TFASTQ" "${TBAM%.bam}.map2.bam" "$POL1" > "$POL2"

  # Sanity: length ≥ 100 bp
  LEN=$(awk '/^>/ {next} {l+=length($0)} END{print l+0}' "$POL2")
  if [[ "${LEN:-0}" -ge 100 ]]; then
    echo "[OK] $t → $POL2 (len=$LEN)" | tee -a "$LOG"
  else
    echo "[WARN] $t consensus too short/empty" | tee -a "$LOG"
  fi
done

echo "[OK] consensus step done." | tee -a "$LOG"
