# BTN2 eDNA Pipeline (ONT long reads)

Robust, alignment-uniqueness–aware detection of **BTN2 (Bivalve Transmissible Neoplasia, lineage 2)**
from ONT whole-genome/environmental DNA. The pipeline uses:
- **Composite mapping**: BTN2 marker panel in competition with non-cancer Mytilus mitogenomes.
- **Marker-only competitive mapping**: to stress-test specificity and expose multi-mappers.
- Strict post-map gates: primary-only, MAPQ, locus-specific aligned-length, and identity from NM/CIGAR.
- Optional **consensus polishing** per target when coverage allows.

> **Call policy**: BTN2 is **called negative** when **zero** BTN2-labelled alignments remain after all
> gates (primary-only, MAPQ≥threshold, locus-specific aligned length, identity≥threshold).
> Apparent BTN2 hits that vanish under these gates are considered **multi-mappers / non-evidence**.

## Inputs

- ONT reads: `fastq/ONT_BMK230829-BO313-006N0002_clean.fq.gz` (example; set in `config.env`)
- BTN2 marker panel (standardized FASTA): `refs/btn_markers.standardized.fasta`
- Non-cancer mitogenomes: `refs/mito_non_cancer_refs.fasta`
- (Optional) Broader mtCOI/other markers: `refs/btn_markers.extra.fasta`

## Requirements

- `bash`, `awk`
- `minimap2 >= 2.24`
- `samtools >= 1.15`
- (optional for consensus) `racon >= 1.5`

## Quick start

```bash
# 0) Configure paths & thresholds (or accept the defaults)
cp config.env config.local && edit config.local
# The scripts will auto-source config.local if present.

# 1) Composite mapping (panel + non-cancer mitogenomes)
bash scripts/01_btn2_map_composite.sh

# 2) Marker-only competitive mapping + audit
bash scripts/02_btn2_marker_only.sh

# 3) Optional per-target consensus (if coverage exists)
bash scripts/03_btn2_consensus.sh

*The main source of the data is the paper: Marisa A Yonemitsu Rachael M Giersch Maria Polo-Prieto Maurine Hammel Alexis Simon Florencia Cremonte Fernando T Avilés Nicolás Merino-Véliz Erika AV Burioli Annette F Muttray James Sherry Carol Reinisch Susan A Baldwin Stephen P Goff Maryline Houssin Gloria Arriagada Nuria Vázquez Nicolas Bierne Michael J Metzger (2019) A single clonal lineage of transmissible cancer identified in two marine mussel species in South America and Europe eLife 8:e47788.
https://doi.org/10.7554/eLife.47788
