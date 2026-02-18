# RxSexID

**RxSexID** infers biological sex from mapped reads using the **Rx** statistic (Skoglund et al. 2013):  
the ratio of **X-chromosome coverage** to an **autosomal baseline**.

This repository contains:
- **`rx_sex_id.py`** — sex inference from one or many BAMs using `samtools idxstats`
- **`rx_sex_id_plotter.R`** (RxSexIDPlotter) — publication-friendly Rx strip plot from one or more result TSVs
- Convenience wrappers: `run_rx_sex_id.sh`, `run_rx_sex_id_plotter.sh`

---

## Quick start

### 1) Run sex inference (RxSexID)
```bash
./run_rx_sex_id.sh /path/to/*.bam /path/to/reference.fna rxsexid_calls.tsv
```

Optional overrides if contig naming is unusual:
```bash
./run_rx_sex_id.sh --x-contig <X_contig_id> --y-contig <Y_contig_id> /path/to/*.bam reference.fna out.tsv
```

### 2) Plot results (RxSexIDPlotter)
```bash
./run_rx_sex_id_plotter.sh sexing \
  Bear-Arc=/path/BearSexing_Ursus-mapped.tsv \
  Bear-Can=/path/BearSexing_Canis-mapped.tsv \
  Marmot-Ict=/path/MarmotSexing_Ictidomys-mapped.tsv \
  Caribou-Odo=/path/CaribouSexing_Odocoileus-mapped.tsv
```

Outputs:
- `sexing_rx_strip.png`
- `sexing_rx_strip.pdf`
- `sexing_rx_strip.svg` (if `svglite` installed)
- `sexing_combined.tsv`

---

## Installation

### Requirements
- **Python 3.8+** (3.8+ recommended)
- **samtools** in `PATH` (used for `samtools idxstats`)
- For plotting: **R** plus packages listed in *README_RxSexIDPlotter.md*

### Conda (recommended)

**Python**
```bash
conda create -n rxsexid -c conda-forge python samtools
conda activate rxsexid
```

**R (plotting)**
```bash
conda create -n rxsexid-r -c conda-forge \
  r-base r-ggplot2 r-dplyr r-readr r-stringr r-forcats r-scales r-tidyr r-svglite
conda activate rxsexid-r
```

---

## Reference assembly guidance (important)

RxSexID is most reliable when reads are mapped to a **chromosome-level** or near chromosome-level reference where:
- autosomes are represented as chromosomes (often named `1..N` or `chr1..chrN`)
- the **X** chromosome is explicitly present and identifiable  
- ideally, **Y** is present and identifiable (helpful but not required)

### Workarounds for scaffold-heavy assemblies
RxSexID can still work for scaffold/unplaced-heavy references **if X (and ideally Y) are identifiable**.

Autosomal baseline selection has two modes:

1) **Numbered autosomes (default; recommended)**
- Uses contigs classified as autosomes by names like `1..N` / `chr1..chrN`
- Applies a length filter (`RX_MIN_AUTOSOME_LENGTH`, default 10,000,000 bp)

2) **Non-sex contigs (automatic fallback)**
- If no numbered autosomes are detected, RxSexID falls back to using **all non-sex contigs**
  (not X, not Y, not MT) above `RX_FALLBACK_MIN_AUTOSOME_LENGTH` (default 1,000,000 bp)

This fallback enables sex inference for assemblies where autosomes are all “unplaced scaffolds”, but:
- confidence may be lower if the baseline contig set is small or biased
- consider increasing the fallback minimum length to avoid tiny scaffolds

---

## Method summary (Rx)

For any contig *i*:
- `cov_i = mapped_reads_i / length_i`

RxSexID computes:
- `cov_X = reads_X / length_X`
- `cov_auto = median(cov_autosome_i)`
- **Rx = cov_X / cov_auto**
- **Ry = cov_Y / cov_auto** (supplementary; only if a Y contig exists)

Expected values (diploid organism mapped to a haploid reference):
- **XX** → Rx ≈ 1.0
- **XY** → Rx ≈ 0.5

Default interpretation thresholds:
- **Rx ≥ 0.80** → Female
- **Rx ≤ 0.60** → Male
- **0.60 < Rx < 0.80** → Ambiguous

---

## Confidence assignment (Low / Medium / High)

RxSexID reports a **Sex** call plus a heuristic **Confidence** label. Confidence is **not** a formal probability; it is intended as a quick quality flag based on:

- whether an autosomal baseline is reliable,
- how close **Rx** is to expected **XX (~1.0)** or **XY (~0.5)** values,
- whether **Ry** agrees (when a Y contig exists in the reference),
- mapped-read depth (used as a final cap).

### 1) “Insufficient data” overrides everything
RxSexID returns `Sex=Unknown` and `Confidence=Insufficient data` if it cannot build a stable autosomal baseline or cannot find X in the reference. This happens if any of the following are true:

- median autosomal coverage is 0  
- too few autosomes/contigs were available for the baseline  
  - **≥3 autosomes** required when `Autosome_Mode=numbered`  
  - **≥1 contig** required when using the non-sex fallback mode  
- total length of baseline autosomal contigs is **< 50,000,000 bp**  
- no X chromosome/contig is detected in the reference  

### 2) Sex call thresholds
- **Rx ≥ 0.80 → Female**
- **Rx ≤ 0.60 → Male**
- **0.60 < Rx < 0.80 → Ambiguous**

### 3) Confidence rules (summary)

**Female (Rx ≥ 0.80)**
- **High** if `0.85 ≤ Rx ≤ 1.15` and `Ry < 0.10`
- **Medium** if `0.80 ≤ Rx ≤ 1.25` and `Ry < 0.20`
- **Low** otherwise  
If a Y exists in the reference and `Ry ≥ 0.15`, RxSexID appends a note (possible male contamination or cross-species mismapping), but does not automatically flip the call.

**Male (Rx ≤ 0.60)**
- If a Y exists in the reference:
  - **High** if `Ry ≥ 0.20` and Rx is near the expected male range (`0.35–0.65`)
  - **Medium** if `Ry ≥ 0.05` (or any Y reads) and `0.40 ≤ Rx ≤ 0.60`
  - **Low** if no Y reads are detected (possible Y degradation or unmappable Y reference)
- If no Y exists in the reference:
  - **Medium** if `0.40 ≤ Rx ≤ 0.60`, else **Low**

**Ambiguous (0.60 < Rx < 0.80)**
- Default: `Sex=Ambiguous`, `Confidence=Low`
- If a Y exists and `Ry ≥ 0.15`, the label is nudged to **Male** (still low confidence)
- If a Y exists and `Y_Reads == 0` and `Ry < 0.02`, the label is nudged to **Female** (still low confidence)

### 4) Read-depth cap (defaults)
Finally, confidence is capped by `Total_Mapped_Reads`:

- `< 100` mapped reads → **Insufficient data**
- `100–999` mapped reads → **High** is capped to **Medium**
- `≥ 1000` mapped reads → no cap


## Output

Produces a tab-separated file (TSV) with one row per BAM. Typical columns include:

- sample / BAM basename
- Rx estimate
- counts / coverage stats used to compute Rx
- optional QC fields (depending on script configuration)

---

## Recommended BAM preprocessing (aDNA)

Your mapping and filters matter a lot. Typical aDNA-ready preprocessing includes:

- remove PCR duplicates
- minimum mapping quality filter (e.g., MQ ≥ 30)
- optional minimum read length filter (e.g., ≥ 24 bp)
- consider excluding problematic contigs (unplaced/alt/decoy) if they inflate noise

---


## Documentation
- **`README_rx_sex_id.md`** — details for the Python sex inference script
- **`README_RxSexIDPlotter.md`** — details for the plotting script

---

## License
MIT — see `LICENSE`.
