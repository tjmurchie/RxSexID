# rx_sex_id.py (RxSexID)

> Documentation updated for **RxSexID v1.0.4**.

## What it does
`rx_sex_id.py` infers biological sex from one or more BAM files by computing:

- **Rx = (X coverage) / (autosomal baseline coverage)**  
- **Ry = (Y coverage) / (autosomal baseline coverage)** (supplementary only)

Coverage is approximated by `mapped_reads / contig_length`, where mapped reads come from `samtools idxstats`.

## Requirements
- Python 3.8+
- `samtools` in PATH

## Usage
```bash
python3 rx_sex_id.py <BAM_file_or_folder_or_glob> [more inputs ...] <reference.fasta> [output.tsv]
```

Examples:
```bash
# Many BAMs via shell glob
python3 rx_sex_id.py /path/to/*.bam /path/to/reference.fna out.tsv

# A folder of BAMs (recursively finds *.bam)
python3 rx_sex_id.py /path/to/bams/ /path/to/reference.fna out.tsv
```

## Contig identification

RxSexID detects X/Y/MT using both FASTA deflines (e.g., `chromosome X`) and common contig-name patterns (e.g., `chrX`, `SUPER_X`, `LGX`). NCBI accession version suffixes like `.1` are handled.

## Contig overrides
If the FASTA headers / contig names don’t clearly identify X/Y/MT, override with:
```bash
python3 rx_sex_id.py --x-contig <ID> --y-contig <ID> --mt-contig <ID> *.bam ref.fna out.tsv
```
Or environment variables:
- `RX_X_CONTIG`, `RX_Y_CONTIG`, `RX_MT_CONTIG`

## Autosomal baseline logic
RxSexID first tries **numbered autosomes** (recommended):
- contigs like `1..N` or `chr1..chrN`
- filtered by `RX_MIN_AUTOSOME_LENGTH` (default 10,000,000 bp)

If *no* numbered autosomes are detected, it automatically falls back to **non-sex contigs**:
- any contig that is not X, Y, or MT
- filtered by `RX_FALLBACK_MIN_AUTOSOME_LENGTH` (default 1,000,000 bp)

You can force baseline mode via:
- `RX_AUTOSOME_MODE=numbered` or `RX_AUTOSOME_MODE=nonsex`

## Interpreting results depending on reference
**Chromosome-level reference**  
Best-case: many autosomes available → stable median autosomal baseline.

**Scaffold/unplaced-heavy reference**  
Works if X (and ideally Y) are identifiable, but the autosomal baseline may be noisier.
Consider tuning the fallback minimum length upward to avoid tiny scaffolds.

**Reference lacking Y**  
Ry cannot be computed; calls rely on Rx.

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


## Output TSV
The output includes:
- `Sample`, `Sex`, `Confidence`, `Rx`, `Ry`
- `Total_Mapped_Reads`, `X_Reads`, `Y_Reads`
- autosome diagnostics (`N_Autosomes`, `Autosome_Mode`, `Min_Autosome_Length`, `Median_Auto_Cov`, etc.)
- reference fields (`Reference`, `Reference_Basename`) used by the plotter
