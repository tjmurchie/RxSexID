# rx_sex_id.py (RxSexID)

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

## Output TSV
The output includes:
- `Sample`, `Sex`, `Confidence`, `Rx`, `Ry`
- `Total_Mapped_Reads`, `X_Reads`, `Y_Reads`
- autosome diagnostics (`N_Autosomes`, `Autosome_Mode`, `Min_Autosome_Length`, `Median_Auto_Cov`, etc.)
- reference fields (`Reference`, `Reference_Basename`) used by the plotter
