#!/usr/bin/env python3
"""
RxSexID: Rx-based Sex Inference for Ancient and Degraded DNA
==================================================

Determines biological sex from BAM files by calculating the Rx statistic:
the ratio of X chromosome coverage to median autosomal coverage.

Method
------
The Rx statistic is defined as:

    Rx = (reads_X / length_X) / median(reads_auto_i / length_auto_i)

For a diploid organism mapped to a haploid reference genome, both X and
autosomal reads pile onto single reference copies. The resulting ratio
reflects the copy-number difference:

    Female (XX):  2 X copies / 2 autosomal copies  →  Rx ≈ 1.0
    Male   (XY):  1 X copy  / 2 autosomal copies   →  Rx ≈ 0.5

Classification thresholds (following Skoglund et al. 2013):

    Rx >= 0.80  →  Female
    Rx <= 0.60  →  Male
    0.60 < Rx < 0.80  →  Ambiguous

Y chromosome data (Ry = Y_coverage / Autosomal_coverage) is used as
supplementary evidence when available. Y chromosomes are frequently
absent from reference assemblies, poorly assembled, or dominated by
repetitive sequence, so Rx is always the primary determinant.

Confidence levels incorporate both the Rx/Ry values and the total number
of mapped reads, since low read counts produce unreliable ratios.

Input Handling
--------------
The script handles NCBI-style accession-based sequence names by parsing
the reference FASTA headers to identify chromosome assignments:

    >CM125500.1 Ictidomys tridecemlineatus ... chromosome X, ...
    >CM125490.1 Ictidomys tridecemlineatus ... chromosome 1, ...

Unlocalized scaffolds (e.g., ``chromosome X_unloc``) are aggregated with
their parent chromosome. Unplaced scaffolds are excluded from the autosomal baseline by default.
If a reference lacks numbered autosomes (common in some NCBI FASTA deflines),
the script automatically falls back to using all non-sex contigs as autosomes.

Requirements
------------
Reference assemblies and contig naming
--------------------------------------
RxSexID is most reliable when reads are mapped to a *chromosome-level* or near chromosome-level
assembly where autosomes are represented as full chromosomes (e.g., 1..N) and the X chromosome is
explicitly identified (and ideally Y as well). This allows the autosomal baseline to be computed
from true autosomes and reduces bias from fragmented scaffold sets.

However, RxSexID can still work with scaffold-heavy / unplaced assemblies provided that:
  - The X contig is identifiable (by name or FASTA header annotation), and
  - A reasonable set of non-sex contigs can be used as an autosomal baseline.

Autosomal baseline selection modes
----------------------------------
The script chooses which contigs to treat as "autosomes" in order to estimate autosomal coverage.

1) Numbered autosomes (default; recommended)
   Uses only contigs classified as autosomes based on names like:
     1..N, chr1..chrN (case-insensitive)
   A minimum length filter is applied (env: RX_MIN_AUTOSOME_LENGTH).

2) Non-sex contigs (fallback for scaffold-heavy references)
   If *no* numbered autosomes are detected (common when autosomes are all "unplaced"),
   RxSexID automatically falls back to using all contigs that are not X, Y, or MT as autosomes,
   again applying a minimum length filter (env: RX_FALLBACK_MIN_AUTOSOME_LENGTH).

You can override contig identification with:
  --x-contig <contig_id>   or env RX_X_CONTIG
  --y-contig <contig_id>   or env RX_Y_CONTIG
  --mt-contig <contig_id>  or env RX_MT_CONTIG

You can influence autosome selection with:
  env RX_AUTOSOME_MODE=numbered|nonsex
  env RX_MIN_AUTOSOME_LENGTH=<bp>              (default tuned for chromosome-level assemblies)
  env RX_FALLBACK_MIN_AUTOSOME_LENGTH=<bp>     (default 1,000,000 bp)

Interpreting results by reference type
--------------------------------------
Chromosome-level reference
  - Autosomal baseline is computed from multiple true autosomes.
  - Rx tends to cluster near ~1.0 (XX) or ~0.5 (XY) with high confidence if coverage is sufficient.

Scaffold/unplaced-heavy reference (fallback mode)
  - Autosomal baseline is computed from large non-sex contigs.
  - Rx can still be informative, but confidence may be lower if the baseline set is small or biased.
  - Increase RX_FALLBACK_MIN_AUTOSOME_LENGTH to avoid tiny scaffolds, or decrease it if no large
    contigs exist (at the cost of potential noise).

Reference without a Y chromosome
  - Ry cannot be computed (Y length = 0) and is omitted/0.
  - Sex calls rely primarily on Rx.


- Python 3.6+
- samtools (must be in PATH; used for ``samtools idxstats``)

Usage
-----
::

    python3 rx_sex_id.py <BAM_file_or_folder> <reference.fasta> [output.tsv]

Examples
--------
::

    # Process all BAMs in a directory
    python3 rx_sex_id.py ./bams/ /path/to/reference.fna results.tsv

    # Process a single BAM file
    python3 rx_sex_id.py sample.bam /path/to/reference.fna sample_sex.tsv

Output
------
1. A tab-separated file (TSV) with per-sample results including Rx, Ry,
   read counts, confidence levels, and explanatory notes.
2. A formatted summary table printed to stderr for quick assessment.

References
----------
Skoglund, P., Stora, J., Gotherstrom, A., & Jakobsson, M. (2013).
    Accurate sex identification of ancient human remains using DNA
    shotgun sequencing. Journal of Archaeological Science, 40(12),
    4477-4482. https://doi.org/10.1016/j.jas.2013.07.004

Lamnidis, T. C., Majander, K., Jeong, C., et al. (2018).
    Ancient Fennoscandian genomes reveal origin and spread of Siberian
    ancestry in Europe. Nature Communications, 9(1), 5018.
    https://doi.org/10.1038/s41467-018-07483-5

License
-------
MIT License. See repository root for details.

Author
------
Tyler Murchie
"""

import os
import sys
import subprocess
import re
import gzip
from pathlib import Path
import statistics
from collections import defaultdict
from datetime import datetime

__version__ = '1.0.0'

# ─── Classification thresholds ──────────────────────────────────────────────
# Based on Skoglund et al. (2013): Female Rx ~ 1.0, Male Rx ~ 0.5
RX_FEMALE_THRESHOLD = 0.80  # Rx >= this → Female
RX_MALE_THRESHOLD = 0.60    # Rx <= this → Male
# Values between these thresholds are classified as Ambiguous.

# Minimum mapped reads for confidence assessment
MIN_READS_HIGH_CONF = 1000   # Below this, confidence capped at Medium
MIN_READS_MEDIUM_CONF = 100  # Below this, confidence set to Insufficient data


# Defaults tuned for chromosome-level assemblies; scaffold-heavy references often need smaller thresholds
DEFAULT_MIN_AUTOSOME_LENGTH = int(os.environ.get('RX_MIN_AUTOSOME_LENGTH', 10_000_000))
FALLBACK_MIN_AUTOSOME_LENGTH = int(os.environ.get('RX_FALLBACK_MIN_AUTOSOME_LENGTH', 1_000_000))

# ─── FASTA Parsing ──────────────────────────────────────────────────────────

def parse_fasta_headers(fasta_path, x_contig=None, y_contig=None, mt_contig=None):
    """
    Parse FASTA headers to build an accession-to-chromosome mapping.

    Handles NCBI-style headers where chromosome identity is encoded in
    the description line, e.g.::

        >CM125500.1 Ictidomys tridecemlineatus isolate ... chromosome X, ...

    Parameters
    ----------
    fasta_path : str or Path
        Path to reference FASTA file (plain text or gzipped).

    Returns
    -------
    dict
        ``{accession: {'chr': chromosome_name, 'type': 'main'|'unloc'|'scaffold'}}``
    """
    fasta_path = Path(fasta_path)
    mapping = {}

    # Optional explicit contig IDs for sex chromosomes / mitochondrion.
    # These should match BAM @SQ SN: entries exactly (e.g., 'NC_079873.1').
    x_contig = (x_contig or os.environ.get('RX_X_CONTIG') or '').strip() or None
    y_contig = (y_contig or os.environ.get('RX_Y_CONTIG') or '').strip() or None
    mt_contig = (mt_contig or os.environ.get('RX_MT_CONTIG') or '').strip() or None

    # Robust header patterns across common reference FASTA styles
    x_pat = re.compile(r'(\bchromosome\s*x\b|\bchrx\b)', re.IGNORECASE)
    y_pat = re.compile(r'(\bchromosome\s*y\b|\bchry\b)', re.IGNORECASE)
    mt_pat = re.compile(r'(mitochond\w*|\bmtDNA\b|\bchrM\b|\bchrMT\b|\bMT\b)', re.IGNORECASE)

    opener = gzip.open if str(fasta_path).endswith('.gz') else open
    mode = 'rt' if str(fasta_path).endswith('.gz') else 'r'

    print(f"Parsing reference FASTA headers from {fasta_path}...", file=sys.stderr)

    with opener(fasta_path, mode) as f:
        for line in f:
            if not line.startswith('>'):
                continue

            header = line[1:].strip()
            accession = header.split()[0]

            # Explicit overrides first (highest priority)
            if x_contig and accession == x_contig:
                mapping[accession] = {'chr': 'X', 'type': 'main'}
                continue
            if y_contig and accession == y_contig:
                mapping[accession] = {'chr': 'Y', 'type': 'main'}
                continue
            if mt_contig and accession == mt_contig:
                mapping[accession] = {'chr': 'MT', 'type': 'main'}
                continue

            unloc_match = re.search(r'_unloc', header, re.IGNORECASE)

            # Robust detection of sex chromosomes / mitochondrion across header styles
            acc_up = accession.upper()
            if acc_up in ('X', 'CHRX'):
                mapping[accession] = {'chr': 'X', 'type': 'unloc' if unloc_match else 'main'}
                continue
            if acc_up in ('Y', 'CHRY'):
                mapping[accession] = {'chr': 'Y', 'type': 'unloc' if unloc_match else 'main'}
                continue
            if acc_up in ('MT', 'M', 'CHRM', 'CHRMT'):
                mapping[accession] = {'chr': 'MT', 'type': 'main'}
                continue

            hay = f"{accession} {header}"
            if x_pat.search(hay):
                mapping[accession] = {'chr': 'X', 'type': 'unloc' if unloc_match else 'main'}
                continue
            if y_pat.search(hay):
                mapping[accession] = {'chr': 'Y', 'type': 'unloc' if unloc_match else 'main'}
                continue
            if mt_pat.search(hay):
                mapping[accession] = {'chr': 'MT', 'type': 'main'}
                continue

            # Generic chromosome parsing (e.g., 'chromosome 1', 'chromosome 12', 'chromosome 3A', etc.)
            chr_match = re.search(r'chromosome\s+([^,\s]+)', header, re.IGNORECASE)
            if chr_match:
                chr_name = chr_match.group(1).strip().rstrip(',').upper()
                chr_type = 'unloc' if unloc_match else 'main'
                mapping[accession] = {'chr': chr_name, 'type': chr_type}
            else:
                # No chromosome annotation in header -> treat as unplaced/scaffold
                mapping[accession] = {'chr': 'unplaced', 'type': 'scaffold'}

    chr_counts = defaultdict(int)
    for info in mapping.values():
        chr_counts[info['chr']] += 1

    print(f"  Found {len(mapping)} sequences", file=sys.stderr)
    print(f"  Chromosomes: {dict(chr_counts)}", file=sys.stderr)

    return mapping


# ─── Chromosome Classification ──────────────────────────────────────────────

def classify_chromosome(chr_name):
    """Classify a parsed chromosome name as X, Y, autosome, MT, or other."""
    c = chr_name.upper()
    if c == 'X':
        return 'X'
    if c == 'Y':
        return 'Y'
    if c in ('MT', 'M', 'MITO'):
        return 'MT'
    if c == 'UNPLACED':
        return 'unplaced'
    if c.isdigit() or re.match(r'^[0-9]+[A-Z]?$', c):
        return 'autosome'
    return 'other'


# ─── BAM Statistics ─────────────────────────────────────────────────────────

def get_idxstats(bam_file):
    """
    Run ``samtools idxstats`` and return parsed output.

    Returns
    -------
    list of dict
        ``[{'accession': str, 'length': int, 'mapped': int, 'unmapped': int}, ...]``
    """
    result = subprocess.run(
        ['samtools', 'idxstats', str(bam_file)],
        capture_output=True, text=True, check=True
    )
    stats = []
    for line in result.stdout.strip().split('\n'):
        if line.startswith('*'):
            continue
        parts = line.split('\t')
        if len(parts) >= 4:
            stats.append({
                'accession': parts[0],
                'length': int(parts[1]),
                'mapped': int(parts[2]),
                'unmapped': int(parts[3])
            })
    return stats


def aggregate_by_chromosome(idxstats, chr_mapping):
    """
    Aggregate reads and lengths by chromosome, combining main sequences
    with their unlocalized scaffolds.

    Returns
    -------
    tuple
        ``(chr_stats_dict, unmapped_stats_dict)``
    """
    chr_stats = defaultdict(lambda: {'length': 0, 'mapped': 0, 'type': None})
    unmapped_stats = {'length': 0, 'mapped': 0}

    for entry in idxstats:
        accession = entry['accession']
        if accession in chr_mapping:
            info = chr_mapping[accession]
            name = info['chr']
            chr_stats[name]['length'] += entry['length']
            chr_stats[name]['mapped'] += entry['mapped']
            chr_stats[name]['type'] = classify_chromosome(name)
        else:
            unmapped_stats['length'] += entry['length']
            unmapped_stats['mapped'] += entry['mapped']

    return dict(chr_stats), unmapped_stats


# ─── Coverage Calculation ───────────────────────────────────────────────────

def calculate_coverage_stats(chr_stats, min_autosome_length=DEFAULT_MIN_AUTOSOME_LENGTH, autosome_mode='numbered'):
    """
    Calculate per-chromosome coverage and autosomal baseline.

    Parameters
    ----------
    chr_stats : dict
        Output of ``aggregate_by_chromosome()``.
    min_autosome_length : int
        Minimum chromosome length (bp) to include in the autosomal
        baseline.  Filters out small or fragmented chromosomes that
        could skew the median.

    Returns
    -------
    dict
        Keys: ``n_autosomes``, ``median_auto_cov``, ``mean_auto_cov``,
        ``auto_std``, ``autosome_info``, ``x_stats``, ``y_stats``.
    """
    autosome_coverages = []
    autosome_info = []
    x_stats = None
    y_stats = None

    for chr_name, stats in chr_stats.items():
        if stats['length'] == 0:
            continue

        coverage = stats['mapped'] / stats['length']
        chr_type = stats['type']

        # Autosomal baseline selection:
        # - numbered: only contigs labeled as numeric chromosomes (1..N)
        # - nonsex: any contig that is not X, Y, or MT (useful for scaffold-heavy references where autosomes are 'unplaced')
        mode = (autosome_mode or 'numbered').lower()
        if mode not in ('numbered', 'nonsex'):
            mode = 'numbered'
        if mode == 'numbered':
            include_as_auto = (chr_type == 'autosome')
        else:  # nonsex
            include_as_auto = (chr_type not in ('X', 'Y', 'MT'))

        if include_as_auto and stats['length'] >= min_autosome_length:
            autosome_coverages.append(coverage)
            autosome_info.append({
                'chr': chr_name,
                'length': stats['length'],
                'reads': stats['mapped'],
                'coverage': coverage
            })
        elif chr_type == 'X':
            x_stats = {
                'length': stats['length'],
                'reads': stats['mapped'],
                'coverage': coverage
            }
        elif chr_type == 'Y':
            y_stats = {
                'length': stats['length'],
                'reads': stats['mapped'],
                'coverage': coverage
            }

    if autosome_coverages:
        median_auto_cov = statistics.median(autosome_coverages)
        mean_auto_cov = statistics.mean(autosome_coverages)
        auto_std = (statistics.stdev(autosome_coverages)
                    if len(autosome_coverages) > 1 else 0)
    else:
        median_auto_cov = 0
        mean_auto_cov = 0
        auto_std = 0

    return {
        'n_autosomes': len(autosome_coverages),
        'median_auto_cov': median_auto_cov,
        'mean_auto_cov': mean_auto_cov,
        'auto_std': auto_std,
        'autosome_info': autosome_info,
        'x_stats': x_stats,
        'y_stats': y_stats,
    }


# ─── Sex Determination ──────────────────────────────────────────────────────

def determine_sex(cov_stats, total_mapped_reads):
    """
    Determine biological sex using the Rx statistic.

    Classification logic
    --------------------
    1. **Primary call** is based on Rx (X/Autosome coverage ratio):

       - Rx >= 0.80 → Female  (expected ~1.0 for XX)
       - Rx <= 0.60 → Male    (expected ~0.5 for XY)
       - Between    → Ambiguous

    2. **Ry** (Y/Autosome ratio) is used as supplementary evidence.
       Y chromosomes are frequently absent, incomplete, or dominated by
       repetitive sequence in reference assemblies, so Ry is never the
       sole basis for a call.

    3. **Confidence** is adjusted based on proximity to theoretical
       values, agreement between Rx and Ry, and total read depth.

    Parameters
    ----------
    cov_stats : dict
        Output of ``calculate_coverage_stats()``.
    total_mapped_reads : int
        Total mapped reads in the BAM (for read-depth quality check).

    Returns
    -------
    tuple
        ``(sex_call, confidence, rx, ry, explanation)``
    """
    median_auto = cov_stats['median_auto_cov']
    x_stats = cov_stats['x_stats']
    y_stats = cov_stats['y_stats']

    # ── Insufficient data checks ────────────────────────────────────────
    autosome_mode_used = cov_stats.get('autosome_mode_used', 'numbered').lower()
    min_autosomes_required = 3 if autosome_mode_used == 'numbered' else 1
    autosome_total_len = sum(a['length'] for a in cov_stats.get('autosome_info', []))

    if median_auto == 0 or cov_stats['n_autosomes'] < min_autosomes_required or autosome_total_len < 50_000_000:
        return ('Unknown', 'Insufficient data', 0, 0,
                'Insufficient autosomal data for a reliable baseline')

    if x_stats is None:
        return ('Unknown', 'Insufficient data', 0, 0,
                'No X chromosome found in reference')

    # ── Calculate Rx and Ry ─────────────────────────────────────────────
    rx = x_stats['coverage'] / median_auto
    ry = (y_stats['coverage'] / median_auto) if y_stats else 0

    # Read-depth confidence cap
    if total_mapped_reads < MIN_READS_MEDIUM_CONF:
        max_conf = 'Insufficient data'
    elif total_mapped_reads < MIN_READS_HIGH_CONF:
        max_conf = 'Medium'
    else:
        max_conf = 'High'

    y_in_ref = y_stats is not None and y_stats['length'] > 0
    y_reads = y_stats['reads'] if y_stats else 0

    sex = 'Unknown'
    confidence = 'Low'
    explanation = ''

    # ── Female call (Rx >= 0.80) ────────────────────────────────────────
    if rx >= RX_FEMALE_THRESHOLD:
        sex = 'Female'

        if 0.85 <= rx <= 1.15 and ry < 0.10:
            confidence = 'High'
        elif 0.80 <= rx <= 1.25 and ry < 0.20:
            confidence = 'Medium'
        else:
            confidence = 'Low'

        explanation = f'Rx={rx:.3f} (female expected ~1.0)'

        if y_in_ref and ry >= 0.15:
            explanation += (f'; NOTE: elevated Ry={ry:.3f} '
                           f'({y_reads:,} Y reads) — possible male '
                           f'contamination or cross-species mismapping')

    # ── Male call (Rx <= 0.60) ──────────────────────────────────────────
    elif rx <= RX_MALE_THRESHOLD:
        sex = 'Male'

        if y_in_ref:
            if ry >= 0.20:
                confidence = 'High' if 0.35 <= rx <= 0.65 else 'Medium'
                explanation = (f'Rx={rx:.3f} (male expected ~0.5), '
                               f'Ry={ry:.3f} confirms Y presence')
            elif ry >= 0.05 or y_reads > 0:
                confidence = 'Medium' if 0.40 <= rx <= 0.60 else 'Low'
                explanation = (f'Rx={rx:.3f} (male expected ~0.5); '
                               f'Y coverage low (Ry={ry:.3f}, '
                               f'{y_reads:,} Y reads)')
            else:
                confidence = 'Low'
                explanation = (f'Rx={rx:.3f} suggests male, but no Y '
                               f'reads detected — possible Y degradation '
                               f'or unmappable Y reference')
        else:
            confidence = 'Medium' if 0.40 <= rx <= 0.60 else 'Low'
            explanation = (f'Rx={rx:.3f} (male expected ~0.5); '
                           f'no Y chromosome in reference for confirmation')

    # ── Ambiguous (0.60 < Rx < 0.80) ───────────────────────────────────
    else:
        sex = 'Ambiguous'
        confidence = 'Low'
        explanation = (f'Rx={rx:.3f} falls between male (~0.5) '
                       f'and female (~1.0) expectations')

        if y_in_ref and ry >= 0.15:
            sex = 'Male'
            explanation += f'; Ry={ry:.3f} supports male'
        elif y_in_ref and y_reads == 0 and ry < 0.02:
            sex = 'Female'
            explanation += '; absence of Y reads favours female'

    # ── Apply read-depth cap ────────────────────────────────────────────
    if max_conf == 'Insufficient data':
        confidence = 'Insufficient data'
        explanation += (f' [WARNING: only {total_mapped_reads:,} mapped '
                        f'reads — estimate unreliable]')
    elif max_conf == 'Medium' and confidence == 'High':
        confidence = 'Medium'
        explanation += (f' [capped: {total_mapped_reads:,} reads — '
                        f'moderate depth]')

    return sex, confidence, rx, ry, explanation


# ─── BAM Processing ─────────────────────────────────────────────────────────

def process_bam(bam_file, chr_mapping, reference_fasta):
    """Process a single BAM file and return sex determination results."""
    bam_path = Path(bam_file)
    sample = bam_path.stem

    # Ensure BAM index exists
    bai = bam_path.with_suffix('.bam.bai')
    if not bai.exists():
        bai_alt = Path(str(bam_path) + '.bai')
        if not bai_alt.exists():
            print(f"  Indexing {sample}...", file=sys.stderr)
            subprocess.run(['samtools', 'index', str(bam_path)], check=True)

    print(f"Processing {sample}...", file=sys.stderr)

    idxstats = get_idxstats(bam_path)
    chr_stats, unmapped = aggregate_by_chromosome(idxstats, chr_mapping)
    autosome_mode = os.environ.get('RX_AUTOSOME_MODE', 'numbered').lower()
    cov_stats = calculate_coverage_stats(chr_stats, autosome_mode=autosome_mode)
    # If no autosomes were detected (common when a reference has X/Y labeled but autosomes are 'unplaced'),
    # fall back to using all non-sex contigs as the autosomal baseline.
    if cov_stats['n_autosomes'] == 0 and autosome_mode == 'numbered':
        print('  NOTE: No numbered autosomes detected in reference; using non-sex contigs as autosomal baseline '
              f"(min length {FALLBACK_MIN_AUTOSOME_LENGTH:,} bp).", file=sys.stderr)
        cov_stats = calculate_coverage_stats(
            chr_stats,
            min_autosome_length=FALLBACK_MIN_AUTOSOME_LENGTH,
            autosome_mode='nonsex'
        )
        cov_stats['autosome_mode_used'] = 'nonsex_fallback'
        cov_stats['min_autosome_length_used'] = FALLBACK_MIN_AUTOSOME_LENGTH
    else:
        cov_stats['autosome_mode_used'] = autosome_mode
        cov_stats['min_autosome_length_used'] = DEFAULT_MIN_AUTOSOME_LENGTH

    total_reads = sum(entry['mapped'] for entry in idxstats)

    sex, confidence, rx, ry, explanation = determine_sex(cov_stats,
                                                          total_reads)

    x_reads = cov_stats['x_stats']['reads'] if cov_stats['x_stats'] else 0
    y_reads = cov_stats['y_stats']['reads'] if cov_stats['y_stats'] else 0
    x_length = cov_stats['x_stats']['length'] if cov_stats['x_stats'] else 0
    y_length = cov_stats['y_stats']['length'] if cov_stats['y_stats'] else 0

    return {
        'Sample': sample,
        'Reference': reference_fasta,
        'Reference_Basename': os.path.basename(reference_fasta),
        'Sex': sex,
        'Confidence': confidence,
        'Rx': rx,
        'Ry': ry,
        'Total_Mapped_Reads': total_reads,
        'X_Reads': x_reads,
        'Y_Reads': y_reads,
        'N_Autosomes': cov_stats['n_autosomes'],
        'Autosome_Mode': cov_stats.get('autosome_mode_used', ''),
        'Min_Autosome_Length': cov_stats.get('min_autosome_length_used', DEFAULT_MIN_AUTOSOME_LENGTH),
        'Median_Auto_Cov': cov_stats['median_auto_cov'],
        'X_Coverage': cov_stats['x_stats']['coverage'] if cov_stats['x_stats'] else 0,
        'Y_Coverage': cov_stats['y_stats']['coverage'] if cov_stats['y_stats'] else 0,
        'X_Length': x_length,
        'Y_Length': y_length,
        'Explanation': explanation,
    }


# ─── Summary Table ──────────────────────────────────────────────────────────

def print_summary_table(results, reference_name=None):
    """
    Print a formatted summary table to stderr.

    The table includes Rx, Ry, mapped read counts, and sex calls so that
    reliability can be assessed at a glance.
    """
    w = 120  # table width

    print("\n" + "=" * w, file=sys.stderr)
    print("  SEX DETERMINATION SUMMARY"
          "  (Rx method; Skoglund et al. 2013)", file=sys.stderr)
    print("=" * w, file=sys.stderr)
    print(f"  Thresholds : Female Rx >= {RX_FEMALE_THRESHOLD:.2f}"
          f"  |  Male Rx <= {RX_MALE_THRESHOLD:.2f}"
          f"  |  Ambiguous {RX_MALE_THRESHOLD:.2f} < Rx < "
          f"{RX_FEMALE_THRESHOLD:.2f}", file=sys.stderr)
    if reference_name:
        print(f"  Reference  : {reference_name}", file=sys.stderr)
    print(f"  Date       : {datetime.now().strftime('%Y-%m-%d %H:%M')}",
          file=sys.stderr)
    print(f"  Version    : RxSexID v{__version__}",
          file=sys.stderr)
    print("=" * w, file=sys.stderr)

    # Check if any sample has Y data in the reference
    any_y = any(r['Y_Length'] > 0 for r in results)

    # Column headers
    hdr = (f"{'Sample':<45s} {'Sex':>8s} {'Confidence':>18s} "
           f"{'Rx':>7s} {'Ry':>7s} {'Mapped Reads':>14s} "
           f"{'X Reads':>10s} {'Y Reads':>9s}")
    print(hdr, file=sys.stderr)
    print("-" * w, file=sys.stderr)

    for r in results:
        name = r['Sample']
        if len(name) > 44:
            name = name[:41] + "..."

        rx_s = f"{r['Rx']:.3f}" if r['Rx'] > 0 else "N/A"

        if r['Y_Length'] == 0:
            ry_s = "-"
            yr_s = "-"
        else:
            ry_s = f"{r['Ry']:.3f}" if r['Ry'] > 0 else "0.000"
            yr_s = f"{r['Y_Reads']:,}"

        row = (f"{name:<45s} {r['Sex']:>8s} {r['Confidence']:>18s} "
               f"{rx_s:>7s} {ry_s:>7s} "
               f"{r['Total_Mapped_Reads']:>14,} "
               f"{r['X_Reads']:>10,} {yr_s:>9s}")
        print(row, file=sys.stderr)

    print("-" * w, file=sys.stderr)

    # Footer notes
    print(f"  Rx = X_cov / Autosomal_cov  |"
          f"  Female (XX) ~ 1.0  |  Male (XY) ~ 0.5", file=sys.stderr)
    if any_y:
        print(f"  Ry = Y_cov / Autosomal_cov  |"
              f"  Supplementary; Y often poorly assembled",
              file=sys.stderr)
    else:
        print(f"  Ry = -  (no Y chromosome in reference)",
              file=sys.stderr)

    low = [r for r in results
           if r['Total_Mapped_Reads'] < MIN_READS_MEDIUM_CONF]
    if low:
        print(f"  WARNING: {len(low)} sample(s) with "
              f"< {MIN_READS_MEDIUM_CONF} mapped reads — "
              f"sex calls unreliable for these samples",
              file=sys.stderr)

    print("=" * w + "\n", file=sys.stderr)


# ─── Main ───────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print("\nError: Missing required arguments", file=sys.stderr)
        print("Usage: python3 rx_sex_id.py "
              "<BAM_file_or_folder> <reference.fasta> [output.tsv]",
              file=sys.stderr)
        sys.exit(1)
    # ------------------------------------------------------------------
    # Input parsing
    # ------------------------------------------------------------------
    # New calling pattern (backwards-compatible):
    #   python3 rx_sex_id.py <bam_or_dir_or_glob> [more ...] <reference.fasta> [output.tsv]
    #
    # If the final argument ends with .tsv/.txt/.csv, treat it as an output path.
    # Otherwise, write to rx_sex_id_results.tsv in the current directory.

    args = sys.argv[1:]

    # Optional flags (kept lightweight to preserve backwards compatibility without argparse):
    #   --x-contig <ID>  --y-contig <ID>  --mt-contig <ID>
    x_contig = os.environ.get('RX_X_CONTIG')
    y_contig = os.environ.get('RX_Y_CONTIG')
    mt_contig = os.environ.get('RX_MT_CONTIG')

    cleaned = []
    i = 0
    while i < len(args):
        a = args[i]
        if a in ('--x-contig', '--x') and i + 1 < len(args):
            x_contig = args[i + 1]
            i += 2
            continue
        if a in ('--y-contig', '--y') and i + 1 < len(args):
            y_contig = args[i + 1]
            i += 2
            continue
        if a in ('--mt-contig', '--mt') and i + 1 < len(args):
            mt_contig = args[i + 1]
            i += 2
            continue
        cleaned.append(a)
        i += 1

    args = cleaned
    if len(args) < 2:
        print("\nError: Missing required arguments", file=sys.stderr)
        print("Usage: python3 rx_sex_id.py \
              <BAM_file_or_folder_or_glob> [more inputs ...] <reference.fasta> [output.tsv]",
              file=sys.stderr)
        sys.exit(1)

    # Output inference
    output_exts = {'.tsv', '.txt', '.csv'}
    maybe_output = Path(args[-1])
    if maybe_output.suffix.lower() in output_exts and len(args) >= 3:
        output = str(maybe_output)
        fasta_path = Path(args[-2])
        input_specs = args[:-2]
    else:
        output = 'rx_sex_id_results.tsv'
        fasta_path = Path(args[-1])
        input_specs = args[:-1]

    if not fasta_path.exists():
        print(f"Error: Reference FASTA not found: {fasta_path}", file=sys.stderr)
        sys.exit(1)

    # Parse chromosome mapping from reference
    chr_mapping = parse_fasta_headers(fasta_path, x_contig=x_contig, y_contig=y_contig, mt_contig=mt_contig)

    # Reference name for summary header
    ref_name = fasta_path.name

    def _expand_inputs(spec: str):
        # If shell didn't expand a glob (e.g. quoted), expand it here.
        if any(ch in spec for ch in ['*', '?', '[']):
            hits = glob.glob(spec)
            return hits if hits else [spec]
        return [spec]

    def collect_bam_files(specs):
        bam_files = []
        for spec in specs:
            for expanded in _expand_inputs(spec):
                p = Path(expanded)
                if p.is_file():
                    if p.suffix.lower() != '.bam':
                        raise ValueError(f"Not a BAM file: {p}")
                    bam_files.append(p)
                elif p.is_dir():
                    bam_files.extend(sorted(p.glob('*.bam')))
                else:
                    raise ValueError(f"Input not found: {p}")

        # Deduplicate while preserving order
        seen = set()
        unique = []
        for b in bam_files:
            s = str(b.resolve())
            if s not in seen:
                seen.add(s)
                unique.append(Path(s))
        return unique

    try:
        bam_files = collect_bam_files(input_specs)
    except Exception as e:
        print(f"Error resolving BAM inputs: {e}", file=sys.stderr)
        sys.exit(1)

    if not bam_files:
        print(f"No BAM files found for inputs: {', '.join(input_specs)}", file=sys.stderr)
    if not bam_files:
        print(f"No BAM files found in inputs: {', '.join(input_specs)}", file=sys.stderr)
        sys.exit(1)

    print(f"\nFound {len(bam_files)} BAM file(s) to process\n",
          file=sys.stderr)

    # Process all BAM files
    results = []
    for bam_file in bam_files:
        try:
            reference_fasta = str(fasta_path)
            result = process_bam(bam_file, chr_mapping, reference_fasta)
            results.append(result)
            print(f"  -> {result['Sex']} ({result['Confidence']}): "
                  f"Rx={result['Rx']:.3f}, Ry={result['Ry']:.3f}, "
                  f"reads={result['Total_Mapped_Reads']:,}",
                  file=sys.stderr)
            print(f"     {result['Explanation']}\n", file=sys.stderr)

        except Exception as e:
            print(f"Error processing {bam_file}: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc(file=sys.stderr)

    # Write detailed TSV output
    if results:
        header = list(results[0].keys())
        with open(output, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for result in results:
                vals = []
                for k in header:
                    v = result[k]
                    if isinstance(v, float):
                        v = f'{v:.6f}'
                    vals.append(str(v))
                f.write('\t'.join(vals) + '\n')

        print(f"Detailed results written to: {output}", file=sys.stderr)

        # Print summary table
        print_summary_table(results, reference_name=ref_name)
    else:
        print("No results to write.", file=sys.stderr)


if __name__ == '__main__':
    main()
