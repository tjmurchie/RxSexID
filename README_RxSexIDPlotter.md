# rx_sex_id_plotter.R (RxSexIDPlotter)

## What it does
Creates a compact strip plot of Rx values from one or more RxSexID result TSVs.

- Fill colors encode **Sex + Confidence** (Female = magenta shades; Male = teal shades)
- Ambiguous is a **single category** (white fill)
- Point size encodes **Total_Mapped_Reads** on a **log10 scale**

The plot also includes dashed guide lines at Rx=0.60 and Rx=0.80.

## Requirements
R packages:
- ggplot2, dplyr, readr, stringr, forcats, scales, tidyr
Optional:
- svglite (for SVG output)

Install in R:
```r
install.packages(c("ggplot2","dplyr","readr","stringr","forcats","scales","tidyr"))
install.packages("svglite") # optional
```

## Usage
```bash
Rscript rx_sex_id_plotter.R --out-prefix sexing --inputs \
  Bear-Arc=/path/BearSexing_Ursus-mapped.tsv \
  Bear-Can=/path/BearSexing_Canis-mapped.tsv
```

## Output
- `<prefix>_rx_strip.png`
- `<prefix>_rx_strip.pdf`
- `<prefix>_rx_strip.svg` (if svglite installed)
- `<prefix>_combined.tsv`

## Read-depth (mapped reads) legend
The mapped-read size breaks are computed **from the input data** and then rounded to “nice” values.
If there are too few samples to estimate breaks robustly, the script falls back to a small set of
sensible defaults.
