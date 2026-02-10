#!/usr/bin/env Rscript
# RxSexIDPlotter (v1.0.0)
# =======================
#
# Creates a compact, publication-friendly Rx strip plot from one or more RxSexID TSV outputs.
#
# What it shows
# -------------
# Each sample is shown as a point at its Rx value:
#   Rx = (X coverage) / (autosomal baseline coverage)
#
# Threshold guides (Skoglund et al. 2013-style defaults):
#   Rx >= 0.80  -> Female
#   Rx <= 0.60  -> Male
#   0.60 < Rx < 0.80 -> Ambiguous
#
# Point aesthetics:
#   - Fill color encodes Sex + Confidence (Female=magenta shades; Male=teal shades)
#   - Ambiguous is a single category (white fill, black outline)
#   - Point size encodes Total_Mapped_Reads on a log10 scale
#
# Inputs
# ------
# One or more TSV files produced by rx_sex_id.py (or older compatible formats).
# Provide inputs as:  Label=/path/to/results.tsv
#
# Output
# ------
# <prefix>_rx_strip.png
# <prefix>_rx_strip.pdf
# <prefix>_rx_strip.svg  (if svglite is installed)
# <prefix>_combined.tsv  (inputs concatenated + normalized columns)
#
# Dependencies
# ------------
# R packages: ggplot2, dplyr, readr, stringr, forcats, scales, tidyr
# Optional:   svglite (for SVG output)
#
# Install example (conda):
#   conda create -n rxsexid-r -c conda-forge r-base r-ggplot2 r-dplyr r-readr r-stringr r-forcats r-scales r-tidyr r-svglite
#
# Usage
# -----
#   Rscript rx_sex_id_plotter.R --out-prefix sexing --inputs \
#     Bear-Arc=/path/BearSexing_Ursus-mapped.tsv \
#     Bear-Can=/path/BearSexing_Canis-mapped.tsv
#
# Notes on references
# -------------------
# RxSexIDPlotter will attempt to display the NCBI assembly accession (GCF_*/GCA_*) beneath each sample
# if it is present in the TSV (columns: Reference or Reference_Basename). For best results, use rx_sex_id.py
# v1.0.0+ which always writes these fields.
#
# License: MIT (see LICENSE)
# Author: Tyler Murchie

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(scales)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  if (!(flag %in% args)) return(default)
  i <- match(flag, args)
  if (i == length(args)) return(default)
  args[[i + 1]]
}

out_prefix <- get_arg("--out-prefix", NULL)
if (is.null(out_prefix)) stop("Missing --out-prefix", call. = FALSE)

if (!("--inputs" %in% args)) stop("Missing --inputs", call. = FALSE)
i_inputs <- match("--inputs", args)
input_specs <- args[(i_inputs + 1):length(args)]
if (length(input_specs) < 1) stop("No inputs after --inputs", call. = FALSE)

parse_input <- function(spec) {
  if (str_detect(spec, "=")) {
    parts <- str_split_fixed(spec, "=", 2)
    list(label = parts[,1], path = parts[,2])
  } else {
    bn <- basename(spec)
    label <- str_replace(bn, "\\.tsv$", "")
    list(label = label, path = spec)
  }
}

dfs <- lapply(input_specs, function(s) {
  inp <- parse_input(s)
  d <- read_tsv(inp$path, show_col_types = FALSE, progress = FALSE)
  d$Dataset <- inp$label
  d
})

df <- bind_rows(dfs)

# Harmonize columns from older formats
if (!("Sex" %in% names(df)) && ("Sex_Call" %in% names(df))) df$Sex <- df$Sex_Call
if (!("Rx" %in% names(df)) && ("X_to_Auto_Ratio" %in% names(df))) df$Rx <- df$X_to_Auto_Ratio
if (!("Ry" %in% names(df)) && ("Y_to_Auto_Ratio" %in% names(df))) df$Ry <- df$Y_to_Auto_Ratio
if (!("Total_Mapped_Reads" %in% names(df))) df$Total_Mapped_Reads <- NA_real_
if (!("Confidence" %in% names(df))) df$Confidence <- NA_character_

df <- df %>%
  mutate(
    Total_Mapped_Reads = parse_number(as.character(Total_Mapped_Reads)),
    Rx = suppressWarnings(as.numeric(Rx)),
    Ry = suppressWarnings(as.numeric(Ry))
  )

# Extract NCBI assembly accession
if ("Reference_Basename" %in% names(df)) {
  df <- df %>% mutate(
    RefAcc = str_extract(as.character(Reference_Basename), "(GCF_\\d+\\.\\d+|GCA_\\d+\\.\\d+)")
  )
} else {
  df <- df %>% mutate(RefAcc = NA_character_)
}
if ("Reference" %in% names(df)) {
  df <- df %>% mutate(
    RefAcc = if_else(is.na(RefAcc) | RefAcc == "",
                     str_extract(as.character(Reference), "(GCF_\\d+\\.\\d+|GCA_\\d+\\.\\d+)"),
                     RefAcc)
  )
}

# Truncate sample names to 25 chars
df <- df %>%
  mutate(
    SampleTrunc = if_else(nchar(Sample) <= 25, Sample, paste0(substr(Sample, 1, 22), "...")),
    LabelPlot = if_else(!is.na(RefAcc) & RefAcc != "", paste0(SampleTrunc, "\n", RefAcc), SampleTrunc)
  )

# Sex / confidence cleanup
df <- df %>%
  mutate(
    Sex = str_to_title(as.character(Sex)),
    Sex = case_when(
      str_detect(Sex, "Unknown") ~ "Ambiguous",
      str_detect(Sex, "Ambiguous") ~ "Ambiguous",
      TRUE ~ Sex
    ),
    Sex = factor(Sex, levels = c("Female", "Male", "Ambiguous")),
    Confidence = str_to_title(as.character(Confidence)),
    Confidence = if_else(Confidence %in% c("Low","Medium","High"), Confidence, "Medium"),
    Confidence = factor(Confidence, levels = c("Low","Medium","High"))
  )

# Legend categories:
df <- df %>%
  mutate(
    SexConf = case_when(
      Sex == "Ambiguous" ~ "Ambiguous",
      TRUE ~ paste0(as.character(Sex), " (", as.character(Confidence), ")")
    )
  )

df$SexConf <- factor(
  df$SexConf,
  levels = c("Female (Low)","Female (Medium)","Female (High)",
             "Male (Low)","Male (Medium)","Male (High)",
             "Ambiguous")
)

female_med <- "#d64aa8"
male_med   <- "#2aa9a0"

fill_map <- c(
  "Female (Low)"="#f4b6d9",
  "Female (Medium)"=female_med,
  "Female (High)"="#8b005f",
  "Male (Low)"="#b7ebe6",
  "Male (Medium)"=male_med,
  "Male (High)"="#006a63",
  "Ambiguous"="#ffffff"
)

# Order samples within each dataset (by Rx)
df <- df %>%
  group_by(Dataset) %>%
  arrange(Rx, SampleTrunc, .by_group = TRUE) %>%
  mutate(i_in_dataset = row_number()) %>%
  ungroup()

# x positions with gap between datasets
dataset_order <- unique(df$Dataset)
dataset_sizes <- df %>% count(Dataset, name = "n")
gap <- 0.35

offset_tbl <- dataset_sizes %>%
  mutate(
    Dataset = factor(Dataset, levels = dataset_order),
    cum_n = cumsum(lag(n, default = 0)),
    offset = cum_n + (seq_along(n) - 1) * gap
  )

df <- df %>%
  left_join(offset_tbl %>% select(Dataset, offset, n), by = "Dataset") %>%
  mutate(x = i_in_dataset + offset)

ranges <- df %>%
  group_by(Dataset) %>%
  summarise(x_min = min(x), x_max = max(x), .groups = "drop") %>%
  mutate(
    x_left = x_min - 0.6,
    x_right = x_max + 0.6,
    x_center = (x_min + x_max)/2
  )

# Size mapping helper
df <- df %>%
  mutate(Total_Mapped_Reads_size = if_else(is.na(Total_Mapped_Reads) | Total_Mapped_Reads <= 0, 1, Total_Mapped_Reads))

# Data-driven "nice" breaks (about 5)
nice_round <- function(v) {
  if (!is.finite(v) || v <= 0) return(v)
  p <- floor(log10(v))
  base <- 10^p
  mant <- v / base
  opts <- c(0.5, 1, 1.5, 2, 3, 5, 7.5, 10)
  mant2 <- opts[which.min(abs(opts - mant))]
  mant2 * base
}

pos_reads <- df$Total_Mapped_Reads_size[df$Total_Mapped_Reads_size > 1]
if (length(pos_reads) >= 2) {
  lo <- log10(min(pos_reads, na.rm = TRUE))
  hi <- log10(max(pos_reads, na.rm = TRUE))
  raw <- 10^(seq(lo, hi, length.out = 5))
  br <- sapply(raw, nice_round)
  br <- sort(unique(br))
  if (length(br) < 4) {
    med <- median(pos_reads, na.rm = TRUE)
    extra <- sapply(c(med/3, med, med*3), nice_round)
    br <- sort(unique(c(br, extra)))
  }
  if (length(br) > 6) br <- br[round(seq(1, length(br), length.out = 6))]
  size_breaks <- br
} else {
  size_breaks <- c(5e5, 1e6, 1.5e6, 3e6, 1e7)
}

legend_df <- tibble(
  SexConf = factor(levels(df$SexConf), levels = levels(df$SexConf)),
  x = min(df$x, na.rm = TRUE),
  Rx = min(df$Rx, na.rm = TRUE),
  Total_Mapped_Reads_size = 1
)

# Layout values
y_top <- 1.25
y_bracket <- y_top - 0.005
y_title <- y_top + 0.07     # ABOVE the top border
bracket_drop <- 0.03

# Region label positions
x_left_text <- min(df$x, na.rm = TRUE) + 0.25

sample_bars <- tibble(x = df$x)

p <- ggplot(df, aes(x = x, y = Rx)) +
  geom_hline(yintercept = 0.60, linetype = "dashed") +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  geom_vline(data = sample_bars, aes(xintercept = x),
             inherit.aes = FALSE, colour = "grey85", linewidth = 0.35) +
  # brackets per dataset
  geom_segment(data = ranges, aes(x = x_left, xend = x_right, y = y_bracket, yend = y_bracket),
               inherit.aes = FALSE, linewidth = 0.7, colour = "black") +
  geom_segment(data = ranges, aes(x = x_left, xend = x_left, y = y_bracket, yend = y_bracket - bracket_drop),
               inherit.aes = FALSE, linewidth = 0.7, colour = "black") +
  geom_segment(data = ranges, aes(x = x_right, xend = x_right, y = y_bracket, yend = y_bracket - bracket_drop),
               inherit.aes = FALSE, linewidth = 0.7, colour = "black") +
  # dataset titles ABOVE the bounding box
  geom_text(data = ranges, aes(x = x_center, y = y_title, label = Dataset),
            inherit.aes = FALSE, size = 3.0, fontface = "bold", vjust = 0) +
  # dummy legend keys
  geom_point(data = legend_df, aes(x = x, y = Rx, fill = SexConf),
             shape = 21, color = "black", stroke = 0.25,
             size = 4, alpha = 0.0, inherit.aes = FALSE, show.legend = TRUE) +
  geom_point(aes(fill = SexConf, size = Total_Mapped_Reads_size),
             shape = 21, color = "black", stroke = 0.25, alpha = 1.0) +
  annotate("text", x = x_left_text, y = 1.07, label = "Female",
           hjust = 0, size = 3.8, colour = female_med, fontface = "italic") +
  annotate("text", x = x_left_text, y = 0.70, label = "Ambiguous",
           hjust = 0, size = 3.8, colour = "#4f4f4f", fontface = "italic") +
  annotate("text", x = x_left_text, y = 0.28, label = "Male",
           hjust = 0, size = 3.8, colour = male_med, fontface = "italic") +
  scale_fill_manual(values = fill_map, name = "Sex (confidence)", drop = FALSE, breaks = names(fill_map)) +
  scale_size_continuous(
    trans = "log10",
    range = c(2.7, 11),
    breaks = size_breaks,
    labels = comma,
    name = "Mapped reads"
  ) +
  scale_x_continuous(
    breaks = df$x,
    labels = df$LabelPlot,
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  coord_cartesian(ylim = c(0, y_top), clip = "off") +
  labs(x = NULL, y = "Rx (X coverage / autosomal baseline)") +
  guides(
    fill = guide_legend(order = 1, override.aes = list(alpha = 1, size = 4, shape = 21, color = "black", stroke = 0.25)),
    size = guide_legend(order = 2)
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 65, vjust = 1, hjust = 1, size = 7.5,
      lineheight = 0.85,
      margin = margin(t = 0)
    ),
    plot.margin = margin(t = 34, r = 5.5, b = 1, l = 14),
    legend.position = "right",
    legend.box = "vertical"
  )

n_total <- nrow(df)
plot_width <- max(6.5, min(9.0, 0.17 * n_total + 4.2))
plot_height <- 6

ggsave(paste0(out_prefix, "_rx_strip.png"), p, width = plot_width, height = plot_height, dpi = 300)
ggsave(paste0(out_prefix, "_rx_strip.pdf"), p, width = plot_width, height = plot_height)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(paste0(out_prefix, "_rx_strip.svg"), p, width = plot_width, height = plot_height)
} else {
  message("NOTE: svglite not installed; skipping SVG output.")
}

df_export <- df %>%
  mutate(Label = if_else(!is.na(RefAcc) & RefAcc != "", paste0(SampleTrunc, " | ", RefAcc), SampleTrunc)) %>%
  select(-LabelPlot)

write_tsv(df_export, paste0(out_prefix, "_combined.tsv"))

message("Wrote: ",
        paste0(out_prefix, "_rx_strip.png, "),
        paste0(out_prefix, "_rx_strip.pdf, "),
        paste0(out_prefix, "_rx_strip.svg (if enabled), "),
        paste0(out_prefix, "_combined.tsv"))
