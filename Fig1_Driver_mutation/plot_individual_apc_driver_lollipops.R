#!/usr/bin/env Rscript

# Per-patient APC driver lollipop plots.
#
# This intentionally reuses the original Fig1_Driver_mutation/driver_plots.Rmd
# visualization approach: trackViewer::lolliplot + grid overlays.

required_pkgs <- c("trackViewer", "GenomicRanges", "grid")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R package(s): ", paste(missing_pkgs, collapse = ", "),
    "\nThis script reuses the old Fig1 visualization and requires trackViewer.\n",
    "Load the same R/Bioc environment used for Fig1_Driver_mutation/driver_plots.Rmd, then rerun."
  )
}

suppressPackageStartupMessages({
  library(trackViewer)
  library(GenomicRanges)
  library(grid)
})

lower_trackviewer_lollipop_labels <- function(y_offset = 0.25) {
  plot_lollipops <- getFromNamespace("plotLollipops", "trackViewer")
  plot_lollipops_txt <- deparse(plot_lollipops)
  old <- "y = this.height + "
  hit <- grep(old, plot_lollipops_txt, fixed = TRUE)
  hit <- hit[grepl("grid.text", plot_lollipops_txt[hit], fixed = TRUE)]
  if (length(hit) != 1) {
    stop("Could not patch trackViewer label y-position; internal plotLollipops layout changed.")
  }
  plot_lollipops_txt[hit] <- sub(
    "y = this.height \\+ $",
    sprintf("y = this.height + "),
    plot_lollipops_txt[hit]
  )
  plot_lollipops_txt[hit + 1] <- sub(
    "feature.height,",
    sprintf("feature.height - %.4f,", y_offset),
    plot_lollipops_txt[hit + 1],
    fixed = TRUE
  )
  patched <- eval(parse(text = paste(plot_lollipops_txt, collapse = "\n")))
  environment(patched) <- asNamespace("trackViewer")
  assignInNamespace("plotLollipops", patched, ns = "trackViewer")
}

lower_trackviewer_lollipop_labels()

`%||%` <- function(x, y) {
  if (length(x) == 0 || is.na(x) || is.null(x)) y else x
}

script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1] %||% "")
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()
repo_dir <- normalizePath(file.path(script_dir, ".."))

metadata_file <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/metadata/fap_database_withexposure.csv"
cnv_file <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/CNV_SV_RT/CNV.csv"
out_dir <- file.path(script_dir, "individual_apc_driver_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_combined_pdf <- file.path(out_dir, "individual_APC_driver_lollipops.pdf")
out_summary_csv <- file.path(out_dir, "individual_APC_driver_lollipops_data.csv")

apc_length <- 2843
apc_genomic_start <- 112707498
apc_genomic_end <- 112846239

# APC MANE/canonical transcript ENST00000257430 on GRCh38, forward strand.
# CDS intervals are intersected with the translation span; final stop codon is
# retained here only for coordinate conversion and capped to apc_length below.
apc_cds <- data.frame(
  start = c(
    112754891, 112766326, 112767189, 112775629, 112780790,
    112792446, 112801279, 112815495, 112818966, 112821896,
    112827108, 112827929, 112828856, 112834951, 112837553
  ),
  end = c(
    112755025, 112766410, 112767390, 112775737, 112780903,
    112792529, 112801383, 112815593, 112819344, 112821991,
    112827247, 112828006, 112828972, 112835165, 112844126
  )
)
apc_cds$cds_start <- cumsum(c(1, head(apc_cds$end - apc_cds$start + 1, -1)))
apc_cds$cds_end <- apc_cds$cds_start + apc_cds$end - apc_cds$start

group_colors <- c(
  "polyp" = "#DB7575",
  "adenoma" = "#F5A623",
  "ACA" = "#377EB8",
  "normal" = "grey",
  "germline APC" = "black"
)

apc_genomic_to_aa <- function(pos) {
  pos <- suppressWarnings(as.integer(pos))
  out <- rep(NA_integer_, length(pos))

  for (i in seq_along(pos)) {
    p <- pos[i]
    if (is.na(p) || p < apc_genomic_start || p > apc_genomic_end) {
      next
    }

    hit <- which(apc_cds$start <= p & apc_cds$end >= p)
    if (length(hit) == 0) {
      # If the breakpoint is intronic but inside APC, place it at the next
      # affected coding amino acid. If it is after the final CDS base, cap to
      # the last amino acid.
      hit <- which(apc_cds$start > p)[1]
      if (is.na(hit)) {
        out[i] <- apc_length
        next
      }
      cds_offset <- apc_cds$cds_start[hit] - 1
    } else {
      cds_offset <- apc_cds$cds_start[hit[1]] + (p - apc_cds$start[hit[1]]) - 1
    }

    out[i] <- min(apc_length, floor(cds_offset / 3) + 1)
  }

  out
}

extract_apc_aa_position <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""

  pos <- rep(NA_integer_, length(x))

  # Common APC protein notations: E1309Dfs*4, Y159*, X438_splice, *177*.
  m1 <- regexec("^[A-Za-z\\*X]([0-9]+)", x, perl = TRUE)
  r1 <- regmatches(x, m1)
  hit1 <- lengths(r1) >= 2
  pos[hit1] <- as.integer(vapply(r1[hit1], `[`, character(1), 2))

  m2 <- regexec("^\\*([0-9]+)\\*", x, perl = TRUE)
  r2 <- regmatches(x, m2)
  hit2 <- is.na(pos) & lengths(r2) >= 2
  pos[hit2] <- as.integer(vapply(r2[hit2], `[`, character(1), 2))

  m3 <- regexec("chr5 deletion:([0-9]+)-([0-9]+)", x, perl = TRUE, ignore.case = TRUE)
  r3 <- regmatches(x, m3)
  hit3 <- is.na(pos) & lengths(r3) >= 3
  pos[hit3] <- apc_genomic_to_aa(as.integer(vapply(r3[hit3], `[`, character(1), 2)))

  pos
}

standardise_sample_group <- function(type) {
  type <- tolower(trimws(as.character(type)))
  out <- rep("ACA", length(type))
  out[type %in% c("polyp", "diminutive_polyp", "large polyp", "large_polyp", "lp")] <- "polyp"
  out[type %in% c("monocryptal_adenoma", "oligocryptal_adenoma")] <- "adenoma"
  out[type %in% c("acf", "aca", "rc", "bc")] <- "ACA"
  out[type %in% c("normal", "crypt", "control")] <- "normal"
  out
}

safe_patient_filename <- function(patient) {
  gsub("[^A-Za-z0-9_.-]+", "_", patient)
}

make_apc_features <- function() {
  features <- GRanges(
    "chr5",
    IRanges(
      c(4, 453, 498, 543, 588, 633, 678, 723, 1020, 1136, 1155, 1173,
        1259, 1372, 1486, 1637, 1841, 1950, 2008, 2224, 2670),
      width = c(52, 42, 42, 42, 42, 42, 42, 42, 15, 15, 15, 15,
                20, 20, 20, 20, 20, 20, 20, 352, 174),
      names = c(
        "Oligomerization",
        rep("Armadillo repeats", 7),
        rep("15aa repeats", 4),
        rep("20aa repeats", 7),
        "Basic domain",
        "EB1-binding domain"
      )
    )
  )
  features$fill <- c(
    "#F5847A",
    rep("#EABE63", 7),
    rep("#F4F3B9", 4),
    rep("#8DD3C7", 7),
    "#8AB1C9",
    "#BD98A2"
  )
  features$color <- features$fill
  features$height <- 0.05
  features
}

draw_trackviewer_legend <- function(max_y, germline_label = NA_character_) {
  legend_labels <- c("Polyps", "Adenoma", "ACA/RC/BC", "Normal", "Germline APC")
  legend_colors <- c(group_colors["polyp"], group_colors["adenoma"], group_colors["ACA"], group_colors["normal"], group_colors["germline APC"])
  y_positions <- unit(max_y + c(0.92, 0.82, 0.72, 0.62, 0.52), "native")

  for (i in seq_along(legend_labels)) {
    grid.circle(
      x = unit(2200, "native"),
      y = y_positions[i],
      r = unit(1.9, "mm"),
      gp = gpar(col = legend_colors[i], fill = legend_colors[i])
    )
    grid.text(
      label = legend_labels[i],
      x = unit(2250, "native"),
      y = y_positions[i],
      just = "left",
      gp = gpar(fontsize = 10)
    )
  }

  if (!is.na(germline_label) && nzchar(germline_label)) {
    grid.text(
      label = germline_label,
      x = unit(2250, "native"),
      y = unit(max_y + 0.44, "native"),
      just = "left",
      gp = gpar(fontsize = 9)
    )
  }
}

plot_patient <- function(patient_id, plot_data) {
  patient_data <- plot_data[plot_data$patient == patient_id, , drop = FALSE]

  positioned <- patient_data[!is.na(patient_data$position), , drop = FALSE]
  germline_label <- unique(patient_data$mutation[patient_data$source == "germline"])
  germline_label <- germline_label[!is.na(germline_label) & germline_label != ""]
  germline_label <- if (length(germline_label) > 0) paste(germline_label, collapse = "; ") else NA_character_
  unpositioned_germline <- unique(patient_data[
    patient_data$source == "germline" & is.na(patient_data$position),
    c("mutation", "position")
  ])

  if (nrow(positioned) > 0) {
    mutation_counts <- aggregate(
      count_group_id ~ position + mutation + group,
      data = positioned,
      FUN = function(x) length(unique(x[!is.na(x) & x != ""]))
    )
    names(mutation_counts)[names(mutation_counts) == "count_group_id"] <- "count"
    mutation_counts$count[mutation_counts$group == "germline APC"] <- 1
  } else {
    mutation_counts <- data.frame(
      position = integer(),
      mutation = character(),
      group = character(),
      count = integer()
    )
  }

  max_y <- max(c(1, mutation_counts$count), na.rm = TRUE)
  y_limit <- max_y + 1
  xaxis <- seq(0, apc_length, by = 400)
  yaxis <- seq(0, ceiling(y_limit), by = 1)

  if (nrow(mutation_counts) > 0) {
    mutation_counts$plot_label <- mutation_counts$mutation
    mutation_counts$plot_label[mutation_counts$group == "germline APC"] <- ""
    sample_gr <- GRanges(
      "chr5",
      IRanges(mutation_counts$position, width = 1, names = mutation_counts$plot_label)
    )
    sample_gr$color <- group_colors[mutation_counts$group]
    sample_gr$border <- group_colors[mutation_counts$group]
    sample_gr$score <- mutation_counts$count
  } else {
    # lolliplot expects at least one point. Use an invisible dummy point if the
    # patient only has a germline hit.
    sample_gr <- GRanges("chr5", IRanges(1, width = 1, names = ""))
    sample_gr$color <- "#FFFFFF00"
    sample_gr$border <- "#FFFFFF00"
    sample_gr$score <- 0
  }

  lolliplot(
    sample_gr,
    make_apc_features(),
    xaxis = xaxis,
    yaxis = yaxis,
    lollipop_style_switch_limit = 1,
    ylab = "Count",
    cex.axis = 0.5,
    cex.lab = 0.5,
    cex.main = 0.5,
    dashline.col = NA
  )

  grid.text(
    paste0(patient_id, " APC driver mutations"),
    x = 0.5,
    y = 0.98,
    just = "top",
    gp = gpar(cex = 1.35, fontface = "bold")
  )

  pushViewport(viewport(
    width = 1,
    height = 1,
    xscale = c(0, apc_length),
    yscale = c(0, y_limit)
  ))

  draw_trackviewer_legend(max_y, germline_label)
  popViewport()

  if (nrow(unpositioned_germline) > 0) {
    grid.text(
      paste0("Unpositioned germline APC event: ", paste(unpositioned_germline$mutation, collapse = "; ")),
      x = 0.02,
      y = 0.04,
      just = "left",
      gp = gpar(fontsize = 8, col = "black")
    )
  }
}

metadata <- read.csv(metadata_file, check.names = FALSE, na.strings = c("", "NA"))

# 2026-08-20: use subgroup_id, not group_id, as the lesion identifier. subgroup_id
# is lineage-resolved (PD44721I_PLP_02_1 vs _2, PD42778D_PLP_03_1 vs _2), and the
# downstream timing/20aa analyses key on that granularity -- APC_onset_and_trunk_
# length.csv is built per lineage. Reading group_id collapses independent
# APC-mutant lineages within one lesion into a single row, which silently halved
# the number of timing-constrained lesions (9 -> 5) when this script was rerun.
# Falls back to group_id where subgroup_id is blank.
if ("subgroup_id" %in% names(metadata)) {
  .sg <- as.character(metadata$subgroup_id)
  .blank <- is.na(.sg) | trimws(.sg) == ""
  .sg[.blank] <- as.character(metadata$group_id)[.blank]
  metadata$group_id <- .sg
}
cnv <- read.csv(cnv_file, check.names = FALSE, na.strings = c("", "NA"))

required_cols <- c("sample", "patient", "type", "group_id", "germlinehit", "somatic1", "somatic2")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) {
  stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
}

if ("label_id" %in% names(metadata)) {
  missing_group <- is.na(metadata$group_id) | trimws(as.character(metadata$group_id)) == ""
  use_label <- missing_group & !is.na(metadata$label_id) & trimws(as.character(metadata$label_id)) != ""
  metadata$group_id[use_label] <- metadata$label_id[use_label]
}
missing_group <- is.na(metadata$group_id) | trimws(as.character(metadata$group_id)) == ""
metadata$group_id[missing_group] <- metadata$sample[missing_group]

somatic_cols <- intersect(c("somatic1", "somatic2", "somatic3"), names(metadata))

somatic_data <- do.call(rbind, lapply(somatic_cols, function(col) {
  data.frame(
    sample = metadata$sample,
    patient = metadata$patient,
    type = metadata$type,
    group_id = metadata$group_id,
    source_column = col,
    mutation = metadata[[col]],
    source = "somatic",
    stringsAsFactors = FALSE
  )
}))

somatic_data <- somatic_data[
  !is.na(somatic_data$mutation) &
    somatic_data$mutation != "" &
    !(somatic_data$mutation %in% c("APC_intragenic_deletion")),
  ,
  drop = FALSE
]
somatic_data$group <- standardise_sample_group(somatic_data$type)
somatic_data$position <- extract_apc_aa_position(somatic_data$mutation)

cnv_apc_start <- cnv[
  cnv$chr == "chr5" &
    !is.na(cnv$start) &
    cnv$start >= apc_genomic_start &
    cnv$start <= apc_genomic_end &
    cnv$end >= apc_genomic_start,
  c("sample", "start")
]
if (nrow(cnv_apc_start) > 0) {
  cnv_apc_start <- cnv_apc_start[order(cnv_apc_start$sample, cnv_apc_start$start), ]
  cnv_apc_start <- cnv_apc_start[!duplicated(cnv_apc_start$sample), ]
  somatic_data <- merge(
    somatic_data,
    cnv_apc_start,
    by = "sample",
    all.x = TRUE,
    sort = FALSE
  )
  loh_rows <- somatic_data$mutation %in% c("Chr5LOH", "Chr5CNNLOH", "APC_gene_deletion")
  somatic_data$position[loh_rows] <- apc_genomic_to_aa(somatic_data$start[loh_rows])
  somatic_data$start <- NULL
}

germline_data <- unique(metadata[, c("patient", "group_id", "germlinehit")])
germline_data <- data.frame(
  sample = NA_character_,
  patient = germline_data$patient,
  type = NA_character_,
  group_id = germline_data$group_id,
  source_column = "germlinehit",
  mutation = germline_data$germlinehit,
  source = "germline",
  group = "germline APC",
  position = extract_apc_aa_position(germline_data$germlinehit),
  stringsAsFactors = FALSE
)
germline_data <- germline_data[!is.na(germline_data$mutation) & germline_data$mutation != "", , drop = FALSE]

plot_data <- rbind(somatic_data, germline_data)
plot_data <- plot_data[!is.na(plot_data$patient) & plot_data$patient != "", , drop = FALSE]
plot_data$count_group_id <- plot_data$group_id
plot_data$count_group_id[
  plot_data$patient == "PD42778" &
    plot_data$group_id == "PD42778D_PLP_002"
] <- "PD42778D_PLP_01"
write.csv(plot_data, out_summary_csv, row.names = FALSE)

patients <- sort(unique(plot_data$patient))

# Match old Fig1-style dimensions more closely than the previous ggplot output.
pdf(out_combined_pdf, width = 7.2, height = 4.2, onefile = TRUE, useDingbats = FALSE)
for (patient_id in patients) {
  plot_patient(patient_id, plot_data)
}
dev.off()

for (patient_id in patients) {
  pdf(
    file.path(out_dir, paste0(safe_patient_filename(patient_id), "_APC_driver_lollipop.pdf")),
    width = 7.2,
    height = 4.2,
    onefile = FALSE,
    useDingbats = FALSE
  )
  plot_patient(patient_id, plot_data)
  dev.off()
}

message("Wrote combined PDF: ", out_combined_pdf)
message("Wrote per-patient PDFs under: ", out_dir)
message("Wrote plotting data: ", out_summary_csv)
