library(RColorBrewer)
library(ape)
library(sigfit)

# --- Path configuration ---
REPO          <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/github_repo"
HDP_DIR       <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/HDP/snv"
SBS_INPUT     <- file.path(REPO, "data/somatic_mutation/sbs_fap_control_perbranch.hdp.txt")
TREE_DIR      <- file.path(REPO, "data/phylogenetic_tree/tree_structure")
DATABASE      <- file.path(REPO, "data/fap_database_withexposure.csv")
FAP_PATIENTS  <- file.path(REPO, "data/phylogenetic_tree/fap_patients.txt")
NORM_PATIENTS <- file.path(REPO, "data/phylogenetic_tree/ctrl_patients.txt")

dir.create(file.path(HDP_DIR, "tree_plots", "fap"),     showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(HDP_DIR, "tree_plots", "control"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(HDP_DIR, "refit"),                 showWarnings=FALSE, recursive=TRUE)

data("cosmic_signatures_v2")

sig_order <- 1:8
names(sig_order) <- c("SBS1","SBS5","SBS18","SBS88","SBS89","SBSC","SBS17b","SBS93")

getPalette <- colorRampPalette(brewer.pal(9, "Set3"))
all_cols <- getPalette(9)[c(1:3,6:8)]
all_cols <- c(all_cols, "magenta", "firebrick")
names(all_cols) <- names(sig_order)

fap_patients  <- read.table(FAP_PATIENTS,  header=FALSE)$V1
norm_patients <- read.table(NORM_PATIENTS, header=FALSE)$V1
database      <- read.csv(DATABASE, stringsAsFactors=FALSE)

get_group <- function(sample_names, patient_list) {
  sample_names[sapply(strsplit(sample_names, '_'), function(x) x[1]) %in% patient_list]
}

build_exposure_matrix <- function(samples, sig_order) {
  mat <- data.frame(matrix(ncol=length(samples), nrow=length(sig_order)))
  rownames(mat) <- names(sig_order)
  colnames(mat) <- samples
  for (sample in samples) {
    f <- file.path(HDP_DIR, "refit", paste0("fit_branch_", sample, ".txt"))
    if (!file.exists(f)) { mat[, sample] <- 0; next }
    exposures <- read.table(f, check.names=FALSE)
    for (sig in names(sig_order)) {
      mat[sig, sample] <- if (sig %in% rownames(exposures)) exposures[sig, sample] else 0
    }
  }
  mat[rowSums(mat) > 0, , drop=FALSE]
}


# Style follows SuppFig_remaining_patient_trees.pdf (plot_remaining_trees_supp.R)
PT     <- 6                      # pdf pointsize; text sizes are cex = pt / PT
FS_TIP <- 5 ; FS_AXIS <- 5 ; FS_TITLE <- 7 ; FS_LEGEND <- 6
LWD    <- 0.2                    # every stroke: edges, bar outlines, axis
# Panel width: 60 mm, one panel of the supplementary figure (180 mm / 3).
# The five Figure 4 patients have the largest trees (46-127 crypts), so they
# get a single Nature column, 89 mm, to leave room for the branches.
FIG4_PATIENTS <- c("PD40734", "PD42778", "PD45542", "PD44721", "PD44720")
PANEL_W_MM    <- 60
FIG4_W_MM     <- 89
MIN_H_MM   <- 34                 # the supplementary figure's row height, for small trees
MM_PER_TIP <- 2.2                # larger trees: room for 5 pt tip labels
LEG_H_IN   <- 0.25               # legend strip at the bottom
MAI        <- c(0.20, 0.03, 0.09, 0.03)

SIG_DISPLAY <- c(SBSC = "SBS112")
display_sig <- function(x) ifelse(x %in% names(SIG_DISPLAY), SIG_DISPLAY[x], x)

plot_tree_sigs <- function(tree, exposure_mat, cols, patient, out_pdf) {
  ref_sigs <- rownames(exposure_mat)
  branches <- sapply(strsplit(colnames(exposure_mat), '_'), function(x) x[[2]])
  labs     <- tree$tip.label
  n_tips   <- length(labs)
  fig_w    <- (if (patient %in% FIG4_PATIENTS) FIG4_W_MM else PANEL_W_MM) / 25.4
  fig_h    <- max(MIN_H_MM, n_tips * MM_PER_TIP) / 25.4 + LEG_H_IN

  pdf(out_pdf, width = fig_w, height = fig_h, pointsize = PT)
  layout(matrix(1:2, ncol = 1), heights = c(fig_h - LEG_H_IN, LEG_H_IN))
  par(mai = MAI, xaxs = "i", lend = 1)

  # room on the right for the sample-id tip labels
  lab_w <- max(strwidth(labs, units = "inches", cex = FS_TIP / PT)) * 1.10
  xmax  <- max(node.depth.edgelength(tree))
  u     <- xmax / (par("pin")[1] - lab_w)        # user units per inch

  # Half a tip of clearance top and bottom so the axis does not run through the
  # last bar; x.lim starts just below 0 so the root connector is not clipped;
  # underscore = TRUE keeps the "_" in sample ids.
  plot(tree, show.tip.label = TRUE, cex = FS_TIP / PT, label.offset = 0.01 * xmax,
       x.lim = c(-0.004 * xmax, xmax + lab_w * u), y.lim = c(0.5, n_tips + 0.5),
       edge.width = LWD, no.margin = FALSE, underscore = TRUE)
  lp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  xx <- lp$xx; yy <- lp$yy

  # signature bars along each branch: segments unstroked, then one outline
  half <- min(0.02 * n_tips, 0.4)
  for (k in seq_along(colnames(exposure_mat))) {
    n <- as.numeric(branches[k])
    if (is.na(n) || n > length(xx)) next
    x_start <- xx[lp$edge[match(n, lp$edge[, 2]), 1]]
    if (length(x_start) == 0 || is.na(x_start)) next
    x_intv <- xx[n] - x_start
    x_bar0 <- x_start
    for (s in ref_sigs) {
      x_end <- x_start + exposure_mat[s, k] * x_intv
      rect(x_start, yy[n] - half, x_end, yy[n] + half, col = cols[s], border = NA)
      x_start <- x_end
    }
    rect(x_bar0, yy[n] - half, x_start, yy[n] + half, col = NA, border = "black", lwd = LWD)
  }

  axisPhylo(side = 1, backward = FALSE, cex.axis = FS_AXIS / PT, lwd = LWD,
            tcl = -0.15, mgp = c(3, 0.15, 0))
  mtext("Substitutions", side = 1, line = 0.9, cex = FS_AXIS / PT)
  mtext(patient, side = 3, line = -0.2, adj = 0, cex = FS_TITLE / PT, font = 2)

  par(mai = c(0.01, 0.02, 0.01, 0.02))
  plot.new()
  legend("center", legend = display_sig(ref_sigs), fill = cols[ref_sigs], border = NA,
         bty = "n", ncol = 4, cex = FS_LEGEND / PT, x.intersp = 0.5)
  invisible(dev.off())
}

# -----------------Plot signatures onto phylogenetic trees-----------------
hdp_counts <- read.table(SBS_INPUT, check.names=FALSE)
hdp_counts <- hdp_counts[rowSums(hdp_counts) > 50, ]
colnames(hdp_counts) <- colnames(cosmic_signatures_v2)

# FAP
cat("== FAP trees ==\n")
fap_counts <- hdp_counts[get_group(rownames(hdp_counts), fap_patients), , drop=FALSE]
fap_pts    <- fap_patients[fap_patients %in% unique(sapply(strsplit(rownames(fap_counts),'_'), function(x) x[1]))]

for (patient in fap_pts) {
  cat("Processing", patient, "\n")
  tree <- try(read.tree(file.path(TREE_DIR,
                                  paste0(patient, "_snp_tree_with_branch_length.tree"))))
  if (inherits(tree, "try-error")) { message("  Tree not found, skipping"); next }

  samples <- grep(paste0(patient, "_"), rownames(fap_counts), value=TRUE)
  if (length(samples) == 0) { message("  No samples in HDP input, skipping"); next }

  expo <- build_exposure_matrix(samples, sig_order)
  if (nrow(expo) == 0) { message("  No exposure data, skipping"); next }

  write.table(expo,
              file.path(HDP_DIR, "refit", paste0("fit_branch_", patient, ".txt")),
              sep="\t", quote=FALSE)

  plot_tree_sigs(tree, expo, all_cols, patient,
                 file.path(HDP_DIR, "tree_plots", "fap",
                           paste0(patient, "_tree_with_sigs.pdf")))
}

# Control
cat("\n== Control trees ==\n")
norm_counts <- hdp_counts[get_group(rownames(hdp_counts), norm_patients), , drop=FALSE]
norm_pts    <- norm_patients[norm_patients %in% unique(sapply(strsplit(rownames(norm_counts),'_'), function(x) x[1]))]

for (patient in norm_pts) {
  cat("Processing", patient, "\n")
  tree <- try(read.tree(file.path(TREE_DIR,
                                  paste0(patient, "_snp_tree_with_branch_length.tree"))))
  if (inherits(tree, "try-error")) { message("  Tree not found, skipping"); next }

  samples <- grep(paste0(patient, "_"), rownames(norm_counts), value=TRUE)
  if (length(samples) == 0) { message("  No samples in HDP input, skipping"); next }

  expo <- build_exposure_matrix(samples, sig_order)
  if (nrow(expo) == 0) { message("  No exposure data, skipping"); next }

  write.table(expo,
              file.path(HDP_DIR, "refit", paste0("fit_branch_", patient, "_control.txt")),
              sep="\t", quote=FALSE)

  plot_tree_sigs(tree, expo, all_cols, patient,
                 file.path(HDP_DIR, "tree_plots", "control",
                           paste0("control_", patient, "_tree_with_ref.pdf")))
}

# --------------- Build combined exposure matrix ---------------
cat("\n== Building exposure matrix ==\n")
exposure_matrix_all <- data.frame(matrix(nrow=0, ncol=length(sig_order)))
colnames(exposure_matrix_all) <- names(sig_order)

for (patient in c(fap_pts, norm_pts)) {
  cohort_counts <- if (patient %in% fap_pts) fap_counts else norm_counts
  samples <- grep(paste0(patient, "_"), rownames(cohort_counts), value=TRUE)
  if (length(samples) == 0) next
  mat <- data.frame(matrix(nrow=length(samples), ncol=length(sig_order)))
  colnames(mat) <- names(sig_order)
  rownames(mat) <- samples
  for (sample in samples) {
    f <- file.path(HDP_DIR, "refit", paste0("fit_branch_", sample, ".txt"))
    if (!file.exists(f)) { mat[sample, ] <- 0; next }
    exposures <- read.table(f, check.names=FALSE)
    for (sig in names(sig_order)) {
      mat[sample, sig] <- if (sig %in% rownames(exposures)) exposures[sig, sample] else 0
    }
  }
  exposure_matrix_all <- rbind(exposure_matrix_all, mat)
}

write.table(exposure_matrix_all,
            file.path(HDP_DIR, "refit", "sig_exposure.txt"),
            sep="\t", quote=FALSE)
cat("Written:", file.path(HDP_DIR, "refit", "sig_exposure.txt"), "\n")
cat("Done.\n")
