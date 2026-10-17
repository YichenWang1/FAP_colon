## Fig4: signature-decomposed phylogenies for the five most-sampled FAP patients
##
## Derived from hdp_plot_trees.R in this directory (itself copied from
## script/hdp_plot_trees.R). That script draws one PDF per patient at whatever
## size the tree needs; this one lays the five largest trees out on a single
## Nature two-column page.
##
## Differences from the parent script, all deliberate:
##   * five patients only, chosen by tip count, on one page
##   * tips carry no text; each lineage is labelled with its subgroup_id, and
##     the lineages of one lesion are joined by a link line -- see below
##   * driver mutations from fap_database are drawn in red, placed by whether
##     they are shared by the whole lineage or subclonal within it
##   * SBSC is displayed as SBS112 (the HDP refit files on disk still use the
##     old name; the databases were renamed, these outputs were not)

suppressPackageStartupMessages({
  library(RColorBrewer)
  library(ape)
  library(sigfit)
})

# --- Path configuration (as in hdp_plot_trees.R) ---
ROOT          <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon"
HDP_DIR       <- file.path(ROOT, "HDP/snv")
SBS_INPUT     <- file.path(ROOT, "LCMB_analysis/sbs_fap_control_perbranch.hdp.txt")
FAP_SNP_DIR   <- file.path(ROOT, "LCMB_analysis/FAP/snp")
DATABASE      <- file.path(ROOT, "metadata/fap_database_withexposure.csv")
OUT_DIR       <- file.path(ROOT, "github_repo/Fig4_Phylogeny")

N_PATIENTS <- 5   # the five most sampled trees

# --- Nature style -------------------------------------------------------
# Two-column width, full page height; all text 5-7 pt Helvetica.
FIG_W <- 180 / 25.4
FIG_H <- 240 / 25.4

FS_TIP    <- 5   # lesion labels and red driver text
FS_AXIS   <- 5
FS_TITLE  <- 7
FS_LEGEND <- 6

# base R sizes text as ps * cex, so set ps once and express everything as cex
PS <- 5
CEX_TIP    <- FS_TIP    / PS
CEX_AXIS   <- FS_AXIS   / PS
CEX_TITLE  <- FS_TITLE  / PS
CEX_LEGEND <- FS_LEGEND / PS

DRIVER_COL <- "red"

# One weight for every stroke in a panel -- tree edges, bar outlines, axis and
# the grey leader/link lines -- so nothing reads heavier than anything else.
# Previously these were four separate literals (0.4, 0.4, 0.25, 0.35) and the
# figure carried three different stroke widths. R's pdf device puts lwd 1 at
# 1/96 inch, so 0.2 is 0.15 pt.
LWD <- 0.2
# Outline around each branch's signature bar. Thinner than LWD: these panels
# carry 30-90 tips, so at the tree-edge weight the outlines merged into a black
# mesh that swamped the colours. 0.1 renders at 0.075 pt.
OUTLINE_LWD <- 0.1


# Driver text is drawn as plotmath so the gene symbol is italic:
# "APC Chr5CNNLOH; KRAS Y64D" -> italic(APC) Chr5CNNLOH; italic(KRAS) Y64D
# A token starting with "+" is the continuation of a same-gene label wrapped
# after its ";" (see wrap_tokens) and has no gene symbol of its own.
as_plotmath <- function(txt) {
  toks  <- strsplit(txt, "; ", fixed = TRUE)[[1]]
  parts <- vapply(toks, function(t)
    if (startsWith(t, "+")) sprintf('"%s"', substring(t, 2)) else
    # concatenation, not ~: plotmath's ~ draws its space from the Symbol font,
    # which would put a fourth, non-Helvetica face in the PDF
    sprintf('italic("%s")*" %s"', sub(" .*$", "", t), sub("^[^ ]+ ", "", t)), "")
  parse(text = paste(parts, collapse = '*"; "*'))[[1]]
}

# --- Signature palette (unchanged from hdp_plot_trees.R) ----------------
data("cosmic_signatures_v2")

sig_order <- 1:8
names(sig_order) <- c("SBS1","SBS5","SBS18","SBS88","SBS89","SBSC","SBS17b","SBS93")

getPalette <- colorRampPalette(brewer.pal(9, "Set3"))
all_cols <- getPalette(9)[c(1:3,6:8)]
all_cols <- c(all_cols, "magenta", "firebrick")
names(all_cols) <- names(sig_order)

# display name only; the refit files on disk are still keyed on SBSC
SIG_DISPLAY <- c(SBSC = "SBS112")
display_sig <- function(x) ifelse(x %in% names(SIG_DISPLAY), SIG_DISPLAY[x], x)

database <- read.csv(DATABASE, stringsAsFactors = FALSE, na.strings = c("", "NA"))

# Panel title "PD40734, 23, Male": age and sex are per patient, so every row of
# the patient must agree; stop rather than silently pick one.
patient_title <- function(patient) {
  db  <- database[database$patient == patient, ]
  age <- unique(na.omit(db$age))
  sex <- unique(na.omit(db$sex.x))
  if (length(age) != 1 || length(sex) != 1)
    stop("expected one age and one sex for ", patient, "; got age = ",
         paste(age, collapse = "/"), ", sex = ", paste(sex, collapse = "/"))
  paste0(patient, ", ", age, ", ", sex)
}

# --- Drivers ------------------------------------------------------------
# Every somatic driver column in fap_database. germlinehit is excluded: it is
# the germline APC allele, shared by every crypt in the patient, so it carries
# no lineage information here.
# Regenerative crypts are labelled only when the lineage carries a somatic APC
# hit -- across the five patients that is one lineage, PD40734I_RC_08
# (T1556Nfs*3). The other 30 RC lineages are single crypts with no second hit
# and their labels crowded out the adenomatous ones.
NO_LABEL_TYPES <- c("rc")

# Clonal means clonal: a driver goes on the left only if every crypt of its
# lineage carries it.
CLONAL_FRAC <- 1

# ...with one exception. The somatic APC second hit is what defines the
# adenomatous lineage, so the crypts carrying it ARE a lineage even when they
# are a subset of the subgroup -- typically an ACA in which only some crypts
# went on to lose the second allele (APC R876* in 4 of 9 crypts of
# PD42778D_ACA_02, APC I231Yfs*21 in 4 of 6 of PD42778E_ACA_02). These are
# placed on the left, at the MRCA of the crypts that carry them. Drivers that
# arrive later within such a lineage (TP53, PIK3CA, BRAF, ALK, ERBB3, FBXW7)
# are subclonal to it and stay on the right.
LINEAGE_DEFINING <- "APC"

# Putative drivers (Methods) are flagged in somatic_drivers.csv's
# is.putative_driver column by script/annotate_putative_driver.R:
#   (1) a hotspot in a top-20 COSMIC colorectal cancer gene or a putative driver
#       of normal colorectal crypts (Lee-Six et al. 2019), or
#   (2) a biallelic truncation of a tumour suppressor (second hit: LOH/CNNLOH of
#       the gene, a second truncation, or an X-linked gene in a male),
# and a row is flagged if any of its crypts qualifies. Every crypt listed on a
# flagged row carries the mutation, so the label goes on all of them.
# APC is read from fap_database (somatic1/somatic2, incl. chr5 LOH) exactly as
# before: every somatic APC hit sits on the germline APC allele.
# LOH events come from fap_database (the gene columns, e.g. TP53
# "R196*;Chr17LOH", and additional_mutations, e.g. "SOX9 Chr17CNNLOH"). An LOH
# is drawn only in a lineage that also carries a putative driver of the same
# gene, so the 17p LOH that precedes TP53 R196* in PD42778J_PLP_09 keeps its
# own stem label.
DRIVERS_CSV <- file.path(ROOT, "metadata/somatic_drivers.csv")

drv <- read.csv(DRIVERS_CSV, stringsAsFactors = FALSE, na.strings = c("", "NA"))
stopifnot("is.putative_driver" %in% names(drv))
drv <- drv[drv$Gene != "APC" & drv$is.putative_driver == "Y", ]
drv_long <- do.call(rbind, lapply(seq_len(nrow(drv)), function(i)
  data.frame(sample = trimws(strsplit(drv$Sample[i], ";", fixed = TRUE)[[1]]),
             gene = drv$Gene[i], mut = drv$Mutation[i], stringsAsFactors = FALSE)))
DRIVER_GENES <- unique(drv$Gene)

# LOH events per crypt, as "GENE ChrNNLOH" tokens
crypt_loh <- function(row) {
  toks <- c()
  for (gene in intersect(DRIVER_GENES, names(row))) {
    v <- row[[gene]]
    if (!is.na(v) && nzchar(trimws(v))) {
      ev <- trimws(strsplit(v, ";", fixed = TRUE)[[1]])
      loh <- grep("LOH", ev, value = TRUE)
      if (length(loh)) toks <- c(toks, paste0(gene, " ", loh))
    }
  }
  v <- row[["additional_mutations"]]
  if (!is.na(v) && nzchar(trimws(v))) {
    ev <- trimws(strsplit(v, ";", fixed = TRUE)[[1]])
    toks <- c(toks, grep("^[A-Za-z0-9]+ Chr[0-9XY]+[A-Za-z]*LOH$", ev, value = TRUE))
  }
  unique(toks)
}

# Putative-driver point mutations (non-APC) carried by one crypt
crypt_point_drivers <- function(row) drv_long[drv_long$sample == row[["sample"]], ]

# genes with a qualifying mutation anywhere in each lineage (subgroup_id); an
# LOH is shown only in such a lineage
lineage_genes <- local({
  out <- list()
  for (i in seq_len(nrow(database))) {
    g <- crypt_point_drivers(database[i, ])$gene
    sg <- database$subgroup_id[i]
    if (length(g) && !is.na(sg)) {
      out[[sg]] <- union(out[[sg]], g)
    }
  }
  out
})

# One crypt -> character vector of driver tokens, e.g. c("APC V1414Pfs*3", "KRAS G12D")
crypt_drivers <- function(row) {
  toks <- c()
  for (col in c("somatic1", "somatic2")) {
    v <- row[[col]]
    if (!is.na(v) && nzchar(trimws(v))) toks <- c(toks, paste0("APC ", trimws(v)))
  }
  pm <- crypt_point_drivers(row)
  if (nrow(pm)) toks <- c(toks, paste0(pm$gene, " ", sub("^p\\.", "", pm$mut)))
  loh <- crypt_loh(row)
  sg <- row[["subgroup_id"]]
  sg_genes <- if (is.na(sg)) character(0) else lineage_genes[[sg]]
  toks <- c(toks, loh[sub(" .*$", "", loh) %in% sg_genes])
  unique(toks)
}

# Several events of the SAME gene on one label are written with the gene named
# once, as fap_database stores them: "TP53 R196*;Chr17LOH", not "TP53 R196*;
# TP53 Chr17LOH". Placement is unaffected -- tokens are still grouped by which
# crypts carry them, so an LOH whose carrier set differs from its mutation's
# keeps its own label on its own branch.
merge_gene_tokens <- function(toks) {
  gene <- sub(" .*$", "", toks)
  ev   <- sub("^[^ ]+ ", "", toks)
  vapply(unique(gene), function(g) paste0(g, " ", paste(ev[gene == g], collapse = ";")),
         "", USE.NAMES = FALSE)
}

# --- Early-onset lineages ----------------------------------------------
# Upper-limit age at somatic APC inactivation per lineage, from APC_timing.Rmd
# (trunk SBS5 / patient SBS5 rate). Lineages initiated within the first decade
# of life get the estimate next to their label, e.g. "PD40734A_PLP_06 (by age 2.3)".
ONSET_CSV <- file.path(ROOT, "github_repo/data/phylogenetic_tree/APC_onset_and_trunk_length.csv")
EARLY_ONSET_MAX <- 10   # years
onset <- read.csv(ONSET_CSV, stringsAsFactors = FALSE)
# Only trunks carrying a somatic driver (trunk_somatic_driver counts them) date
# an APC inactivation; PD40734P_PLP_03 has none (trunk_somatic_driver = 0).
onset <- onset[onset$trunk_somatic_driver > 0, ]
onset_by_lineage <- tapply(onset$age_onset_upperlim, onset$lineage, function(x) {
  stopifnot(length(unique(x)) == 1)   # one estimate per lineage
  x[1]
})
early_onset <- onset_by_lineage[onset_by_lineage < EARLY_ONSET_MAX]
cat("Early-onset lineages (onset <", EARLY_ONSET_MAX, "years):\n")
print(round(sort(early_onset), 1))

lineage_label <- function(lineage, apc_run) {
  if (apc_run && lineage %in% names(early_onset))
    sprintf("%s (by age %.1f)", lineage, early_onset[[lineage]])
  else lineage
}

# --- Exposures (unchanged from hdp_plot_trees.R) ------------------------
build_exposure_matrix <- function(samples, sig_order) {
  mat <- data.frame(matrix(ncol = length(samples), nrow = length(sig_order)))
  rownames(mat) <- names(sig_order)
  colnames(mat) <- samples
  for (sample in samples) {
    f <- file.path(HDP_DIR, "refit", paste0("fit_branch_", sample, ".txt"))
    if (!file.exists(f)) { mat[, sample] <- 0; next }
    exposures <- read.table(f, check.names = FALSE)
    for (sig in names(sig_order)) {
      mat[sig, sample] <- if (sig %in% rownames(exposures)) exposures[sig, sample] else 0
    }
  }
  mat[rowSums(mat) > 0, , drop = FALSE]
}

# --- Pick the five most sampled patients --------------------------------
hdp_counts <- read.table(SBS_INPUT, check.names = FALSE)
hdp_counts <- hdp_counts[rowSums(hdp_counts) > 50, ]
colnames(hdp_counts) <- colnames(cosmic_signatures_v2)

tree_path <- function(p) file.path(FAP_SNP_DIR, p, paste0(p, "_snp_tree_with_branch_length.tree"))

all_fap <- unique(database$patient)
n_tips  <- sapply(all_fap, function(p) {
  if (!file.exists(tree_path(p))) return(NA_integer_)
  length(read.tree(tree_path(p))$tip.label)
})
n_tips   <- n_tips[!is.na(n_tips)]
patients <- names(sort(n_tips, decreasing = TRUE))[seq_len(N_PATIENTS)]
cat("Five most sampled patients:\n")
print(sort(n_tips, decreasing = TRUE)[seq_len(N_PATIENTS)])

# --- Per-patient annotation -------------------------------------------
# Lineage = subgroup_id, not group_id: a group_id can hold several independent
# APC lineages (PD40734C_PLP_07_1/_2/_3) plus its trapped bystander crypts, and
# a driver shared by one of those is not shared by the others. The label the
# reader sees is the subgroup_id, so each lineage is named separately and a
# lesion's trapped bystander crypts get their own label. The lesion is put back
# together by the LINK lines, which join the APC-mutant subgroup labels sharing
# a group_id -- i.e. they mark polyclonality directly, and unlike the previous
# rule they do not depend on whether a lesion's crypts happen to be contiguous
# in the tree (PD42778D_ACA_07 and PD42778D_PLP_03 are polyclonal but
# contiguous, and were silently unmarked).
annotate <- function(patient) {
  tree <- read.tree(tree_path(patient))
  db   <- database[database$patient == patient, ]
  rownames(db) <- db$sample
  db <- db[tree$tip.label, ]        # tree order

  is_lesion <- !is.na(db$type) & db$type != "normal" &
               !is.na(db$group_id) & db$group_id != "normal"

  # group_id labels: one per group, at the mid-row of its crypts. Regenerative
  # crypts are not labelled -- they are single-crypt lineages, so in PD40734
  # alone they contributed 22 labels that crowded out the adenomatous ones.
  # Normal crypts are already excluded by is_lesion.
  has_apc  <- (!is.na(db$somatic1) & nzchar(trimws(db$somatic1))) |
              (!is.na(db$somatic2) & nzchar(trimws(db$somatic2)))
  lab_key <- db$subgroup_id            # one label per lineage
  apc_key <- unique(lab_key[is_lesion & has_apc])
  is_label <- is_lesion & (!(!is.na(db$type) & db$type %in% NO_LABEL_TYPES) |
                           lab_key %in% apc_key)
  # 21 of the 63 groups are not a contiguous block of the tree (a group's
  # _trapped bystander crypts sit far from its adenomatous ones). Label each
  # contiguous run separately, repeating the group_id, rather than labelling
  # one run and leaving the others as bare ticks.
  tip_y <- node.height(tree)[seq_along(tree$tip.label)]
  labs  <- list()
  for (g in unique(lab_key[is_label])) {
    idx <- which(lab_key == g & is_label)
    idx <- idx[order(tip_y[idx])]
    run <- cumsum(c(1, diff(tip_y[idx]) > 1.5))
    for (r in unique(run)) {
      tips_r <- idx[run == r]
      # apc: does this run's crypts carry a somatic APC mutation? Only such
      # runs are joined by a link line below.
      labs[[length(labs) + 1]] <- list(text = lineage_label(g, any(has_apc[tips_r])),
                                       tips = tips_r,
                                       apc = any(has_apc[tips_r]),
                                       group = db$group_id[tips_r[1]])
    }
  }

  # drivers, keyed on subgroup_id
  drv_shared <- list()   # one entry per subgroup with a fully shared driver set
  drv_tip    <- list()   # one entry per crypt carrying a subclonal driver
  # Lesion crypts are grouped by subgroup_id. Normal crypts are taken one at
  # a time: they all share subgroup_id "normal", so grouping them would merge
  # unrelated crypts into a single lineage, and each is an independent lineage
  # anyway. Before this they were skipped entirely, because is_lesion excludes
  # type == "normal", so a driver in a normal crypt was never labelled
  # (FBXW7 R224* in PD42778g_lo0016, APC R213* in PD44897b_lo0001).
  normal_with_driver <- which(!is_lesion &
                              vapply(seq_len(nrow(db)),
                                     function(i) length(crypt_drivers(db[i, ])) > 0, TRUE))
  idx_groups <- c(lapply(unique(db$subgroup_id[is_lesion]),
                         function(sg) which(db$subgroup_id == sg & is_lesion)),
                  lapply(normal_with_driver, function(i) i))
  for (idx in idx_groups) {
    toks <- lapply(idx, function(i) crypt_drivers(db[i, ]))
    all_toks <- unique(unlist(toks))
    if (length(all_toks) == 0) next
    # clonal for the lineage: carried by CLONAL_FRAC of its crypts, or
    # lineage-defining (the somatic APC hit), which makes its carriers a
    # lineage in their own right. A one-crypt lineage is trivially clonal.
    frac   <- sapply(all_toks, function(t) mean(sapply(toks, function(x) t %in% x)))
    is_ld  <- sub(" .*$", "", all_toks) %in% LINEAGE_DEFINING
    shared <- all_toks[frac >= CLONAL_FRAC | is_ld]
    # Group the clonal drivers by which crypts carry them: two APC hits with
    # different carrier sets are different sub-lineages and must not share a
    # label, while drivers with the same carriers do.
    if (length(shared) > 0) {
      key <- vapply(shared, function(t)
        paste(which(vapply(toks, function(x) t %in% x, TRUE)), collapse = ","), "")
      for (k in unique(key)) {
        carriers <- idx[as.integer(strsplit(k, ",", fixed = TRUE)[[1]])]
        drv_shared[[length(drv_shared) + 1]] <-
          list(text = paste(merge_gene_tokens(shared[key == k]), collapse = "; "), tips = carriers)
      }
    }
    # Subclonal drivers, grouped by carrier set exactly as the clonal ones
    # above: one label per set of crypts, placed at their common ancestor.
    # Grouping on each crypt's joined driver text instead split nested events
    # into two labels -- in PD42778J_PLP_09 the 17p LOH (4 crypts) and TP53
    # R196* (3 of those 4) came out as "TP53 R196*; TP53 Chr17LOH" plus a
    # second "TP53 Chr17LOH", although the LOH sits on the stem of all four.
    sub <- setdiff(all_toks, shared)
    if (length(sub) > 0) {
      key <- vapply(sub, function(t)
        paste(which(vapply(toks, function(x) t %in% x, TRUE)), collapse = ","), "")
      for (k in unique(key)) {
        carriers <- idx[as.integer(strsplit(k, ",", fixed = TRUE)[[1]])]
        drv_tip[[length(drv_tip) + 1]] <-
          list(text = paste(merge_gene_tokens(sub[key == k]), collapse = "; "), tips = carriers)
      }
    }
  }

  samples <- grep(paste0("^", patient, "_"), rownames(hdp_counts), value = TRUE)
  expo    <- build_exposure_matrix(samples, sig_order)

  list(tree = tree, expo = expo, labs = labs,
       drv_shared = drv_shared, drv_tip = drv_tip)
}

ann <- setNames(lapply(patients, annotate), patients)

# every driver label drawn, so the rule's effect can be checked against the table
cat("\ndriver labels drawn (clonal | subclonal):\n")
for (p in patients) {
  for (D in ann[[p]]$drv_shared) cat(sprintf("  %-9s clonal     %s\n", p, D$text))
  for (D in ann[[p]]$drv_tip)    cat(sprintf("  %-9s subclonal  %s\n", p, D$text))
}

# tally, so a missing label is visible here rather than only in the figure
cat("\nlabels / clonal drivers / subclonal drivers per patient:\n")
for (p in patients)
  cat(sprintf("  %-9s %2d labels  %2d clonal  %2d subclonal\n", p,
              length(ann[[p]]$labs), length(ann[[p]]$drv_shared), length(ann[[p]]$drv_tip)))

# --- One panel ----------------------------------------------------------
# Text widths must be known before x.lim is chosen, but plot.new() consumes a
# page on the output device, so measure once up front on a throwaway device.
# Long driver labels (PD42778's "APC E1408*; APC L304Qfs*31; FBXW7 R479Q;
# PIK3CA E542K") are wider than the column they have to live in, and shrinking
# the type below 5 pt is not an option, so pack their tokens onto several lines.
wrap_tokens <- function(txt, max_in) {
  toks <- strsplit(txt, "; ", fixed = TRUE)[[1]]
  # A same-gene token wider than the column on its own ("APC E1408*;L304Qfs*31")
  # is split after each ";": "APC E1408*;" then "+L304Qfs*31" on the next line.
  tok_w <- function(t) strwidth(as_plotmath(t), units = "inches", cex = CEX_TIP)
  toks <- unlist(lapply(toks, function(t) {
    if (!grepl(";", t, fixed = TRUE) || tok_w(t) <= max_in) return(t)
    ev <- strsplit(sub("^[^ ]+ ", "", t), ";", fixed = TRUE)[[1]]
    n <- length(ev)
    c(paste0(sub(" .*$", "", t), " ", ev[1], ";"),
      if (n > 2) paste0("+", ev[2:(n - 1)], ";"),
      paste0("+", ev[n]))
  }))
  lines <- character(0); cur <- character(0)
  for (t in toks) {
    cand <- c(cur, t)
    too_wide <- strwidth(as_plotmath(paste(cand, collapse = "; ")),
                         units = "inches", cex = CEX_TIP) > max_in
    if (length(cur) && too_wide) {
      lines <- c(lines, paste(cur, collapse = "; ")); cur <- t
    } else cur <- cand
  }
  c(lines, paste(cur, collapse = "; "))
}

# The tree must keep at least this share of the panel; the clonal column is
# wrapped harder until it does. Without it PD42778's long clonal labels left
# the tree only 30% of the panel and squashed every branch.
TREE_MIN_FRAC <- 0.62

measure <- function(a, pin_w, shr_cap = Inf) {
  wid <- function(ll) max(vapply(unlist(ll), function(x)
    strwidth(as_plotmath(x), units = "inches", cex = CEX_TIP), 0), 0)
  lab <- max(strwidth(vapply(a$labs, `[[`, "", "text"),
                      units = "inches", cex = CEX_TIP), 0)
  gap <- strwidth("nn", units = "inches", cex = CEX_TIP)
  # subclonal labels are overlaid on the tree, so their width is not budgeted
  drv <- lapply(a$drv_tip, function(D) wrap_tokens(D$text, 0.30 * pin_w))
  for (frac in seq(0.30, 0.06, by = -0.01)) {
    shr  <- lapply(a$drv_shared, function(D) wrap_tokens(D$text, frac * pin_w))
    side <- wid(shr) + lab + 3 * gap
    if (side <= (1 - TREE_MIN_FRAC) * pin_w && wid(shr) <= shr_cap) break
  }
  g  <- vapply(a$labs, `[[`, "", "text")
  gi <- vapply(a$labs, function(l) l$group, "")
  ga <- vapply(a$labs, function(l) isTRUE(l$apc), TRUE)
  n_lk <- sum(vapply(unique(gi[ga]), function(x)
                length(unique(g[gi == x & ga])) > 1, TRUE))
  list(lab = lab, drv = wid(drv), shr = wid(shr), gap = gap,
       n_linked = n_lk,
       shr_lines = shr, drv_lines = drv,
       shr_nl = max(lengths(shr), 1), drv_nl = max(lengths(drv), 1))
}

# Draw one multi-line label, centred on y
text_block <- function(x, y, lines, adj_x, d, col) {
  n <- length(lines)
  for (j in seq_len(n))
    text(x, y + ((n - 1) / 2 - (j - 1)) * d, as_plotmath(lines[j]),
         adj = c(adj_x, 0.5), cex = CEX_TIP, col = col)
}

# Push labels apart just enough that they do not collide, keeping them as close
# to their true row as possible, and never leaving the panel.
spread <- function(y, d, lo, hi) {
  if (length(y) == 0) return(y)
  o <- order(y); v <- y[o]
  for (i in seq_along(v)[-1]) v[i] <- max(v[i], v[i - 1] + d)          # push up
  if (v[length(v)] > hi) {
    v <- v - (v[length(v)] - hi)
    for (i in rev(seq_along(v))[-1]) v[i] <- min(v[i], v[i + 1] - d)   # push back down
    v <- pmax(v, lo)
  }
  out <- numeric(length(y)); out[o] <- v; out
}

# txt is a list, so it can hold plotmath expressions as well as strings
draw_column <- function(txt, y_true, x_text, adj_x, col, d, lo, hi, x_leader) {
  if (length(txt) == 0) return(invisible())
  y <- spread(y_true, d, lo, hi)
  moved <- abs(y - y_true) > d / 4
  if (any(moved))
    segments(x_leader, y_true[moved], x_text, y[moved], col = "black", lwd = LWD)
  for (i in seq_along(txt))
    text(x_text, y[i], txt[[i]], adj = c(adj_x, 0.5), cex = CEX_TIP, col = col)
  y
}


panel <- function(a, patient, w, pin_w) {
  tree <- a$tree; expo <- a$expo

  # Four columns across the panel: clonal drivers to the LEFT of the whole
  # tree, the tree, subclonal drivers immediately right of the tips, then the
  # lesion labels last.
  lab_w <- w$lab; shr_w <- w$shr; gap_w <- w$gap
  # Only the clonal column and the label column are reserved. Subclonal
  # drivers are overlaid on the branches instead of getting a column of their
  # own, so the tree runs on until it almost touches the labels.
  link_w <- if (w$n_linked > 0) min(0.05 * pin_w, (w$n_linked + 0.5) * 0.35 * gap_w) else 0
  side_w <- min(shr_w + lab_w + link_w + 3 * gap_w, 0.70 * pin_w)
  xmax   <- max(node.depth.edgelength(tree))
  u      <- xmax / (pin_w - side_w)          # user units per inch
  gap    <- gap_w * u

  x_shr   <- -gap                             # clonal column, right-justified
  x_left  <- x_shr - shr_w * u
  x_lab   <- xmax + gap                       # label column, left-justified
  x_right <- x_lab + lab_w * u + 0.25 * gap + link_w * u

  # y.lim: plot.phylo defaults to c(1, Ntip), putting the first and last tips
  # exactly on the floor and ceiling of the plot region -- so axis(1) drew its
  # line straight through the bottom crypt's signature bar and that one bar read
  # thicker than the rest. Half a tip's clearance at each end fixes it.
  ntip <- length(tree$tip.label)
  plot(tree, show.tip.label = FALSE, x.lim = c(x_left, x_right),
       y.lim = c(0.5, ntip + 0.5), edge.width = LWD, no.margin = FALSE)
  lp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  xx <- lp$xx; yy <- lp$yy

  # signature bars along each branch (same construction as hdp_plot_trees.R)
  if (nrow(expo) > 0) {
    ref_sigs <- rownames(expo)
    branches <- sapply(strsplit(colnames(expo), "_"), function(x) x[[2]])
    half     <- min(0.02 * length(tree$tip.label), 0.4)
    for (k in seq_along(colnames(expo))) {
      n <- as.numeric(branches[k])
      if (is.na(n) || n > length(xx)) next
      x_start <- xx[lp$edge[match(n, lp$edge[, 2]), 1]]
      if (length(x_start) == 0 || is.na(x_start)) next
      x_intv <- xx[n] - x_start
      # Fill the signature segments unstroked, then outline the branch ONCE at
      # OUTLINE_LWD. Stroking each segment would draw every internal boundary
      # twice -- the right edge of one on the left edge of the next -- and the
      # doubled strokes read as heavier lines than the tree's own edges.
      x_bar0 <- x_start
      for (s in ref_sigs) {
        x_end <- x_start + expo[s, k] * x_intv
        rect(xleft = x_start, xright = x_end,
             ybottom = yy[n] - half, ytop = yy[n] + half,
             col = all_cols[s], border = NA)
        x_start <- x_end
      }
      rect(xleft = x_bar0, xright = x_start,
           ybottom = yy[n] - half, ytop = yy[n] + half,
           col = NA, border = "black", lwd = OUTLINE_LWD)
    }
  }

  d  <- strheight("Mg", cex = CEX_TIP) * 1.25   # minimum row pitch for text
  lo <- 1; hi <- length(tree$tip.label)

  # Clonal drivers, in a column to the left of the whole tree. The label is
  # tied to its lineage by a leader running to the START of that lineage's stem
  # (the parent of the MRCA of its crypts) -- the earliest point the event
  # could have occurred.
  if (length(a$drv_shared) > 0) {
    anch <- t(vapply(a$drv_shared, function(D) {
      nd <- if (length(D$tips) > 1) getMRCA(tree, tree$tip.label[D$tips]) else D$tips
      pa <- lp$edge[match(nd, lp$edge[, 2]), 1]        # NA if nd is the root
      c(if (is.na(pa)) xx[nd] else xx[pa],
        if (length(D$tips) > 1) yy[nd] else yy[D$tips])
    }, numeric(2)))
    y_shr <- spread(anch[, 2], d * w$shr_nl, lo, hi)
    segments(x_shr, y_shr, anch[, 1], anch[, 2], col = "black", lwd = LWD)
    for (i in seq_along(a$drv_shared))
      text_block(x_shr, y_shr[i], w$shr_lines[[i]], 1, d, DRIVER_COL)
  }

  # Subclonal drivers, drawn on top of the branch that carries them: anchored
  # at the start of the stem subtending the mutated crypt(s), the same rule the
  # clonal column's leaders use. A label that would collide with one already
  # drawn slides right along its own branch rather than moving off it.
  if (length(a$drv_tip) > 0) {
    placed <- list()
    ord <- order(vapply(a$drv_tip, function(D) mean(range(yy[D$tips])), 0))
    for (D in a$drv_tip[ord]) {
      nd <- if (length(D$tips) > 1) getMRCA(tree, tree$tip.label[D$tips]) else D$tips
      pa <- lp$edge[match(nd, lp$edge[, 2]), 1]        # NA if nd is the root
      x0 <- if (is.na(pa)) xx[nd] else xx[pa]
      ya <- if (length(D$tips) > 1) yy[nd] else yy[D$tips]

      lines_i <- w$drv_lines[[which(vapply(a$drv_tip, function(q)
        identical(q, D), TRUE))[1]]]
      tw <- max(vapply(lines_i, function(x)
        strwidth(as_plotmath(x), cex = CEX_TIP), 0))
      nl <- length(lines_i)
      xt <- min(x0, max(0, xmax - tw))
      bx <- function(x, up = 0) c(x, x + tw, ya + (0.35 + up) * d, ya + (0.35 + up + nl) * d)
      hits <- function(b) any(vapply(placed, function(q)
        b[1] < q[2] && b[2] > q[1] && b[3] < q[4] && b[4] > q[3], TRUE))
      # On a collision, first stack the label one line up at its own branch
      # start, so nested events on one lineage (a 17p LOH on the stem, then
      # TP53 R196* on the sub-clade) read left to right in tree order; slide
      # right along the branch only if a single lift does not clear it.
      up <- 0
      if (hits(bx(xt)) && !hits(bx(xt, 1))) up <- 1
      step <- 0.015 * xmax
      while (hits(bx(xt, up)) && xt + tw < xmax) xt <- xt + step
      placed[[length(placed) + 1]] <- bx(xt, up)
      for (j in seq_len(nl))
        text(xt, ya + (nl - j + 0.35 + up) * d, as_plotmath(lines_i[j]),
             adj = c(0, 0), cex = CEX_TIP, col = DRIVER_COL)
    }
  }

  # lesion label, once per contiguous run of its crypts
  x_brk  <- x_lab - 0.55 * gap
  anchor <- numeric(length(a$labs))
  span   <- vector("list", length(a$labs))
  for (i in seq_along(a$labs)) {
    ys <- yy[a$labs[[i]]$tips]
    if (length(ys) > 1) segments(x_brk, min(ys), x_brk, max(ys), col = "black", lwd = LWD)
    # a drawn dash, not pch = 45: that is a text glyph, and at 0.7 * 5 pt it
    # would put 3.5 pt type in a figure specified as 5-7 pt
    else segments(x_brk - 0.12 * gap, ys, x_brk + 0.12 * gap, ys, col = "black", lwd = LWD)
    anchor[i] <- mean(range(ys))
    span[[i]] <- range(ys)
  }
  # Bridge the gaps between the runs of one group with a dotted line at the
  # same x, so the blocks of a discontiguous group read as one lineage rather
  # than as unrelated entries that happen to repeat a name.
  y_lab <- draw_column(as.list(vapply(a$labs, `[[`, "", "text")),
                       anchor, x_lab, 0, "black", d, lo, hi, x_brk)

  # Join the blocks of a group that is split across the tree. The joins get a
  # column of their own to the RIGHT of the labels: run in the same narrow
  # gutter as the brackets they would read as one line down the panel, since
  # several split groups interleave. One vertical per group, staggered across
  # the column, with a tick into each of its labels.
  grp     <- vapply(a$labs, `[[`, "", "text")
  gid     <- vapply(a$labs, function(l) l$group, "")
  lab_apc <- vapply(a$labs, function(l) isTRUE(l$apc), TRUE)
  # Join the APC-mutant lineages of one lesion: same group_id, more than one
  # distinct subgroup label. Trapped bystander subgroups carry no APC hit and
  # are excluded, so a link always asserts "these lineages are the same lesion",
  # never "this crypt sits inside that polyp".
  cand    <- unique(gid[lab_apc])
  linked  <- cand[vapply(cand, function(g)
                    length(unique(grp[gid == g & lab_apc])) > 1, TRUE)]
  if (length(linked) > 0) {
    x0   <- x_lab + lab_w * u + 0.25 * gap
    step <- (x_right - x0) / (length(linked) + 0.5)
    for (n in seq_along(linked)) {
      k    <- which(gid == linked[n] & lab_apc)
      x_ln <- x0 + n * step
      segments(x_ln, min(y_lab[k]), x_ln, max(y_lab[k]), col = "black", lwd = LWD)
      segments(x0 - 0.15 * gap, y_lab[k], x_ln, y_lab[k], col = "black", lwd = LWD)
    }
  }

  axis(1, at = pretty(c(0, xmax), 4), cex.axis = CEX_AXIS, lwd = LWD,
       tcl = -0.15, mgp = c(3, 0.15, 0))
  mtext("Substitutions", side = 1, line = 0.9, cex = CEX_AXIS)
  mtext(patient_title(patient), side = 3, line = -0.2, adj = 0, cex = CEX_TITLE, font = 2)
}

# --- Page layout --------------------------------------------------------
# Two columns, panel heights proportional to tip count so every tree gets the
# same vertical space per crypt. Left column takes the two largest trees, right
# column the remaining three; the totals happen to balance closely.
nt   <- sapply(ann, function(a) length(a$tree$tip.label))
col1 <- patients[1:2]
col2 <- patients[3:5]

LEG_H <- 0.030   # bottom strip for the shared signature legend
TOP_H <- 0.012   # top strip so the top-row panel titles clear the page edge

col_fig <- function(pts, x0, x1) {
  h <- nt[pts] / sum(nt[pts]) * (1 - LEG_H - TOP_H)
  top <- 1 - TOP_H
  out <- list()
  for (p in pts) {
    out[[p]] <- c(x0, x1, top - h[[p]], top)
    top <- top - h[[p]]
  }
  out
}
# gutter between the columns so the left column's driver text does not run
# up against the right column's tree
GUTTER <- 0.014
figs <- c(col_fig(col1, 0.00, 0.50 - GUTTER / 2),
          col_fig(col2, 0.50 + GUTTER / 2, 1.00))

# bottom just clears the axis ticks, tick labels and "Substitutions";
# top just clears the panel title
MAI <- c(0.16, 0.02, 0.07, 0.02)

# measure on a scratch device with the same font metrics as the output
local({
  pdf(tempfile(), width = FIG_W, height = FIG_H, family = "Helvetica", pointsize = PS)
  par(ps = PS)
  plot.new()
  widths <<- lapply(patients, function(p)
    measure(ann[[p]], (figs[[p]][2] - figs[[p]][1]) * FIG_W - (MAI[2] + MAI[4])))
  names(widths) <<- patients
  # Trees in one column start at the same x: wrap the clonal-driver column of
  # every lower panel to fit the top panel's width where the tokens allow, then
  # give the whole column the widest remaining width.
  for (pts in list(col1, col2)) {
    cap <- widths[[pts[1]]]$shr
    for (p in pts[-1])
      widths[[p]] <<- measure(ann[[p]], (figs[[p]][2] - figs[[p]][1]) * FIG_W - (MAI[2] + MAI[4]),
                              shr_cap = cap)
    col_shr <- max(vapply(pts, function(p) widths[[p]]$shr, 0))
    for (p in pts) widths[[p]]$shr <<- col_shr
    cat("clonal-driver column width (in):", paste(pts, collapse = "/"), "=", round(col_shr, 3), "\n")
  }
  invisible(dev.off())
})

out_pdf <- file.path(OUT_DIR, "Fig4_trees_with_signatures.pdf")
pdf(out_pdf, width = FIG_W, height = FIG_H, family = "Helvetica", pointsize = PS)
par(ps = PS, mai = MAI, xaxs = "i", lend = 1)

first <- TRUE
for (p in patients) {
  par(fig = figs[[p]], new = !first)
  first <- FALSE
  pin_w <- (figs[[p]][2] - figs[[p]][1]) * FIG_W - (MAI[2] + MAI[4])
  panel(ann[[p]], p, widths[[p]], pin_w)
}

# shared signature legend
par(fig = c(0, 1, 0, LEG_H), new = TRUE, mai = c(0.01, 0.02, 0.01, 0.02))
plot.new()
present <- names(sig_order)[names(sig_order) %in%
                            unique(unlist(lapply(ann, function(a) rownames(a$expo))))]
legend("center", legend = display_sig(present), fill = all_cols[present],
       border = NA, bty = "n", horiz = TRUE, cex = CEX_LEGEND,
       x.intersp = 0.5, text.width = NA)
invisible(dev.off())
cat("Wrote", out_pdf, "\n")
