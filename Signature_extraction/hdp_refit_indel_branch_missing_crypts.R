# Attribute whole-crypt indel signature exposure for FAP/control crypts whose
# entire tip-to-root path has no individually HDP-fit branch (every branch on
# the path falls below the sum>30 indel threshold used by hdp_refit_indel.R /
# merge_control_expsoure.R)
#
# For each such crypt, its own 83-channel spectrum (crypt_hdp_matrix_indel_rounded.txt)
# is refit with sigfit against the signatures active (fraction > 0.05) in the
# same patient's fitted branches (HDP/indel/refit/fit_<branch>.txt), or against
# ID1/ID2/ID5 if the patient has none. Results go to the ID columns of the
# database and to HDP/indel/refit/fit_<sample>.txt.

library(sigfit)

ROOT              <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon"
HDP_DIR           <- file.path(ROOT, "HDP/indel")
COUNT_MATRIX_PATH <- file.path(HDP_DIR, "crypt_hdp_matrix_indel_rounded.txt")
REFIT_DIR         <- file.path(HDP_DIR, "refit")

# FAP database has CRLF line endings, the control database LF; both are kept
DATABASES <- list(
  list(path = file.path(ROOT, "metadata/fap_database_withexposure.csv"),  eol = '\r\n'),
  list(path = file.path(ROOT, "metadata/ctrl_database_withexposure.csv"), eol = '\n')
)

ID_COLS            <- c("ID1", "ID2", "ID5", "ID18", "ID14", "IDA", "IDB")
ACTIVE_FRAC_MIN    <- 0.05                     # a branch's signature counts as active above this
DEFAULT_SIGS       <- c("ID1", "ID2", "ID5")   # for patients with no fitted branch

# Write a csv quoting only cells that need it, with NA written as a blank cell,
# so cells this script does not touch keep their original formatting.
write_csv_minimal <- function(df, file, eol = '\r\n') {
  quote_if_needed <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ''
    needs_quote <- grepl('[",\n]', x)
    x[needs_quote] <- paste0('"', gsub('"', '""', x[needs_quote]), '"')
    x
  }
  quoted <- data.frame(lapply(df, quote_if_needed), stringsAsFactors = FALSE, check.names = FALSE)
  con <- file(file, 'wb'); on.exit(close(con))
  lines <- c(paste(quote_if_needed(colnames(df)), collapse = ','),
             do.call(paste, c(quoted, sep = ',')))
  writeLines(lines, con, sep = eol)
}

# Read a database with every column as text; only blank cells become NA.
read_database <- function(path) {
  read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,
           colClasses = 'character', na.strings = '')
}

# TRUE for crypts with no indel signature attribution (every ID column empty or "NA")
has_no_id_attribution <- function(db) {
  id_values <- as.matrix(as.data.frame(lapply(db[, ID_COLS], function(x) suppressWarnings(as.numeric(x)))))
  apply(id_values, 1, function(r) all(is.na(r)))
}

## ---------------- reference signatures ----------------
# final_sigs: one row per signature, 83 indel channels as columns
final_sigs   <- data.frame(t(read.table(file.path(HDP_DIR, "final_sigs.txt"), check.names = FALSE)),
                           check.names = FALSE)
sig_rank     <- setNames(seq_len(nrow(final_sigs)), rownames(final_sigs))   # canonical order
default_sigs <- intersect(DEFAULT_SIGS, rownames(final_sigs))

## ---------------- each crypt's own 83-channel spectrum (FAP + control) ----------------
crypt_counts <- read.table(COUNT_MATRIX_PATH, check.names = FALSE)
stopifnot(all(colnames(crypt_counts) == colnames(final_sigs)))

## ---------------- signatures active in each patient's fitted branches ----------------
# One refit directory covers FAP and control branches; files are fit_<patient>_<node>.txt
refit_files  <- list.files(REFIT_DIR, pattern = "^fit_.*\\.txt$", full.names = FALSE)
branch_ids   <- sub("^fit_", "", sub("\\.txt$", "", refit_files))
branch_patient <- sub("_[^_]+$", "", branch_ids)

active_sigs_by_patient <- list()
for (k in seq_along(refit_files)) {
  branch  <- branch_ids[k]
  patient <- branch_patient[k]
  branch_fit <- try(read.table(file.path(REFIT_DIR, refit_files[k]), check.names = FALSE), silent = TRUE)
  # skip unreadable files and files that are not single-branch fits
  if (inherits(branch_fit, "try-error") || !(branch %in% colnames(branch_fit))) next
  sig_frac <- setNames(as.numeric(branch_fit[[branch]]), rownames(branch_fit))
  active   <- names(sig_frac)[!is.na(sig_frac) & sig_frac > ACTIVE_FRAC_MIN]
  active_sigs_by_patient[[patient]] <- union(active_sigs_by_patient[[patient]], active)
}

## ---------------- refit the unattributed crypts of each database ----------------
for (db_spec in DATABASES) {
  db_path <- db_spec$path
  cat("\n==== ", db_path, " ====\n", sep = "")

  db <- read_database(db_path)
  is_missing <- has_no_id_attribution(db)
  missing_samples <- db$sample[is_missing]
  cat("Missing crypts (no branch-level indel attribution found):", length(missing_samples), "\n")
  print(table(db$patient[is_missing]))

  if (length(missing_samples) == 0) { cat("Nothing to do for", db_path, "\n"); next }

  stopifnot(all(missing_samples %in% rownames(crypt_counts)))
  patient_of_sample <- setNames(db$patient, db$sample)

  # absolute exposures per refit crypt: one row per crypt, one column per signature
  new_exposures <- data.frame(matrix(nrow = 0, ncol = length(sig_rank)))
  colnames(new_exposures) <- names(sig_rank)

  for (sample in missing_samples) {
    patient  <- patient_of_sample[[sample]]
    ref_sigs <- active_sigs_by_patient[[patient]]
    ref_sigs <- ref_sigs[ref_sigs %in% rownames(final_sigs)]

    if (length(ref_sigs) == 0) {
      ref_sigs <- default_sigs
      cat(sample, "(patient", patient, "): no other attributed crypts for this patient -> using default",
          paste(default_sigs, collapse = ","), "\n")
    } else {
      ref_sigs <- names(sort(sig_rank[ref_sigs]))
      cat(sample, "(patient", patient, "): borrowing", paste(ref_sigs, collapse = ","), "\n")
    }

    sample_counts <- crypt_counts[sample, , drop = FALSE]
    crypt_total   <- sum(as.numeric(sample_counts[1, ]))

    fit <- fit_signatures(counts = sample_counts, signatures = final_sigs[ref_sigs, , drop = FALSE],
                          iter = 5000, warmup = 2000, model = "poisson", chains = 5)

    # posterior mean fraction per signature; a single signature takes everything
    if (length(ref_sigs) > 1) {
      posterior <- retrieve_pars(fit, par = "exposures", hpd_prob = 0.95)
      sig_frac  <- setNames(as.numeric(posterior$mean[1, ref_sigs]), ref_sigs)
    } else {
      sig_frac <- setNames(1, ref_sigs)
    }

    exposure_row <- setNames(rep(0, length(sig_rank)), names(sig_rank))
    exposure_row[names(sig_frac)] <- sig_frac * crypt_total
    new_exposures[sample, ] <- exposure_row

    # Also write the fractions to HDP/indel/refit/fit_<sample>.txt, in the same
    # format as the branch fits (rows = signatures, one column named after the
    # sample), so that directory holds every indel attribution.
    frac_out <- setNames(rep(0, length(sig_rank)), names(sig_rank))
    frac_out[names(sig_frac)] <- sig_frac
    frac_df <- data.frame(x = frac_out, check.names = FALSE)
    colnames(frac_df) <- sample
    write.table(frac_df, file.path(REFIT_DIR, paste0("fit_", sample, ".txt")),
                sep = "\t", quote = FALSE)
  }

  cat("Refit", nrow(new_exposures), "/", length(missing_samples), "missing crypts.\n")

  ## ---------------- write the exposures into the database ----------------
  db_rows <- match(rownames(new_exposures), db$sample)
  stopifnot(!any(is.na(db_rows)))
  for (col in colnames(new_exposures)) {
    if (!(col %in% colnames(db))) next
    db[[col]][db_rows] <- as.character(new_exposures[[col]])
  }
  write_csv_minimal(db, db_path, eol = db_spec$eol)
  cat("Wrote", db_path, "\n")

  ## ---------------- check: refit exposures sum to each crypt's indel burden ----------------
  db_after <- read_database(db_path)
  rownames(db_after) <- db_after$sample
  refit_rows <- db_after[rownames(new_exposures), ]
  burden_gap <- as.numeric(refit_rows$indeltotal_corr) - rowSums(new_exposures)
  cat("Max |indeltotal_corr - sum(newly-attributed ID cols)| across", length(burden_gap), "crypts:",
      max(abs(burden_gap), na.rm = TRUE), "\n")
  remaining_na <- sum(has_no_id_attribution(db_after))
  cat("Remaining crypts with no indel signature attribution at all:", remaining_na, "\n")
}
