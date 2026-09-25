# Optional input preparation: install.packages("haven"), then Rscript jtpa/prepare.R.
# An optional argument supplies a local copy of the same checksum-pinned Stata file.
cky_prepare_jtpa <- function(root, source_file = NULL) {
  if (!requireNamespace("haven", quietly = TRUE)) {
    stop("Install the optional Stata reader first: install.packages('haven')")
  }
  source_url <- "http://fmwww.bc.edu/repec/bocode/j/jtpa.dta"
  expected_sha256 <- "b2e1a9ef217fa79919f5e6c8fd9bc85fea912ec80d04a071e3936605a7fe8ec6"
  if (is.null(source_file)) {
    source_file <- tempfile("jtpa-", fileext = ".dta")
    on.exit(unlink(source_file), add = TRUE)
    utils::download.file(source_url, source_file, mode = "wb", quiet = TRUE)
  }
  source_sha256 <- unname(tools::sha256sum(source_file))
  if (!identical(source_sha256, expected_sha256)) {
    stop("Unexpected source SHA-256; no input CSVs were written")
  }
  original <- haven::read_dta(source_file)
  fields <- c("earnings", "assignmt", "prevearn", "age")
  if (!identical(dim(original), c(11204L, 26L)) || !all(fields %in% names(original))) {
    stop("Unexpected public dataset schema; no input CSVs were written")
  }
  source_data <- as.data.frame(lapply(original[fields], as.numeric))
  if (any(!is.finite(as.matrix(source_data)))) stop("Unexpected missing or nonfinite source data")
  # Strip Stata labels and use the source's three-decimal prior-earnings precision.
  source_data$prevearn <- round(source_data$prevearn, 3)
  for (name in c("earnings", "assignmt", "age")) {
    source_data[[name]] <- as.integer(source_data[[name]])
  }
  n_input <- nrow(source_data)
  output_files <- character(2L)

  for (p in 1:2) {
    # Use a separate, reproducible seed-12 sequence for each input.
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(12)
    raw <- source_data
    required <- c("earnings", "assignmt", "prevearn", if (p == 2L) "age")
    if (!all(required %in% names(raw))) stop("Missing columns: ",
      paste(setdiff(required, names(raw)), collapse = ", "))
    if (!all(vapply(raw[required], is.numeric, logical(1)))) {
      stop("Required columns must be numeric")
    }
    raw <- raw[complete.cases(raw[, required]), , drop = FALSE]
    n_complete <- nrow(raw)
    if (!all(is.finite(as.matrix(raw[required])))) stop("Nonfinite input values")
    if (!all(raw$assignmt %in% c(0, 1))) stop("Assignment must be coded 0 or 1")
    if (p == 1L) {
      raw <- raw[raw$prevearn >= 0 & raw$prevearn <= 6000, , drop = FALSE]
      fraction <- 0.5
      cap <- 30L
      restrictions <- "0 <= prevearn <= 6000"
      thinning <- "singleton unchanged; otherwise min(max(1, floor(0.5 * group_n)), 30)"
      filename <- "jtpa_prevearn_0_5000_subsample.csv"
    } else {
      raw <- raw[raw$prevearn >= 0 & raw$prevearn <= 5000 & raw$age <= 50, , drop = FALSE]
      fraction <- 0.2
      cap <- 200L
      restrictions <- "0 <= prevearn <= 5000; age <= 50"
      thinning <- "singleton unchanged; otherwise min(max(1, floor(0.2 * group_n)), 200)"
      filename <- "jtpa_prevearn_0_5000_subsample_p2.csv"
    }
    if (!nrow(raw)) stop("No observations satisfy preprocessing restrictions")
    raw$.row_id <- seq_len(nrow(raw))
    prevearn_groups <- split(raw, raw$prevearn)
    sample_group <- function(df) {
      n_group <- nrow(df)
      if (n_group <= 1L) return(df)
      keep_n <- min(max(1L, floor(fraction * n_group)), cap)
      df[sample.int(n_group, size = keep_n, replace = FALSE), , drop = FALSE]
    }
    subsampled <- do.call(rbind, lapply(prevearn_groups, sample_group))
    subsampled <- subsampled[order(subsampled$prevearn, subsampled$.row_id), , drop = FALSE]
    subsampled$.row_id <- NULL
    names(subsampled)[names(subsampled) == "earnings"] <- "y"
    names(subsampled)[names(subsampled) == "assignmt"] <- "d"
    names(subsampled)[names(subsampled) == "prevearn"] <- "X1"
    if (p == 2L) names(subsampled)[names(subsampled) == "age"] <- "X2"
    subsampled <- subsampled[, c("y", "d", "X1", if (p == 2L) "X2"), drop = FALSE]

    output_csv <- file.path(root, "jtpa", filename)
    out_dir <- dirname(output_csv)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    write.csv(subsampled, file = output_csv, row.names = FALSE)
    group_sizes <- vapply(prevearn_groups, nrow, integer(1))
    metadata <- list(
      source_url = source_url, source_sha256 = source_sha256,
      source_columns = ncol(original), haven_version = as.character(utils::packageVersion("haven")),
      output_file = basename(output_csv), output_md5 = unname(tools::md5sum(output_csv)),
      n_input = n_input, n_complete = n_complete, n_restricted = nrow(raw),
      n_output = nrow(subsampled), prevearn_values = length(prevearn_groups),
      mass_points = sum(group_sizes > 1L), seed = 12L, rng_kind = RNGkind(),
      restrictions = restrictions, thinning = thinning,
      ordering = "prevearn, then original row position within the restricted sample",
      r_version = R.version.string
    )
    dput(metadata, file = paste0(output_csv, ".metadata.R"))
    cat("Input rows:", n_input, "\n")
    cat("Complete rows:", n_complete, "\n")
    cat("Restricted rows:", nrow(raw), "\n")
    cat("Output rows:", nrow(subsampled), "\n")
    cat("Saved:", output_csv, "\n")
    output_files[p] <- output_csv
  }
  invisible(output_files)
}

if (sys.nframe() == 0L) {
  script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script) != 1L) stop("Run this script with Rscript")
  root <- dirname(dirname(normalizePath(sub("^--file=", "", script))))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("Use Rscript jtpa/prepare.R [source.dta]")
  cky_prepare_jtpa(root, if (length(args)) args[1] else NULL)
}
