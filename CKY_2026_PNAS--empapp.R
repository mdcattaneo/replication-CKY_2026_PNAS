#!/usr/bin/env Rscript
# Accuracy Limits of Causal Trees for Individualized Treatment Effects
# Matias D. Cattaneo, Jason M. Klusowski, and Ruiqi (Rae) Yu

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run this file using Rscript")
script_file <- normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE)
root <- dirname(script_file)
source(file.path(root, "src", "estimators.R"))
cky_load_engine(root)

get_env_int <- function(name, default = 2000L) {
  value <- Sys.getenv(name, as.character(default))
  parsed <- suppressWarnings(as.integer(value))
  if (!grepl("^[0-9]+$", value) || is.na(parsed) || parsed < 1L) {
    stop(name, " must be a positive integer")
  }
  parsed
}

dimension_setting <- Sys.getenv("CKY_DIMENSIONS", "1,2")
dimensions <- trimws(strsplit(dimension_setting, ",", fixed = TRUE)[[1]])
if (!grepl("^[[:space:]]*[12]([[:space:]]*,[[:space:]]*[12])?[[:space:]]*$", dimension_setting) ||
    !length(dimensions) || any(!dimensions %in% c("1", "2")) || anyDuplicated(dimensions)) {
  stop("CKY_DIMENSIONS must be 1, 2, or 1,2")
}
dimensions <- sort(as.integer(dimensions))

output_dir <- path.expand(Sys.getenv("CKY_OUTPUT_DIR", "output"))
if (!grepl("^/", output_dir)) output_dir <- file.path(root, output_dir)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(output_dir) || file.access(output_dir, 2L) != 0L) {
  stop("Output directory is not writable: ", output_dir)
}
output_dir <- normalizePath(output_dir, mustWork = TRUE)

all_files <- character()
for (p in dimensions) {
  m.reps <- get_env_int(paste0("CKY_M_REPS_JTPA_P", p))
  seed <- 12L
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  K.list <- 1:5
  xi <- 2/3 # Known randomized offer probability in the National JTPA Study.
  started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)


  input_csv <- file.path(root, "jtpa", if (p == 1L)
    "jtpa_prevearn_0_5000_subsample.csv" else "jtpa_prevearn_0_5000_subsample_p2.csv")
  raw <- read.csv(input_csv, stringsAsFactors = FALSE)
  n_input <- nrow(raw)
  covariates <- paste0("X", seq_len(p))
  required <- c("y", "d", covariates)
  if (!all(required %in% names(raw))) stop("Missing required columns in ", input_csv)
  if (!all(vapply(raw[required], is.numeric, logical(1)))) {
    stop("JTPA inputs must be numeric: ", input_csv)
  }
  keep <- complete.cases(raw[, required])
  raw <- raw[keep, , drop = FALSE]
  if (!nrow(raw) || !all(is.finite(as.matrix(raw[required])))) {
    stop("Empty or nonfinite JTPA input after complete-case selection")
  }
  if (!all(raw$d %in% c(0, 1)) || !all(c(0, 1) %in% raw$d)) {
    stop("JTPA input must contain both treatment arms coded 0 and 1")
  }

  tauTRUE <- mean(raw$y[raw$d == 1]) - mean(raw$y[raw$d == 0])

  if (p == 1L) {
    grid_vals <- sort(unique(raw$X1))
    if (length(grid_vals) > 201) {
      keep_idx <- unique(round(seq(1, length(grid_vals), length.out = 201)))
      grid_vals <- grid_vals[keep_idx]
    }

    grid <- data.frame(X1 = grid_vals)
  } else {
    x1.seq <- seq(quantile(raw$X1, 0.01, na.rm = TRUE),
                  quantile(raw$X1, 0.99, na.rm = TRUE),
                  length.out = 41)
    x2.seq <- seq(quantile(raw$X2, 0.01, na.rm = TRUE),
                  quantile(raw$X2, 0.99, na.rm = TRUE),
                  length.out = 31)
    grid <- expand.grid(X1 = x1.seq, X2 = x2.seq)
  }

  n <- nrow(raw)
  half <- floor(n / 2)
  train_idx <- seq_len(half)
  est_idx <- seq.int(half + 1L, n)
  stopifnot(length(intersect(train_idx, est_idx)) == 0L,
            identical(c(train_idx, est_idx), seq_len(n)))
  if (!all(c(0, 1) %in% raw$d[train_idx]) ||
      !all(c(0, 1) %in% raw$d[est_idx])) {
    stop("Both honest folds must contain treatment and control observations")
  }

  print(paste("n =", n))

  y <- raw$y
  d <- raw$d
  X1 <- raw$X1
  if (p == 2L) X2 <- raw$X2

  if (p == 1L) {
    gen_sample <- function() {
      x_perm <- X1
      treat_idx <- which(d == 1)
      control_idx <- which(d == 0)
      x_perm[treat_idx] <- sample(X1[treat_idx], replace = FALSE)
      x_perm[control_idx] <- sample(X1[control_idx], replace = FALSE)
      data.frame(
        y = y,
        d = d,
        X1 = x_perm
      )
    }
  } else {
    gen_sample <- function() {
      x1_perm <- X1
      x2_perm <- X2
      treat_idx <- which(d == 1)
      control_idx <- which(d == 0)

      # Permute the covariate pairs jointly within treatment arms.
      x1x2_treat <- cbind(X1[treat_idx], X2[treat_idx])
      x1x2_ctrl <- cbind(X1[control_idx], X2[control_idx])
      treat_ord <- sample.int(nrow(x1x2_treat), replace = FALSE)
      ctrl_ord <- sample.int(nrow(x1x2_ctrl), replace = FALSE)

      x1_perm[treat_idx] <- x1x2_treat[treat_ord, 1]
      x2_perm[treat_idx] <- x1x2_treat[treat_ord, 2]
      x1_perm[control_idx] <- x1x2_ctrl[ctrl_ord, 1]
      x2_perm[control_idx] <- x1x2_ctrl[ctrl_ord, 2]

      data.frame(y = y, d = d, X1 = x1_perm, X2 = x2_perm)
    }
  }


  methods <- c("NSS-DIM", "NSS-IPW", "NSS-SSE",
               "HON-DIM", "HON-IPW", "HON-SSE")

  rmse <- array(
    0,
    dim = c(nrow(grid), length(K.list), length(methods)),
    dimnames = list(NULL, paste0("K", K.list), methods)
  )

  diagnostics <- NULL
  for (rep in seq_len(m.reps)) {
    samp <- gen_sample()
    fit <- tryCatch(cky_experiment(samp, covariates, grid, K.list, xi,
                                   train_idx, est_idx),
                    error=function(e) stop("JTPA p=", p, " rep ", rep, ": ", conditionMessage(e)))
    diagnostics <- if(is.null(diagnostics)) fit$diagnostics else diagnostics + fit$diagnostics

    for (k in seq_along(K.list)) {
      K <- K.list[k]
      for (j in seq_along(methods)) {
        tauhat <- fit$prediction[,k,j]
        if (length(tauhat) != nrow(grid) || any(!is.finite(tauhat))) {
          stop(sprintf("Invalid prediction: rep=%d, K=%d, method=%s",
                       rep, K, methods[j]))
        }
        rmse[, k, j] <- rmse[, k, j] + (tauhat - tauTRUE)^2
      }
    }

    if (rep %% 50 == 0) cat("rep", rep, "finished\n")
  }

  rmse <- sqrt(rmse / m.reps)

  if (any(!is.finite(rmse))) stop("Nonfinite RMSE; no result file written")

  rmse_csv <- do.call(rbind, lapply(seq_along(methods), function(j) {
    data.frame(
      lapply(grid, base::rep, times = length(K.list)),
      K = rep(K.list, each = nrow(grid)),
      method = methods[j],
      RMSE = c(rmse[, , j]),
      n = n, m_reps = m.reps, tau = tauTRUE, seed = seed, n_valid = m.reps
    )
  }))

  write.csv(
    rmse_csv,
    file = file.path(output_dir, if (p == 1L) "rmse_K1to5.csv" else "rmse_K1to5_p2.csv"),
    row.names = FALSE
  )

  write.csv(cky_diagnostics(diagnostics, methods, K.list, m.reps),
            file.path(output_dir, sprintf("jtpa_diagnostics-p%d.csv", p)), row.names=FALSE)

  # Machine-readable provenance accompanies every successful run.
  metadata <- list(
    experiment = "JTPA covariate permutation", p = p,
    input_file = basename(input_csv), input_md5 = unname(tools::md5sum(input_csv)),
    n_input = n_input, n_complete = n, n_dropped_missing = n_input - n,
    n_treated = sum(d == 1), n_control = sum(d == 0),
    training_n = length(train_idx), estimation_n = length(est_idx),
    training_treated = sum(d[train_idx] == 1),
    estimation_treated = sum(d[est_idx] == 1),
    fold_rule = "first floor(n/2) rows train; all remaining rows estimate",
    m_reps = m.reps, n_valid = m.reps, seed = seed, rng_kind = RNGkind(),
    tau = tauTRUE, depths = K.list, methods = methods, grid_points = nrow(grid),
    started_utc = started_at, finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    r_version = R.version.string,
    implementation = cky_provenance(),
    treatment_probability = xi,
    propensity_source = "National JTPA Study randomized offer: two thirds; fixed for both folds",
    script_md5 = unname(tools::md5sum(script_file))
  )
  dput(metadata, file = file.path(output_dir, sprintf("jtpa_p%d_metadata.R", p)))
  writeLines(capture.output(sessionInfo()),
             file.path(output_dir, sprintf("jtpa_p%d_sessionInfo.txt", p)))
  cat("Saved JTPA p=", p, " results and metadata to ", output_dir, "\n", sep = "")
  all_files <- c(all_files, file.path(output_dir, c(
    if (p == 1L) "rmse_K1to5.csv" else "rmse_K1to5_p2.csv",
    sprintf("jtpa_diagnostics-p%d.csv", p), sprintf("jtpa_p%d_metadata.R", p),
    sprintf("jtpa_p%d_sessionInfo.txt", p))))
}
