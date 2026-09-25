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

data_dir <- path.expand(Sys.getenv("CKY_OUTPUT_DIR", "output"))
if (!grepl("^/", data_dir)) data_dir <- file.path(root, data_dir)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(data_dir) || file.access(data_dir, 2L) != 0L) {
  stop("Output directory is not writable: ", data_dir)
}
data_dir <- normalizePath(data_dir, mustWork = TRUE)

methods <- c("NSS-DIM", "NSS-IPW", "NSS-SSE", "HON-DIM", "HON-IPW", "HON-SSE")
simulation_configs <- lapply(dimensions, function(p) {
  list(
    label = paste0("p", p),
    p = p,
    seed = 123L,
    n = 1000L,
    m_reps = get_env_int(paste0("CKY_M_REPS_P", p)),
    K_list = 1:5,
    grid = if (p == 1L) data.frame(X1 = seq(0, 1, length.out = 201)) else
      expand.grid(X1 = seq(0, 1, length.out = 51), X2 = seq(0, 1, length.out = 51))
  )
})

safe_name <- function(x) gsub("[^[:alnum:]_]+", "_", x)

run_simulation <- function(config) {
  # Explicit RNG settings make independent p=1 and p=2 runs reproducible even
  # when an interactive startup file has changed the R RNG defaults.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(config$seed)
  p <- config$p
  n <- config$n
  m_reps <- config$m_reps
  K_list <- config$K_list
  grid <- config$grid
  tauTRUE <- 1
  half <- floor(n / 2)
  covariates <- paste0("X", seq_len(p))

  gen_sample <- function() {
    X <- as.data.frame(replicate(p, runif(n)))
    names(X) <- covariates
    d <- rbinom(n, 1, 0.5)
    y0 <- rnorm(n)
    y1 <- 1 + rnorm(n)
    y <- d * y1 + (1 - d) * y0
    data.frame(y = y, d = d, X)
  }

  rmse_ss <- array(0, dim = c(nrow(grid), length(K_list), length(methods)),
                   dimnames = list(NULL, paste0("K", K_list), methods))
  rmse_count <- rmse_ss
  diagnostics <- NULL
  cat("Starting ", config$label, ": n = ", n, ", reps = ", m_reps,
      ", tau = ", tauTRUE, ", seed = ", config$seed,
      ", grid rows = ", nrow(grid), "\n", sep = "")

  for (rep_id in seq_len(m_reps)) {
    samp <- gen_sample()
    fit <- tryCatch(cky_experiment(samp, covariates, grid, K_list, xi=0.5,
                                   train_idx=seq_len(half), est_idx=seq.int(half+1L,n)),
                    error=function(e) stop(config$label, " rep ", rep_id, ": ", conditionMessage(e)))
    diagnostics <- if (is.null(diagnostics)) fit$diagnostics else diagnostics + fit$diagnostics
    for (k in seq_along(K_list)) {
      K <- K_list[k]
      for (j in seq_along(methods)) {
        context <- paste(config$label, "rep", rep_id, "K", K, methods[j])
        tauhat <- fit$prediction[,k,j]
        if (length(tauhat) != nrow(grid) || any(!is.finite(tauhat))) {
          stop(context, " produced missing, nonfinite, or incorrectly sized predictions.",
               call. = FALSE)
        }
        err2 <- (tauhat - tauTRUE)^2
        if (any(!is.finite(err2))) {
          stop(context, " produced nonfinite squared errors.", call. = FALSE)
        }
        rmse_ss[, k, j] <- rmse_ss[, k, j] + err2
        rmse_count[, k, j] <- rmse_count[, k, j] + 1L
      }
    }
    if (rep_id %% 50L == 0L || rep_id == m_reps) {
      cat(config$label, "rep", rep_id, "finished\n")
    }
  }
  stopifnot(all(rmse_count == m_reps))
  rmse <- sqrt(rmse_ss / rmse_count)
  if (any(!is.finite(rmse))) stop("Nonfinite aggregate RMSE for ", config$label, ".")
  list(config = config, tau = tauTRUE, rmse = rmse, count = rmse_count, diagnostics = diagnostics)
}

write_simulation_csvs <- function(result) {
  config <- result$config
  K_cols <- paste0("K", config$K_list)
  count_cols <- paste0("N_", K_cols)
  files <- character(length(methods))
  for (j in seq_along(methods)) {
    rmse <- as.data.frame(result$rmse[, , j])
    counts <- as.data.frame(result$count[, , j])
    names(rmse) <- K_cols
    names(counts) <- count_cols
    out <- cbind(data.frame(p = config$p, method = methods[j], n = config$n,
                            m_reps = config$m_reps, tau = result$tau, seed = config$seed),
                 config$grid, rmse, counts)
    files[j] <- file.path(data_dir, sprintf("rmse_%s-%s.csv", safe_name(methods[j]), config$label))
    write.csv(out, files[j], row.names = FALSE)
  }
  diagnostic_file <- file.path(data_dir, paste0("simulation_diagnostics-", config$label, ".csv"))
  write.csv(cky_diagnostics(result$diagnostics, methods, config$K_list, config$m_reps),
            diagnostic_file, row.names=FALSE)
  files <- c(files, diagnostic_file)
  # Saved only after every requested replication and all six method files finish.
  # dget() reads this self-contained record; no additional package is needed.
  metadata <- list(
    status = "complete",
    completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    design = list(mu0 = 0, mu1 = 1, tau = result$tau, potential_outcome_variance = 1,
                  treatment_probability = 0.5, covariates = "independent Uniform[0,1]"),
    config = config,
    methods = methods,
    counts_per_grid_point_and_depth = config$m_reps,
    rng_kind = RNGkind(),
    script = basename(script_file),
    script_md5 = unname(tools::md5sum(script_file)),
    implementation = cky_provenance(),
    output_files = basename(files),
    session_info = utils::sessionInfo()
  )
  metadata_file <- file.path(data_dir, paste0("simulation_metadata-", config$label, ".R"))
  dput(metadata, file = metadata_file)
  c(files, metadata_file)
}

all_files <- character()
for (config in simulation_configs) {
  result <- run_simulation(config)
  all_files <- c(all_files, write_simulation_csvs(result))
}
cat("Saved simulation files:\n", paste(" -", all_files, collapse = "\n"), "\n")
