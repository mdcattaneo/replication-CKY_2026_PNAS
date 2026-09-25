# Shared canonical CKY causal-tree engine. Requires base R and a C++ compiler.
# Source this file, then call cky_load_engine(root) once per R process.
.cky_engine <- new.env(parent = emptyenv())

cky_load_engine <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  cpp <- file.path(root, "src", "cky_theory.cpp")
  checksum <- unname(tools::md5sum(cpp))
  if (is.na(checksum)) stop("Missing canonical tree source: ", cpp)
  build <- file.path(root, ".build", paste0(substr(checksum, 1, 16), "-", getRversion()))
  dir.create(build, recursive = TRUE, showWarnings = FALSE)
  dll <- file.path(build, paste0("cky_theory", .Platform$dynlib.ext))
  if (!file.exists(dll)) {
    # Compile within an isolated directory so .o/.so files never enter src/.
    build_source <- file.path(build, "cky_theory.cpp")
    file.copy(cpp, build_source, overwrite = TRUE)
    oldwd <- setwd(build)
    on.exit(setwd(oldwd), add = TRUE)
    output <- system2(file.path(R.home("bin"), "R"),
                      c("CMD", "SHLIB", "cky_theory.cpp", "-o", shQuote(basename(dll))),
                      stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(output, "status")) || !file.exists(dll)) {
      stop("Canonical tree build failed:\n", paste(output, collapse = "\n"))
    }
  }
  loaded <- dyn.load(dll)
  .cky_engine$fit <- getNativeSymbolInfo("cky_fit", PACKAGE = loaded)
  .cky_engine$predict <- getNativeSymbolInfo("cky_predict", PACKAGE = loaded)
  .cky_engine$root <- root
  invisible(dll)
}

cky_numeric_matrix <- function(x, label) {
  if (is.data.frame(x) && !all(vapply(x, is.numeric, logical(1)))) {
    stop(label, " must have numeric columns")
  }
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L || !ncol(x) || any(!is.finite(x))) {
    stop(label, " must be a finite numeric matrix with at least one column")
  }
  storage.mode(x) <- "double"
  x
}

cky_fit <- function(x, y, d, rule, xi, maxdepth,
                    est_x = x, est_y = y, est_d = d) {
  if (is.null(.cky_engine$fit)) stop("Call cky_load_engine() first")
  rule <- match.arg(rule, c("IPW", "DIM", "SSE"))
  x <- cky_numeric_matrix(x, "x")
  est_x <- cky_numeric_matrix(est_x, "est_x")
  if (!nrow(x) || ncol(x) != ncol(est_x)) stop("Incompatible construction/estimation covariates")
  check_response <- function(y, d, n, label) {
    if (!is.numeric(y) || length(y) != n || any(!is.finite(y)) ||
        !is.numeric(d) || length(d) != n || anyNA(d) || any(!d %in% c(0, 1))) {
      stop("Invalid ", label, " outcomes/assignments")
    }
  }
  check_response(y, d, nrow(x), "construction")
  check_response(est_y, est_d, nrow(est_x), "estimation")
  if (length(xi) != 1L || !is.finite(xi) || xi <= 0 || xi >= 1) stop("xi must lie in (0,1)")
  if (length(maxdepth) != 1L || !is.finite(maxdepth) || maxdepth != floor(maxdepth) ||
      maxdepth < 0 || maxdepth > 30) stop("maxdepth must be an integer in [0,30]")
  nodes <- .Call(.cky_engine$fit, x, as.double(y), as.integer(d),
                 est_x, as.double(est_y), as.integer(est_d),
                 match(rule, c("IPW", "DIM", "SSE")) - 1L, as.double(xi), as.integer(maxdepth))
  colnames(nodes) <- c("left", "right", "feature", "cut", "depth", "estimate",
                       "n_train", "n1_train", "n_est", "n1_est", "gain", "n0_train", "n0_est")
  structure(list(nodes = nodes, p = ncol(x), covariates = colnames(x), rule = rule,
                 xi = xi, maxdepth = as.integer(maxdepth)), class = "cky_tree")
}

predict.cky_tree <- function(object, newdata, depths = object$maxdepth, ...) {
  x <- cky_numeric_matrix(newdata, "newdata")
  if (ncol(x) != object$p) stop("Prediction covariate count differs from construction")
  if (!is.null(object$covariates) && !is.null(colnames(x)) &&
      !identical(colnames(x), object$covariates)) stop("Prediction covariate order differs")
  if (!length(depths) || any(!is.finite(depths)) || any(depths != floor(depths)) ||
      any(depths < 0 | depths > object$maxdepth) || any(diff(depths) <= 0)) {
    stop("depths must be strictly increasing integers between 0 and maxdepth")
  }
  result <- .Call(.cky_engine$predict, object$nodes, x, as.integer(depths))
  colnames(result) <- paste0("K", depths)
  result
}

cky_provenance <- function() {
  files <- c("src/estimators.R", "src/cky_theory.cpp")
  list(backend = "CKY canonical causal trees", version = 1L,
       source_md5 = setNames(unname(tools::md5sum(file.path(.cky_engine$root, files))), files),
       threshold = "largest observed value in left child; left <= cut, right > cut",
       ties = "first coordinate, then smallest observed cutoff; numerical ties within 64 double eps * max(absolute scores)",
       stopping = "maximum depth or no admissible split; zero-gain maximizers allowed",
       honesty = "construction sample only selects splits; estimation sample only supplies node means",
       empty_leaves = "IPW zero iff empty; DIM/SSE zero if either arm absent",
       r_version = R.version.string)
}


# Fit once to maximum depth, then evaluate every requested truncation. This is
# equivalent to separate fits because the split rule does not depend on K.
cky_experiment <- function(sample, covariates, grid, depths, xi, train_idx, est_idx) {
  methods <- c("NSS-DIM", "NSS-IPW", "NSS-SSE", "HON-DIM", "HON-IPW", "HON-SSE")
  x <- cky_numeric_matrix(sample[covariates], "sample covariates")
  grid <- cky_numeric_matrix(grid[covariates], "grid")
  prediction <- array(0, c(nrow(grid), length(depths), length(methods)),
                      dimnames = list(NULL, paste0("K", depths), methods))
  diagnostics <- matrix(0, length(methods)*length(depths), 6,
                         dimnames = list(NULL, c("leaves", "empty_est_leaves", "missing_arm_est_leaves",
                                                "zero_convention_leaves", "singleton_train_leaves", "minimum_train_size")))
  for (j in seq_along(methods)) {
    part <- strsplit(methods[j], "-", fixed = TRUE)[[1]]
    construct <- if (part[1] == "HON") train_idx else seq_len(nrow(sample))
    estimate <- if (part[1] == "HON") est_idx else seq_len(nrow(sample))
    tree <- cky_fit(x[construct,,drop=FALSE], sample$y[construct], sample$d[construct],
                    part[2], xi, max(depths), est_x=x[estimate,,drop=FALSE],
                    est_y=sample$y[estimate], est_d=sample$d[estimate])
    prediction[,,j] <- predict(tree, grid, depths=depths)
    for (k in seq_along(depths)) {
      nodes <- tree$nodes
      leaf <- nodes[nodes[,"depth"] == depths[k] |
                    (nodes[,"depth"] < depths[k] & nodes[,"left"] == 0),,drop=FALSE]
      empty <- leaf[,"n_est"] == 0
      missing_arm <- leaf[,"n1_est"] == 0 | leaf[,"n0_est"] == 0
      diagnostics[(j-1L)*length(depths)+k,] <- c(nrow(leaf), sum(empty), sum(missing_arm),
              sum(if(part[2] == "IPW") empty else missing_arm),
              sum(leaf[,"n_train"] == 1), min(leaf[,"n_train"]))
    }
  }
  list(prediction=prediction, diagnostics=diagnostics)
}

cky_diagnostics <- function(totals, methods, depths, replications) {
  data.frame(method=rep(methods,each=length(depths)), K=rep(depths,length(methods)),
             replications=replications, totals, row.names=NULL)
}
