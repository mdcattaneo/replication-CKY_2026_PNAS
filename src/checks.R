# Validate numerical outputs and test the shared estimator against direct formulas.
cky_check_outputs <- function(root, output_dir, expected = 2000L) {
  if (length(expected) != 1L || !is.finite(expected) || expected < 1 ||
      expected != as.integer(expected)) stop("Expected a positive replication count")
  require_fields <- function(x, fields, label) {
    missing <- setdiff(fields, names(x))
    if (length(missing)) stop(label, " is missing required fields: ", paste(missing, collapse = ", "))
  }
  require_scalars <- function(x, fields, label) {
    require_fields(x, fields, label)
    valid <- vapply(x[fields], function(value) is.atomic(value) && length(value) == 1L && !is.na(value), logical(1))
    if (any(!valid)) stop(label, " has invalid scalar fields: ", paste(fields[!valid], collapse = ", "))
  }
  check_implementation <- function(meta, label) {
    require_fields(meta, "implementation", label)
    info <- meta$implementation
    require_fields(info, c("backend", "version", "source_md5"), label)
    if (!identical(info$backend, "CKY canonical causal trees") || info$version != 1L) {
      stop(label, " does not use the canonical estimator backend")
    }
    files <- c("src/estimators.R", "src/cky_theory.cpp")
    actual <- setNames(unname(tools::md5sum(file.path(root, files))), files)
    if (!identical(info$source_md5, actual)) stop(label, " estimator source checksum mismatch")
  }
  check_diagnostics <- function(f, expected) {
    z <- read.csv(f, check.names=FALSE)
    fields <- c("leaves","empty_est_leaves","missing_arm_est_leaves",
                "zero_convention_leaves","singleton_train_leaves","minimum_train_size")
    require_fields(z,c("method","K","replications",fields),f)
    stopifnot(nrow(z)==30L, !anyDuplicated(z[c("method","K")]),
              setequal(z$method, methods), setequal(z$K,1:5), all(z$replications==expected),
              all(is.finite(as.matrix(z[fields]))), all(as.matrix(z[fields])>=0),
              all(z$leaves>=expected), all(z$leaves<=expected*2^z$K),
              all(z$minimum_train_size>=expected),
              all(z$empty_est_leaves<=z$missing_arm_est_leaves),
              all(z$missing_arm_est_leaves<=z$leaves))
    ipw <- grepl("IPW$",z$method)
    stopifnot(all(z$zero_convention_leaves[ipw]==z$empty_est_leaves[ipw]),
              all(z$zero_convention_leaves[!ipw]==z$missing_arm_est_leaves[!ipw]),
              all(z$singleton_train_leaves[!ipw]==0),
              all(z$zero_convention_leaves[grepl("^NSS",z$method)]==0))
  }
  sim_dir <- output_dir
  emp_dir <- output_dir
  methods <- c("NSS-DIM", "NSS-IPW", "NSS-SSE", "HON-DIM", "HON-IPW", "HON-SSE")
  for (p in 1:2) {
    for (method in methods) {
      f <- file.path(sim_dir, sprintf("rmse_%s-p%d.csv", gsub("-", "_", method), p))
      x <- read.csv(f, check.names = FALSE)
      coords <- if (p == 1) "X1" else c("X1", "X2")
      require_fields(x, c("p", "method", "n", "m_reps", "tau", "seed", coords,
                          paste0("K", 1:5), paste0("N_K", 1:5)), f)
      stopifnot(nrow(x) == if (p == 1) 201L else 2601L,
                all(x$p == p), all(x$method == method), all(x$n == 1000L),
                all(x$m_reps == expected), all(x$tau == 1), all(x$seed == 123L),
                all(is.finite(as.matrix(x[paste0("K", 1:5)]))),
                all(as.matrix(x[paste0("K", 1:5)]) >= 0),
                all(as.matrix(x[paste0("N_K", 1:5)]) == expected))
      stopifnot(!anyDuplicated(x[coords]), all(as.matrix(x[coords]) >= 0), all(as.matrix(x[coords]) <= 1))
    }
    meta_file <- file.path(sim_dir, sprintf("simulation_metadata-p%d.R", p))
    meta <- dget(meta_file)
    check_implementation(meta, meta_file)
    expected_script_md5 <- unname(tools::md5sum(file.path(root,"CKY_2026_PNAS--simuls.R")))
    stopifnot(identical(meta$script_md5, expected_script_md5))
    require_fields(meta, c("status", "config", "design"), meta_file)
    require_scalars(meta, "status", meta_file)
    require_scalars(meta$config, c("p", "n", "m_reps", "seed"), paste(meta_file, "config"))
    require_scalars(meta$design, c("tau", "treatment_probability"), paste(meta_file, "design"))
    stopifnot(meta$status == "complete", meta$config$p == p, meta$config$n == 1000L,
              meta$config$m_reps == expected, meta$config$seed == 123L, meta$design$tau == 1,
              meta$design$treatment_probability == .5)
    check_diagnostics(file.path(sim_dir, sprintf("simulation_diagnostics-p%d.csv",p)), expected)
  }
  for (p in 1:2) {
    f <- file.path(emp_dir, if (p == 1) "rmse_K1to5.csv" else "rmse_K1to5_p2.csv")
    x <- read.csv(f, check.names = FALSE)
    coords <- if (p == 1) c("method", "K", "X1") else c("method", "K", "X1", "X2")
    require_fields(x, c(coords, "RMSE", "n", "m_reps", "tau", "seed", "n_valid"), f)
    meta_file <- file.path(emp_dir, sprintf("jtpa_p%d_metadata.R", p))
    meta <- dget(meta_file)
    check_implementation(meta, meta_file)
    expected_script_md5 <- unname(tools::md5sum(file.path(root,"CKY_2026_PNAS--empapp.R")))
    stopifnot(identical(meta$script_md5, expected_script_md5))
    require_scalars(meta, c("p", "n_complete", "training_n", "estimation_n", "m_reps",
                            "n_valid", "seed", "tau", "grid_points"), meta_file)
    stopifnot(setequal(unique(x$method), methods), setequal(unique(x$K), 1:5),
              all(is.finite(x$RMSE)), all(x$RMSE >= 0), all(x$m_reps == expected),
              all(x$n_valid == expected), all(x$seed == 12L),
              meta$p == p, meta$seed == 12L, meta$m_reps == expected,
              all(x$n == meta$n_complete),
              meta$training_n + meta$estimation_n == meta$n_complete,
              meta$n_valid == expected, nrow(x) == 30L * meta$grid_points)
    if (p == 1L) stopifnot(meta$n_complete == 2475L, meta$training_n == 1237L, meta$estimation_n == 1238L)
    stopifnot(!anyDuplicated(x[coords]), isTRUE(all.equal(meta$treatment_probability, 2/3, tolerance=1e-14)))
    check_diagnostics(file.path(emp_dir, sprintf("jtpa_diagnostics-p%d.csv",p)), expected)
  }
  cat("Validated all four experiments, all six methods, five depths, and", expected, "replications per cell.\n")
  invisible(TRUE)
}

# Independent oracle: enumerate candidate partitions; SSE uses OLS residual sums.
cky_test_engine <- function(root) {
  cky_load_engine(root)
  checks <- 0L
  assert <- function(ok, description) {
    if (!isTRUE(ok)) stop(description, call. = FALSE)
    checks <<- checks + 1L
  }
  close <- function(a, b, description) {
    assert(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-10)), description)
  }
  oracle <- function(x, y, d, rule, xi, K, ex, ey, ed, grid) {
    leaf <- function(rows) {
      if (rule == "IPW") return(if (length(rows)) mean(ey[rows] * (ed[rows] - xi) / (xi * (1-xi))) else 0)
      t <- rows[ed[rows] == 1]; c <- rows[ed[rows] == 0]
      if (!length(t) || !length(c)) 0 else mean(ey[t]) - mean(ey[c])
    }
    fitmean <- function(rows) {
      if (rule == "IPW") mean(y[rows] * (d[rows] - xi) / (xi*(1-xi))) else
        mean(y[rows[d[rows] == 1]]) - mean(y[rows[d[rows] == 0]])
    }
    rss <- function(rows) {
      # Direct residual fit of y on intercept and treatment, independently in a child.
      residual <- lm.fit(cbind(1, d[rows]), y[rows])$residuals
      sum(residual^2)
    }
    grow <- function(rows, est, points, depth) {
      best <- -Inf; chosen <- NULL
      if (depth < K && length(rows) > 1L) {
        for (j in seq_len(ncol(x))) {
          cuts <- head(sort(unique(x[rows,j])), -1L)
          for (cut in cuts) {
            L <- rows[x[rows,j] <= cut]; R <- rows[x[rows,j] > cut]
            if (rule != "IPW" && (!all(c(0,1) %in% d[L]) || !all(c(0,1) %in% d[R]))) next
            score <- if (rule == "SSE") rss(rows)-rss(L)-rss(R) else
              length(L)*length(R)/length(rows)*(fitmean(L)-fitmean(R))^2
            if (!is.finite(best) || score > best + 64*.Machine$double.eps*max(abs(best),abs(score))) {best <- score; chosen <- list(j=j,cut=cut,L=L,R=R)}
          }
        }
      }
      if (is.null(chosen)) return(rep(leaf(est), length(points)))
      left <- grid[points,chosen$j] <= chosen$cut
      est_left <- ex[est,chosen$j] <= chosen$cut
      answer <- numeric(length(points))
      answer[left] <- grow(chosen$L, est[est_left], points[left], depth+1L)
      answer[!left] <- grow(chosen$R, est[!est_left], points[!left], depth+1L)
      answer
    }
    grow(seq_len(nrow(x)), seq_len(nrow(ex)), seq_len(nrow(grid)), 0L)
  }

  # IPW is a transformed-outcome mean even when local arm proportions differ.
  x <- matrix(1:4, ncol=1); y <- c(1,1,1,0); d <- c(1,1,1,0)
  ipw <- cky_fit(x,y,d,"IPW",.5,0)
  dim <- cky_fit(x,y,d,"DIM",.5,0)
  close(predict(ipw,x), rep(1.5,4), "IPW leaf mean must be 1.5")
  close(predict(dim,x), rep(1,4), "DIM leaf mean must be 1")

  # One-arm IPW construction leaves are admissible; DIM/SSE omit those splits.
  tiny_x <- matrix(1:2,ncol=1)
  for (rule in c("IPW","DIM","SSE")) {
    tr <- cky_fit(tiny_x,c(1,0),c(1,0),rule,.5,1)
    assert(nrow(tr$nodes) == if(rule == "IPW") 3L else 1L, paste(rule,"admissibility"))
  }

  # Empty and missing-arm honest estimation leaves return the specified values.
  tr <- cky_fit(tiny_x,c(1,0),c(1,0),"IPW",.5,1,
                est_x=matrix(1,ncol=1),est_y=3,est_d=1)
  close(predict(tr,tiny_x),c(6,0),"IPW uses one-arm mean and zero only when empty")
  for (rule in c("DIM","SSE")) {
    tx <- matrix(1:4,ncol=1); ty <- c(2,0,8,1); td <- c(1,0,1,0)
    tr <- cky_fit(tx,ty,td,rule,.5,1,est_x=matrix(c(1,3,4),ncol=1),
                  est_y=c(20,10,2),est_d=c(1,1,0))
    close(predict(tr,tx),c(0,0,8,8),paste(rule,"must not inherit the parent estimate"))
  }
  for (rule in c("IPW","DIM","SSE")) {
    tr <- cky_fit(x,y,d,rule,.5,3,est_x=matrix(numeric(),ncol=1),est_y=numeric(),est_d=numeric())
    close(predict(tr,x),rep(0,4),paste(rule,"entire estimation sample empty"))
  }

  # Observed-value thresholds, <= convention, coordinate ties, zero-gain splits.
  tr <- cky_fit(cbind(c(0,2),c(0,2)),c(3,0),c(1,0),"IPW",.5,1)
  assert(tr$nodes[1,"feature"]==1 && tr$nodes[1,"cut"]==0,"First coordinate and observed-value cutoff")
  close(predict(tr,cbind(c(0,1,2),c(0,1,2))),c(6,0,0),"Prediction between observations goes right of observed cutoff")
  for (rule in c("IPW","DIM","SSE")) {
    tr <- cky_fit(matrix(1:8,ncol=1),rep(0,8),rep(0:1,4),rule,.5,2)
    assert(nrow(tr$nodes)>1 && tr$nodes[1,"gain"]==0,paste(rule,"no positive-gain stopping restriction"))
    same_x <- cky_fit(matrix(1,8,2),1:8,rep(0:1,4),rule,.5,3)
    assert(nrow(same_x$nodes)==1,paste(rule,"no splitting equal covariates"))
  }

  # Hundreds of small comparisons against independently enumerated direct formulas.
  set.seed(97431)
  for (trial in 1:70) {
    n <- sample(8:30,1); p <- sample(1:3,1); ne <- sample(0:25,1)
    x <- matrix(if (trial %% 2) runif(n*p) else sample(0:9,n*p,TRUE),n,p)
    y <- rnorm(n,mean=2); d <- rbinom(n,1,.4)
    ex <- matrix(if (trial %% 2) runif(ne*p) else sample(0:9,ne*p,TRUE),ne,p)
    ey <- rnorm(ne); ed <- rbinom(ne,1,.6)
    grid <- rbind(x,ex,matrix(rep(c(-1,11),p),2,p))
    for (rule in c("IPW","DIM","SSE")) {
      tr <- cky_fit(x,y,d,rule,.4,3,est_x=ex,est_y=ey,est_d=ed)
      pred <- predict(tr,grid,depths=0:3)
      for (K in 0:3) {
        ref <- oracle(x,y,d,rule,.4,K,ex,ey,ed,grid)
        close(pred[,K+1],ref,paste("Oracle mismatch",trial,rule,K))
        separate <- cky_fit(x,y,d,rule,.4,K,est_x=ex,est_y=ey,est_d=ed)
        close(pred[,K+1],predict(separate,grid),paste("Depth truncation mismatch",trial,rule,K))
      }
      changed <- cky_fit(x,y,d,rule,.4,3,est_x=ex,est_y=ey+100,est_d=1-ed)
      columns <- c("left","right","feature","cut","depth","n_train","n1_train","gain","n0_train")
      assert(identical(tr$nodes[,columns],changed$nodes[,columns]),paste(rule,"honest outcomes leaked into construction"))
    }
  }

  # Fixed DIM objective fixture: the maximizing cutoff is 6.
  fixture_y <- c(-2,4,0,-3,-2,0,3,0,0,4,1,0)
  fixture_d <- c(0,0,0,0,1,1,0,1,1,0,1,1)
  fixture_x <- matrix(1:12,ncol=1)
  fixture <- cky_fit(fixture_x,fixture_y,fixture_d,"DIM",.5,1)
  assert(fixture$nodes[1,"cut"]==6,"DIM must select the maximizing child-contrast cutoff")
  for (scale in c(1e-9,1e9)) {
    scaled <- cky_fit(fixture_x,fixture_y*scale,fixture_d,"DIM",.5,1)
    assert(scaled$nodes[1,"cut"]==6,"Numerical tie rule must preserve outcome scaling")
  }

  # No stochastic fits: calls must not consume the simulation RNG stream.
  before <- .Random.seed
  invisible(cky_fit(x,y,d,"IPW",.4,3))
  assert(identical(before,.Random.seed),"Tree fitting consumed RNG")
  assert(inherits(try(cky_fit(x,y,d,"IPW",0,1),silent=TRUE),"try-error"),"Reject invalid propensity")
  assert(inherits(try(cky_fit(x,y,d,"DIM",.4,-1),silent=TRUE),"try-error"),"Reject invalid depth")
  cat("Passed", checks,"checks against direct formulas, edge cases, honest isolation, and depth equivalence.\n")
  invisible(TRUE)
}

# The manifest uses repository-relative paths, including the bundled outputs.
cky_manifest <- function(root, write = FALSE) {
  oldwd <- setwd(root)
  on.exit(setwd(oldwd), add = TRUE)
  manifest <- file.path("output", "SHA256SUMS")
  if (write) {
    files <- c("CKY_2026_PNAS--simuls.R", "CKY_2026_PNAS--empapp.R",
               "CKY_2026_PNAS--plots.R", "src/estimators.R", "src/checks.R",
               "src/cky_theory.cpp",
               "jtpa/prepare.R",
               "jtpa/jtpa_prevearn_0_5000_subsample.csv",
               "jtpa/jtpa_prevearn_0_5000_subsample_p2.csv",
               list.files("output", pattern = "\\.(csv|R|txt)$", full.names = TRUE))
    files <- sort(files)
    sums <- system2("sha256sum", shQuote(files), stdout = TRUE)
    if (!is.null(attr(sums, "status")) || length(sums) != length(files)) {
      stop("Could not calculate package checksums")
    }
    writeLines(sums, manifest)
    cat("Saved", length(files), "checksums to output/SHA256SUMS\n")
  } else {
    result <- system2("sha256sum", c("--check", "--status", shQuote(manifest)))
    if (result != 0L) stop("Package checksum verification failed")
    cat("Package checksums verified.\n")
  }
  invisible(TRUE)
}

# Optional verification utility; the three main scripts remain the public workflow.
if (sys.nframe() == 0L) {
  script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script) != 1L) stop("Run with Rscript src/checks.R [check|test]")
  root <- dirname(dirname(normalizePath(sub("^--file=", "", script))))
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args)) args[1] else "check"
  if (length(args) > 1L || !mode %in% c("check", "test")) {
    stop("Use Rscript src/checks.R [check|test]")
  }
  source(file.path(root, "src", "estimators.R"))
  if (mode == "test") {
    cky_test_engine(root)
  } else {
    cky_check_outputs(root, file.path(root, "output"))
    cky_manifest(root)
  }
}
