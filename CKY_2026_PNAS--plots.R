# Rebuild all 24 manuscript panels from the numerical CSV outputs.
# Paths may be absolute or relative to the repository root.
cky_plots <- function(root, output_dir, plot_dir) {
  # Shared base-R plotting functions for all 24 manuscript panels.
  # The manuscript includes panels at 0.31 * textwidth = about 2.17 inches.
  # At a PDF width of 3.2 inches, 11-point text remains about 7.45 points.
  # Keep all text at this size or larger; do not scale it with cex < 1.

  cky_methods <- c("NSS-IPW", "NSS-DIM", "NSS-SSE", "HON-IPW", "HON-DIM", "HON-SSE")
  cky_depth_colors <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#E69F00")

  cky_path <- function(path, root) {
    if (grepl("^(/|~)", path)) path.expand(path) else file.path(root, path)
  }

  cky_require_columns <- function(df, columns, source) {
    missing <- setdiff(columns, names(df))
    if (length(missing)) stop(source, ": missing columns: ", paste(missing, collapse = ", "))
  }

  cky_depths <- function(depths) {
    if (!is.numeric(depths) || anyNA(depths) || any(depths != as.integer(depths)) ||
        any(!depths %in% seq_along(cky_depth_colors))) {
      stop("Tree depths must be integers from 1 through 5.")
    }
    sort(unique(as.integer(depths)))
  }

  cky_limit <- function(values) {
    if (!length(values) || any(!is.finite(values)) || any(values < 0)) {
      stop("Every plotted RMSE must be finite and nonnegative.")
    }
    upper <- max(values)
    if (upper == 0) upper <- 1
    c(0, max(pretty(c(0, upper), n = 4)))
  }

  cky_surface <- function(df, value_column = "RMSE") {
    x <- sort(unique(df$X1))
    y <- sort(unique(df$X2))
    if (length(x) < 2L || length(y) < 2L ||
        nrow(df) != length(x) * length(y) || anyDuplicated(df[c("X1", "X2")])) {
      stop("Each surface must contain one row for every point on a rectangular grid.")
    }
    z <- matrix(NA_real_, nrow = length(x), ncol = length(y))
    z[cbind(match(df$X1, x), match(df$X2, y))] <- df[[value_column]]
    if (any(!is.finite(z))) stop("Surface contains missing or nonfinite RMSE values.")
    list(x = x, y = y, z = z)
  }

  cky_legend <- function(depths, surface = FALSE) {
    par(fig = c(0, 1, 0, 0.20), new = TRUE, mar = rep(0, 4), xpd = NA)
    plot.new()
    order <- as.vector(matrix(seq_len(ceiling(length(depths) / 3) * 3), ncol = 3, byrow = TRUE))
    depths <- depths[order[order <= length(depths)]]
    cols <- cky_depth_colors[depths]
    # Three columns keep all five depth labels readable at the final panel size.
    args <- list(x = "center", legend = paste("K =", depths), ncol = 3,
                 bty = "n", cex = 1, x.intersp = 0.5, y.intersp = 1.15,
                 seg.len = 1.1)
    if (surface) {
      args$fill <- cols
      args$border <- NA
    } else {
      args$col <- cols
      args$lty <- depths
      args$lwd <- 1.4
    }
    do.call(legend, args)
  }

  cky_pdf <- function(file) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    if (!capabilities("cairo")) stop("R must support cairo to embed publication fonts.")
    cairo_pdf(file, width = 3.2, height = 2.9, pointsize = 11,
              family = "Liberation Sans", onefile = TRUE)
    par(cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1,
        las = 1, mgp = c(1.7, 0.45, 0), tcl = -0.22)
  }

  cky_plot_lines <- function(df, file, ylim, xlab = "x", ylab = "RMSE",
                             xscale = 1, zscale = 1) {
    depths <- cky_depths(df$K)
    cky_pdf(file)
    on.exit(dev.off(), add = TRUE)
    par(fig = c(0, 1, 0.20, 1), mar = c(2.5, 3.2, 0.4, 0.5), xpd = FALSE)
    plot(range(df$X1) / xscale, ylim / zscale, type = "n",
         xlab = xlab, ylab = ylab, yaxs = "i")
    for (k in depths) {
      d <- df[df$K == k, , drop = FALSE]
      if (anyDuplicated(d$X1)) stop("Each depth needs exactly one RMSE per X1 value.")
      d <- d[order(d$X1), , drop = FALSE]
      # Join the supplied evaluation points directly: no interpolation or trimming.
      lines(d$X1 / xscale, d$RMSE / zscale,
            col = cky_depth_colors[k], lty = k, lwd = 1.4)
    }
    cky_legend(depths)
    invisible(file)
  }

  cky_plot_surfaces <- function(df, file, zlim, xlab = "X1", ylab = "X2",
                                zlab = "RMSE", xscale = 1, zscale = 1) {
    depths <- cky_depths(df$K)
    surfaces <- lapply(depths, function(k) cky_surface(df[df$K == k, , drop = FALSE]))
    for (s in surfaces[-1L]) {
      if (!identical(s$x, surfaces[[1]]$x) || !identical(s$y, surfaces[[1]]$y)) {
        stop("All depths must use the same evaluation grid.")
      }
    }
    cky_pdf(file)
    on.exit(dev.off(), add = TRUE)
    par(fig = c(0, 1, 0.20, 1), mar = c(1.7, 1.2, 0.3, 1.5),
        mgp = c(1.4, 0.25, 0), xpd = NA)
    for (i in seq_along(depths)) {
      s <- surfaces[[i]]
      if (i > 1L) par(new = TRUE)
      projection <- persp(s$x / xscale, s$y, s$z / zscale,
            theta = 35, phi = 25, expand = 0.7,
            col = adjustcolor(cky_depth_colors[depths[i]], alpha.f = 0.45),
            shade = 0.20, border = NA,
            ticktype = "detailed", nticks = 3,
            xlab = "", ylab = "", zlab = "",
            zlim = zlim / zscale, axes = FALSE, box = i == 1L)
    }
    # Draw axes in the projected coordinates. Default persp labels overlap the
    # tick labels at the small physical width used in the manuscript.
    xr <- range(surfaces[[1]]$x) / xscale
    yr <- range(surfaces[[1]]$y)
    zr <- zlim / zscale
    usr <- par("usr")
    dx <- diff(usr[1:2]); dy <- diff(usr[3:4])
    ticks <- function(r) {
      t <- pretty(r, n = 2)
      t[t >= r[1] & t <= r[2]]
    }
    labels <- function(t) format(t, trim = TRUE, scientific = FALSE)
    projected <- function(x, y, z) trans3d(x, y, z, projection)
    xt <- ticks(xr); yt <- ticks(yr); zt <- ticks(zr)
    a <- projected(xt, rep(yr[1], length(xt)), rep(zr[1], length(xt)))
    segments(a$x, a$y, a$x - 0.010 * dx, a$y - 0.016 * dy)
    text(a$x - 0.040 * dx, a$y - 0.057 * dy, labels(xt), cex = 1)
    b <- projected(rep(xr[2], length(yt)), yt, rep(zr[1], length(yt)))
    segments(b$x, b$y, b$x + 0.014 * dx, b$y - 0.010 * dy)
    text(b$x + 0.080 * dx, b$y - 0.024 * dy, labels(yt), cex = 1)
    # The zero vertical tick coincides with the horizontal-axis origin; omit its
    # duplicate label while retaining the shared zero lower bound.
    zt <- zt[zt > zr[1]]
    c <- projected(rep(xr[1], length(zt)), rep(yr[1], length(zt)), zt)
    segments(c$x, c$y, c$x - 0.015 * dx, c$y)
    text(c$x - 0.032 * dx, c$y, labels(zt), adj = 1, cex = 1)
    a <- projected(mean(xr), yr[1], zr[1])
    text(a$x - 0.055 * dx, a$y - 0.133 * dy, xlab, srt = -30, cex = 1)
    b <- projected(xr[2], mean(yr), zr[1])
    text(b$x + 0.215 * dx, b$y - 0.022 * dy, ylab, srt = 52, cex = 1)
    mtext(zlab, side = 2, line = 0.0, las = 0, cex = 1)
    cky_legend(depths, surface = TRUE)
    invisible(file)
  }

  output_dir <- cky_path(output_dir, root)
  plot_dir <- cky_path(plot_dir, root)

  inputs <- list()
  for (p in 1:2) for (method in cky_methods) {
    stem <- sprintf("rmse_%s-p%d", gsub("-", "_", method, fixed = TRUE), p)
    path <- file.path(output_dir, paste0(stem, ".csv"))
    if (!file.exists(path)) stop("Missing simulation output: ", path)
    d <- read.csv(path, check.names = FALSE)
    k_cols <- grep("^K[0-9]+$", names(d), value = TRUE)
    if (!length(k_cols)) stop("No depth columns in ", path)
    cky_require_columns(d, c("p", "method", "X1", if (p == 2L) "X2"), path)
    if (anyNA(d$p) || any(d$p != p) || anyNA(d$method) || any(d$method != method)) {
      stop("Dimension or method disagrees with the filename: ", path)
    }
    depths <- cky_depths(as.numeric(sub("^K", "", k_cols)))
    long <- do.call(rbind, lapply(depths, function(k) {
      out <- d[c("X1", if (p == 2L) "X2")]
      out$K <- k
      out$RMSE <- d[[paste0("K", k)]]
      out
    }))
    inputs[[stem]] <- list(p = p, df = long)
  }
  # Use common vertical limits across all six methods within each dimension.
  limits <- lapply(1:2, function(p) {
    cky_limit(unlist(lapply(inputs, function(x) if (x$p == p) x$df$RMSE)))
  })
  for (stem in names(inputs)) {
    x <- inputs[[stem]]
    out <- file.path(plot_dir, paste0(stem, ".pdf"))
    if (x$p == 1L) cky_plot_lines(x$df, out, limits[[1L]])
    else cky_plot_surfaces(x$df, out, limits[[2L]])
    message("Saved ", out)
  }

  paths <- file.path(output_dir, c("rmse_K1to5.csv", "rmse_K1to5_p2.csv"))
  inputs <- lapply(seq_along(paths), function(p) {
    path <- paths[p]
    if (!file.exists(path)) stop("Missing JTPA resampling output: ", path)
    d <- read.csv(path, check.names = FALSE)
    cky_require_columns(d, c("X1", if (p == 2L) "X2", "K", "method", "RMSE"), path)
    cky_depths(d$K)
    if (anyNA(d$method) || !setequal(unique(d$method), cky_methods)) {
      stop("Expected all six NSS/HON methods in ", path)
    }
    d
  })
  limits <- lapply(inputs, function(d) cky_limit(d$RMSE))
  for (p in 1:2) for (method in cky_methods) {
    d <- inputs[[p]][inputs[[p]]$method == method, , drop = FALSE]
    out <- file.path(plot_dir, sprintf("jtpa_p%d_%s.pdf", p, method))
    # Dollar axes are expressed in thousands only for legibility; CSVs stay in dollars.
    if (p == 1L) {
      cky_plot_lines(d, out, limits[[p]], xlab = "X1 ($1000)", ylab = "RMSE ($1000)",
                     xscale = 1000, zscale = 1000)
    } else {
      cky_plot_surfaces(d, out, limits[[p]], xlab = "X1 ($1000)", ylab = "X2 (years)",
                        zlab = "RMSE ($1000)", xscale = 1000, zscale = 1000)
    }
    message("Saved ", out)
  }

  invisible(plot_dir)
}

script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script) != 1L) stop("Run this script with Rscript")
root <- dirname(normalizePath(sub("^--file=", "", script), mustWork = TRUE))
path_setting <- function(name, fallback) {
  path <- path.expand(Sys.getenv(name, fallback))
  if (grepl("^/", path)) path else file.path(root, path)
}
output_dir <- path_setting("CKY_OUTPUT_DIR", "output")
plot_dir <- path_setting("CKY_PLOT_DIR", "plots")
# Validate bundled results before plotting; custom outputs can use other replication counts.
main_output <- normalizePath(output_dir) == normalizePath(file.path(root, "output"))
if (main_output) {
  source(file.path(root, "src", "checks.R"))
  cky_check_outputs(root, output_dir)
}
cky_plots(root, output_dir, plot_dir)
if (main_output) cky_manifest(root, write = TRUE)
