#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The data.table R package is required. Install it with install.packages('data.table').")
  }
})

parse_args <- function(args) {
  out <- list(chromosomes = "1-22", denominator = "n_minus_1", strip_prefix = character())
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--a-dir", "--b-dir", "--out", "--chromosomes", "--denominator", "--strip-prefix")) {
      if (i == length(args)) stop("Missing value for ", key)
      value <- args[[i + 1]]
      if (key == "--a-dir") out$a_dir <- value
      if (key == "--b-dir") out$b_dir <- value
      if (key == "--out") out$out <- value
      if (key == "--chromosomes") out$chromosomes <- value
      if (key == "--denominator") out$denominator <- value
      if (key == "--strip-prefix") out$strip_prefix <- c(out$strip_prefix, value)
      i <- i + 2
    } else if (key %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript make_adjusted_ldscore_correlation.R ",
        "--a-dir matrices/A --b-dir matrices/B --out R.txt ",
        "[--chromosomes 1-22] [--denominator n_minus_1|n] [--strip-prefix PREFIX]\n",
        sep = ""
      )
      quit(status = 0)
    } else {
      stop("Unknown option: ", key)
    }
  }
  required <- c("a_dir", "b_dir", "out")
  missing <- required[!required %in% names(out)]
  if (length(missing) > 0) stop("Missing required option(s): ", paste(missing, collapse = ", "))
  if (!out$denominator %in% c("n_minus_1", "n")) stop("--denominator must be n_minus_1 or n")
  out
}

parse_chromosomes <- function(value) {
  if (grepl("-", value) && !grepl("[, ]", value)) {
    bounds <- strsplit(value, "-", fixed = TRUE)[[1]]
    return(as.character(seq.int(as.integer(bounds[[1]]), as.integer(bounds[[2]]))))
  }
  strsplit(value, "[, ]+")[[1]]
}

read_matrix_dir <- function(dir, chromosomes) {
  chunks <- list()
  for (chr in chromosomes) {
    path_gz <- file.path(dir, paste0("L2.", chr, ".txt.gz"))
    path_txt <- file.path(dir, paste0("L2.", chr, ".txt"))
    path <- if (file.exists(path_gz)) path_gz else path_txt
    if (!file.exists(path)) stop("Missing matrix file: ", path)
    message("[make_adjusted_ldscore_correlation] reading ", path)
    chunks[[chr]] <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  }
  data <- do.call(rbind, chunks)
  as.matrix(data)
}

strip_prefixes <- function(names, prefixes) {
  out <- names
  for (prefix in prefixes) {
    out <- sub(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", prefix)), "", out)
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
chromosomes <- parse_chromosomes(args$chromosomes)

A <- read_matrix_dir(args$a_dir, chromosomes)
B <- read_matrix_dir(args$b_dir, chromosomes)

if (nrow(A) != nrow(B)) {
  stop("A and B have different numbers of SNP rows: ", nrow(A), " vs ", nrow(B))
}

storage.mode(A) <- "double"
storage.mode(B) <- "double"

message("[make_adjusted_ldscore_correlation] residualizing ", ncol(A), " annotations on ", ncol(B), " background annotations")
X <- cbind(intercept = 1, B)
qr_x <- qr(X)
coef <- qr.coef(qr_x, A)
A_resid <- A - X %*% coef

sd_resid <- apply(A_resid, 2, stats::sd)
keep <- is.finite(sd_resid) & sd_resid > 0
if (!all(keep)) {
  warning("Dropping ", sum(!keep), " annotations with zero or non-finite residual variance")
  A_resid <- A_resid[, keep, drop = FALSE]
}

A_scaled <- scale(A_resid)
denominator <- if (args$denominator == "n") nrow(A_scaled) else nrow(A_scaled) - 1
R <- crossprod(A_scaled) / denominator
diag(R) <- 1

rn <- strip_prefixes(colnames(A_scaled), args$strip_prefix)
rownames(R) <- rn
colnames(R) <- rn

dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)
write.table(R, file = args$out, sep = "\t", quote = FALSE, col.names = NA)
message("[make_adjusted_ldscore_correlation] wrote ", args$out)
