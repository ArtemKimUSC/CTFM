#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(susieR))

parse_args <- function(args) {
  out <- list(
    wd = getwd(),
    custom_dir = "data/SLDSC_custom",
    sldsc_out = "out/SLDSC_custom",
    out_dir = "out/ctfm_custom",
    L = 10,
    n = 1190321,
    coverage = 0.96,
    min_abs_corr = 0.3,
    alpha_s_threshold = 3,
    alpha_coverage = 0.9
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--wd", "--trait", "--custom-dir", "--sldsc-out", "--out-dir",
                   "--L", "--n", "--coverage", "--min-abs-corr",
                   "--alpha-s-threshold", "--alpha-coverage")) {
      if (i == length(args)) stop("Missing value for ", key)
      value <- args[[i + 1]]
      if (key == "--wd") out$wd <- value
      if (key == "--trait") out$trait <- value
      if (key == "--custom-dir") out$custom_dir <- value
      if (key == "--sldsc-out") out$sldsc_out <- value
      if (key == "--out-dir") out$out_dir <- value
      if (key == "--L") out$L <- as.integer(value)
      if (key == "--n") out$n <- as.numeric(value)
      if (key == "--coverage") out$coverage <- as.numeric(value)
      if (key == "--min-abs-corr") out$min_abs_corr <- as.numeric(value)
      if (key == "--alpha-s-threshold") out$alpha_s_threshold <- as.numeric(value)
      if (key == "--alpha-coverage") out$alpha_coverage <- as.numeric(value)
      i <- i + 2
    } else if (key %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript launch_susie_custom.R --trait TRAIT ",
        "[--wd /path/to/CTFM] [--custom-dir data/SLDSC_custom] ",
        "[--sldsc-out out/SLDSC_custom] [--out-dir out/ctfm_custom]\n",
        sep = ""
      )
      quit(status = 0)
    } else {
      stop("Unknown option: ", key)
    }
  }
  if (is.null(out$trait)) stop("Missing required option: --trait")
  out
}

read_annotation_names <- function(path) {
  first_line <- readLines(path, n = 1)
  has_header <- grepl("(^|\\s)Name(\\s|$)", first_line)
  tab <- read.table(path, header = has_header, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if ("Name" %in% colnames(tab)) return(tab$Name)
  if (ncol(tab) == 1) return(tab[[1]])
  stop("Could not identify annotation names in ", path, "; expected a Name column")
}

read_r_matrix <- function(path) {
  mat <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  mat
}

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

path_from_wd <- function(path, wd) {
  if (grepl("^/", path)) return(path)
  file.path(wd, path)
}

read_sldsc_results <- function(trait_dir, trait) {
  files <- list.files(
    trait_dir,
    pattern = paste0("^", escape_regex(trait), "\\.[^.]+\\.cell_type_results\\.txt$"),
    full.names = TRUE
  )
  if (length(files) == 0) {
    stop("No S-LDSC cell_type_results files found in ", trait_dir)
  }
  files <- sort(files)
  data <- do.call(rbind, lapply(files, function(path) {
    read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }))
  data
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

wd <- normalizePath(args$wd, mustWork = TRUE)
custom_dir <- path_from_wd(args$custom_dir, wd)
trait_dir <- file.path(path_from_wd(args$sldsc_out, wd), args$trait)
out_dir <- path_from_wd(args$out_dir, wd)

annots <- read_annotation_names(file.path(custom_dir, "annots.txt"))
cor_res <- read_r_matrix(file.path(custom_dir, "R.txt"))

common_annots <- annots[annots %in% rownames(cor_res) & annots %in% colnames(cor_res)]
cor_res <- cor_res[common_annots, common_annots, drop = FALSE]

if (nrow(cor_res) == 0) {
  stop("No overlap between annots.txt and R.txt annotation names")
}
missing_from_r <- setdiff(annots, rownames(cor_res))
if (length(missing_from_r) > 0) {
  warning(
    "The following annotations are in annots.txt but missing from R.txt rows/columns: ",
    paste(head(missing_from_r, 20), collapse = ", "),
    ifelse(length(missing_from_r) > 20, " ...", "")
  )
}

data <- read_sldsc_results(trait_dir, args$trait)
data <- data[!duplicated(data$Name), , drop = FALSE]
data$Name <- sub("L2$", "", data$Name)

missing_from_sldsc <- setdiff(rownames(cor_res), data$Name)
if (length(missing_from_sldsc) > 0) {
  stop(
    "The following annotations are in R.txt/annots.txt but missing from S-LDSC output: ",
    paste(head(missing_from_sldsc, 20), collapse = ", "),
    ifelse(length(missing_from_sldsc) > 20, " ...", "")
  )
}

data <- data[data$Name %in% rownames(cor_res), , drop = FALSE]
data <- data[match(rownames(cor_res), data$Name), , drop = FALSE]

if (!all(c("Coefficient", "Coefficient_std_error") %in% colnames(data))) {
  stop("S-LDSC output must contain Coefficient and Coefficient_std_error columns")
}

data$Zscore <- data$Coefficient / data$Coefficient_std_error
data$TRAIT <- args$trait

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.table(
  data,
  file = file.path(out_dir, paste0(args$trait, "_CTFM_sldsc.txt")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

positive <- data[data$Zscore > 0, , drop = FALSE]
if (nrow(positive) < 2) {
  warning("Fewer than two annotations have positive S-LDSC Z-scores; writing empty CT-FM PIP and CS outputs")
  write.table(data.frame(), file = file.path(out_dir, paste0(args$trait, "_CTFM_pips.txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(), file = file.path(out_dir, paste0(args$trait, "_CTFM_CS.txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  quit(status = 0)
}

cor_temp <- cor_res[rownames(cor_res) %in% positive$Name, colnames(cor_res) %in% positive$Name, drop = FALSE]
positive <- positive[match(rownames(cor_temp), positive$Name), , drop = FALSE]

rss <- susie_rss(
  R = cor_temp,
  z = positive$Zscore,
  estimate_residual_variance = FALSE,
  L = args$L,
  n = args$n
)
colnames(rss$alpha) <- positive$Name

pip <- as.data.frame(rss$pip)
pip$Name <- positive$Name
pip$Zscore <- positive$Zscore
pip$Coeff <- positive$Coefficient
pip$TRAIT <- args$trait
pip$INDEX <- seq_len(nrow(pip))
colnames(pip)[1] <- "PIP"

all_cs <- data.frame()
cs_obj <- susie_get_cs(rss, Xcorr = cor_temp, coverage = args$coverage, min_abs_corr = args$min_abs_corr)

if (length(cs_obj$cs) > 0) {
  for (j in seq_along(cs_obj$cs)) {
    component_name <- names(cs_obj$cs)[j]
    component_idx <- suppressWarnings(as.integer(sub("^L", "", component_name)))
    if (is.na(component_idx)) component_idx <- j
    cs <- data.frame(
      CS_Num = component_name,
      Size = length(cs_obj$cs[[j]]),
      annots = cs_obj$cs[[j]],
      TRAIT = args$trait
    )
    cs$annot_name <- pip$Name[match(cs$annots, pip$INDEX)]
    cs$INDEX <- pip$INDEX[match(cs$annot_name, pip$Name)]
    cs$PIP <- pip$PIP[match(cs$annots, pip$INDEX)]
    cs$Zscore <- pip$Zscore[match(cs$annots, pip$INDEX)]
    cs$Alpha <- rss$alpha[component_idx, cs$INDEX]
    cs$Type <- "SuSiE"
    all_cs <- rbind(all_cs, cs)
  }
}

if (nrow(all_cs) >= 1) all_cs$S <- 0

L_keep <- data.frame()
S_all <- data.frame()
if (nrow(rss$alpha) > 0) {
  for (i in seq_len(nrow(rss$alpha))) {
    alpha_i <- rss$alpha[i, ]
    alpha_i <- alpha_i[alpha_i > 0]
    S <- sum(alpha_i * log(alpha_i * ncol(rss$alpha)))
    if (is.finite(S) && S > args$alpha_s_threshold) {
      L_keep <- rbind(L_keep, data.frame(L = paste0("L", i)))
      S_all <- rbind(S_all, data.frame(L = paste0("L", i), S = S))
    }
  }
}

cs_alpha <- data.frame()
if (nrow(L_keep) > 0) {
  for (i in seq_len(nrow(L_keep))) {
    if (!(L_keep$L[i] %in% names(cs_obj$cs))) {
      b <- as.data.frame(rss$alpha[i, ])
      b$annot_name <- rownames(b)
      b <- b[order(b[, 1], decreasing = TRUE), , drop = FALSE]
      total <- 0
      j <- 1
      cs <- data.frame()
      while (total < args$alpha_coverage && j <= nrow(b)) {
        cs <- rbind(cs, b[j, , drop = FALSE])
        total <- total + b[j, 1]
        j <- j + 1
      }
      cs$CS_Num <- L_keep$L[i]
      cs$Size <- nrow(cs)
      cs$PIP <- pip$PIP[match(cs$annot_name, pip$Name)]
      cs$TRAIT <- args$trait
      cs$Zscore <- pip$Zscore[match(cs$annot_name, pip$Name)]
      cs$Type <- "Alpha"
      colnames(cs)[1] <- "Alpha"
      cs$S <- S_all$S[match(cs$CS_Num, S_all$L)]
      cs$INDEX <- pip$INDEX[match(cs$annot_name, pip$Name)]
      cs$annots <- pip$INDEX[match(cs$annot_name, pip$Name)]
      if (nrow(all_cs) > 0) cs <- cs[, colnames(all_cs), drop = FALSE]
      cs_alpha <- rbind(cs_alpha, cs)
    }
  }
}

all_cs <- rbind(all_cs, cs_alpha)

write.table(
  pip,
  file = file.path(out_dir, paste0(args$trait, "_CTFM_pips.txt")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  all_cs,
  file = file.path(out_dir, paste0(args$trait, "_CTFM_CS.txt")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("DONE ", args$trait)
message("Number of annotations in R: ", nrow(cor_res))
message("Number of positive-Z annotations analyzed by SuSiE: ", nrow(positive))
n_cs <- if ("CS_Num" %in% colnames(all_cs)) length(unique(all_cs$CS_Num)) else 0
message("Number of total CS: ", n_cs)


