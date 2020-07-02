#!/usr/bin/env Rscript

options(stringsAsFactors = F)

library(argparse)
library(data.table)
library(dplyr)
library(stringr)

logsum <- function(x, base = exp(1), na.rm = FALSE) {
  xm <- max(x)
  return(xm + log(sum(exp(x - xm), na.rm = na.rm), base = base))
}

# cf. https://github.com/annahutch/corrcoverage/blob/master/R/posterior_prob_functions.R
# Note: W is variance here, not sd
ABF <- function(beta, stderr, W = 0.04) {
  z <- beta / stderr
  V <- stderr^2
  r <- W / (W + V)
  lbf <- 0.5 * (log(1 - r) + (r * z^2))
  denom <- logsum(lbf)
  prob <- exp(lbf - denom)
  return(data.frame(lbf, prob))
}

get_cs <- function(variant, prob, coverage = 0.95) {
  ordering <- order(prob, decreasing = T)
  idx <- which(cumsum(prob[ordering]) > coverage)[1]
  cs <- variant[ordering][1:idx]
  return(cs)
}

main <- function(args) {
  df <- data.table::fread(args$z, data.table = F) %>%
    dplyr::mutate(cs = -1, cs_99 = -1)
  cols <- colnames(df)

  if (!("v" %in% cols)) {
    df <- dplyr::mutate(df, v = stringr::str_c(chromosome, position, allele1, allele2, sep = ":"))
  }

  ret <- ABF(df$beta, df$se, args$prior_variance) %>%
    dplyr::mutate(v = df$v)
  cs <- get_cs(ret$v, ret$prob, coverage = 0.95)
  cs_99 <- get_cs(ret$v, ret$prob, coverage = 0.99)

  df <- dplyr::left_join(df, ret)
  df$cs[df$v %in% cs] <- 1
  df$cs_99[df$v %in% cs_99] <- 1
  df <- dplyr::select(df, -v)

  cred <- data.frame(
    cs = 1,
    cs_log10bf = logsum(ret$lbf, 10),
    cs_size = length(cs)
  )

  write.table(df, args$snp, sep = "\t", row.names = F, quote = F)
  write.table(cred, args$cred, sep = "\t", row.names = F, quote = F)
}

parser <- ArgumentParser()

parser$add_argument("--prefix", type = "character")
parser$add_argument("--z", "-z", type = "character", required = TRUE)
parser$add_argument("--out", type = "character")
parser$add_argument("--snp", type = "character")
parser$add_argument("--cred", type = "character")
parser$add_argument("--log", type = "character")
parser$add_argument("--prior-variance", "-W", type = "double", default = 0.04)

args <- parser$parse_args()

if (is.null(args$prefix)) {
  args$prefix <- tools::file_path_sans_ext(args$z)
}

if (is.null(args$out)) {
  args$out <- "tmp"
}
if (is.null(args$snp)) {
  args$snp <- paste0(args$out, ".abf.snp")
}
if (is.null(args$cred)) {
  args$cred <- paste0(args$out, ".abf.cred")
}
if (is.null(args$log)) {
  args$log <- paste0(args$out, ".abf.log")
}

logfile <- file(args$log, open = "w")
sink(logfile, type = "output", split = TRUE)
sink(logfile, type = "message")
print(args)

if (is.null(args$out) & any(sapply(list(args$snp, args$cred, args$log), is.null))) {
  stop("Either --out or all of --snp, --cred, and --log should be specified.")
}

print("Analysis started")
tryCatch(
  {
    main(args)
  },
  error = function(e) {
    sink()
    message(as.character(e))
    sink(type = "message")
    stop(e)
  }
)

print("Finished!")
