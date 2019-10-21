setwd("/home/panwei/wuxx0845/TWAS_casual/Realdata")

library(Rcpp)
library(RcppArmadillo)
library(bigmemory)
library(mvtnorm)
sourceCpp("GMaSPU_support.cpp")
source("aSPU.R")
source("dist_support2.R")
source("JointRidge.R")
source("SMI.R")
source("FOCUS_support.R")
suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))
library(data.table)
suppressMessages(library("dplyr"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("mvnfast"))

args <- commandArgs(TRUE)
job <- (eval(parse(text = args[[1]])))
pfx <- (as.character(args[[2]]))
phenotype <- (as.character(args[[3]]))

#job = 62
#pfx = "LHSsnps_9"
#phenotype = "scz3"

conditional.method <- "ridgeCOJO"
ld.reference <- "LHS"
fine.model <- "reduced"
weight.input <- "net"

loci.indx <- job %% 150 + 1


if (pfx == "LHSsnps_9") {
  prune.dir <- "/home/panwei/wuxx0845/TWAS_casual/Realdata/LHSsnps_9"
  ld.reference <- "LHS"
}

if (pfx == "LHSsnps_99") {
  prune.dir <- "/home/panwei/wuxx0845/TWAS_casual/Realdata/LHSsnps_99"
  ld.reference <- "LHS"
}

if (pfx == "1000Gsnps_9") {
  prune.dir <- "/home/panwei/wuxx0845/TWAS_casual/Realdata/1000Gsnps_9"
  ld.reference <- "1000G"
}

if (pfx == "1000Gsnps_99") {
  prune.dir <- "/home/panwei/wuxx0845/TWAS_casual/Realdata/1000Gsnps_99"
  ld.reference <- "1000G"
}

if (ld.reference == "1000G") {
  ref_ld_chr <- "/home/panwei/wuxx0845/TWAS_casual/LDREF/1000G.EUR."
}

if (ld.reference == "LHS") {
  ref_ld_chr <- "/home/panwei/wuxx0845/lunghealth/LHSref_1000G_chr"
}

weights <- "/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/CMC.twas.com.net.weight.pos"

opt <- list(
    weights = weights,
    weights_dir = "/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/",
    ref_ld_chr = ref_ld_chr,
    fine.model = fine.model
)

outd <- paste("/home/panwei/wuxx0845/TWAS_casual/Realdata/", pfx, "_", phenotype, "_rev", sep = "")
system(paste("mkdir -p ", outd, sep = ""))


tmp.loci <- readRDS("/home/panwei/wuxx0845/TWAS_casual/Realdata/sig_loci_scz3.rds")
tmp.loci[, 1] <- as.numeric(gsub("chr", "", tmp.loci[, 1]))


wgtlist <- read.table(weights, head = T, as.is = T)

################################################
#check Can the authors report the reference panel used for the genes reported by FOGS? A CMC-prioritized approach
#was taken, but it isn 't clear how many of the genes identified were initially taken from CMC or "filled in"
# using the other available tissues/panels.
tmp = wgtlist[, "WGT"]
sum(grepl("CMC", tmp))
################################################


risk.region.used <- t(as.matrix(tmp.loci[loci.indx,]))
wgtlist0 <- wgtlist[wgtlist[, 3] == risk.region.used[1, 1],]
wgtlist0 <- wgtlist0[(wgtlist0$P0 < risk.region.used[1, 3] & wgtlist0$P0 > risk.region.used[1, 2]),]

if (dim(wgtlist0)[1] == 0) {
  out.res <- as.data.frame(matrix(NA, 1, 37))
  out.res[1,] <- -1
  out.file <- paste(outd, "/out_", loci.indx, ".rds", sep = "")
  saveRDS(out.res, out.file)
}

SNP <- NULL
for (w in 1:dim(wgtlist0)[1]) {
  wgt.file <- paste("/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/", wgtlist0$WGT[w], sep = "")
  load(wgt.file)
  wgt.matrix[is.na(wgt.matrix)] <- 0
  tmp <- wgt.matrix[, "enet"]
  tmp <- tmp[tmp != 0]
  SNP <- c(SNP, names(tmp))
}

start.p0 <- min(wgtlist0$P0) - 500 * 1000

if (start.p0 < 0) {
  start.p0 <- 0
}
start.p1 <- max(wgtlist0$P1) + 500 * 1000
start.chr <- wgtlist0$CHR[1]

chr.id <- risk.region.used[1, 1]
system(paste0("plink --bfile ", opt$ref_ld_chr, chr.id, " --chr ", chr.id, " --from-kb ", start.p0 / 1e3, " --to-kb ", start.p1 / 1e3, " --maf 0.01 --make-bed --out ", outd, "/LDref_tmp_", job))

genos <- read_plink(paste0(outd, "/LDref_tmp_", job), impute = "avg")

tmp.bim <- genos$bim
snp.inf <- tmp.bim[, 2]

# load summary statistics
if (phenotype == "scz3") {
  sumstat.orgin <- readRDS("/home/panwei/wuxx0845/TWAS_casual/Realdata/scz3.rds")
  # number of varinats that passed the filter critera
  # sumstat.orgin = sumstat.orgin[sumstat.orgin[,2]>0.01,]

  # a1 = sumstat.orgin[,5]

  # a2 = sumstat.orgin[,6]
  # sumstat.orgin = sumstat.orgin[!((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C")) & (a1=="A" | a1=="T" | a1=="G" | a1=="C") & (a2=="A" | a2=="T" | a2=="G" | a2=="C") ,]

  sumstat.orgin <- sumstat.orgin[, c(1, 5, 6, 7, 8, 11, 10)]
}

colnames(sumstat.orgin) <- c("SNP", "A1", "A2", "beta", "se", "N", "Z")

wgtlist <- read.table(weights, head = T, as.is = T)

tmp.bim <- genos$bim
snp.inf <- tmp.bim[, 2]

m <- match(snp.inf, sumstat.orgin$SNP)
## For each wgt file:
sumstat.orgin <- sumstat.orgin[m,]

# rm(genos)
genos2 <- genos
out.res <- as.data.frame(matrix(NA, nrow(wgtlist0), 38))

if (dim(wgtlist0)[1] < 2) {
  out.res[1,] <- -1
  out.file <- paste(outd, "/out_", loci.indx, ".rds", sep = "")
  saveRDS(out.res, out.file)
}

causal.gene.id <- 99 # dummy variable to make sure the code run smoothly

for (w in 1:nrow(wgtlist0)) {
  # out.fun <- function(w) {
  tryCatch({
    sumstat <- sumstat.orgin
    start.time = proc.time()[3]
    genos <- genos2 # read_plink(paste0(outd,"/LDref_tmp_",job),impute="avg")

    save.name <- paste(prune.dir, "/snp_", loci.indx, "_", w, ".rds", sep = "")

    out <- readRDS(save.name)

    prune.snp.length.preproc <- out$prune.snp.length.preproc
    prune.snp <- out$prune.snp
    twas.weight.snp.set <- out$twas.weight.snp.set


    # quality control
    m <- match(genos$bim[, 2], sumstat$SNP)
    sum.missing <- is.na(m)
    sumstat <- sumstat[m,]
    sumstat$SNP <- genos$bim[, 2]
    sumstat$A1[sum.missing] <- genos$bim[sum.missing, 5]
    sumstat$A2[sum.missing] <- genos$bim[sum.missing, 6]

    # QC / allele-flip the input and output
    qc <- allele.qc(sumstat$A1, sumstat$A2, genos$bim[, 5], genos$bim[, 6])

    # Flip Z-scores for mismatching alleles
    sumstat$Z[qc$flip] <- -1 * sumstat$Z[qc$flip]
    sumstat$beta[qc$flip] <- -1 * sumstat$beta[qc$flip]

    sumstat$A1[qc$flip] <- genos$bim[qc$flip, 5]
    sumstat$A2[qc$flip] <- genos$bim[qc$flip, 6]

    # Remove strand ambiguous SNPs (if any)
    if (sum(!qc$keep) > 0) {
      genos$bim <- genos$bim[qc$keep,]
      genos$bed <- genos$bed[, qc$keep]
      sumstat <- sumstat[qc$keep,]
    }

    # make sure prune SNPs is in the summary data and corresponding reference
    cur.genos <- genos$bed
    prune.snp <- prune.snp[prune.snp %in% colnames(cur.genos)]
    twas.weight.snp.set <- twas.weight.snp.set[twas.weight.snp.set %in% colnames(cur.genos)]

    prune.snp <- prune.snp[prune.snp %in% sumstat$SNP]
    twas.weight.snp.set <- twas.weight.snp.set[twas.weight.snp.set %in% sumstat$SNP]

    # make sure prune SNPs is in the summary data and corresponding reference
    cur.genos <- genos$bed
    prune.snp <- prune.snp[prune.snp %in% colnames(cur.genos)]
    twas.weight.snp.set <- twas.weight.snp.set[twas.weight.snp.set %in% colnames(cur.genos)]


    # calcualte the conditional score for each SNPs
    sumstat.final <- NULL
    for (i in 1:length(twas.weight.snp.set)) {
      snp.used <- c(twas.weight.snp.set[i], prune.snp)
      cur.genos <- genos$bed
      cur.genos <- cur.genos[, snp.used]

      sumstat.tmp2 <- sumstat[sumstat[, 1] %in% snp.used,]

      rownames(sumstat.tmp2) <- sumstat.tmp2[, 1]
      sumstat.tmp2 <- sumstat.tmp2[snp.used,]

      B <- sumstat.tmp2[, "beta"]
      S <- sumstat.tmp2[, "se"]
      N <- sumstat.tmp2[, "N"]

      cur.genos <- cur.genos - matrix(rep(colMeans(cur.genos), each = dim(cur.genos)[1]), dim(cur.genos)[1], dim(cur.genos)[2])

      XX <- cov(cur.genos)

      if (conditional.method == "regCOJO") {
        XX <- XX + diag(0.1, dim(XX)[1], dim(XX)[2])
        tmp <- JointSum(B, S, N, XX)
      }

      if (conditional.method == "ridgeCOJO") {
        lambda <- 0.1
        tmp <- JointRidge(B, S, N, XX, lambda)
      }

      beta.tmp <- tmp$beta
      se.tmp <- sqrt(diag(tmp$cov))
      z.tmp <- beta.tmp / se.tmp

      sumstat.tmp2$cond.z <- z.tmp
      sumstat.final <- rbind(sumstat.final, sumstat.tmp2[1,])
    }

    # sumstat.final2 = sumstat.final
    sumstat.tmp <- sumstat.final

    tmp.index <- genos$bim[, 2] %in% sumstat.tmp[, 1]
    genos$bim <- genos$bim[tmp.index,]
    genos$bed <- genos$bed[, tmp.index]

    wgt.file <- paste(opt$weights_dir, "/", wgtlist0$WGT[w], sep = "")
    load(wgt.file)

    # Remove NAs (these should not be here)
    wgt.matrix[is.na(wgt.matrix)] <- 0

    # Match up the SNPs and weights
    m <- match(snps[, 2], genos$bim[, 2])
    m.keep <- !is.na(m)
    snps <- snps[m.keep,]
    wgt.matrix <- wgt.matrix[m.keep,]
    cur.genos <- scale(genos$bed[, m[m.keep]])
    cur.bim <- genos$bim[m[m.keep],]
    # Flip WEIGHTS for mismatching alleles
    qc <- allele.qc(snps[, 5], snps[, 6], cur.bim[, 5], cur.bim[, 6])
    wgt.matrix[qc$flip,] <- -1 * wgt.matrix[qc$flip,]
    rm(snps)

    if (class(wgt.matrix) != "matrix") {
      wgt.matrix <- as.matrix(wgt.matrix)
      wgt.matrix <- t(wgt.matrix)
    }

    cur.FAIL <- FALSE

    # Match up the SNPs and the summary stats
    m <- match(cur.bim[, 2], sumstat.tmp$SNP)

    cur.Z <- sumstat.tmp$cond.z[m]
    Z <- sumstat.tmp$Z[m]

    # Identify the best model -- use enet
    mod.best <- (which.max(cv.performance[1,]))

    if (names(mod.best) == "top1") {
      # cat( "WARNING: top eQTL is the best predictor for this gene, continuing with 2nd-best model\n" )
      mod.best <- names(which.max(cv.performance[1, colnames(cv.performance) != "top1"]))
      mod.best <- which(colnames(cv.performance) == mod.best)
    }

    if (weight.input == "net") {
      mod.best <- "enet"
    }

    if (sum(abs(wgt.matrix[, mod.best])) == 0) {
      tmp.res <- c(wgtlist0$CHR[w], as.character(wgtlist0$ID[w]), wgtlist0$P0[w], wgtlist0$P1[w])
      out.res[w, 1:length(tmp.res)] <- tmp.res
      next
    }

    fine.p <- calc.pvalue(cur.Z)
    org.p <- calc.pvalue(Z)

    tmp.res <- c(wgtlist0$CHR[w], as.character(wgtlist0$ID[w]), wgtlist0$P0[w], wgtlist0$P1[w], length(cur.Z), length(prune.snp), causal.gene.id)
    out.res[w, 1:length(tmp.res)] <- tmp.res
    out.res[w, (length(tmp.res) + 1):(length(tmp.res) + length(fine.p))] <- fine.p
    out.res[w, (length(tmp.res) + length(fine.p) + 1):(length(tmp.res) + length(fine.p) + length(org.p))] <- org.p
    out.res[w, 38] = proc.time()[3] - start.time

    rm(genos)
    rm(cur.genos)
    cat(w)
    # return(out.res)
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })
}

colnames(out.res) <- c("CHR", "ID", "P0", "P1", "n.SNP", "n.condSNP", "casual.id", "afSPU", "fSUM", "fSSU", "afSSU", "fSPU(1)", "fSPU(2)", "fSPU(3)", "fSPU(4)", "fSPU(5)", "fSPU(6)", "fSPU(Inf)", "afSPU", "ZfSum", "aSPU", "SUM", "SSU", "aSSU", "SPU(1)", "SPU(2)", "SPU(3)", "SPU(4)", "SPU(5)", "SPU(6)", "SPU(Inf)", "aSPU", "ZSum", "Focus-Sum", "Focus-SSU", "Focus-aSPU", "Focus-aSSU", "runing_time")

out.file <- paste(outd, "/res_", loci.indx, ".rds", sep = "")
saveRDS(out.res, out.file)


max_iter <- 10

prior_chisq <- 40
prb <- 1e-3
tol <- 2.220446e-14

wgtlist <- wgtlist0
tmp.na <- is.na(out.res[, 33])
causal.gene.id <- causal.gene.id - sum(tmp.na[1:(causal.gene.id - 1)])

out.res <- out.res[!is.na(out.res[, 33]),]
out.res[, "casual.id"] <- causal.gene.id
# cascual id should be channged

chr <- unique(out.res[, 1])
chr <- chr[1]

rownames(out.res) <- paste(out.res[, 2], out.res[, 3], out.res[, 4], sep = "-")
rownames(wgtlist) <- paste(wgtlist[, 2], wgtlist[, 4], wgtlist[, 5], sep = "-")

wgtlist <- wgtlist[rownames(wgtlist) %in% rownames(out.res),]

tryCatch({
  tmp <- FOCUS(out.res[, "SUM"], out.res[, "ZSum"])
  rownames(tmp) <- paste(tmp[, "ID"], tmp[, "P0"], tmp[, "P1"], sep = "-")
  tmp <- tmp[rownames(out.res),]
  out.res[, "Focus-Sum"] <- tmp[, "PIP"]
  out.file <- paste(outd, "/res_", loci.indx, ".rds", sep = "")
  saveRDS(out.res, out.file)
}, error = function(e) {
  cat("ERROR :", conditionMessage(e), "\n")
})
