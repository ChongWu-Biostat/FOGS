setwd("/home/panwei/wuxx0845/TWAS_casual/Realdata")

library(Rcpp)
library(RcppArmadillo)
library(bigmemory)
library(mvtnorm)
sourceCpp("GMaSPU_support.cpp")
source("aSPU.R")
source("dist_support.R")
source("JointSum.R")
source("SMI.R")
suppressMessages(library('plink2R'))
suppressMessages(library("optparse"))
library(data.table)

args = commandArgs(TRUE)
job = (eval(parse(text = args[[1]])))
pfx = (as.character(args[[2]]))
#conditional.method = (as.character(args[[3]]))
#prune.snp.method = (as.character(args[[8]]))


source("finmap_support.R") # cutoff 0.9
ld.reference = "LHS"

fine.model = "reduced"
weight.input = "net"
conditional.method = "yangCOJO"
causal.prop = 0.05 # input parameter
causal.effct = 0.05 #input parameter

loci.indx = job %% 150 + 1

#loci.indx = job%%24 + 1
outd = paste("/home/panwei/wuxx0845/TWAS_casual/Realdata/", pfx, sep = "")

system(paste("mkdir -p ", outd, sep = ""))


tmp.loci = readRDS("/home/panwei/wuxx0845/TWAS_casual/Realdata/sig_loci_scz3.rds")

tmp.loci[, 1] = as.numeric(gsub("chr", "", tmp.loci[, 1]))

weights = "/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/CMC.twas.com.net.weight.pos"
wgtlist = read.table(weights, head = T, as.is = T)


if (ld.reference == "1000G") {
  ref_ld_chr = "/home/panwei/wuxx0845/TWAS_casual/LDREF/1000G.EUR."
}

if (ld.reference == "LHS") {
  ref_ld_chr = "/home/panwei/wuxx0845/lunghealth/LHSref_1000G_chr"
}

opt = list(
weights = weights,
weights_dir = "/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/",
ref_ld_chr = ref_ld_chr,
fine.model = fine.model)

#res.tmp = matrix(NA,24,2)
#for(loci.indx in 1:24) {
risk.region.used = t(as.matrix(tmp.loci[loci.indx,]))
wgtlist0 = wgtlist[wgtlist[, 3] == risk.region.used[1, 1],]
wgtlist0 = wgtlist0[(wgtlist0$P0 < risk.region.used[1, 3] & wgtlist0$P0 > risk.region.used[1, 2]),]
#    res.tmp[loci.indx,] = c(loci.indx, dim(wgtlist0)[1])
#}


SNP = NULL
for (w in 1:dim(wgtlist0)[1]) {
  wgt.file = paste("/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/", wgtlist0$WGT[w], sep = '')
  load(wgt.file)
  wgt.matrix[is.na(wgt.matrix)] = 0
  tmp = wgt.matrix[, "enet"]
  tmp = tmp[tmp != 0]
  SNP = c(SNP, names(tmp))
}

start.p0 = min(wgtlist0$P0) - 500 * 1000

if (start.p0 < 0) {
  start.p0 = 0
}
start.p1 = max(wgtlist0$P1) + 500 * 1000
start.chr = wgtlist0$CHR[1]

chr.id = risk.region.used[1, 1]
system(paste0("plink --bfile ", opt$ref_ld_chr, chr.id, " --chr ", chr.id, " --from-kb ", start.p0 / 1e3, " --to-kb ", start.p1 / 1e3, " --maf 0.01 --make-bed --out ", outd, "/LDref_tmp_", job))

#system(paste0("plink --bfile /home/panwei/wuxx0845/TWAS_casual/simulation/finemap_chr22  --chr 22 --from-kb ",start.p0/1e3," --to-kb ", start.p1/1e3," --make-bed --out ",outd,"/LDref_tmp_",job))

# Load in reference data
genos = read_plink(paste0(outd, "/LDref_tmp_", job), impute = "avg")

tmp.bim = genos$bim
snp.inf = tmp.bim[, 2]

# randomly slect one gene to be causal

# load summary statistics

sumstat.orgin = readRDS("/home/panwei/wuxx0845/TWAS_casual/Realdata/scz3.rds")
sumstat.orgin = sumstat.orgin[, c(1, 5, 6, 7, 8, 11, 10)]
colnames(sumstat.orgin) = c("SNP", "A1", "A2", "beta", "se", "N", "Z")
# create summary statistics

# Load in list of weights
# TODO : TEST FOR NO HEADER HERE
wgtlist = read.table(weights, head = T, as.is = T)


# Load in reference data
#system(paste0("plink --bfile /home/panwei/wuxx0845/TWAS_casual/simulation/finemap_chr22  --chr 22 --from-kb ",start.p0/1e3," --to-kb ", start.p1/1e3," --make-bed --out ",outd,"/LDref_tmp_",job))


#system(paste0("plink --bfile ", opt$ref_ld_chr,"22 --chr 22 --from-kb ",start.p0/1e3," --to-kb ", start.p1/1e3," --maf 0.01 --make-bed --out ",outd,"/LDref_tmp_",job))

genos = read_plink(paste0(outd, "/LDref_tmp_", job), impute = "avg")

tmp.bim = genos$bim
snp.inf = tmp.bim[, 2]

m = match(snp.inf, sumstat.orgin$SNP)
## For each wgt file:
sumstat.orgin = sumstat.orgin[m,]

#rm(genos)

genos2 = genos
out.res = as.data.frame(matrix(NA, nrow(wgtlist0), 32))

if (dim(wgtlist0)[1] < 2) {
  out.res[1, 24] = -1
  out.file = paste(outd, "/out_", job, ".rds", sep = "")
  saveRDS(out.res, out.file)
}

quantile(sumstat.orgin$Z, na.rm = T)
# input sumstat.orgin
# input
for (w in 1:nrow(wgtlist0)) {
  #out.fun <- function(w) {
  tryCatch({
    sumstat = sumstat.orgin
    genos = genos2 #read_plink(paste0(outd,"/LDref_tmp_",job),impute="avg")

    save.name = paste(outd, "/snp_", loci.indx, "_", w, ".rds", sep = "")
    out = pre.process()

    prune.snp.length.preproc = out$prune.snp.length.preproc
    prune.snp = out$prune.snp
    twas.weight.snp.set = out$twas.weight.snp.set

    if (length(prune.snp) >= 2) {
      saveRDS(out, save.name)
    }
  }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n") })
}




