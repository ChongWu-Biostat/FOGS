pre.process <- function() {
  # Match summary data to input, record NA where summary data is missing
  m = match(genos$bim[, 2], sumstat$SNP)
  sum.missing = is.na(m)
  sumstat = sumstat[m,]
  sumstat$SNP = genos$bim[, 2]
  sumstat$A1[sum.missing] = genos$bim[sum.missing, 5]
  sumstat$A2[sum.missing] = genos$bim[sum.missing, 6]

  # QC / allele-flip the input and output
  qc = allele.qc(sumstat$A1, sumstat$A2, genos$bim[, 5], genos$bim[, 6])

  # Flip Z-scores for mismatching alleles
  sumstat$Z[qc$flip] = -1 * sumstat$Z[qc$flip]
  sumstat$beta[qc$flip] = -1 * sumstat$beta[qc$flip]

  sumstat$A1[qc$flip] = genos$bim[qc$flip, 5]
  sumstat$A2[qc$flip] = genos$bim[qc$flip, 6]

  # Remove strand ambiguous SNPs (if any)
  if (sum(!qc$keep) > 0) {
    genos$bim = genos$bim[qc$keep,]
    genos$bed = genos$bed[, qc$keep]
    sumstat = sumstat[qc$keep,]
  }

  ## get the pruned SNP set
  prune.cutoff = 0.9
  if (opt$fine.model == "full") {
    tmp = genos$bim[, 4]
    tmp.index = !is.na(sumstat$Z)
    sumstat.tmp = sumstat[tmp.index,]

    snp.list.name = paste(outd, "/", wgtlist0$ID[w], "_", job, "_SNPlist.txt", sep = "")
    write.table(sumstat.tmp[, 1], snp.list.name, quote = FALSE, row.names = FALSE, col.names = FALSE)

    system(paste0("plink --bfile ", outd, "/LDref_tmp_", job, " --extract ", snp.list.name, " --make-bed --out ", outd, "/", wgtlist0$ID[w], "_", job))
    system(paste0("plink --bfile ", outd, "/", wgtlist0$ID[w], "_", job, " --indep-pairwise 1000 1 ", prune.cutoff, " --out ", outd, "/", wgtlist0$ID[w]))
    prune.snp = read.table(paste(outd, "/", wgtlist0$ID[w], ".prune.in", sep = ""), stringsAsFactors = FALSE)
    prune.snp = prune.snp[, 1]
    system(paste0("rm ", outd, "/", wgtlist0$ID[w], "*"))
  }

  if (opt$fine.model == "reduced") {
    prune.snp = NULL
    index.tmp = 1:dim(wgtlist0)[1]
    index.tmp = index.tmp[-w]

    for (fine.index in index.tmp) {
      wgt.file = paste(opt$weights_dir, "/", wgtlist0$WGT[fine.index], sep = '')
      load(wgt.file)
      wgt.matrix[is.na(wgt.matrix)] = 0

      if (weight.input == "optimal") {
        mod.best = (which.max(cv.performance[1,]))

        if (names(mod.best) == "top1") {
          # cat( "WARNING: top eQTL is the best predictor for this gene, continuing with 2nd-best model\n" )
          mod.best = names(which.max(cv.performance[1, colnames(cv.performance) != "top1"]))
          mod.best = which(colnames(cv.performance) == mod.best)
        }
      }

      if (weight.input == "net") {
        mod.best = "enet"
      }

      tmp.weight = wgt.matrix[, mod.best]
      tmp.weight = tmp.weight[tmp.weight != 0]

      snp.set = tmp.weight
      snp.set.tmp = names(snp.set)
      prune.snp = c(prune.snp, snp.set.tmp)
    }

    prune.snp = unique(prune.snp)

    wgt.file = paste(opt$weights_dir, "/", wgtlist0$WGT[w], sep = '')
    load(wgt.file)
    wgt.matrix[is.na(wgt.matrix)] = 0

    if (weight.input == "optimal") {

      mod.best = (which.max(cv.performance[1,]))

      if (names(mod.best) == "top1") {
        # cat( "WARNING: top eQTL is the best predictor for this gene, continuing with 2nd-best model\n" )
        mod.best = names(which.max(cv.performance[1, colnames(cv.performance) != "top1"]))
        mod.best = which(colnames(cv.performance) == mod.best)
      }
    }

    if (weight.input == "net") {
      mod.best = "enet"
    }

    tmp.weight = wgt.matrix[, mod.best]
    tmp.weight = tmp.weight[tmp.weight != 0]

    prune.snp = prune.snp[!prune.snp %in% names(tmp.weight)]
    prune.snp = prune.snp[prune.snp %in% sumstat[, 1]]

    snp.list.name = paste(outd, "/", wgtlist0$ID[w], "_", job, "_SNPlist.txt", sep = "")
    write.table(prune.snp, snp.list.name, quote = FALSE, row.names = FALSE, col.names = FALSE)

    system(paste0("plink --bfile ", outd, "/LDref_tmp_", job, " --extract ", snp.list.name, " --make-bed --out ", outd, "/", wgtlist0$ID[w], "_", job))
    system(paste0("plink --bfile ", outd, "/", wgtlist0$ID[w], "_", job, " --indep-pairwise 1000 1 ", prune.cutoff, " --out ", outd, "/", wgtlist0$ID[w]))

    prune.snp = read.table(paste(outd, "/", wgtlist0$ID[w], ".prune.in", sep = ""), stringsAsFactors = FALSE)

    prune.snp = prune.snp[, 1]
    system(paste0("rm ", outd, "/", wgtlist0$ID[w], "*"))
  }

  prune.snp.length.preproc = length(prune.snp)

  wgt.file = paste(opt$weights_dir, "/", wgtlist0$WGT[w], sep = '')
  load(wgt.file)
  wgt.matrix[is.na(wgt.matrix)] = 0

  if (weight.input == "optimal") {

    mod.best = (which.max(cv.performance[1,]))

    if (names(mod.best) == "top1") {
      mod.best = names(which.max(cv.performance[1, colnames(cv.performance) != "top1"]))
      mod.best = which(colnames(cv.performance) == mod.best)
    }
  }

  if (weight.input == "net") {
    mod.best = "enet"
  }

  tmp.weight = wgt.matrix[, mod.best]
  tmp.weight = tmp.weight[tmp.weight != 0]

  snp.set = tmp.weight
  snp.set = snp.set[names(snp.set) %in% sumstat[, 1]]
  twas.weight.snp.set = names(snp.set)

  prune.snp = prune.snp[!prune.snp %in% names(snp.set)]

  prune.snp = unique(prune.snp)

  # remove the SNPs with NAN, Inf, etc...
  qual.sumstat = sumstat[sumstat[, 1] %in% prune.snp,]
  qual.sumstat = qual.sumstat[!duplicated(qual.sumstat[, 1]),]
  rownames(qual.sumstat) = qual.sumstat[, 1]
  qual.sumstat = qual.sumstat[prune.snp,]
  prune.snp = prune.snp[rowSums(is.na(qual.sumstat)) == 0]

  qual.sumstat = sumstat[sumstat[, 1] %in% twas.weight.snp.set,]
  qual.sumstat = qual.sumstat[!duplicated(qual.sumstat[, 1]),]
  rownames(qual.sumstat) = qual.sumstat[, 1]
  qual.sumstat = qual.sumstat[twas.weight.snp.set,]

  twas.weight.snp.set = twas.weight.snp.set[rowSums(is.na(qual.sumstat)) == 0]

  if (length(snp.set) <= 1) {
    tmp.res = c(wgtlist0$CHR[w], as.character(wgtlist0$ID[w]), wgtlist0$P0[w], wgtlist0$P1[w])
    out.res[w, 1:length(tmp.res)] = tmp.res
    next # if no twas snp lies in the lipidis data, test next gene
  }

  # remove pruned SNPs set that are highly correlated with the twas.weight.snp.set
  out = list(prune.snp.length.preproc = prune.snp.length.preproc, prune.snp = prune.snp, twas.weight.snp.set = twas.weight.snp.set)
  return(out)
}
