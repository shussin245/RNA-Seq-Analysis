BiocManager::install("Rsubread")

targets <- read.delim("targets.txt", header=TRUE)
targets
group <- factor(paste0(targets$CellType, ".", targets$Status))

output.files <- sub("^","GSE60450/",sub(".gz", ".subread-align.bam", targets$FileName))
library(Rsubread)
fc <- featureCounts(output.files,annot.ext = "Mus_musculus.GRCm39/Mus_musculus.GRCm39.113.gtf.gz",isGTFAnnotationFile = T)
colnames(fc$counts) <- 1:12

library(edgeR)
y <- DGEList(fc$counts, group=group)
colnames(y) <- targets$GEO

require(org.Mm.eg.db)
Symbol <- mapIds(org.Mm.eg.db, keys=rownames(y), keytype="ENSEMBL",column="SYMBOL")
Symbol[is.na(Symbol)] <- names(Symbol[is.na(Symbol)])
y$genes <- data.frame(Symbol=Symbol)

keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y$samples

y.calc <- calcNormFactors(y) ## same as normLibSizes for this data
y.calc$samples$norm.factors

sapply(1:12,function(x){plotMD(cpm(y.calc, log=TRUE), column=x);abline(h=0, col="red", lty=2, lwd=2)})

points <- c(0,1,2,15,16,17)
colors <- rep(c("blue", "darkgreen", "red"), 2)
library(limma)
plotMDS(y, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design

contrast.matrix <- makeContrasts(
  Lact_vs_Preg_in_B=B.lactate-B.pregnant,
  Lact_vs_Preg_in_L=L.lactate-L.pregnant,     
  Diff_LP=(L.lactate-L.pregnant)-(B.lactate-B.pregnant),
  Preg_vs_Virg_in_B=B.pregnant-B.virgin,
  Preg_vs_Virg_in_L=L.pregnant-L.virgin,     
  Diff_PV=(L.pregnant-L.virgin)-(B.pregnant-B.virgin),
  levels=design)

voom.fit <- voom(y.calc,design,plot=TRUE)
limma.lm.fit <- lmFit(voom.fit, design)
limma.contrasts.fit <- contrasts.fit(limma.lm.fit, contrast.matrix)
limma.ebayes.fit    <- eBayes(limma.contrasts.fit,robust=T)

limma.tt <- topTable(limma.ebayes.fit, coef = 2, number = 100, sort.by = "B") # M == LogFC
library(DT)
datatable(limma.tt)

edgeR.estimateDisp <- estimateDisp(y.calc, design, robust=TRUE)
plotBCV(edgeR.estimateDisp)

edgeR.glmQLF.fit <- glmQLFit(edgeR.estimateDisp, design, robust=TRUE)
plotQLDisp(edgeR.glmQLF.fit)
edgeR.glmQLF.estimate <- glmQLFTest(edgeR.glmQLF.fit, contrast = contrast.matrix)

edgeR.glmQLF.estimate.coef2 <- glmQLFTest(edgeR.glmQLF.fit, contrast=contrast.matrix[,2])
edgeR.tt <- topTags(edgeR.glmQLF.estimate.coef2, n = 100) # by default sorted by Pvalue
datatable(edgeR.tt$table)

limma.tt <- topTable(limma.ebayes.fit, coef = 2, number = Inf)
edgeR.tt <- topTags(edgeR.glmQLF.estimate.coef2, n = Inf)

summary(decideTests(limma.ebayes.fit,p.value=0.1))

summary(decideTests(edgeR.glmQLF.estimate.coef2,p.value=0.1))

hist(limma.tt$P.Value)
hist(edgeR.tt$table$PValue) 

fullcounts <- read.delim("GSE60450/GSE60450_Lactation-GenewiseCounts.txt")

sample_names <- c("MCL1.LA_BC2CTUACXX_GATCAG_L001_R1","MCL1.LB_BC2CTUACXX_TGACCA_L001_R1","MCL1.LC_BC2CTUACXX_GCCAAT_L001_R1","MCL1.LD_BC2CTUACXX_GGCTAC_L001_R1","MCL1.LE_BC2CTUACXX_TAGCTT_L001_R1","MCL1.LF_BC2CTUACXX_CTTGTA_L001_R1","MCL1.DG_BC2CTUACXX_ACTTGA_L002_R1","MCL1.DH_BC2CTUACXX_CAGATC_L002_R1","MCL1.DI_BC2CTUACXX_ACAGTG_L002_R1","MCL1.DJ_BC2CTUACXX_CGATGT_L002_R1","MCL1.DK_BC2CTUACXX_TTAGGC_L002_R1","MCL1.DL_BC2CTUACXX_ATCACG_L002_R1")

ref_ids <- c("GSM1480291","GSM1480292","GSM1480293","GSM1480294","GSM1480295","GSM1480296","GSM1480297","GSM1480298","GSM1480299","GSM1480300","GSM1480301","GSM1480302")

names(ref_ids) <- sample_names

names(fullcounts)[3:14] <- ref_ids[names(fullcounts)[3:14]]

y.full <- DGEList(fullcounts[,3:14], group=group)
library(org.Mm.eg.db)
Symbol <- mapIds(org.Mm.eg.db, keys=as.character(fullcounts$EntrezGeneID), keytype="ENTREZID",column="SYMBOL")
Symbol[is.na(Symbol)] <- names(Symbol[is.na(Symbol)])
y.full$genes <- data.frame(Symbol=Symbol)
keep <- filterByExpr(y.full)
y.full <- y.full[keep, , keep.lib.sizes=FALSE]
y.full$samples
y.fullnorm <- calcNormFactors(y.full) ## same as normLibSizes for this data

sapply(1:2,function(x){plotMD(cpm(y.fullnorm, log=TRUE), column=x);abline(h=0, col="red", lty=2, lwd=2)})
points <- c(0,1,2,15,16,17)
colors <- rep(c("blue", "darkgreen", "red"), 2)
plotMDS(y.fullnorm, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)

voom.fit.full <- voom(y.fullnorm,design,plot=TRUE)
limma.lm.fit.full <- lmFit(voom.fit.full, design)
limma.contrasts.fit.full <- contrasts.fit(limma.lm.fit.full, contrast.matrix)
limma.ebayes.fit.full    <- eBayes(limma.contrasts.fit.full)
limma.tt.full <- topTable(limma.ebayes.fit.full, coef = 2, number = Inf)

edgeR.estimateDisp.full <- estimateDisp(y.fullnorm, design)
plotBCV(edgeR.estimateDisp.full)
edgeR.glmQLF.fit.full <- glmQLFit(edgeR.estimateDisp.full, design)
plotQLDisp(edgeR.glmQLF.fit.full)
edgeR.glmQLF.estimate.coef2.full <- glmQLFTest(edgeR.glmQLF.fit.full, contrast=contrast.matrix[,2])
edgeR.tt.full <- topTags(edgeR.glmQLF.estimate.coef2.full, n = Inf) # by default sorted by Pvalue

summary(decideTests(limma.ebayes.fit.full,p.value=0.1))

summary(decideTests(edgeR.glmQLF.estimate.coef2.full,p.value=0.1))

hist(limma.tt.full$P.Value)
hist(edgeR.tt.full$table$PValue) 

edgeR.glmQLF.estimate.coef1.full <- glmQLFTest(edgeR.glmQLF.fit.full, contrast=contrast.matrix[,1])
edgeR.glmQLF.estimate.coef2.full <- glmQLFTest(edgeR.glmQLF.fit.full, contrast=contrast.matrix[,2])
edgeR.glmQLF.estimate.coef3.full <- glmQLFTest(edgeR.glmQLF.fit.full, contrast=contrast.matrix[,3])
edgeR.glmQLF.estimate.coef4.full <- glmQLFTest(edgeR.glmQLF.fit.full, contrast=contrast.matrix[,4])
edgeR.glmQLF.estimate.coef5.full <- glmQLFTest(edgeR.glmQLF.fit.full, contrast=contrast.matrix[,5])
edgeR.glmQLF.estimate.coef6.full <- glmQLFTest(edgeR.glmQLF.fit.full, contrast=contrast.matrix[,6])

edgeR.tt.full.coef1 <- topTags(edgeR.glmQLF.estimate.coef1.full, n = Inf)
edgeR.tt.full.coef2 <- topTags(edgeR.glmQLF.estimate.coef2.full, n = Inf)
edgeR.tt.full.coef3 <- topTags(edgeR.glmQLF.estimate.coef3.full, n = Inf)
edgeR.tt.full.coef4 <- topTags(edgeR.glmQLF.estimate.coef4.full, n = Inf)
edgeR.tt.full.coef5 <- topTags(edgeR.glmQLF.estimate.coef5.full, n = Inf)
edgeR.tt.full.coef6 <- topTags(edgeR.glmQLF.estimate.coef6.full, n = Inf)

library(FDRestimation)
pi0.coef1 <- get.pi0(edgeR.tt.full.coef1$table$PValue)
pi0.coef1
pi0.coef2 <- get.pi0(edgeR.tt.full.coef2$table$PValue)
pi0.coef2
pi0.coef3 <- get.pi0(edgeR.tt.full.coef3$table$PValue)
pi0.coef3
pi0.coef4 <- get.pi0(edgeR.tt.full.coef4$table$PValue)
pi0.coef4
pi0.coef5 <- get.pi0(edgeR.tt.full.coef5$table$PValue)
pi0.coef5
pi0.coef6 <- get.pi0(edgeR.tt.full.coef6$table$PValue)
pi0.coef6

p.fdr.obj.coef1 <- p.fdr(p=edgeR.tt.full.coef1$table$PValue,set.pi0 = pi0.coef1)
p.fdr.obj.coef2 <- p.fdr(p=edgeR.tt.full.coef2$table$PValue,set.pi0 = pi0.coef2)
p.fdr.obj.coef3 <- p.fdr(p=edgeR.tt.full.coef3$table$PValue,set.pi0 = pi0.coef3)
p.fdr.obj.coef4 <- p.fdr(p=edgeR.tt.full.coef4$table$PValue,set.pi0 = pi0.coef4)
p.fdr.obj.coef5 <- p.fdr(p=edgeR.tt.full.coef5$table$PValue,set.pi0 = pi0.coef5)
p.fdr.obj.coef6 <- p.fdr(p=edgeR.tt.full.coef6$table$PValue,set.pi0 = pi0.coef6)

p.fdr.obj.coef1 <- p.fdr(p=edgeR.tt.full.coef1$table$PValue,set.pi0 = pi0.coef1)
p.fdr.obj.coef2 <- p.fdr(p=edgeR.tt.full.coef2$table$PValue,set.pi0 = pi0.coef2)
p.fdr.obj.coef3 <- p.fdr(p=edgeR.tt.full.coef3$table$PValue,set.pi0 = pi0.coef3)
p.fdr.obj.coef4 <- p.fdr(p=edgeR.tt.full.coef4$table$PValue,set.pi0 = pi0.coef4)
p.fdr.obj.coef5 <- p.fdr(p=edgeR.tt.full.coef5$table$PValue,set.pi0 = pi0.coef5)
p.fdr.obj.coef6 <- p.fdr(p=edgeR.tt.full.coef6$table$PValue,set.pi0 = pi0.coef6)

summary(decideTests(edgeR.glmQLF.estimate.coef1.full,p.value=0.01))
summary(decideTests(edgeR.glmQLF.estimate.coef2.full,p.value=0.01))
summary(decideTests(edgeR.glmQLF.estimate.coef3.full,p.value=0.01))
summary(decideTests(edgeR.glmQLF.estimate.coef1.full,p.value=0.01))
summary(decideTests(edgeR.glmQLF.estimate.coef2.full,p.value=0.01))
summary(decideTests(edgeR.glmQLF.estimate.coef3.full,p.value=0.01))

summary(p.fdr.obj.coef1$fdrs <= 0.01)
summary(p.fdr.obj.coef2$fdrs <= 0.01)
summary(p.fdr.obj.coef3$fdrs <= 0.01)
summary(p.fdr.obj.coef4$fdrs <= 0.01)
summary(p.fdr.obj.coef5$fdrs <= 0.01)
summary(p.fdr.obj.coef6$fdrs <= 0.01)

library(wordcloud2)
termdata <- edgeR.tt.full.coef1$table$Symbol[p.fdr.obj.coef1$fdrs <= 0.01 & edgeR.tt.full.coef1$table$logFC > 0]
freqdata <- edgeR.tt.full.coef1$table$logFC[p.fdr.obj.coef1$fdrs <= 0.01 & edgeR.tt.full.coef1$table$logFC > 0]

termdataFC <- edgeR.tt.full.coef1$table$logFC[p.fdr.obj.coef1$fdrs <= 0.01 & edgeR.tt.full.coef1$table$logFC > 0]
termdata <- termdata[order(termdataFC,decreasing=TRUE)]
freqdata <- freqdata[order(termdataFC,decreasing=TRUE)]

# limit to top 500 terms
termdata <- termdata[1:500]
freqdata <- freqdata[1:500]

# add a title to wordcloud2, from
# https://stackoverflow.com/questions/66957909/add-a-title-and-remove-whitespace-from-wordcloud2
layout(matrix(c(1, 2), nrow = 2), heights = c(1, 1))
par(mar = rep(0, 4))
plot.new()
text(
  x = 0.5,
  y = 0.5,
  "Overexpressed in B cells in Lactating vs Pregnant Mice",
  cex = 1.5,
  font = 2
)

wordcloud2(data = data.frame(word=termdata,freq=freqdata),size=0.2)

library(wordcloud2)
termdata <- edgeR.tt.full.coef2$table$Symbol[p.fdr.obj.coef2$fdrs <= 0.01 & edgeR.tt.full.coef2$table$logFC > 0]
freqdata <- edgeR.tt.full.coef2$table$logFC[p.fdr.obj.coef2$fdrs <= 0.01 & edgeR.tt.full.coef2$table$logFC > 0]

termdataFC <- edgeR.tt.full.coef2$table$logFC[p.fdr.obj.coef2$fdrs <= 0.01 & edgeR.tt.full.coef2$table$logFC > 0]
termdata <- termdata[order(termdataFC,decreasing=TRUE)]
freqdata <- freqdata[order(termdataFC,decreasing=TRUE)]

# limit to top 500 terms
termdata <- termdata[1:500]
freqdata <- freqdata[1:500]

# add a title to wordcloud2, from
# https://stackoverflow.com/questions/66957909/add-a-title-and-remove-whitespace-from-wordcloud2
layout(matrix(c(1, 2), nrow = 2), heights = c(1, 1))
par(mar = rep(0, 4))
plot.new()
text(
  x = 0.5,
  y = 0.5,
  "Overexpressed in L cells in Lactating vs Pregnant Mice",
  cex = 1.5,
  font = 2
)

wordcloud2(data = data.frame(word=termdata,freq=freqdata),size=0.2)

library(wordcloud2)
termdata <- edgeR.tt.full.coef3$table$Symbol[p.fdr.obj.coef3$fdrs <= 0.01 & edgeR.tt.full.coef3$table$logFC > 0]
freqdata <- edgeR.tt.full.coef3$table$logFC[p.fdr.obj.coef3$fdrs <= 0.01 & edgeR.tt.full.coef3$table$logFC > 0]

termdataFC <- edgeR.tt.full.coef3$table$logFC[p.fdr.obj.coef3$fdrs <= 0.01 & edgeR.tt.full.coef3$table$logFC > 0]
termdata <- termdata[order(termdataFC,decreasing=TRUE)]
freqdata <- freqdata[order(termdataFC,decreasing=TRUE)]

# limit to top 500 terms
termdata <- termdata[1:500]
freqdata <- freqdata[1:500]

# add a title to wordcloud2, from
# https://stackoverflow.com/questions/66957909/add-a-title-and-remove-whitespace-from-wordcloud2
layout(matrix(c(1, 2), nrow = 2), heights = c(1, 1))
par(mar = rep(0, 4))
plot.new()
text(
  x = 0.5,
  y = 0.5,
  "Excess Overexpression in L cells Relative to B cells in Lactating vs Pregnant Mice",
  cex = 1.5,
  font = 2
)

wordcloud2(data = data.frame(word=termdata,freq=freqdata),size=0.2)