library(fgsea)
library(data.table)
library(ggplot2)
library(io)

source("../R/common.R")

set.seed(1234);

# default coef = 2 is because coef = 1 is the intercept
mast_zstat <- function(zfit, coef=2) {
	# calculate the signed Wald statistic for the continuous component
	# variables: intercept, coefficient
	d.mean.c <- zfit@coefC[, coef];
	d.var.c <- zfit@vcovC[coef, coef, ];
	w.c <- d.mean.c / sqrt(d.var.c)	;

	# repeat same thing for the discrete component
	d.mean.d <- zfit@coefD[, coef];
	d.var.d <- zfit@vcovD[coef, coef, ];
	w.d <- d.mean.d / sqrt(d.var.d)	;

	# combine the continuous and discrete Wald statistics for the hurdle
	w.h <- w.c + w.d;

	data.frame(
		discrete = w.d,
		continuous = w.c,
		hurdle = w.h
	)
}

wd_to_stats <- function(wd, var="hurdle_swald") {
	stats <- wd[[var]];
	if (!is.null(wd$gene)) {
		names(stats) <- wd$gene;
	} else {
		names(stats) <- rownames(wd);
	}
	stats[!is.na(stats)]
}

# ----

h.all <- qread("../msigdb/h.all.v7.5.1.symbols.gmt")$data;
names(h.all) <- rename_hallmarks(names(h.all));

# ----

in.fn <- as.filename("v160_wd_cd8_durable-vs-transient.rds");
wd <- qread(in.fn);

out.fn <- filename(in.fn$fstem, tag=c(in.fn$tag, "gsea"));
pdf.fn <- insert(out.fn, ext="pdf");

gsea.hurdle <- fgsea(h.all, wd_to_stats(wd, "hurdle_swald"), minSize=15, maxSize=500, eps=0.0);
gsea.hurdle <- gsea.hurdle[order(gsea.hurdle$pval), ];
gsea.hurdle[gsea.hurdle$pval < 0.05, ]

gsea.cont <- fgsea(h.all, wd_to_stats(wd, "continuous_swald"), minSize=15, maxSize=500, eps=0.0);
gsea.cont <- gsea.cont[order(gsea.cont$pval), ];
gsea.cont[gsea.cont$pval < 0.05, ]

gsea.disc <- fgsea(h.all, wd_to_stats(wd, "discrete_swald"), minSize=15, maxSize=500, eps=0.0);
gsea.disc <- gsea.disc[order(gsea.disc$pval), ];
gsea.disc[gsea.cont$pval < 0.05, ]

stats.hurdle.w <- sqrt(wd_to_stats(wd, "hurdle_wald"));
gsea.hurdle.w <- fgsea(h.all, stats.hurdle.w, minSize=15, maxSize=500, eps=0.0);
gsea.hurdle.w <- gsea.hurdle.w[order(gsea.hurdle.w$pval), ];
gsea.hurdle.w[gsea.hurdle.w$padj < 0.1, ]

pathways <- c(
	tnfa = "TNF-alpha signaling via NFKB",
	ifng = "Interferon gamma response",
	myc = "MYC targets v1",
	mtorc = "mTORC1 signaling",
	il6 = "IL6 Jak STAT3 signaling",
	il2 = "IL2 STAT5 signaling",
	glycolysis = "Glycolysis",
	e2f = "E2F targets",
	pi3k = "PI3K Akt mTOR signaling",
	ifna = "Interferon alpha response"
);

for (i in 1:length(pathways)) {
	pathway <- pathways[i];
	q <- gsea.hurdle.w$padj[match(pathway, gsea.hurdle.w$pathway)];
	qdraw(
		plotEnrichment(h.all[[pathway]], stats.hurdle.w) + 
			labs(title = sprintf("%s (q = %s)", pathway, format(q, digits=2)))
		,
		width = 4, height = 2,
		file = insert(pdf.fn, tag=c("hurdle-wald", names(pathways)[i]))
	)
}

# ----

in.fn <- as.filename("v160_wd_cd8_durable_c-vs-b.rds");
wd <- qread(in.fn);

out.fn <- filename(in.fn$fstem, tag=c(in.fn$tag, "gsea"));
pdf.fn <- insert(out.fn, ext="pdf");

wstats.hurdle <- wd_to_stats(wd, "hurdle_swald");
gsea.hurdle.w <- fgsea(h.all, wstats.hurdle, minSize=15, maxSize=500, eps=0.0);
gsea.hurdle.w <- gsea.hurdle.w[order(gsea.hurdle.w$pval), ];
gsea.hurdle.w[gsea.hurdle.w$padj < 0.1, ]

pathways <- c(
	myc1 = "MYC targets v1",
	myc2 = "MYC targets v2",
	tnfa = "TNF-alpha signaling via NFKB",
	mtorc = "mTORC1 signaling",
	il2 = "IL2 STAT5 signaling",
	e2f = "E2F targets",
	ifna = "Interferon alpha response",
	ifng = "Interferon gamma response"
);

for (i in 1:length(pathways)) {
	pathway <- pathways[i];
	q <- gsea.hurdle$padj[match(pathway, gsea.hurdle.w$pathway)];
	qdraw(
		plotEnrichment(h.all[[pathway]], wstats.hurdle) + 
			labs(title = sprintf("%s (q = %s)", pathway, format(q, digits=2)))
		,
		width = 4, height = 2,
		file = insert(pdf.fn, tag=c("hurdle-wald", names(pathways)[i]))
	)
}

# ----

zfit.fn <- as.filename("v160_cd8_durable_mast-zlm_c-vs-b.rds");
zfit <- qread(zfit.fn);

out.fn <- filename(in.fn$fstem, tag=c(in.fn$tag, "gsea"));
pdf.fn <- insert(out.fn, ext="pdf");

zstat <- mast_zstat(zfit);

hist(zstat$hurdle, breaks=100)
hist(zstat$continuous, breaks=100)
hist(zstat$discrete, breaks=100)

zstats.hurdle <- wd_to_stats(zstat, "hurdle");
gsea.hurdle.z <- fgsea(h.all, zstats.hurdle, minSize=15, maxSize=500, eps=0.0);
gsea.hurdle.z <- gsea.hurdle.z[order(gsea.hurdle.z$pval), ];
gsea.hurdle.z[gsea.hurdle.z$padj < 0.1, ]

gsea.continuous <- fgsea(h.all, wd_to_stats(zstat, "continuous"), minSize=15, maxSize=500, eps=0.0);
gsea.continuous <- gsea.continuous[order(gsea.continuous$pval), ];
gsea.continuous[gsea.continuous$padj < 0.1, ]

	q <- gsea.hurdle$padj[match(pathway, gsea.hurdle.w$pathway)];
	qdraw(
		plotEnrichment(h.all[[pathway]], zstats.hurdle) + 
			labs(title = sprintf("%s (q = %s)", pathway, format(q, digits=2)))
		,
		width = 4, height = 2,
		file = insert(pdf.fn, tag=c("hurdle-z", names(pathways)[i]))
	)

# ----

zfit.fn <- as.filename("v160_cd8_mast-zlm.rds");
zfit <- qread(zfit.fn);

out.fn <- filename(in.fn$fstem, tag=c(in.fn$tag, "gsea"));
pdf.fn <- insert(out.fn, ext="pdf");

zstat <- mast_zstat(zfit);

zstats.hurdle <- wd_to_stats(zstat, "hurdle");
gsea.hurdle.z <- fgsea(h.all, zstats.hurdle, minSize=15, maxSize=500, eps=0.0);
gsea.hurdle.z <- gsea.hurdle.z[order(gsea.hurdle.z$pval), ];
gsea.hurdle.z[gsea.hurdle.z$padj < 0.25, ]

