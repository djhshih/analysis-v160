library(fgsea)
library(data.table)
library(ggplot2)
library(io)
library(RColorBrewer)
library(dplyr)

source("../R/common.R")

set.seed(1337);

mast_zstat <- function(zfit, coef) {
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

# contrast coef1 vs. coef0
mast_zstat_diff <- function(zfit, coef0, coef1) {
	# calculate the signed Wald statistic for the continuous component
	# variables: intercept, coefficient
	components <- c(C="C", D="D");
	w <- lapply(components,
		function(comp) {
			d.mean <- coef(zfit, comp)[, coef1] - coef(zfit, comp)[, coef0];
			v <- vcov(zfit, comp);
			d.var <- v[coef0, coef0, ] + v[coef1, coef1, ] - 2 * v[coef0, coef1, ];
			d.mean / sqrt(d.var)
		}
	);

	# combine the continuous and discrete Wald statistics for the hurdle
	w$H <- w$C + w$D;

	data.frame(
		discrete = w$D,
		continuous = w$C,
		hurdle = w$H
	)
}

make_stats <- function(wd, var) {
	stats <- wd[[var]];
	if (!is.null(wd$gene)) {
		names(stats) <- wd$gene;
	} else {
		names(stats) <- rownames(wd);
	}
	stats[!is.na(stats)]
}

plot_fgsea_col <- function(d) {
	# order pathway by FDR of hurdle component
	d$pathway <- factor(d$pathway, levels=c(unique(d$pathway[order(d$NES)])));

	cols <- brewer.pal(9, "RdYlBu");
	col.high <- cols[1];
	col.low <- cols[9];
	col.mid <- cols[5];

	q.limits <- c(0.25, 1e-3);

	ggplot(d, aes(x=pathway, y=NES, fill=(NES > 0), alpha=bound(padj, q.limits))) + 
		theme_classic() +
		geom_col() + 
		scale_fill_manual(name="", values=c(col.low, col.high), guide=FALSE) +
		scale_alpha_continuous(trans=revlog_trans(), name="FDR", limits=q.limits) +
		theme(legend.position="bottom") +
		coord_flip() +
		xlab("") + ylab("normalized enrichment score")
}


plot_fgseas <- function(gseas, gsets, target=NULL, fdr.cut=0.1, z.cut=1) {
	d <- do.call(rbind, mapply(
		function(gsea, g) {
			select(gsea, pathway, NES, pval, padj) %>%
				mutate(group = g)
		}
		,
		gseas,
		names(gseas),
		SIMPLIFY = FALSE
	));
	d$group <- factor(d$group, levels=names(gseas));

	d.f <- filter(d, pathway %in% unlist(gsets));

	if (is.null(target)) {
		# use the last group as the target
		target <- length(gseas);
	}

	# order set by effect size of target group
	d.f$pathway <- factor(d.f$pathway, levels=c(unique(gseas[[target]]$pathway[order(gseas[[target]]$NES)])));

	cols <- brewer.pal(9, "RdYlBu");
	col.high <- cols[1];
	col.low <- cols[9];
	col.mid <- cols[5];

	z.limits <- c(-3, 3);
	q.limits <- c(0.25, 1e-6);

	ggplot(d.f, aes(y=pathway, x=group, colour=bound(NES, z.limits), size=bound(padj, q.limits))) + 
		theme_classic() +
		geom_point() + 
		scale_colour_gradient2(name="", low=col.low, mid=col.mid, high=col.high, limits=z.limits) +
		scale_size_continuous(trans=revlog_trans(), range=c(0, 4), name="FDR", limits=q.limits) +
		theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
		theme(legend.position="right") +
		coord_fixed() +
		xlab("") + ylab("")
}



# ---

fdr.cut <- 0.1;

my_fgsea <- function(gset, stat, min.size=15, max.size=500) {
	d <- fgsea(gset, stat, minSize=min.size, maxSize=max.size, eps=0.0)
	d[order(d$pval), ];
}

h.all <- qread("../msigdb/h.all.v7.5.1.symbols.gmt")$data;
names(h.all) <- rename_hallmarks(names(h.all));

pathways.f <- c(
	myc1 = "MYC targets v1",
	myc2 = "MYC targets v2",
	tnfa = "TNF-alpha signaling via NFKB",
	mtorc = "mTORC1 signaling",
	il2 = "IL2 STAT5 signaling",
	pi3k = "PI3K Akt mTOR signaling",
	glycolysis = "Glycolysis",
	e2f = "E2F targets",
	ifna = "Interferon alpha response",
	ifng = "Interferon gamma response"
);

gset.name <- "hallmark";
gset <- h.all;

# ---

zfit.fn <- as.filename("v160_cd8_mast-zlm.rds");
plot.height <- 5;

out.fn <- filename(zfit.fn$fstem, tag=c(setdiff(zfit.fn$tag, "mast-zlm"), "fgsea"));
pdf.fn <- insert(out.fn, ext="pdf");

zfit <- qread(zfit.fn);

zstat.durable <- mast_zstat(zfit, "responsedurable");
stats.durable.h <- make_stats(zstat.durable, "hurdle");
gsea.durable.h.z <- my_fgsea(gset, stats.durable.h);
gsea.durable.h.z[gsea.durable.h.z$padj < fdr.cut, ]

zstat.transient <- mast_zstat(zfit, "responsetransient");
stats.transient.h <- make_stats(zstat.transient, "hurdle");
gsea.transient.h.z <- my_fgsea(gset, stats.transient.h);
gsea.transient.h.z[gsea.transient.h.z$padj < fdr.cut, ]

zstat.diff <- mast_zstat_diff(zfit, "responsetransient", "responsedurable");
stats.diff.h <- make_stats(zstat.diff, "hurdle");
gsea.diff.h.z <- my_fgsea(gset, zstats.diff.h);
gsea.diff.h.z[gsea.diff.h.z$padj < fdr.cut, ]

qdraw(
	plot_fgsea_col(gsea.diff.h.z[gsea.diff.h.z$padj < fdr.cut, ]),
	width = 5, height = plot.height,
	file = insert(pdf.fn, tag=c("hurdle-z", "col", "durable-vs-transient"))
)

pathways.sig <- gsea.diff.h.z$pathway[gsea.diff.h.z$padj < fdr.cut];
pathways <- intersect(pathways.f, pathways.sig);
names(pathways) <- names(pathways.f)[match(pathways, pathways.f)];

for (i in 1:length(pathways)) {
	pathway <- pathways[i];
	q <- gsea.diff.h.z$padj[match(pathway, gsea.diff.h.z$pathway)];
	qdraw(
		plotEnrichment(gset[[pathway]], stats.diff.h) + 
			labs(title = sprintf("%s (q = %s)", pathway, format(q, digits=2)))
		,
		width = 4, height = 2,
		file = insert(pdf.fn, tag=c("hurdle-z", names(pathways)[i], "durable-vs-transient"))
	)
}

gseas.response <- list(
	"transient" = gsea.transient.h.z,
	"durable" = gsea.durable.h.z,
	"durable vs. transient" = gsea.diff.h.z
);

qdraw(
	plot_fgseas(gseas.response, pathways.sig),
	width=6, height=4,
	file=insert(pdf.fn, c("durable-vs-transient", gset.name))
);


# ---

zfit.fn <- as.filename("v160_cd8_durable_mast-zlm_c-vs-b.rds");
plot.height <- 4;

out.fn <- filename(zfit.fn$fstem, tag=c(setdiff(zfit.fn$tag, "mast-zlm"), "fgsea"));
pdf.fn <- insert(out.fn, ext="pdf");

zfit <- qread(zfit.fn);

zstat <- mast_zstat(zfit, "cluster2C");

zstats.h <- make_stats(zstat, "hurdle");
gsea.h.z <- my_fgsea(gset, zstats.h);
gsea.h.z[gsea.h.z$padj < fdr.cut, ]

zstats.c <- make_stats(zstat, "continuous");
gsea.c.z <- my_fgsea(gset, zstats.c);
gsea.c.z[gsea.c.z$padj < fdr.cut, ]

zstats.d <- make_stats(zstat, "discrete");
gsea.d.z <- my_fgsea(gset, zstats.d);
gsea.d.z[gsea.d.z$padj < fdr.cut, ]

qdraw(
	plot_fgsea_col(gsea.h.z[gsea.h.z$padj < fdr.cut, ]),
	width = 5, height = plot.height,
	file = insert(pdf.fn, tag=c("hurdle-z", "col", "c-vs-b"))
)

pathways.sig <- gsea.h.z$pathway[gsea.h.z$padj < fdr.cut];
pathways <- intersect(pathways.f, pathways.sig);
names(pathways) <- names(pathways.f)[match(pathways, pathways.f)];

for (i in 1:length(pathways)) {
	pathway <- pathways[i];
	q <- gsea.h.z$padj[match(pathway, gsea.h.z$pathway)];
	qdraw(
		plotEnrichment(gset[[pathway]], zstats.h) + 
			labs(title = sprintf("%s (q = %s)", pathway, format(q, digits=2)))
		,
		width = 4, height = 2,
		file = insert(pdf.fn, tag=c("hurdle-z", names(pathways)[i], "c-vs-b"))
	)
}

zfit.b <- qread("v160_cd8_mast-zlm_durable-b.rds");
zfit.c <- qread("v160_cd8_mast-zlm_durable-c.rds");

zstat.b <- mast_zstat(zfit.b, "durableB");
zstat.c <- mast_zstat(zfit.c, "durableC");

zstats.b.h <- make_stats(zstat.b, "hurdle");
gsea.b.h.z <- my_fgsea(gset, zstats.b.h);
gsea.b.h.z[gsea.b.h.z$padj < fdr.cut, ]

zstats.c.h <- make_stats(zstat.c, "hurdle");
gsea.c.h.z <- my_fgsea(gset, zstats.c.h);
gsea.c.h.z[gsea.h.z$padj < fdr.cut, ]

gseas.bc <- list(
	"durable B" = gsea.b.h.z,
	"durable C" = gsea.c.h.z,
	"durable C vs. durable B" = gsea.h.z
);

qdraw(
	plot_fgseas(gseas.bc, pathways.sig),
	width=6, height=4,
	file=insert(pdf.fn, c("c-vs-b", gset.name))
);


# ---

sce <- qread("v160_sce.rds");

out.fn <- filename("v160", tag=c("expr", "response2"));
pdf.fn <- insert(out.fn, ext="pdf");

idx <- colData(sce)$t_cell == "CD8" & colData(sce)$response %in% c("transient", "durable") &
	(!is.na(colData(sce)$cluster2) | colData(sce)$response == "transient");
sce.sel <- sce[, idx];

with(colData(sce.sel), table(response, cluster2, useNA="always"))
colData(sce.sel)$group <- factor(NA, levels=c("transient", "durable B", "durable C"));
colData(sce.sel)$group[colData(sce.sel)$response == "transient"] <- "transient";
colData(sce.sel)$group[colData(sce.sel)$response == "durable" & colData(sce.sel)$cluster2 == "B"] <- "durable B";
colData(sce.sel)$group[colData(sce.sel)$response == "durable" & colData(sce.sel)$cluster2 == "C"] <- "durable C";

cols.response2 <- c(nonresponsive="grey60", transient="#DF8F44", "durable B"="#00C2CC" , "durable C"="#0073D1");

genes <- c(
	"CCL4", "CCL5", "XCL1", "XCL2", "IFNG", "GZMB", "PGAM1",
	"MAP2K3", "SLC7A5", "NAMPT", "M6PR",
	"IRF1", "RIPK2", "GBP4", "CASP8"
);

for (gene in genes) {
	qdraw(
		gene_expr_plot(sce.sel, gene, "group", cols = cols.response2),
		width = 1.25, height = 4,
		file = insert(pdf.fn, c("gene", tolower(gene)))
	)
}

