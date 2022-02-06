library(DropletUtils)
library(io)
library(scuttle)  # logNormCounts
library(MAST)
library(dplyr)
library(randomForest)
library(scran)
library(miQC)  # requires Bioconductor 3.14 and R 4.1
library(scater)
library(binom)
library(ggplot2)
library(ggsci)
library(bluster)


source("R/common.R")


options(mc.cores=100);

set.seed(1337);


m_fix_table <- function(d) {
	d$label <- factor(d$label, levels=sort(as.integer(unique(d$label))));

	d$lower <- pmax(0, d$lower);
	d$upper <- pmin(1, d$upper);

	d$response <- factor(d$response, levels=levels.response);

	d
}

mod_rd_plot <- function(g) {
	g + theme(legend.position="bottom", legend.box = "horizontal") + coord_fixed()
}


indir <- "../aggr/5p/pt26/outs/count"

out.fn <- filename("v160");
rds.fn <- insert(out.fn, ext="rds");
pdf.fn <- insert(out.fn, ext="pdf");

barcode.d <- qread("../../tcr-profiling/tcr/merged/tcr_aggr_responsive_barcodes.csv");

x <- read10xCounts(file.path(indir, "filtered_feature_bc_matrix"));
mcold <- DataFrame(rowData(x));

# sort by ensembl ID
# so that the earlier entry will be used when we match on gene symbol
#ens.order.idx <- order(rowData(x)$ID);

# deal with duplicate gene symbol entries
dup.idx <- which(duplicated(rowData(x)$Symbol));
cat("Duplicated gene symbols:\n")
print(rowData(x)$Symbol[dup.idx])
rowData(x)$Symbol[dup.idx] <- paste(rowData(x)$Symbol[dup.idx], rowData(x)$ID[dup.idx], sep="_")

mcold <- DataFrame(rowData(x));
rownames(x) <- mcold$Symbol;

base.gset <- read.table("../annot/gset/gobp_immune_response.grp", skip=2)[,1];
tact.gset <- read.table("../annot/gset/gobp_t_cell_activation.grp", skip=2)[,1];
pbmc.markers.d <- qread("../annot/pbmc_markers.csv");
tmem.markers.d <- qread("../annot/t-mem_markers.tsv");

tcell.gset <- unique(sort(c(
	tact.gset,
	pbmc.markers.d$gene_symbol[pbmc.markers.d$cell_type == "T cell"],
	tmem.markers.d$gene_symbol
)));

base.gset.idx <- na.omit(match(base.gset, mcold$Symbol));
tcell.gset.idx <- na.omit(match(tcell.gset, mcold$Symbol));
tmem.gset.idx <- na.omit(match(tmem.markers.d$gene_symbol, mcold$Symbol));
tmem.core.gset.idx <- na.omit(match(tmem.markers.d$gene_symbol[tmem.markers.d$varying == 1], mcold$Symbol));

# There are no spike-ins
cat("Spike-ins:\n")
print(grep("^ERCC-", rowData(x)$ID, value=TRUE))

ctrl.features <- list(
	mito = rownames(x)[grep("^MT-", rowData(x)$Symbol)]
);
x <- addPerCellQC(x, subsets = ctrl.features);

lib.factors <- librarySizeFactors(x);

qdraw(
	with(colData(x),
		hist(total, breaks=100)
	),
	file = insert(pdf.fn, c("qc", "total"))
);

qdraw(
	hist(lib.factors, breaks=200, xlim=c(0, 2))
	,
	file = insert(pdf.fn, c("qc", "lib-factors"))
);

qdraw(
	with(colData(x),
		smoothScatter(detected, subsets_mito_percent, ylim=c(0, 20))
	),
	file = insert(pdf.fn, c("qc", "detected", "mito-pct"))
);

qdraw(
	plotMetrics(x),
	file = insert(pdf.fn, c("qc", "metrics"))
)

# only one mixture component found under all models
# (a linear, spline, and polynomial, and one_dimensional)
#qc.mod <- mixtureModel(x, model_type = "one_dimensional");
#summary(qc.mod)

# remove poor quality barcodes
x.f <- with(colData(x), x[, lib.factors > 0.2 & detected > 600 & subsets_mito_percent < 8]);


####

br <- barcodeRanks(counts(x.f));

qdraw(
	{
		with(br, plot(rank, total, log="xy", xlab="Rank", ylab="Total", pch='.'))
		idx <- order(br$rank);
		with(br, lines(rank[idx], fitted[idx], col="red"))
		abline(h=metadata(br)$knee, col="dodgerblue", lty=2)
		abline(h=metadata(br)$inflection, col="forestgreen", lty=2)
		legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
			legend=c("knee", "inflection"))
	},
	file = insert(pdf.fn, "barcode-rank-plot")
);

# Remove empty droplets
# However, empty droplets have already been removed
#empty <- emptyDrops(counts(x));

# Remove entries that are due to barcode swapping
# However, it is unclear whether the aggregated sample have barcode swapping effects
# removed.
#info.paths <- c(
#	"../count/Ctrl_5p/outs/molecule_info.h5",
#	"../count/IE1_5p/outs/molecule_info.h5"
#);
#swapped <- swappedDrops(info.paths, min.frac=0.9, get.swapped=TRUE, get.diagnostics=TRUE);

####

x.f <- logNormCounts(x.f);

dec <- modelGeneVar(x.f)
qdraw(
	{
		plot(dec$mean, dec$total, xlab="mean log-expression", ylab="variance");
		curve(metadata(dec)$trend(x), col="blue", add=TRUE);
	},
	file = insert(pdf.fn, c("gene-var"))
)

set.seed(1000);
x.f <- denoisePCA(x.f, dec, subset.row=base.gset.idx);
cl.nn <- clusterCells(x.f, use.dimred = "PCA", BLUSPARAM = bluster::NNGraphParam(k=10));
colLabels(x.f) <- factor(cl.nn);
x.f <- runTSNE(x.f, dimred="PCA");

qdraw(
	mod_rd_plot(plotTSNE(x.f, colour_by="label", text_by="label")),
	width = 6, height = 6,
	file = insert(pdf.fn, c("tsne", "clusters-nn"))
);

markers <- scoreMarkers(x.f);

markers.sel <- lapply(markers,
	function(d) {
		idx <- rownames(d) %in% pbmc.markers.d$gene_symbol;
		d.f <- d[idx, ];
		d.f <- d.f[d.f$mean.AUC > 0.5, ];
		d.f[order(d.f$mean.logFC.cohen, decreasing=TRUE), ]
	}
);

print(lapply(markers.sel, function(x) as.data.frame(x[, 1:4])))

# prune B cells and monocytes
cl.b.cells <- 16;
cl.dendritic <- c(4, 15);
cl.remove <- c(cl.b.cells, cl.dendritic);
x.p <- x.f[, ! cl.nn %in% cl.remove];

frac.expressed <- apply(logcounts(x.p), 1, function(z) mean(z > 1));
summary(frac.expressed)
features.expressed <- rownames(x.p)[frac.expressed > 0.05];
length(features.expressed)
qwrite(features.expressed, insert(rds.fn, "genes-expressed"));

set.seed(2000);
dec.p <- modelGeneVar(x.p);
x.p <- denoisePCA(x.p, dec.p, subset.row=tcell.gset.idx);
#x.p <- fixedPCA(x.p, subset.row=tcell.gset.idx, rank=10);


cl.nn <- clusterCells(x.p, use.dimred = "PCA", BLUSPARAM = bluster::NNGraphParam(cluster.fun="louvain", k=15));
#cl.nn <- clusterCells(x.p, use.dimred = "PCA", BLUSPARAM = bluster::NNGraphParam());
colLabels(x.p) <- factor(cl.nn);
cold$cluster <- cold$label;
cold$cluster_sel <- cold$label;
cl.sel <- c("5", "9", "7", "11");
cold$cluster_sel[! cold$cluster %in% cl.sel] <- NA;

x.p <- runTSNE(x.p, dimred="PCA");

table(cl.nn)

table(cl.nn) / length(cl.nn)

qdraw(
	mod_rd_plot(plotTSNE(x.p, colour_by="label", text_by="label")) +
		scale_colour_manual(values=scater:::.get_palette("tableau20"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("pruned", "tsne", "clusters-nn"))
);

qdraw(
	mod_rd_plot(plotTSNE(x.p, colour_by="subsets_mito_percent")) + guides(colour=guide_legend("% mito")),
	width = 6, height = 6,
	file = insert(pdf.fn, c("pruned", "tsne", "mito-pct"))
);

markers.p <- scoreMarkers(x.p);

markers.p.sel <- lapply(markers.p,
	function(d) {
		idx <- rownames(d) %in% tcell.gset;
		d.f <- d[idx, ];
		#d.f <- d;
		d.f <- d.f[d.f$mean.AUC > 0.5, ];
		d.f[order(d.f$mean.logFC.cohen, decreasing=TRUE), ]
	}
);

print(lapply(markers.p.sel, function(x) as.data.frame(x[1:10, 1:4])))
# cluster 1 is NKG2D+ CD8+ T cells
# cluster 4 contains FOXP3+ Treg cells
# cluster 5 and 8 are Th2 cells: CD4+ CXCR4+ CCR4+ GATA3+,
#                         also expresses of ANXA1
# cluster 9 is Th1 cells: CD4+ IRF1+ IL27RA+ TBX21+
#                         also expresses ANXA1
# cluster 7 is effector T cells expressing TNF, IFNG, GZMB, GZMH, GZMA, PRF1
# cluster 11 is CD8+ effector T cells expressing GZMB, GZMA, PRF1, TNF, IFNG
# all clusters express TCR, mostly alpha-beta

as.data.frame(markers.p.sel[[7]][1:100, 1:4])

as.data.frame(markers.p.sel[[11]][1:100, 1:4])

#genes.to.plot <- c("CCL5", "NKG7", "PRF1", "GZMB", "XCL1", "CRTAM", "VSIR", "PDCD1");
#genes.to.plot <- c("FOXP3", "CD4", "CD8A", "TRAC", "TRDC", "TNF", "IFNG", "LAMP1", "ANXA1");
#genes.to.plot <- c("CXCR4", "CCR4", "GATA3");
#genes.to.plot <- c("IRF1", "IL27RA");
#genes.to.plot <- c("IRF1", "TBX21");
genes.to.plot <- c("CCR7");

genes.to.plot <- c(
	"FOXP3", "CD4", "CD8A", "TRAC", "TRDC", "TNF", "IFNG", "LAMP1", "ANXA1",
	"CXCR4", "CCR4", "GATA3", "IRF1", "IL27RA", "IRF1", "TBX21",
	"CCL5", "NKG7", "PRF1", "GZMB", "XCL1", "CRTAM", "VSIR", "PDCD1"
);

for (gene in genes.to.plot) {
	qdraw(
		mod_rd_plot(plotTSNE(x.p, colour_by=gene)),
		width = 6, height = 6,
		file = insert(pdf.fn, c("pruned", "tsne", tolower(gene)))
	);
}

qdraw(
	plotExpression(x.p, x="label", colour_by="label", features=tmem.markers.d$gene_symbol[tmem.markers.d$category == "core_marker"]),
	width = 6, height = 5,
	file = insert(pdf.fn, c("pruned", "expr", "t-mem-markers", "core"))
);

qdraw(
	plotExpression(x.p, x="label", colour_by="label", features=c("CD3G", "CD4", "CD8A", "CD8B")),
	width = 6, height = 5,
	file = insert(pdf.fn, c("pruned", "expr", "cd4-cd8"))
);


#PD-1	PDCD1	suppression
#KLRG-1	KLRG1	senescence
#CD57	B3GAT1	senescence
#CD161	KLRB1	miscellaneous
qdraw(
	plotExpression(x.p, x="label", colour_by="label", features=c("PDCD1", "KLRG1", "B3GAT1", "KLRB1")),
	width = 6, height = 5,
	file = insert(pdf.fn, c("pruned", "expr", "t-te"))
);
# clusters 6 and 7 may be T_EMRA based on 
# high KLRG-1, high CD57, and high CD161


#GranzymeA	GZMA	cytolysis
#GranzymeB	GZMB	cytolysis
#Perforin	PRF1	cytolysis
qdraw(
	plotExpression(x.p, x="label", colour_by="label", features=c("GZMA", "GZMB", "GZMH", "PRF1")),
	width = 6, height = 5,
	file = insert(pdf.fn, c("pruned", "expr", "gzm"))
);

# best marker for T_EM are CD45RO+ CCR7-, but but we don't have CD45RO vs. CD45RA

tcr.gset <- c("TRAC", "TRBC1", "TRBC2", "TRGC1", "TRGC2", "TRDC");
qdraw(
	plotExpression(x.p, x="label", colour_by="label", features=tcr.gset),
	width = 6, height = 5,
	file = insert(pdf.fn, c("pruned", "expr", "tcr"))
);

cytokines.gset <- c("TNF", "IL2", "IFNG", "LAMP1");
qdraw(
	plotExpression(x.p, x="label", colour_by="label", features=cytokines.gset),
	width = 6, height = 5,
	file = insert(pdf.fn, c("pruned", "expr", "cytokines"))
);

# join barcode data
cold <- left_join(as.data.frame(colData(x.p)), barcode.d, by=c("Barcode"="barcode"));
cold$response[is.na(cold$response)] <- "nonresponsive";
cold$response <- factor(cold$response, levels.response);

gset_score2 <- function(sce, gset1, gset2) {
	a <- t(logcounts(sce)[gset1, , drop=FALSE]);
	b <- t(logcounts(sce)[gset2, , drop=FALSE]);
	a <- apply(a, 1, prod);
	s <- numeric(ncol(sce));
	for (j in 1:ncol(b)) {
		s <- a * b[, j];
	}
	s / ncol(b)
}

gset_score <- function(sce, gset) {
	a <- t(logcounts(sce)[gset, , drop=FALSE]);
	s <- apply(a, 1, prod);
	s^(1/ncol(a))
}

th1.gset <- c("CD4", "TBX21", "IRF1", "IL27RA");
th2.gset <- c("CD4", "GATA3", "CXCR4", "CCR4");
#th2.gset <- c("CD4", "GATA3", "CXCR4");

gset_score(x.p, th1.gset)

#cold$th1 <- gset_score2(x.p, th1.gset[1:2], th1.gset[-c(1:2)]);
#cold$th2 <- gset_score2(x.p, th2.gset[1:2], th2.gset[-c(1:2)]);

cold$th1 <- gset_score(x.p, th1.gset);
cold$th2 <- gset_score(x.p, th2.gset);

colData(x.p) <- DataFrame(cold);

# Given cateogorical variables x and y,
# estimate proportions of each level of y within each level of x
proportions <- function(x, y) {
	levels.y <- unique(y);
	ys <- split(y, x);
	d <- do.call(rbind,
		mapply(
			function(y, x) {
				counts <- table(y);
				levels.missing <- as.character(setdiff(levels.y, names(counts)));
				if (length(levels.missing) > 0) {
					counts[levels.missing] <- 0;
				}
				data.frame(
					binom.confint(counts, sum(counts), method="agresti-coull"),
					group = x,
					value = names(counts)
				)
			},
			ys,
			names(ys),
			SIMPLIFY = FALSE
		)
	);
	rownames(d) <- NULL;

	d$lower <- pmax(0, d$lower);
	d$upper <- pmin(1, d$upper);

	d
}


enrich_test <- function(d, alternative="greater") {
	c01 <- sum(d$x);
	c00 <- sum(d$n) - c01;

	ds <- split(d, factor(d$group, levels=unique(d$group)));
	hs <- lapply(ds,
		function(d) {
			c11 <- sum(d$x);
			c10 <- sum(d$n) - c11;
			m <- matrix(c(c00 - c10, c01 - c11, c10, c11), nrow=2, byrow=TRUE);
			fisher.test(m, alternative=alternative)
		}
	);

	stopifnot(d$group == names(hs))

	d$group <- factor(d$group, levels=unique(d$group));
	
	d2 <- data.frame(
		d,
		odds_ratio = unlist(lapply(hs, function(h) h$estimate)),
		odds_ratio_lower = unlist(lapply(hs, function(h) h$conf.int[1])),
		odds_ratio_upper = unlist(lapply(hs, function(h) h$conf.int[2])),
		p = unlist(lapply(hs, function(h) h$p.value))
	);
	rownames(d2) <- NULL;
	d2$q <- p.adjust(d2$p, method="BH");

	d2
}

cl.th1 <- proportions(cold$label, cold$th1 > 0);
cl.th1.h <- enrich_test(cl.th1[cl.th1$value == TRUE, ]);

cl.th2 <- proportions(cold$label, cold$th2 > 0);
cl.th2.h <- enrich_test(cl.th2[cl.th2$value == TRUE, ]);

cl.th <- rbind(
	data.frame(cl.th1.h, th_type = "Th1"),
	data.frame(cl.th2.h, th_type = "Th2")
);

qdraw(
	ggplot(cl.th, aes(x=group, y=mean, ymin=lower, ymax=upper, alpha=q)) +
		theme_classic() + theme(strip.background = element_blank()) +
		geom_col(fill=cols.response["transient"]) + 
		geom_errorbar(width=0.3, show.legend=FALSE) +
		scale_alpha_continuous(trans=revlog_trans(), breaks=c(0.01, 0.05, 0.25, 0.5), name="FDR", limits=c(1, 0.01)) +
		facet_wrap(~ th_type, ncol=1) +
		scale_y_continuous(n.breaks=3) +
		theme(legend.position="bottom") +
		xlab("cluster") + ylab("proportion")
	,
	width = 3.5, height = 3,
	file = insert(pdf.fn, c("cl", "prop", "th"))
);


qdraw(
	mod_rd_plot(plotTSNE(x.p, colour_by="sizeFactor")) + guides(colour=guide_legend("size factor")),
	width = 6, height = 6,
	file = insert(pdf.fn, c("pruned", "tsne", "size-factor"))
);

qdraw(
	mod_rd_plot(plotTSNE(x.p, colour_by="detected")),
	width = 6, height = 6,
	file = insert(pdf.fn, c("pruned", "tsne", "detected"))
);

qdraw(
	mod_rd_plot(plotTSNE(x.p[, order(colData(x.p)$response)], colour_by="response")) + 
		scale_colour_manual(values=cols.response),
	width = 6, height = 6,
	file = insert(pdf.fn, c("pruned", "tsne", "response"))
);

qdraw(
	mod_rd_plot(plotTSNE(x.p[, order(colData(x.p)$th1)], colour_by="th1")),
	width = 6, height = 6,
	file = insert(pdf.fn, c("pruned", "tsne", "th1"))
);

qdraw(
	mod_rd_plot(plotTSNE(x.p[, order(colData(x.p)$th2)], colour_by="th2")),
	width = 6, height = 6,
	file = insert(pdf.fn, c("pruned", "tsne", "th2"))
);

qdraw(
	plotColData(x.p, y="response", x="label", colour_by="response"),
	file = insert(pdf.fn, c("pruned", "cl", "response"))
);

qdraw(
	plotColData(x.p, y="th1", x="label", colour_by="label"),
	file = insert(pdf.fn, c("pruned", "cl", "th1"))
);

qdraw(
	plotColData(x.p, y="th2", x="label", colour_by="label"),
	file = insert(pdf.fn, c("pruned", "cl", "th2"))
);

cold.cl.s <- split(cold, cold$label);
cl.response.d <- do.call(rbind,
	mapply(
	function(d, label) {
		tt <- table(d$response);
		data.frame(
			binom.confint(tt, sum(tt), method="agresti-coull"),
			response = names(tt),
			label = label
		)
	},
	cold.cl.s,
	names(cold.cl.s),
	SIMPLIFY = FALSE
));
cl.response.d <- m_fix_table(cl.response.d);
cl.response.d <- cl.response.d[cl.response.d$response != "nonresponsive", ];

qdraw(
	ggplot(cl.response.d, aes(x=response, y=mean, ymin=lower, ymax=upper, fill=label)) +
		theme_classic() + theme(strip.background = element_blank()) +
		geom_col() + geom_errorbar(width=0.3) +
		guides(fill="none") +
		facet_wrap(~ label, ncol=4, scales="free_y") +
		xlab("") + ylab("proportion") +
		scale_fill_manual(values=scater:::.get_palette("tableau20")) +
		theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
	,
	width = 4.5, height = 5,
	file = insert(pdf.fn, c("cl", "prop", "response"))
);

cold.response.s <- split(cold, cold$response);
response.cl.d <- do.call(rbind,
	mapply(
	function(d, response) {
		tt <- table(d$label);
		data.frame(
			binom.confint(tt, sum(tt), method="agresti-coull"),
			label = names(tt),
			response = response
		)
	},
	cold.response.s,
	names(cold.response.s),
	SIMPLIFY = FALSE
));
response.cl.d <- m_fix_table(response.cl.d);
response.cl.d$group <- response.cl.d$label;
response.cl.d$value <- response.cl.d$response;

response.cl.h <- rbind(
	enrich_test(filter(response.cl.d, response == "transient")),
	enrich_test(filter(response.cl.d, response == "durable"))
);

qdraw(
	ggplot(response.cl.h, aes(x=label, y=mean, ymin=lower, ymax=upper, fill=response, alpha=q)) +
		theme_classic() + theme(strip.background = element_blank()) +
		geom_col() + geom_errorbar(width=0.3, show.legend=FALSE) +
		scale_fill_manual(values=cols.response) +
		guides(fill="none") +
		scale_alpha_continuous(trans=revlog_trans(), breaks=c(0.01, 0.05, 0.25, 0.5), name="FDR", limits=c(1, 0.01)) +
		facet_wrap(~ response, ncol=1) +
		theme(legend.position="bottom") +
		xlab("cluster") + ylab("proportion")
	,
	width = 3.5, height = 3,
	file = insert(pdf.fn, c("response", "prop", "cl"))
);

####

cold0 <- cold;

# convert to SingellCellAssay for use with MAST
sca <- SceToSingleCellAssay(x.p);
cold <- colData(sca);
mcold <- data.table::as.data.table(mcols(sca));


# --- T cells

#ens.cd4 <- mcold$ID[mcold$Symbol == "CD4"]
#ens.cd8a <- mcold$ID[mcold$Symbol == "CD8A"]
#ens.cd8b <- mcold$ID[mcold$Symbol == "CD8B"]
ens.cd4 <- "CD4";
ens.cd8a <- "CD8A";
ens.cd8b <- "CD8B";

qdraw(
	{
		smoothScatter(as.numeric(log(counts(sca[ens.cd8a, ]) + 1)), as.numeric(log(counts(sca[ens.cd8b, ]) + 1)))
	}
	,
	file = insert(pdf.fn, tag=c("cd8a", "cd8b"))
)

qdraw(
	{
		smoothScatter(
			as.numeric(log(counts(sca[ens.cd4, ]) + 1)),
			as.numeric(log(counts(sca[ens.cd8a, ]) + counts(sca[ens.cd8b, ]) + 1))
		)
	}
	,
	file = insert(pdf.fn, tag=c("cd4", "cd8"))
)

t.cell.markers <- c(ens.cd4, ens.cd8a, ens.cd8b);
sca.t.cells <- counts(sca[t.cell.markers, ]);
table(as.numeric(sca.t.cells))

d.t.cells <- as.data.frame(as.matrix(t(sca.t.cells > 0)));

with(d.t.cells, table(CD8A, CD8B))
with(d.t.cells, cor(CD8A, CD8B))
with(d.t.cells, fisher.test(table(CD8A, CD8B)))

with(d.t.cells, table(CD4, CD8A | CD8B))
with(d.t.cells, prop.table(table(CD4, CD8A | CD8B), 1))
with(d.t.cells, prop.table(table(CD4, CD8A | CD8B), 2))
with(d.t.cells, cor(CD4, CD8A | CD8B))
with(d.t.cells, fisher.test(table(CD4, CD8A | CD8B)))

cold$t_cell <- "Unknown";
cold$t_cell[d.t.cells$CD4] <- "CD4";
cold$t_cell[d.t.cells$CD8A | d.t.cells$CD8B] <- "CD8";
cold$t_cell0 <- cold$t_cell;

with(cold, table(t_cell))
print(with(cold, table(t_cell, response)))

expressed <- counts(sca) > 0;
#expr.cd4 <- apply(expressed[, cold$t_cell == "CD4"], 1, function(xs) mean(xs) > 0.5);
#expr.cd8 <- apply(expressed[, cold$t_cell == "CD8"], 1, function(xs) mean(xs) > 0.5);
#pid.tc.inform.idx <- which(xor(expr.cd4, expr.cd8));
expr.cd4 <- apply(expressed[, cold$t_cell == "CD4"], 1, any);
expr.cd8 <- apply(expressed[, cold$t_cell == "CD8"], 1, any);
pid.tc.inform.idx <- which(expr.cd4 | expr.cd8);
length(pid.tc.inform.idx)

known.cells <- cold$t_cell != "Unknown";
expr.sel <- expressed[pid.tc.inform.idx, known.cells];
label <- factor(cold$t_cell[known.cells]);

prop.table(table(
	expr.sel[t.cell.markers[1], ], label
), 2)

prop.table(table(
	expr.sel[t.cell.markers[2], ], label
), 2)

prop.table(table(
	expr.sel[t.cell.markers[3], ], label
), 2)

rf.fn <- insert(rds.fn, "rforest");

rforest <- qcache(
	{
		#rforest <- randomForest(as.matrix(t(expr.sel)), label);
		#expr.sel2 <- expr.sel[! rownames(expr.sel) %in% t.cell.markers, ];
		#rforest2 <- randomForest(as.matrix(t(expr.sel2)), label);

		# parallel randomForest
		library(doMC)
		library(randomForest)
		mc.cores <- getOption("mc.cores");
		registerDoMC(mc.cores);

		rforest <- foreach(
			ntree=rep(ceiling(500/mc.cores), mc.cores),
			.combine=randomForest::combine,
			.multicombine=TRUE,
			.packages="randomForest") %dopar% {
			randomForest(as.matrix(t(expr.sel)), label, ntree=ntree)
		}

		train.pred <- predict(rforest, data=t(expr.sel));
		prop.table(table(train.pred, label), 2)

		rforest
	},
	file = as.character(rf.fn, simplify=TRUE)
);

# use trained classifier to predict unknown T cells
pred <- predict(rforest, t(expressed[pid.tc.inform.idx, !known.cells]));

cold$t_cell[!known.cells] <- as.character(pred);

print(with(cold, table(t_cell0, response)))
print(with(cold, table(t_cell, response)))
print(with(cold, prop.table(table(t_cell, response), 1)))

print(with(cold, table(t_cell, response)[, -1]))
print(with(cold, prop.table(table(t_cell, response)[, -1], 1)))
print(with(cold, fisher.test(table(t_cell, response)[, -1])))
# odds ratio = 50.07384
# p-value = 4.937e-15

response.p <- with(cold, fisher.test(table(t_cell, response != "nonresponsive"))$p.value);

cold.s <- split(cold, cold$t_cell);
response.d <- do.call(rbind, mapply(
	function(d, cell) {
		tt <- table(d$response != "nonresponsive");
		names(tt) <- c("nonresponsive", "responsive");
		data.frame(
			binom.confint(tt, sum(tt), method="agresti-coull"),
			response = names(tt),
			t_cell = cell
		)
	},
	cold.s,
	names(cold.s),
	SIMPLIFY = FALSE
));

qdraw(
	ggplot(filter(response.d, response == "responsive"), aes(x=t_cell, y=mean * 100, ymin=lower * 100, ymax=upper * 100, colour=t_cell)) +
		theme_classic() +
		geom_point() +
		geom_errorbar(width=0.3) +
		scale_y_log10() +
		scale_colour_npg() +
		xlab("") + ylab("% responsive cells") + 
		guides(colour="none") +
		ggtitle(sprintf("p = %s", format(response.p, digits=2)))
	,
	width = 2, height = 3,
	file = insert(pdf.fn, c("prop", "t-cell", "response"))
)


responset.p <- with(cold, fisher.test(table(t_cell, response)[, -1])$p.value);

responset.d <- do.call(rbind, mapply(
	function(d, cell) {
		tt <- table(d$response)[-1];
		data.frame(
			binom.confint(tt, sum(tt), method="agresti-coull"),
			response = names(tt),
			t_cell = cell
		)
	},
	cold.s,
	names(cold.s),
	SIMPLIFY = FALSE
));

qdraw(
	ggplot(responset.d, aes(x=response, y=mean, ymin=lower, ymax=upper, fill=t_cell)) +
		theme_classic() +
		geom_col(position=position_dodge()) + 
		geom_errorbar(position=position_dodge(width=0.9), aes(group=t_cell), width=0.3) +
		scale_fill_npg(name="") +
		xlab("") + ylab("proportion") + 
		ggtitle(sprintf("p = %s", format(responset.p, digits=2)))
	,
	width = 3, height = 3,
	file = insert(pdf.fn, c("prop", "t-cell", "response-type"))
)

# ---

colData(sca) <- DataFrame(cold);
colData(x.p) <- DataFrame(cold);


# --- Analyze the responsives clones

with(cold, table(response, t_cell))

with(cold[cold$cluster %in% c("5", "9", "7", "11"), ], table(response, t_cell, as.character(cluster)))


x.responsive <- x.p[, cold$response %in% c("responsive", "transient")];

set.seed(1000);
dec.responsive <- modelGeneVar(x.responsive);
x.responsive <- denoisePCA(x.responsive, dec.responsive, subset.row=tcell.gset.idx);
cl.responsive.nn <- clusterCells(x.responsive, use.dimred = "PCA", BLUSPARAM = bluster::NNGraphParam(cluster.fun="louvain", k=5));
colLabels(x.responsive) <- factor(cl.responsive.nn);
x.responsive <- runTSNE(x.responsive, dimred="PCA");

qdraw(
	mod_rd_plot(plotTSNE(x.responsive, colour_by="response")) +
		scale_colour_manual(values=cols.response),
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("responsive", "tsne", "response"))
);

qdraw(
	mod_rd_plot(plotTSNE(x.responsive, colour_by="cluster", text_by="cluster")),
	width = 6, height = 6,
	file = insert(pdf.fn, c("responsive", "tsne", "clusters-nn1"))
);

qdraw(
	mod_rd_plot(plotTSNE(x.responsive, colour_by="t_cell")),
	width = 6, height = 6,
	file = insert(pdf.fn, c("responsive", "tsne", "cd4-cd8"))
);

qdraw(
	mod_rd_plot(plotTSNE(x.responsive, colour_by="label", text_by="label")),
	width = 6, height = 6,
	file = insert(pdf.fn, c("responsive", "tsne", "clusters-nn2"))
);


x.durable <- x.p[, cold$response == "durable"];

set.seed(1000);
dec.durable <- modelGeneVar(x.durable);
x.durable <- denoisePCA(x.durable, dec.durable, subset.row=tcell.gset.idx);
cl.durable.nn <- clusterCells(x.durable, use.dimred = "PCA", BLUSPARAM = bluster::NNGraphParam(cluster.fun="louvain", k=5));
colLabels(x.durable) <- factor(cl.durable.nn);
x.durable <- runTSNE(x.durable, dimred="PCA");

qdraw(
	mod_rd_plot(plotTSNE(x.durable, colour_by="label", text_by="label")),
	width = 6, height = 6,
	file = insert(pdf.fn, c("durable", "tsne", "clusters-nn"))
);


x.transient <- x.p[, cold$response == "transient"];

set.seed(1000);
dec.transient <- modelGeneVar(x.transient);
x.transient <- denoisePCA(x.transient, dec.transient, subset.row=tcell.gset.idx);
cl.transient.nn <- clusterCells(x.transient, use.dimred = "PCA", BLUSPARAM = bluster::NNGraphParam(cluster.fun="louvain", k=5));
colLabels(x.transient) <- factor(cl.transient.nn);
x.transient <- runTSNE(x.transient, dimred="PCA", perplexity=5);

qdraw(
	mod_rd_plot(plotTSNE(x.transient, colour_by="label", text_by="label")),
	width = 6, height = 6,
	file = insert(pdf.fn, c("transient", "tsne", "clusters-nn"))
);


x.rsel <- x.p[, cold$cluster %in% cl.sel & cold$response != "nonresponsive"];
dec.rsel <- modelGeneVar(x.rsel);
x.rsel <- denoisePCA(x.rsel, dec.rsel, subset.row=tcell.gset.idx);
x.rsel <- runUMAP(x.rsel, dimred="PCA");

qdraw(
	mod_rd_plot(plotUMAP(x.rsel, colour_by="response", shape_by="t_cell")) +
		scale_colour_manual(values=cols.response)
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("responsive-sel", "umap", "response"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.rsel, colour_by="cluster", shape_by="t_cell")) +
		scale_colour_manual(values=scater:::.get_palette("tableau20")[as.integer(cl.sel)])
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("responsive-sel", "umap", "cluster"))
);


set.seed(1000);
x.sel <- x.p[, colData(x.p)$cluster %in% cl.sel];
dec.sel <- modelGeneVar(x.sel);
x.sel <- denoisePCA(x.sel, dec.rsel, subset.row=tcell.gset.idx);
#x.sel <- fixedPCA(x.sel, rank=50, subset.row=tcell.gset.idx);
x.sel <- runDiffusionMap(x.sel, dimred="PCA");
x.sel <- runUMAP(x.sel, dimred="PCA");
cl2.n <- clusterCells(x.sel, use.dimred="UMAP",
	BLUSPARAM = bluster::KmeansParam(centers=3));
table(cl2.n)
#colLabels(x.sel) <- factor(cl2.n);
cl2 <- factor(cl2.n, levels=c(2, 1, 3), labels=LETTERS[1:3]);
colLabels(x.sel) <- cl2;

cold$cluster2 <- factor(NA, levels=levels(cl2));
cold$cluster2[match(names(cl2), rownames(cold))] <- cl2;
table(cold$cluster2)

x.sel.reordered <- x.sel[, order(colData(x.sel)$response)];

qdraw(
	mod_rd_plot(plotDiffusionMap(x.sel.reordered, colour_by="response", shape_by="t_cell")) +
		scale_colour_manual(values=cols.response)
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("selected", "diffusion", "response"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.reordered,
		colour_by="response", shape_by="t_cell")) +
		scale_colour_manual(values=cols.response)
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("selected", "umap", "response"))
);

cols.cl <- scater:::.get_palette("tableau20");
cols.cl.sel <- cols.cl[sort(as.integer(cl.sel))];

rgb2col <- function(x) {
	rgb(x[1, ], x[2, ], x[3, ], maxColorValue=255)
}

cols.cl2 <- c(
	A = rgb2col(col2rgb(cols.cl.sel[1])*0.6 + col2rgb(cols.cl.sel[3])*0.4),
	B = cols.cl.sel[2],
	C = cols.cl.sel[4]
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel, colour_by="cluster", shape_by="t_cell")) +
		scale_colour_manual(values=cols.cl.sel)
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("selected", "umap", "cluster"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel, colour_by="label", shape_by="t_cell"))
		#scale_colour_manual(values=cols.cl2)
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("selected", "umap", "label"))
);

qdraw(
	mod_rd_plot(plotPCA(x.sel, colour_by="label", shape_by="t_cell")) +
		scale_colour_manual(values=cols.cl2)
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("selected", "pca", "label"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel, colour_by="sizeFactor", shape_by="t_cell"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("selected", "umap", "size-factor"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel, colour_by="t_cell0", shape_by="t_cell")) +
		scale_colour_npg()
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("selected", "umap", "t-cell0"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel, colour_by="t_cell", shape_by="t_cell")) +
		scale_colour_npg()
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("selected", "umap", "t-cell"))
);


set.seed(1000);
x.sel.cd8 <- x.p[, colData(x.p)$cluster %in% cl.sel & colData(x.p)$t_cell == "CD8"];
#dec.sel.cd8 <- modelGeneVar(x.sel.cd8);
#x.sel.cd8 <- denoisePCA(x.sel.cd8, dec.sel.cd8, subset.row=tmem.core.gset.idx);
x.sel.cd8 <- fixedPCA(x.sel.cd8, rank=20, subset.row=tmem.core.gset.idx);
#x.sel.cd8 <- runDiffusionMap(x.sel.cd8, dimred="PCA");
x.sel.cd8 <- runUMAP(x.sel.cd8, dimred="PCA", n_neighbors=30);
x.sel.cd8 <- runTSNE(x.sel.cd8, dimred="PCA");
cl.cd8.n <- clusterCells(x.sel.cd8, use.dimred="UMAP",
	BLUSPARAM = bluster::KmeansParam(centers=3));
table(cl.cd8.n)
colLabels(x.sel.cd8) <- factor(cl.cd8.n);
#cl2 <- factor(cl2.n, levels=c(2, 1, 3), labels=LETTERS[1:3]);
#colLabels(x.sel) <- cl2;

#qdraw(
#	mod_rd_plot(plotDiffusionMap(x.sel.reordered, colour_by="response", shape_by="t_cell")) +
#		scale_colour_manual(values=cols.response)
#	,
#	width = 6, height = 6,
#	file = insert(pdf.fn, c("cd8-sel", "diffusion", "response"))
#);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="cluster")) +
		scale_colour_manual(values=cols.cl.sel)
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "cluster"))
);

qdraw(
	mod_rd_plot(plotUMAP(
		x.sel.cd8[, order(colData(x.sel.cd8)$t_cell)], colour_by="response")) +
		scale_colour_manual(values=cols.response)
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "response"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="CCR7"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "ccr7"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="CD28"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "cd28"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="FAS"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "fas-cd95"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="CD27"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "cd27"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="IL7R"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "il17r-cd127"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="IL2RB"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "il2rb-cd122"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="KLRG1"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "klrg1"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="KLRB1"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "klrb1-cd161"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="SELL"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "sell"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="CXCR4"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "cxcr4"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="GZMB"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "gzmb"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="PRF1"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "prf1"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="IL6ST"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "il6st-cd130"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="KIT"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "kit"))
);

qdraw(
	mod_rd_plot(plotUMAP(x.sel.cd8, colour_by="CCR4"))
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "umap", "ccr4"))
);



qdraw(
	mod_rd_plot(plotTSNE(
		x.sel.cd8[, order(colData(x.sel.cd8)$t_cell)], colour_by="cluster")) +
		scale_colour_manual(values=cols.cl.sel)
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "tsne", "cluster"))
);

qdraw(
	mod_rd_plot(plotTSNE(
		x.sel.cd8[, order(colData(x.sel.cd8)$t_cell)], colour_by="response")) +
		scale_colour_manual(values=cols.response)
	,
	width = 6, height = 6,
	file = insert(pdf.fn, c("cd8-sel", "tsne", "response"))
);


tmem.clf <- qread("../data/etabm40_pamr.rds");

x.sel.cd8.tmem.cl <- tmem.clf$predict(tmem.clf, logcounts(x.sel.cd8));

x.sel.cd8.tmem.post <- tmem.clf$predict(tmem.clf, logcounts(x.sel.cd8), type="posterior");
x.sel.cd8.tmem.post

table(x.sel.cd8.tmem.cl)
table(colData(x.sel.cd8)$response, x.sel.cd8.tmem.cl)

prop.table(table(colData(x.sel.cd8)$response, x.sel.cd8.tmem.cl), margin=1)


x.p.cd8 <- x.p[, colData(x.p)$t_cell == "CD8"];

x.p.cd8.tmem.cl <- tmem.clf$predict(tmem.clf, logcounts(x.p.cd8));

table(x.p.cd8.tmem.cl)
table(colData(x.p.cd8)$response, x.p.cd8.tmem.cl)

prop.table(table(colData(x.p.cd8)$response, x.p.cd8.tmem.cl), margin=1)


# --- End T cells

colData(sca) <- DataFrame(cold);
colData(x.p) <- DataFrame(cold);

qwrite(cold, insert(rds.fn, "cells"));
qwrite(mcold, insert(rds.fn, "features"));

qwrite(x.p, insert(rds.fn, "sce"));

sca.cd4 <- sca[, cold$t_cell == "CD4"];
sca.cd8 <- sca[, cold$t_cell == "CD8"];

qwrite(sca, insert(rds.fn, "sca"));
qwrite(sca.cd4, insert(rds.fn, c("sca", "cd4")));
qwrite(sca.cd8, insert(rds.fn, c("sca", "cd8")));

