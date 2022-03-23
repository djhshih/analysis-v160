library(io)
library(MAST)
library(ggplot2)
library(ggrepel)
library(dplyr)

source("../R/common.R")

options(mc.cores=50);

set.seed(1337);

csv.fdr.cut <- 0.25;
fdr.cut <- 0.05;

gene_logfc <- function(sca, groups, genes=NULL) {
	group.levels <- levels(groups);
	names(group.levels) <- group.levels;
	means <- lapply(group.levels, function(g) {
		idx <- which(groups == g);
		if (is.null(genes)) {
			rowMeans(logcounts(sca)[, idx, drop=FALSE])
		} else {
			rowMeans(logcounts(sca)[genes, idx, drop=FALSE])
		}
	});
	means[[2]] - means[[1]]
}

wd_summary_table <- function(wd, logfc, genes=NULL) {
	d <- data.frame(
		gene = rownames(wd),
		logfc = logfc,
		discrete_wald = wd[, "disc", "lambda"],
		continuous_wald = wd[, "cont", "lambda"],
		hurdle_wald = wd[, "hurdle", "lambda"],
		discrete_p = wd[, "disc", "Pr(>Chisq)"],
		continuous_p = wd[, "cont", "Pr(>Chisq)"],
		hurdle_p = wd[, "hurdle", "Pr(>Chisq)"]
	);

	# convert Wald statistic to signed Wald statistic
	wald_to_swald <- function(w, logfc) {
		sign(logfc) * sqrt(w)
	}

	d$discrete_swald <- wald_to_swald(d$discrete_wald, logfc);
	d$continuous_swald <- wald_to_swald(d$continuous_wald, logfc);
	d$hurdle_swald <- wald_to_swald(d$hurdle_wald, logfc);

	if (!is.null(genes)) {
		d <- d[match(genes, d$gene), ];
	}

	d$discrete_q <- p.adjust(d$discrete_p, "BH");
	d$continuous_q <- p.adjust(d$continuous_p, "BH");
	d$hurdle_q <- p.adjust(d$hurdle_p, "BH");

	d[order(d$hurdle_p, -abs(d$logfc)), ]
}

volcano_limits <- function() {
	list(
		xlim(-3, 3),
		scale_y_continuous(trans=revlog_trans(), name="FDR", limits=c(1, 1e-15), breaks=c(1, 0.05, 1e-3, 1e-6, 1e-9, 1e-12, 1e-15))
	)
}

wd_results_plot <- function(d, fdr.cut=0.05, n.top=30, n.ns=500) {
	d$group <- ifelse(
		d$hurdle_q < fdr.cut,
		ifelse(d$logfc > 0, "Up", "Down"),
		"NS"
	);

	d <- d[order(d$hurdle_q), ];
	d$label <- NA;
	d$label[1:n.top] <- d$gene[1:n.top];
	d$label[d$group == "NS"] <- NA;

	# downsample NS
	sig.idx <- which(d$group != "NS");
	if (sum(d$group == "NS", na.rm=TRUE) > n.ns) {
		ns.idx <- which(d$group == "NS");
		prob <- -log(pmin(d$hurdle_q[ns.idx], 0.999));
		valid <- !is.na(prob) & prob > 0;

		idx.valid <- ns.idx[valid];
		n.valid <- length(idx.valid);
		if (n.valid >= n.ns) {
			idx <- union(sig.idx, sample(idx.valid, n.ns, prob=prob[valid]));
		} else {
			idx <- union(sig.idx, which(idx.valid));
		}
	}

	ggplot(d[idx, ],
		aes(
			x=discrete_q, y=continuous_q, colour=group
		)
	) +
		theme_classic() +
		geom_point(alpha=0.3) +
		scale_colour_manual(values=c(NS="grey30", Up="firebrick3", Down="royalblue4")) +
		geom_text_repel(aes(label=label), colour="black", show.legend=FALSE) +
		scale_x_continuous(trans=revlog_trans(), name="FDR of discrete component") +
		scale_y_continuous(trans=revlog_trans(), name="FDR of continuous component")
}

volcano_plot <- function(d, fdr.cut=0.05, n.top=30, n.ns=500) {
	d$group <- ifelse(
		d$q < fdr.cut,
		ifelse(d$logfc > 0, "Up", "Down"),
		"NS"
	);

	d <- d[order(d$hurdle_q), ];
	d$label <- NA;
	d$label[1:n.top] <- d$gene[1:n.top];
	d$label[d$group == "NS"] <- NA;
	
	# downsample NS
	sig.idx <- which(d$group != "NS");
	if (sum(d$group == "NS", na.rm=TRUE) > n.ns) {
		ns.idx <- which(d$group == "NS");
		prob <- -log(pmin(d$hurdle_q[ns.idx], 0.999));
		valid <- !is.na(prob) & prob > 0;

		idx.valid <- ns.idx[valid];
		n.valid <- length(idx.valid);
		if (n.valid >= n.ns) {
			idx <- union(sig.idx, sample(idx.valid, n.ns, prob=prob[valid]));
		} else {
			idx <- union(sig.idx, which(idx.valid));
		}
	}

	ggplot(d[idx, ],
		aes(
			x=logfc, y=q, colour=group
		)
	) +
		theme_classic() +
		geom_vline(xintercept=0, colour="grey60") +
		geom_hline(yintercept=0.05, colour="grey60") +
		geom_point(, alpha=0.3) +
		guides(colour="none") +
		scale_colour_manual(values=c(NS="grey30", Up="firebrick3", Down="royalblue4")) +
		geom_text_repel(aes(label=label), colour="black", show.legend=FALSE) +
		scale_alpha_continuous(trans=revlog_trans(), name="FDR") +
		scale_y_continuous(trans=revlog_trans(), name="FDR") +
		scale_x_continuous(name="log FC")
}

sanitize <- function(x) {
	gsub(".", "-", x, fixed=TRUE)
}



in.fn <- as.filename("v160_sca.rds");
mcold <- qread("v160_features.rds");

out.fn <- filename(in.fn$fstem, tag=setdiff(in.fn$tag, "sca"));
rds.fn <- insert(out.fn, ext="rds");
csv.fn <- insert(out.fn, ext="csv");
pdf.fn <- insert(out.fn, ext="pdf");

sca <- qread(in.fn);

# ignore TCR variable genes
genes.tcrv <- grep("^TR(A|B|G)V[0-9]*", rownames(sca), value=TRUE);
genes.ribo <- grep("^RP((LP)|L|S)[0-9]+[A-Z]?", rownames(sca), value=TRUE);
genes.ignore <- c(genes.tcrv, genes.ribo);
genes <- setdiff(rownames(sca), genes.ignore);


zfit <- zlm(~ response, sca);
qwrite(zfit, insert(rds.fn, "mast-zlm"));

# segfaults during calculation of logFC
#sm <- summary(zfit);

# segfaults during LRT
#sm.lr <- summary(zfit, logFC=FALSE, doLRT="responsetransient");

# extract coefficients
# components
# D: discrete, C: continuous, S: summary/combined
coefs <- summary(zfit, logFC=FALSE)$datatable;
coefs <- coefs[contrast != "(Intercept)", ];

# run log-likelihood ratio test
# components
# D: discrete, C: continuous, H: hurdle
lr <- lrTest(zfit, "response");
qwrite(lr, insert(rds.fn, "mast-zlm-lr"));

# ----

wd.durable <- waldTest(zfit, Hypothesis("responsedurable", c("(Intercept)", "responsedurable")));
wd.transient <- waldTest(zfit, Hypothesis("responsetransient", c("(Intercept)", "responsetransient")));
wd.durable.vs.transient <- waldTest(zfit, Hypothesis("responsedurable-responsetransient", c("responsetransient", "responsedurable")));
wd <- list(durable = wd.durable, transient = wd.transient, durable.vs.transient = wd.durable.vs.transient);
qwrite(wd, insert(rds.fn, "mast-zlm-wd"));

# convert from platform primer IDs to genes
pid_to_genes <- function(pids) {
	unique(mcold$Symbol[match(pids, mcold$ID)])
}

# res waldTest or lrTest results
extract_sig_genes <- function(res, component="hurdle", fdr.cut=0.05, convert.symbol=FALSE) {
	ps <- sort(res[, component, "Pr(>Chisq)"]);
	qs <- p.adjust(ps, "BH");

	idx <- qs < fdr.cut;
	sigs <- names(qs)[idx]
	if (convert.symbol) {
		pid_to_genes(sigs)
	} else {
		sigs
	}
}

compare_sets <- function(xs) {
	universe <- unique(unlist(xs));
	idx <- lapply(
		xs,
		function(x) {
			universe %in% x
		}
	);
	table(idx)
}

sgenes.wd <- lapply(wd, extract_sig_genes);
compare_sets(sgenes.wd)

sgenes.wd.disc <- lapply(wd, function(x) extract_sig_genes(x, component="disc"));
compare_sets(sgenes.wd.disc)


response.contrasts <- list(
	durable = c("nonresponsive", "durable"),
	transient = c("nonresponsive", "transient"),
	durable.vs.transient = c("transient", "durable")
);

logfcs <- mapply(
	function(cont) {
		gene_logfc(
			sca,
			factor(colData(sca)$response, cont)
		)
	},
	response.contrasts,
	SIMPLIFY = FALSE
);

wd.ds <- mapply(
	function(d, logfc) wd_summary_table(d, logfc, genes=genes),
	wd,
	logfcs,	
	SIMPLIFY = FALSE
);

components <- c("hurdle", "discrete", "continuous");

for (stratum in names(wd.ds)) {
	for (component in components) {
		d <- wd.ds[[stratum]];
		d$q <- d[[sprintf("%s_q", component)]];

		qdraw(
			volcano_plot(d)
			,
			width = 5,
			file = insert(pdf.fn, c("wd", sanitize(stratum), "volcano", component))
		);
	}
}

for (stratum in names(wd.ds)) {
	qdraw(
		wd_results_plot(wd.ds[[stratum]])
		,
		width = 6,
		file = insert(pdf.fn, c("wd", sanitize(stratum), "disc-cont-fdr"))
	);
}

for (stratum in names(wd.ds)) {
	qwrite(wd.ds[[stratum]], 
		file = insert(rds.fn, c("wd", sanitize(stratum)))
	)
	qwrite(filter(wd.ds[[stratum]], hurdle_q < csv.fdr.cut), 
		file = insert(csv.fn, c("wd", sanitize(stratum)))
	)
}

spids.lr <- extract_sig_genes(lr);
sgenes.lr <- pid_to_genes(spids.lr);
qwrite(sgenes.lr, insert(out.fn, c("lr", "sig-genes"), ext="vtr"));

# ----

with(colData(sca), table(response, t_cell))

with(colData(sca), table(cluster2, response, t_cell))


# durable C vs. nonresponsive C: nothing

sca.c <- sca[, which(colData(sca)$cluster2 == "C")];
zfit.c <- zlm(~ response, sca.c);
wd.c.durable <- waldTest(zfit.c, Hypothesis("responsedurable", c("(Intercept)", "responsedurable")));
logfc.c.durable <- gene_logfc(sca.c, factor(colData(sca.c)$response, c("nonresponsive", "durable")));
d.c.durable <- wd_summary_table(wd.c.durable, logfc.c.durable, genes);

qdraw(
	volcano_plot(mutate(d.c.durable, q=hurdle_q)) + volcano_limits(),
	width = 5,
	file = insert(pdf.fn, c("wd", "cluster2-c", "durable", "volcano", "hurdle"))
);

qwrite(
	filter(d.c.durable, hurdle_q < csv.fdr.cut),
	file = insert(csv.fn, c("wd", "cluster2-c", "durable"))
)

qwrite(
	d.c.durable,
	file = insert(rds.fn, c("wd", "cluster2-c", "durable"))
)


# durable B vs. nonresponsive B: nothing

sca.b <- sca[, which(colData(sca)$cluster2 == "B")];
zfit.b <- zlm(~ response, sca.b);
wd.b.durable <- waldTest(zfit.b, Hypothesis("responsedurable", c("(Intercept)", "responsedurable")));
logfc.b.durable <- gene_logfc(sca.b, factor(colData(sca.b)$response, c("nonresponsive", "durable")));
d.b.durable <- wd_summary_table(wd.b.durable, logfc.b.durable, genes);

qdraw(
	volcano_plot(mutate(d.b.durable, q=hurdle_q)) + volcano_limits(),
	width = 5,
	file = insert(pdf.fn, c("wd", "cluster2-b", "durable", "volcano", "hurdle"))
);

qwrite(
	filter(d.b.durable, hurdle_q < csv.fdr.cut),
	file = insert(csv.fn, c("wd", "cluster2-b", "durable"))
);

qwrite(
	d.b.durable,
	file = insert(rds.fn, c("wd", "cluster2-b", "durable"))
);


# transient A vs. nonresponsive A: nothing

sca.a <- sca[, which(colData(sca)$cluster2 == "A")];
zfit.a <- zlm(~ response, sca.a);
wd.a.transient <- waldTest(zfit.a, Hypothesis("responsetransient", c("(Intercept)", "responsetransient")));
logfc.a.transient <- gene_logfc(sca.a, factor(colData(sca.a)$response, c("nonresponsive", "transient")));
d.a.transient <- wd_summary_table(wd.a.transient, logfc.a.transient, genes);

qdraw(
	volcano_plot(mutate(d.a.transient, q=hurdle_q)) + volcano_limits(),
	width = 5,
	file = insert(pdf.fn, c("wd", "cluster2-a", "transient", "volcano", "hurdle"))
);

qwrite(
	filter(d.a.transient, hurdle_q < csv.fdr.cut),
	file = insert(csv.fn, c("wd", "cluster2-a", "transient"))
);

qwrite(
	d.a.transient,
	file = insert(rds.fn, c("wd", "cluster2-a", "transient"))
);



# CD8+ durable C vs. CD8+ durable B: many differentially expressed genes
# GZMB, XCL1, XCL2, CCL4, CCL4L2
# PGAM1, PKM, GAPDH, 
# HSP90AB1, HSPA8
# GNAS, BCL2A1, EIF4A1

sca.cd8.durable.bc <- sca[,
	which(
		colData(sca)$t_cell == "CD8" &
		colData(sca)$response == "durable" & 
		colData(sca)$cluster2 %in% c("B", "C")
	)
];
with(colData(sca.cd8.durable.bc), table(response, cluster2))
colData(sca.cd8.durable.bc)$cluster2 <- factor(colData(sca.cd8.durable.bc)$cluster2, levels=c("B", "C"));
zfit.durable.bc <- zlm(~ cluster2, sca.cd8.durable.bc);
qwrite(
	zfit.durable.bc,
	file = insert(rds.fn, c("cd8", "durable", "mast-zlm", "c-vs-b"))
)

wd.durable.bc <- waldTest(zfit.durable.bc, Hypothesis("cluster2C", c("(Intercept)", "cluster2C")));
logfc.durable.bc <- gene_logfc(sca.cd8.durable.bc, factor(colData(sca.cd8.durable.bc)$cluster2, c("B", "C")));
d.durable.bc <- wd_summary_table(wd.durable.bc, logfc.durable.bc, genes);

qdraw(
	volcano_plot(mutate(d.durable.bc, q=hurdle_q)),
	width = 5,
	file = insert(pdf.fn, c("wd", "cd8", "durable", "c-vs-b", "volcano", "hurdle"))
);

qdraw(
	wd_results_plot(d.durable.bc)
	,
	width = 6,
	file = insert(pdf.fn, c("wd", "cd8", "durable", "c-vs-b", "disc-cont-fdr"))
);

qwrite(
	filter(d.durable.bc, hurdle_q < csv.fdr.cut),
	file = insert(csv.fn, c("wd", "cd8", "durable", "c-vs-b"))
);

qwrite(
	d.durable.bc,
	file = insert(rds.fn, c("wd", "cd8", "durable", "c-vs-b"))
)


# CD8+ durable vs. CD8+ nonresponsive: meaningful
# PGAM1, GAPDH, PKM
# GZMB
# IFNG, XCL1, XCL2, CCL4L2, CCL4

sca.cd8 <- sca[, which(colData(sca)$t_cell == "CD8")];
zfit.cd8 <- zlm(~ response, sca.cd8);
summary(coef(zfit.cd8, "C"))
summary(coef(zfit.cd8, "D"))
qwrite(
	zfit.cd8,
	file = insert(rds.fn, c("cd8", "mast-zlm"))
)

wd.cd8.durable <- waldTest(zfit.cd8, Hypothesis("responsedurable", c("(Intercept)", "responsedurable")));
logfc.cd8.durable <- gene_logfc(sca.cd8, factor(colData(sca.cd8)$response, c("nonresponsive", "durable")));
d.cd8.durable <- wd_summary_table(wd.cd8.durable, logfc.cd8.durable, genes);

qdraw(
	volcano_plot(mutate(d.cd8.durable, q=hurdle_q)),
	width = 5,
	file = insert(pdf.fn, c("wd", "cd8", "durable", "volcano", "hurdle"))
);

qdraw(
	wd_results_plot(d.cd8.durable)
	,
	width = 6,
	file = insert(pdf.fn, c("wd", "cd8", "durable", "disc-cont-fdr"))
);

qwrite(
	filter(d.cd8.durable, hurdle_q < csv.fdr.cut),
	file = insert(csv.fn, c("wd", "cd8", "durable"))
);

qwrite(
	d.cd8.durable,
	file = insert(rds.fn, c("wd", "cd8", "durable"))
);


# Compare durable cells of cluster B or C to other CD8 T cells

with(colData(sca.cd8), table(cluster, cluster2))
# cluster2 B == cluster 7
# cluster2 C == cluster 11
with(colData(sca.cd8), table(response, cluster))

with(colData(sca), table(response, cluster2, useNA="always"))
with(colData(sca.cd8), table(response, cluster2, useNA="always"))

colData(sca.cd8)$durableB <- 0;
idx <- colData(sca.cd8)$response == "durable" & colData(sca.cd8)$cluster2 == "B";
colData(sca.cd8)$durableB[idx] <- 1;
with(colData(sca.cd8), table(durableB))

colData(sca.cd8)$durableC <- 0;
idx <- colData(sca.cd8)$response == "durable" & colData(sca.cd8)$cluster2 == "C";
colData(sca.cd8)$durableC[idx] <- 1;
with(colData(sca.cd8), table(durableC))

# NB We could not fit a variable that partitions the cells into (B, C, other),
#    likely because the durableB and the durableC are colinear

zfit.cd8.durableb <- zlm(~ durableB, sca.cd8);
summary(coef(zfit.cd8.durableb, "C"))
summary(coef(zfit.cd8.durableb, "D"))
qwrite(
	zfit.cd8.durableb,
	file = insert(rds.fn, c("cd8", "mast-zlm", "durable-b"))
)

zfit.cd8.durablec <- zlm(~ durableC, sca.cd8);
summary(coef(zfit.cd8.durablec, "C"))
summary(coef(zfit.cd8.durablec, "D"))
qwrite(
	zfit.cd8.durablec,
	file = insert(rds.fn, c("cd8", "mast-zlm", "durable-c"))
)


# CD4+ transient vs. CD4+ nonresponsive: nothing interesting
# IL4I1, SRPK2, TNFRSF18

sca.cd4 <- sca[, which(colData(sca)$t_cell == "CD4")];
zfit.cd4 <- zlm(~ response, sca.cd4);
qwrite(
	zfit.cd4,
	file = insert(rds.fn, c("cd4", "mast-zlm"))
)

wd.cd4.transient <- waldTest(zfit.cd4, Hypothesis("responsetransient", c("(Intercept)", "responsetransient")));
logfc.cd4.transient <- gene_logfc(sca.cd4, factor(colData(sca.cd4)$response, c("nonresponsive", "transient")));
d.cd4.transient <- wd_summary_table(wd.cd4.transient, logfc.cd4.transient, genes);

options(ggrepel.max.overlaps=20)
qdraw(
	volcano_plot(mutate(d.cd4.transient, q=hurdle_q)) + 
		scale_y_continuous(trans=revlog_trans(), name="FDR", limits=c(1, 1e-3), breaks=c(1, 0.05, 0.01, 1e-3))
	,
	width = 5,
	file = insert(pdf.fn, c("wd", "cd4", "transient", "volcano", "hurdle"))
);

qdraw(
	wd_results_plot(d.cd4.transient)
	,
	width = 6,
	file = insert(pdf.fn, c("wd", "cd4", "transient", "disc-cont-fdr"))
);

qwrite(
	filter(d.cd4.transient, hurdle_q < csv.fdr.cut),
	file = insert(csv.fn, c("wd", "cd4", "transient"))
);

qwrite(
	d.cd4.transient,
	file = insert(rds.fn, c("wd", "cd4", "transient"))
);


# CD8+ durable vs. CD8+ transient: meaningful
# PGAM, PKM, GAPDH
# NKG7, KLRD1
# CCL4, CCL5
# RACK1-

# already previously run
#sca.cd8 <- sca[, which(colData(sca)$t_cell == "CD8")];
#zfit.cd8 <- zlm(~ response, sca.cd8);

wd.cd8.durable.t <- waldTest(zfit.cd8, Hypothesis("responsedurable-responsetransient", c("responsetransient", "responsedurable")));
logfc.cd8.durable.t <- gene_logfc(sca.cd8, factor(colData(sca.cd8)$response, c("transient", "durable")));
d.cd8.durable.t <- wd_summary_table(wd.cd8.durable.t, logfc.cd8.durable.t, genes);

qdraw(
	volcano_plot(mutate(d.cd8.durable.t, q=hurdle_q)),
	width = 5,
	file = insert(pdf.fn, c("wd", "cd8", "durable-vs-transient", "volcano", "hurdle"))
);

qdraw(
	wd_results_plot(d.cd8.durable.t)
	,
	width = 6,
	file = insert(pdf.fn, c("wd", "cd8", "durable-vs-transient", "disc-cont-fdr"))
);

qwrite(
	filter(d.cd8.durable.t, hurdle_q < csv.fdr.cut),
	file = insert(csv.fn, c("wd", "cd8", "durable-vs-transient"))
);

qwrite(
	d.cd8.durable.t,
	file = insert(rds.fn, c("wd", "cd8", "durable-vs-transient"))
);


# CD8+ durable non-C vs. CD8+ transient non-C: nothing

sca.cd8.nonc <- sca[, which(colData(sca)$t_cell == "CD8" & colData(sca)$cluster2 != "C")];
zfit.cd8.nonc <- zlm(~ response, sca.cd8.nonc);

wd.cd8.nonc.durable <- waldTest(zfit.cd8.nonc, Hypothesis("responsedurable-responsetransient", c("responsetransient", "responsedurable")));
logfc.cd8.nonc.durable <- gene_logfc(sca.cd8.nonc, factor(colData(sca.cd8.nonc)$response, c("transient", "durable")));
d.cd8.nonc.durable <- wd_summary_table(wd.cd8.nonc.durable, logfc.cd8.nonc.durable, genes);

qdraw(
	volcano_plot(mutate(d.cd8.nonc.durable, q=hurdle_q)) + volcano_limits(),
	width = 5,
	file = insert(pdf.fn, c("wd", "cd8", "non-c", "durable-vs-transient", "volcano", "hurdle"))
);

qdraw(
	wd_results_plot(d.cd8.nonc.durable)
	,
	width = 6,
	file = insert(pdf.fn, c("wd", "cd8", "non-c", "durable-vs-transient", "disc-cont-fdr"))
);

qwrite(
	filter(d.cd8.nonc.durable, hurdle_q < csv.fdr.cut),
	file = insert(csv.fn, c("wd", "cd8", "non-c", "durable-vs-transient"))
);

qwrite(
	d.cd8.nonc.durable,
	file = insert(rds.fn, c("wd", "cd8", "non-c", "durable-vs-transient"))
);

