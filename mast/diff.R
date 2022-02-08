library(io)
library(MAST)
library(ggplot2)
library(ggrepel)
library(dplyr)

source("R/common.R")

options(mc.cores=128);

set.seed(1337);

fdr.cut <- 0.05;

in.fn <- as.filename("v160_sca.rds");
mcold <- qread("v160_features.rds");

out.fn <- filename(in.fn$fstem, tag=setdiff(in.fn$tag, "sca"));
rds.fn <- insert(out.fn, ext="rds");
pdf.fn <- insert(out.fn, ext="pdf");

sca <- qread(in.fn);

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
	function(d, logfc) {
		data.frame(
			gene = rownames(d),
			logfc = logfc,
			discrete_q = p.adjust(d[, "disc", "Pr(>Chisq)"], "BH"),
			continuous_q = p.adjust(d[, "cont", "Pr(>Chisq)"], "BH"),
			hurdle_q = p.adjust(d[, "hurdle", "Pr(>Chisq)"], "BH")
		)
	},
	wd,
	logfcs,	
	SIMPLIFY = FALSE
);

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
	prob <- -log(d$hurdle_q);
	prob[is.na(prob)] <- 0;
	idx <- union(sig.idx, sample(1:nrow(d), n.ns, prob=prob));

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
	prob <- -log(d$q);
	prob[is.na(prob)] <- 0;
	idx <- union(sig.idx, sample(1:nrow(d), n.ns, prob=prob));

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


spids.lr <- extract_sig_genes(lr);
sgenes.lr <- pid_to_genes(spids.lr);
qwrite(sgenes.lr, insert(out.fn, c("lr", "sig-genes"), ext="vtr"));

coefs.sig <- coefs[match(spids.lr, coefs$primerid), ];
sca.d.sig <- as(sca[spids.lr, ], "data.table");

