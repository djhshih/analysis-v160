# Perform gene set enrichment anlaysis
# using the blood transcription modules database
# (https://dx.doi.org/10.1038%2Fni.2789)

library(GSEABase)
library(MAST)
library(io)
library(limma)
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(RColorBrewer)
library(fgsea)


source("../R/common.R")

options(mc.cores=50)


#in.fn <- as.filename("v160_mast-zlm.rds");
in.fn <- as.filename("v160_cd8_mast-zlm.rds");
#in.fn <- as.filename("v160_cd8_durable_mast-zlm_c-vs-b.rds");

hypotheses <- list(
	durable = CoefficientHypothesis("responsedurable", c("(Intercept)", "responsedurable")),
	transient = CoefficientHypothesis("responsetransient", c("(Intercept)", "responsetransient"))
);


zfit <- qread(in.fn);
mcold <- qread("v160_features.rds");

out.fn <- filename(in.fn$fstem, tag=setdiff(in.fn$tag, "mast-zlm"));
rds.fn <- insert(out.fn, ext="rds");
pdf.fn <- insert(out.fn, ext="pdf");
boots.fn <- insert(rds.fn, c("mast-zlm", "boots"));

boots <- qcache(
	{
		boots <- bootVcov1(zfit, R=100);
	},
	file = as.character(boots.fn, simplify=TRUE)
);

run_gseas <- function(zfit, gids, hypotheses, min.genes=5) {

	indices <- limma::ids2indices(gids, mcold$Symbol);
	# filter based on number of genes in gene set
	indices <- indices[sapply(indices, length) >= min.genes];

	gseas <- lapply(hypotheses,
		function(hypothesis) {
			gseaAfterBoot(zfit, boots, indices, hypothesis)
		}
	);

	gseas
}

plot_gseas <- function(gseas, fdr.cut=0.1, z.cut=1) {
	sms <- lapply(gseas, function(gsea) summary(gsea, testType="normal"));

	sig.sets <- lapply(sms, function(sm) {
		filter(sm, combined_adj < fdr.cut, abs(combined_Z) > z.cut)$set
	});

	sm <- do.call(rbind, mapply(
		function(sm, response) {
			select(sm, set, disc_effect, combined_Z, combined_P, combined_adj) %>%
				mutate(group = response)
		}
		,
		sms,
		names(sms),
		SIMPLIFY = FALSE
	));

	# summary for sets that are significant in any group
	sig.sm <- filter(sm,
		set %in% unlist(sig.sets)
	);

	# order set by combined_Z
	sig.sm$set <- factor(sig.sm$set, levels=c(unique(sms[[1]]$set[order(sms[[1]]$combined_Z)])));
	#sig.sm$group <- factor(sig.sm$group, levels=unique(sig.sm$group));

	cols <- brewer.pal(9, "RdYlBu");
	col.high <- cols[1];
	col.low <- cols[9];
	col.mid <- cols[5];

	z.limits <- c(-6, 6);
	q.limits <- c(0.25, 1e-6);

	ggplot(sig.sm, aes(y=set, x=group, colour=bound(combined_Z, z.limits), size=bound(combined_adj, q.limits))) + 
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

module.fn <- list.files(system.file("extdata", package="MAST"), pattern = "BTM", full.names = TRUE);
btm <- getGmt(module.fn);
gids <- geneIds(btm);
# remove gene sets that are unknown (TBA)
# remove gene sets pertaining to non-T cells
gids <- gids[
	grep("(B cell)|(plasma cell)|(neutrophil)|(NK cell)|(monocyte)|(myeloid cell)|(dendritic cell)|(DC)|(TBA)",
		names(gids), invert=TRUE, ignore.case=TRUE
	)
];

gseas.btm <- run_gseas(zfit, gids, hypotheses);

qdraw(plot_gseas(gseas.btm, fdr.cut=0.01, z.cut=2),
	width=6, height=12,
	file=insert(pdf.fn, c("gsea", "btm"))
);

# -

h.all <- getGmt("../msigdb/h.all.v7.5.1.symbols.gmt");
gids <- geneIds(h.all);
names(gids) <- rename_hallmarks(names(gids));

gseas.h <- run_gseas(zfit, gids, hypotheses);

qdraw(plot_gseas(gseas.h, fdr.cut=0.1, z.cut=2),
	width=6, height=4,
	file=insert(pdf.fn, c("gsea", "hallmark"))
);


wd.fn <- as.filename("v160_wd_cd8_durable-vs-transient.rds");
wd <- qread(wd.fn);

# positive gene set
genes <- gids[["TNF-alpha signaling via NFKB"]];
wd.sel <- filter(wd[genes, ], hurdle_p < 0.05)
# re-perform multiple hypothesis correction
wd.sel$hurdle_q <- p.adjust(wd.sel$hurdle_p, "BH");
filter(wd.sel, hurdle_q < 0.05, logfc > 0)

# negative gene set
genes <- gids[["Interferon alpha response"]];
wd.sel <- filter(wd[genes, ], hurdle_p < 0.05)
# re-perform multiple hypothesis correction
wd.sel$hurdle_q <- p.adjust(wd.sel$hurdle_p, "BH");
filter(wd.sel, hurdle_q < 0.05, logfc < 0)

# TODO
# plot genes in positive and negative gene sets

