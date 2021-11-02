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

options(mc.cores=120)


zfit <- qread("v160_mast-zlm.rds");
mcold <- qread("v160_features.rds");

out.fn <- filename("v160");
rds.fn <- insert(out.fn, ext="rds");
pdf.fn <- insert(out.fn, ext="pdf");
boots.fn <- insert(rds.fn, c("mast-zlm", "boots"));

if (file.exists(tag(boots.fn)) {
	boots <- qread(boots.fn);
} else {
	boots <- bootVcov1(zfit, R=100);
	qwrite(boots, boots.fn);
}

module <- "BTM";
module.min.genes <- 5;
fdr.cut <- 0.05;

module.fn <- list.files(system.file("extdata", package="MAST"), pattern = module, full.names = TRUE);
gsets <- getGmt(module.fn);
gids <- geneIds(gsets);

# remove gene sets that are unknown (TBA)
# remove gene sets pertaining to non-T cells
gids <- gids[
	grep("(B cell)|(plasma cell)|(neutrophil)|(NK cell)|(monocyte)|(myeloid cell)|(dendritic cell)|(DC)|(TBA)",
		names(gids), invert=TRUE, ignore.case=TRUE
	)
];

indices <- limma::ids2indices(gids, mcold$Symbol);
# filter based on number of genes in modules
indices <- indices[sapply(indices, length) >= module.min.genes];


hypotheses <- list(
	durable = CoefficientHypothesis("responsedurable", c("(Intercept)", "responsedurable")),
	transient = CoefficientHypothesis("responsetransient", c("(Intercept)", "responsetransient"))
	#durable.vs.transient = CoefficientHypothesis("responsedurable-responsetransient", c("responsetransient", "responsedurable"))
);

gseas <- lapply(hypotheses,
	function(hypothesis) {
		gseaAfterBoot(zfit, boots, indices, hypothesis)
	}
);

sms <- lapply(gseas, function(gsea) summary(gsea, testType="normal"));

res.ds <- lapply(sms, function(sm) {
	sigs <- sm[p.adjust(combined_adj, "BH") < fdr.cut];
	data.table::melt(sigs[, .(set, disc_Z, cont_Z, combined_Z)], id.vars="set")
})


bound <- function(xs, lim) {
	ifelse(xs < lim[1], lim[1],
		ifelse(xs > lim[2], lim[2], xs)
	)
}

m_plot_gsea <- function(res.d, ...) {
	res.d$variable <- factor(res.d$variable, levels=c("disc_Z", "cont_Z", "combined_Z"),
		labels=c("discrete", "continuous", "combined"));
	qdraw(
		ggplot(res.d, aes(y=set, x=variable, fill=bound(value, c(-12, 12)))) + theme_classic() +
			geom_tile() + 
			#scale_fill_gradientn(colours=brewer.pal(7, "BrBG"), name="", limits=c(-12, 12)) +
			scale_fill_viridis(name="", option="mako", limits=c(-12, 12)) +
			theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
			theme(legend.position="bottom") +
			xlab("") + ylab("")
		,
		...
	)
}

m_plot_gsea(res.ds$durable, file=insert(pdf.fn, c("gsea", "durable")), width=5.5, height=13)
m_plot_gsea(res.ds$transient, file=insert(pdf.fn, c("gsea", "transient")), width=4.9, height=4.7)

