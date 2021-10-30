# Perform gene set enrichment anlaysis
# using the blood transcription modules database
# (https://dx.doi.org/10.1038%2Fni.2789)

library(GSEABase)
library(MAST)
library(io)
library(limma)
library(viridis)

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
fdr.cut <- 0.10;

module.fn <- list.files(system.file("extdata", package="MAST"), pattern = module, full.names = TRUE);
gsets <- getGmt(module.fn);
gids <- geneIds(gsets);

# remove gene sets that are unknown (TBA)
# remove gene sets pertaining to non-T cells
gids <- gids[
	! names(gids) %like% "B cell" &
	! names(gids) %like% "neutrophil" &
	! names(gids) %like% "NK cell" &
	! names(gids) %like% "monocyte" &
	! names(gids) %like% "myeloid cell" &
	! names(gids) %like% "dendritic cell" &
	! names(gids) %like% "TBA"
];

indices <- limma::ids2indices(gids, mcold$Symbol);
# filter based on number of genes in modules
indices <- indices[sappy(indices, length) >= module.min.genes];


hypotheses <- list(
	durable = Hypothesis("responsedurable", c("(Intercept)", "responsedurable")),
	transient = Hypothesis("responsetransient", c("(Intercept)", "responsetransient")),
	durable.vs.transient = Hypothesis("responsedurable-responsetransient", c("responsetransient", "responsedurable"))
);

gseas <- lapply(hypotheses,
	function(hypothesis) {
		gseaAfterBoot(zfit, boots, indices, hypothesis)
	}
);

sms <- lapply(gseas, function(gsea) summary(gsea, testType="normal"));

res.ds <- lapply(sms, function(sm) {
	sigs <- sm[p.adjust(combined_adj, "BH") < fdr.cut];
	melt(sigs[, .(set, disc_Z, cont_Z, combined_Z)], id.vars="set")
})

m_plot_gsea <- function(res.d, fn) {
	qdraw(
		ggplot(res.d, aes(y=set, x=variable, fill=value)) + theme_classic() +
			geom_raster() + scale_fill_viridis()
		file = fn
	)
}

for (contrast in names(hypotheses)) {
	m_plot_gsea(res.ds[[contrast]], insert(pdf.fn, c("gsea", contrast)))
}

