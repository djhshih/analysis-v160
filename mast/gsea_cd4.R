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

source("R/common.R")

options(mc.cores=120)


in.fn <- as.filename("v160_cd4_mast-zlm.rds");
zfit <- qread(in.fn);
mcold <- qread("v160_features.rds");

out.fn <- filename(in.fn$fstem, tag=setdiff(in.fn$tag, "mast-zlm"));
rds.fn <- insert(out.fn, ext="rds");
pdf.fn <- insert(out.fn, ext="pdf");
boots.fn <- insert(rds.fn, c("mast-zlm", "boots"));

if (file.exists(tag(boots.fn))) {
	boots <- qread(boots.fn);
	message("Using cached data from ", tag(boots.fn));
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
	transient = CoefficientHypothesis("responsetransient", c("(Intercept)", "responsetransient"))
);

gseas <- lapply(hypotheses,
	function(hypothesis) {
		gseaAfterBoot(zfit, boots, indices, hypothesis)
	}
);

sms <- lapply(gseas, function(gsea) summary(gsea, testType="normal"));

sig.sets <- lapply(sms, function(sm) {
	filter(sm, combined_adj < fdr.cut)$set
	#filter(sm, combined_adj < fdr.cut, abs(combined_Z) > 3)$set
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

qdraw(
	ggplot(sig.sm, aes(y=set, x=group, fill=bound(combined_Z, c(-12, 12)))) + theme_classic() +
		geom_tile() + 
		scale_fill_viridis(name="", option="mako", limits=c(-12, 12)) +
		theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
		theme(legend.position="bottom") +
		coord_fixed() +
		xlab("") + ylab("")
	,
	width=6, height=10,
	file=insert(pdf.fn, c("gsea"))
);

