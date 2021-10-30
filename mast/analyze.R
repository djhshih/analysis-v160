library(DropletUtils)
library(io)
library(scuttle)  # logNormCounts
library(MAST)
library(dplyr)

indir <- "../aggr/5p/pt26/outs/count"

options(mc.cores=128);

set.seed(1334);

fdr.cut <- 0.05;


out.fn <- filename("v160");
rds.fn <- insert(out.fn, ext="rds");
pdf.fn <- insert(out.fn, ext="pdf");

barcode.d <- qread("../../tcr-profiling/tcr/merged/tcr_aggr_responsive_barcodes.csv");

x <- read10xCounts(file.path(indir, "filtered_feature_bc_matrix"));

br <- barcodeRanks(counts(x));

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

# convert to SingellCellAssay for use with MAST
sca <- SceToSingleCellAssay(logNormCounts(x));

cold <- left_join(as.data.frame(colData(sca)), barcode.d, by=c("Barcode"="barcode"));
cold$response[is.na(cold$response)] <- "nonresponsive";
cold$response <- factor(cold$response, c("nonresponsive", "transient", "durable"));

colData(sca) <- DataFrame(cold);

mcold <- data.table::as.data.table(mcols(sca));
qwrite(mcold, insert(rds.fn, "features"));

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
extract_sig_genes <- function(res, component="hurdle", fdr.cut=0.05, symbol=TRUE) {
	ps <- sort(res[, component, "Pr(>Chisq)"]);
	qs <- p.adjust(ps, "BH");

	idx <- qs < fdr.cut;
	sigs <- names(hurdle.qs)[idx]
	if (symbol) {
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

sgenes.wd <- list(
	durable = durable <- extract_sig_genes(wd$durable),
	transient = extract_sig_genes(wd$transient),
	durable.vs.transient = extract_sig_genes(wd$durable.vs.transient)
);

compare_sets(sgenes.wd)

spids.lr <- extract_sig_genes(lr, symbol=FALSE);
sgenes.lr <- pid_to_genes(spids.lr);
qwrite(sgenes.lr, insert(out.fn, c("lr", "sig-genes"), ext="vtr"));

coefs.sig <- coefs[match(spids.lr, coefs$primerid), ];
sca.d.sig <- as(sca[spids.lr, ], "data.table");

