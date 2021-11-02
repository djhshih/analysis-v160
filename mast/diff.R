library(io)
library(MAST)


options(mc.cores=128);

set.seed(1334);

fdr.cut <- 0.05;

in.fn <- "v160_sca.rds";

out.fn <- as.filename(in.fn);
out.fn$tag <- setdiff(out.fn$tag, "sca");

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

