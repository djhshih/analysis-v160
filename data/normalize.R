library(io)
library(affy)
library(frma)
library(hgu133plus2frmavecs)
library(hgu133plus2barcodevecs)
library(AnnotationDbi)
library(hgu133plus2.db)
library(pamr)
library(switchBox)
library(scran)


pheno0 <- qread("annot/E-TABM-40.sdrf.txt", type="tsv", comment="");
batch <- ReadAffy(filenames=list.files("cel", full=TRUE));
eset <- frma(batch);

features.use <- qread("../mast/v160_genes-expressed.rds");

out.fn <- filename("etabm40");
rds.fn <- insert(out.fn, ext="rds");
pdf.fn <- insert(out.fn, ext="pdf");

pheno <- data.frame(
	array_id = colnames(eset),
	sample_id = pheno0$"Sample Name"[match(colnames(eset), pheno0[["Array Data File"]])]
);
pheno$group <- factor(gsub("\\d+", "", pheno$sample_id), levels=c("N", "CM", "EM", "EMRA"));


# arrays with median GNUSE >= 1.25 is considered poor quality
GNUSE(eset, type="stats")

qwrite(eset, insert(rds.fn, "eset"));

# get feature from matrix (private)
get_feature <- function(z, gene) {
	d <- data.frame(
		group = pheno$group,
		t(z[gene, , drop=FALSE])
	)
	rownames(d) <- NULL;

	d[order(d$group), ]
}

aggregate_rows <- function(x, rows) {
	x.a <- aggregate(x, by=list(rows), max);
	rownames(x.a) <- x.a[,1];
	as.matrix(x.a[, -1]);
}

# convert to within-sample percentiles
percentile_within <- function(x) {
	apply(x, 2,
		function(z) {
			idx <- which(z > 0);
			z[idx] <- rank(z[idx], ties.method="average") / length(idx);
			z
		}
	)
}



probes <- rownames(eset);
genes <- mapIds(hgu133plus2.db, keys=probes, column="SYMBOL", keytype="PROBEID");

bc <- barcode(eset);

ccr7.probes <- probes[which(genes == "CCR7")];
gzmb.probes <- probes[which(genes == "GZMB")];

get_feature(exprs(eset), ccr7.probes)

get_feature(bc, ccr7.probes)

bc.a <- aggregate_rows(bc, genes)
# discriminative ability of CCR7 is lost after barcoding
bc.a["CCR7", ]

z <- barcode(eset, output="z-score")
z[z < 0] <- 0;

get_feature(z, ccr7.probes)
get_feature(z, gzmb.probes)

z.a <- aggregate_rows(z, genes);

# discriminative ability is retained
gene_expr(z.a, "CCR7")

qdraw(
	{hist(z, breaks=100)}
	,
	file = insert(pdf.fn, c("z", "hist"))
)

qdraw(
	{hist(z.a, breaks=100)}
	,
	file = insert(pdf.fn, c("z-agg", "hist"))
)

z.a.s0 <- z.a[rownames(z.a) %in% features.use, ];

scores <- scoreMarkers(z.a.s0, pheno$group);
marker.sets <- lapply(scores,
	function(s)
		rownames(s)[
			pmax(s$max.AUC, 1 - s$max.AUC) > 0.9 &
			pmax(s$mean.AUC, 1 - s$mean.AUC) > 0.8 &
			abs(s$max.logFC.cohen) > 3 &
			abs(s$mean.logFC.cohen) > 1
		]
);
markers <- Reduce(union, marker.sets)

length(markers)
lapply(marker.sets, length)
lapply(marker.sets, function(s) length(intersect(s, markers)) / length(markers))

# find mean of scores over classes
auc.mat <- do.call(cbind, lapply(scores,
	function(s) {
		pmax(s$mean.AUC, 1 - s$mean.AUC)
	}
));
auc.best <- apply(auc.mat, 1, max);

cohen.mat <- do.call(cbind, lapply(scores,
	function(s) {
		abs(s$mean.logFC.cohen)
	}
));
cohen.best <- apply(cohen.mat, 1, max);

# sort markers by best score
markers <- markers[order(auc.best[markers], cohen.best[markers], decreasing=TRUE)];


z.a.s1 <- z.a.s0[markers, ];

z.a.r1 <- percentile_within(z.a.s1);
summary(z.a.r1)

get_feature(z.a.r1, "CCR7")
get_feature(z.a.r1, "GZMB")

####

# select features that are different across groups

data1 <- list(x = z.a.r1, y = pheno$group);
fit1 <- pamr.train(data1);
fit1

cv1 <- pamr.cv(fit1, data1);
cv1

# minimum threshold with the minimum error
threshold1 <- max(cv1$threshold[cv1$error == min(cv1$error)]);

yhat <- pamr.predict(fit1, data1$x, threshold = threshold1);
table(yhat, y=pheno$group)

idx.nonzero1 <- pamr.predict(fit1, threshold = threshold1, type = "nonzero");
features.nonzero <- rownames(data1$x)[idx.nonzero1];

lapply(marker.sets, function(s) length(intersect(s, features.nonzero)) / length(features.nonzero))

grep("CCR", features.nonzero, value=TRUE)
grep("GZM", features.nonzero, value=TRUE)


z.a.s2 <- z.a.s[idx.nonzero1, ];
z.a.r2 <- percentile_within(z.a.s2);

data2 <- list(x = z.a.r2, y = pheno$group);
fit2 <- pamr.train(data2);
fit2

cv2 <- pamr.cv(fit2, data2);
cv2

threshold2 <- max(cv2$threshold[cv2$error == min(cv2$error)]);

idx.nonzero2 <- pamr.predict(fit2, threshold = threshold2, type = "nonzero");
features.nonzero <- rownames(data2$x)[idx.nonzero2];

grep("CCR", features.nonzero, value=TRUE)
grep("GZM", features.nonzero, value=TRUE)

z.a.s3 <- z.a.s[idx.nonzero2, ];
z.a.r3 <- percentile_within(z.a.s3);

data3 <- list(x = z.a.r3, y = pheno$group);
fit3 <- pamr.train(data3);
fit3

cv3 <- pamr.cv(fit3, data3);
cv3

threshold3 <- max(cv3$threshold[cv3$error == min(cv3$error)]);
idx.nonzero3 <- pamr.predict(fit3, threshold = threshold3, type = "nonzero");
features.nonzero <- rownames(data3$x)[idx.nonzero3];

z.a.s4 <- z.a.s[idx.nonzero3, ];
z.a.r4 <- percentile_within(z.a.s4);

data4 <- list(x = z.a.r4, y = pheno$group);
fit4 <- pamr.train(data4);
fit4

cv4 <- pamr.cv(fit4, data4);
cv4

threshold4 <- max(cv4$threshold[cv4$error == min(cv4$error)]);

idx.nonzero4 <- pamr.predict(fit4, threshold = threshold4, type = "nonzero");
features.nonzero <- rownames(data4$x)[idx.nonzero4];

z.a.s5 <- z.a.s[idx.nonzero4, ];
z.a.r5 <- percentile_within(z.a.s5);

data5 <- list(x = z.a.r5, y = pheno$group);
fit5 <- pamr.train(data5);
fit5

cv5 <- pamr.cv(fit5, data5);
cv5

threshold5 <- max(cv5$threshold[cv5$error == min(cv5$error)]);

idx.nonzero5 <- pamr.predict(fit5, threshold = threshold5, type = "nonzero");
features.nonzero <- rownames(data5$x)[idx.nonzero5];


# select a model
fit <- fit1;
fit$threshold.opt <- threshold1;
fit$predict <- function(fit, x, type="class") {
	library(pamr)

	# extract relevant features from data matrix x
	features <- rownames(fit$centroids);
	idx <- match(features, rownames(x));
	stopifnot(!is.na(idx));
	x <- x[idx, ];

	# apply rank transformation
	x <- apply(x, 2,
		function(z) {
			idx <- which(z > 0);
			out <- numeric(length(z));
			out[idx] <- rank(z[idx], ties.method="average") / length(idx);
			out
		}
	)

	pamr.predict(fit, x, threshold = fit$threshold.opt, type=type)
}

rownames(fit$centroids)
grep("CCR", rownames(fit$centroids), value=TRUE)
grep("GZM", rownames(fit$centroids), value=TRUE)

lapply(marker.sets, function(s) length(intersect(s, rownames(fit$centroids))))

lapply(marker.sets, function(s) length(intersect(s, rownames(fit$centroids))) / length(fit$centroids))


yhat <- fit$predict(fit, z.a);
table(yhat, y=pheno$group)

qwrite(fit, insert(rds.fn, "pamr"));

####

#SWAP.Train.KTSP(data1$x, data1$y)

