library(DropletUtils)
library(io)
library(scuttle)  # logNormCounts
library(MAST)
library(dplyr)
#library(miQC)  # requires Bioconductor 3.14 and R 4.1

indir <- "../aggr/5p/pt26/outs/count"

options(mc.cores=100);

set.seed(1334);

fdr.cut <- 0.05;


out.fn <- filename("v160");
rds.fn <- insert(out.fn, ext="rds");
pdf.fn <- insert(out.fn, ext="pdf");

barcode.d <- qread("../../tcr-profiling/tcr/merged/tcr_aggr_responsive_barcodes.csv");

x <- read10xCounts(file.path(indir, "filtered_feature_bc_matrix"));

# There are no spike-ins
grep("^ERCC-", rowData(x)$ID, value=TRUE) 

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

# remove poor quality barcodes
x.f <- with(colData(x), x[, lib.factors > 0.2 & detected > 600 & subsets_mito_percent < 8]);


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

# convert to SingellCellAssay for use with MAST
sca <- SceToSingleCellAssay(logNormCounts(x.f));

cold <- left_join(as.data.frame(colData(sca)), barcode.d, by=c("Barcode"="barcode"));
cold$response[is.na(cold$response)] <- "nonresponsive";
cold$response <- factor(cold$response, c("nonresponsive", "transient", "durable"));



# --- T cells

ens.cd4 <- mcold$ID[mcold$Symbol == "CD4"]
ens.cd8a <- mcold$ID[mcold$Symbol == "CD8A"]
ens.cd8b <- mcold$ID[mcold$Symbol == "CD8B"]

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
colnames(d.t.cells) <- mcold$Symbol[match(t.cell.markers, mcold$ID)];


with(d.t.cells, table(CD8A, CD8B))
with(d.t.cells, cor(CD8A, CD8B))
with(d.t.cells, fisher.test(table(CD8A, CD8B)))

with(d.t.cells, table(CD4, CD8A | CD8B))
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
	file = rf.fn
);

# use trained classifier to predict unknown T cells
pred <- predict(rforest, t(expressed[pid.tc.inform.idx, !known.cells]));

cold$t_cell[!known.cells] <- as.character(pred);

print(with(cold, table(t_cell0, response)))
print(with(cold, table(t_cell, response)))

print(with(cold, prop.table(table(t_cell, response)[, -1], 2)))
print(with(cold, prop.table(table(t_cell, response)[, -1], 1)))
print(with(cold, fisher.test(table(t_cell, response)[, -1])))
# odds ratio = 50.07384
# p-value = 4.937e-15

# --- End T cells

colData(sca) <- DataFrame(cold);

sca.cd4 <- sca[, cold$t_cell == "CD4"];
sca.cd8 <- sca[, cold$t_cell == "CD8"];

mcold <- data.table::as.data.table(mcols(sca));
qwrite(mcold, insert(rds.fn, "features"));

qwrite(sca, insert(rds.fn, "sca"));
qwrite(sca.cd4, insert(rds.fn, c("sca", "cd4")));
qwrite(sca.cd8, insert(rds.fn, c("sca", "cd8")));

