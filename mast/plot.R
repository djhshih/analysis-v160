library(scater)
library(io)
library(MAST)
library(ggplot2)
library(ggsci)
library(viridis)

options(mc.cores=128);

set.seed(1334);

in.fn <- as.filename("v160_sce.rds");

out.fn <- filename(in.fn$fstem, tag=setdiff(in.fn$tag, "sce"));
rds.fn <- insert(out.fn, ext="rds");
pdf.fn <- insert(out.fn, ext="pdf");

sce <- qread(in.fn);

sce <- logNormCounts(sce);
sce <- runUMAP(sce);
sce <- runTSNE(sce);

sce.r <- sce[, colData(sce)$response != "nonresponsive"];
#colData(sce.r)$response <- droplevels(colData(sce.r)$response);
sce.r <- runUMAP(sce.r);
sce.r <- runTSNE(sce.r);

cols.response <- pal_jama()(3);
names(cols.response) <- levels(colData(sce)$response);

alpha <- 0.5;

qdraw(
	ggcells(sce, mapping=aes(x=UMAP.1, y=UMAP.2, colour=response)) +
		theme_classic() +
		facet_grid(response ~ t_cell) +
		geom_point(size=0.5, alpha=alpha) +
		#stat_density_2d() +
		coord_fixed() +
		scale_colour_manual(values = cols) +
		theme(legend.position="bottom")
	,
	width = 6, height = 12,
	file = insert(pdf.fn, c("umap", "t-cell-type", "response"))
)

qdraw(
	ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=response)) +
		theme_classic() +
		facet_grid(response ~ t_cell) +
		geom_point(size=0.5, alpha=alpha) +
		#stat_density_2d() +
		coord_fixed() +
		scale_colour_manual(values = cols) +
		theme(legend.position="bottom")
	,
	width = 6, height = 12,
	file = insert(pdf.fn, c("tsne", "t-cell-type", "response"))
)


qdraw(
	ggcells(sce.r, mapping=aes(x=UMAP.1, y=UMAP.2, colour=response)) +
		theme_classic() +
		facet_wrap(~ t_cell) +
		geom_point(alpha=alpha) +
		coord_fixed() +
		scale_colour_manual(values = cols) +
		theme(legend.position="bottom")
	,
	width = 6,
	file = insert(pdf.fn, c("umap", "responsive", "t-cell-type"))
)

qdraw(
	ggcells(sce.r, mapping=aes(x=UMAP.1, y=UMAP.2, colour=response)) +
		theme_classic() +
		facet_wrap(~ t_cell0) +
		geom_point(alpha=alpha) +
		coord_fixed() +
		scale_colour_manual(values = cols) +
		theme(legend.position="bottom")
	,
	width = 9,
	file = insert(pdf.fn, c("umap", "responsive", "t-cell-type0"))
)

rowData(sce)$ID[match("IFNG", rowData(sce)$Symbol)]

qdraw(
	ggcells(sce.r, mapping=aes(x=UMAP.1, y=UMAP.2, colour=ENSG00000111537)) +
		theme_classic() +
		facet_grid(response ~ t_cell) +
		geom_point(alpha=alpha) +
		scale_colour_viridis() +
		coord_fixed() +
		theme(legend.position="bottom")
	,
	width = 6,
	file = insert(pdf.fn, c("umap", "responsive", "t-cell-type", "response", "ifng"))
)

rowData(sce)$ID[match("PGAM1", rowData(sce)$Symbol)]

qdraw(
	ggcells(sce.r, mapping=aes(x=UMAP.1, y=UMAP.2, colour=ENSG00000171314)) +
		theme_classic() +
		facet_grid(response ~ t_cell) +
		geom_point(alpha=alpha) +
		scale_colour_viridis() +
		coord_fixed() +
		theme(legend.position="bottom")
	,
	width = 6,
	file = insert(pdf.fn, c("umap", "responsive", "t-cell-type", "response", "pgam1"))
)

rowData(sce)$ID[match("CCL4", rowData(sce)$Symbol)]

qdraw(
	ggcells(sce.r, mapping=aes(x=UMAP.1, y=UMAP.2, colour=ENSG00000275302)) +
		theme_classic() +
		facet_grid(response ~ t_cell) +
		geom_point(alpha=alpha) +
		scale_colour_viridis() +
		coord_fixed() +
		theme(legend.position="bottom")
	,
	width = 6,
	file = insert(pdf.fn, c("umap", "responsive", "t-cell-type", "response", "ccl4"))
)

rowData(sce)$ID[match("GZMB", rowData(sce)$Symbol)]

qdraw(
	ggcells(sce.r, mapping=aes(x=UMAP.1, y=UMAP.2, colour=ENSG00000100453)) +
		theme_classic() +
		facet_grid(response ~ t_cell) +
		geom_point(alpha=alpha) +
		scale_colour_viridis() +
		coord_fixed() +
		theme(legend.position="bottom")
	,
	width = 6,
	file = insert(pdf.fn, c("umap", "responsive", "t-cell-type", "response", "gzmb"))
)

rowData(sce)$ID[match("TNFRSF18", rowData(sce)$Symbol)]

qdraw(
	ggcells(sce.r, mapping=aes(x=UMAP.1, y=UMAP.2, colour=ENSG00000186891)) +
		theme_classic() +
		facet_grid(response ~ t_cell) +
		geom_point(alpha=alpha) +
		scale_colour_viridis() +
		coord_fixed() +
		theme(legend.position="bottom")
	,
	width = 6,
	file = insert(pdf.fn, c("umap", "responsive", "t-cell-type", "response", "tnfrsf18"))
)

rowData(sce)$ID[match("SRPK2", rowData(sce)$Symbol)]

qdraw(
	ggcells(sce.r, mapping=aes(x=UMAP.1, y=UMAP.2, colour=ENSG00000135250)) +
		theme_classic() +
		facet_grid(response ~ t_cell) +
		geom_point(alpha=alpha) +
		scale_colour_viridis() +
		coord_fixed() +
		theme(legend.position="bottom")
	,
	width = 6,
	file = insert(pdf.fn, c("umap", "responsive", "t-cell-type", "response", "srpk2"))
)

qdraw(
	ggcells(sce.r, mapping=aes(x=TSNE.1, y=TSNE.2, colour=response)) +
		theme_classic() +
		facet_wrap(~ t_cell) +
		geom_point(alpha=alpha) +
		coord_fixed() +
		scale_colour_manual(values = cols) +
		theme(legend.position="bottom")
	,
	width = 6,
	file = insert(pdf.fn, c("tsne", "responsive", "t-cell-type"))
)

qdraw(
	ggcells(sce.r, mapping=aes(x=TSNE.1, y=TSNE.2, colour=response)) +
		theme_classic() +
		facet_wrap(~ t_cell0) +
		geom_point(alpha=alpha) +
		coord_fixed() +
		scale_colour_manual(values = cols) +
		theme(legend.position="bottom")
	,
	width = 9,
	file = insert(pdf.fn, c("tsne", "responsive", "t-cell-type0"))
)

