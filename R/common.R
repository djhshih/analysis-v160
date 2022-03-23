bound <- function(xs, lim) {
	lower <- min(lim);
	upper <- max(lim);
	ifelse(xs < lower, lower,
		ifelse(xs > upper, upper, xs)
	)
}

revlog_trans <- function(base = 10) {
	library(scales)
	trans <- function(x) -log(x, base)
	inv <- function(x) base^(-x)
	trans_new(
		paste0("revlog-", format(base)), trans, inv,
		log_breaks(base = base), 
		domain = c(1e-100, Inf)
	)
}

levels.response <- c("nonresponsive", "transient", "durable");

library(ggsci)
cols.response <- pal_jama()(3);
names(cols.response) <- levels.response;

capitalize <- function(s) {
	ifelse(is.na(s),
		NA, 
		paste0(toupper(substring(s, 1, 1)), substring(s, 2))	
	)
}

rename_hallmarks <- function(x, ...) {
	UseMethod("rename_hallmarks", x)
}

rename_hallmarks.default <- function(x, ...) {
	x	<- tolower(gsub("_", " ", sub("HALLMARK_", "", x)));

	library(magrittr)
	x %<>% gsub("kras ", "KRAS ", .) %>%
		gsub("nfkb", "NFKB", .) %>%
		gsub("e2f ", "E2F ", .) %>%
		gsub("tnfa ", "TNF-alpha ", .) %>%
		gsub("myc ", "MYC ", .) %>%
		gsub("il(\\d) ", "IL\\1 ", .) %>%
		gsub("stat(\\d) ", "STAT\\1 ", .) %>%
		gsub("mtorc", "mTORC", .) %>%
		gsub("mtor", "mTOR", .) %>%
		gsub("pi3k ", "PI3K ", .) %>%
		gsub("beta catenin", "beta-catenin", .) %>%
		gsub("tgf beta", "TGF-beta", .) %>%
		gsub("akt ", "Akt ", .) %>%
		gsub("jak ", "Jak ", .) %>%
		gsub("wnt ", "Wnt ", .) %>%
		gsub("notch ", "Notch ", .) %>%
		gsub("uv ", "UV ", .) %>%
		gsub("g2m ", "G2M ", .) %>%
		gsub("dna ", "DNA ", .);
	x <- capitalize(x) %>%
		gsub("P53", "p53", .) %>%
		gsub("MTORC", "mTORC", .);

	x
}

rename_hallmarks.array <- 
rename_hallmarks.matrix <-
rename_hallmarks.data.frame <- function(x, ...) {
	rownames(x) <- rename_hallmarks(rownames(x));
	x
}

library(ggsci)
library(ggtext)
library(ggplot2)
library(ggpubr)
library(cowplot)

tint <- function(x, fac) {
    x <- c(grDevices::col2rgb(x))
    x <- (x + (255 - x) * fac) / 255
    grDevices::rgb(x[1], x[2], x[3])
}

shade <- function(x, fac) {
    x <- c(grDevices::col2rgb(x))
    x <- x * (1 - fac) / 255
    grDevices::rgb(x[1], x[2], x[3])
}


gene_expr_plot <- function(sce, gene, group_by, cols=NULL, rasterize=FALSE) {
	g.prop <- gene_prop_plot(sce, gene, group_by, cols=cols) +
		theme(axis.text.x = element_blank());
	g.vio <- gene_violin_plot(sce, gene, group_by, cols=cols, rasterize=rasterize);

	annotate_figure(
		cowplot::plot_grid(g.prop, g.vio,
			ncol=1, align="v", rel_heights=c(0.6, 1)),
		top = gene
	)
}

gene_violin_plot <- function(sce, gene, group_by, cols=NULL, max.pts=100, rasterize=FALSE) {

	gene.d <- data.frame(
		group = colData(sce)[[group_by]],
		logcount = as.numeric(logcounts(sce[gene, ]))
	);

	if (is.null(cols)) {
		cols <- pal_npg()(length(levels(gene.d$group)));
	}

	# only consider showing non-zero values
	gene.d$show = gene.d$logcount > 0;

	# on average, show a maximum of max.pts points, split among the groups
	for (g in levels(gene.d$group)) {
		show.idx <- gene.d$show & gene.d$group == g;
		n.show <- sum(show.idx);
		if (max.pts < n.show) {
			p.show <- max.pts / n.show / length(levels(gene.d$group));
			gene.d$show[show.idx] <- sample(c(FALSE, TRUE),
				size=sum(show.idx), prob=c(1 - p.show, p.show), replace=TRUE);
		}
	}

	p <- kruskal.test(logcount ~ group, gene.d)$p.value;

	summary.d <- group_by(gene.d, group) %>% summarize(mean = mean(logcount));

	if (is.null(names(cols))) {
		names(cols) <- levels(gene.d$group);
	}
	cols2 <- unlist(lapply(cols, function(cc) shade(cc, 0.3)));

	if (rasterize) {
		rasterize.f <- ggrastr::rasterize;
	} else {
		rasterize.f <- identity;
	}

	ggplot(gene.d, aes(x=group, y=logcount, fill=group)) +
		theme_classic() +
		geom_violin(colour=NA) + 
		rasterize.f(geom_jitter(aes(colour=group), width=0.1,
			data=filter(gene.d, show), alpha=0.6)) +
		geom_violin(fill=NA) + 
		geom_point(aes(x=group, y=mean), data=summary.d, shape=23, size=2, fill="white") +
		scale_fill_manual(values=cols, guide="none") +
		scale_colour_manual(values=cols2, guide="none") +
		xlab("") + ylab("log count") +
		ggtitle(sprintf("*p* = %s", format(p, digits=2))) +
		theme(
			plot.title = element_markdown(size=10),
			axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
		)
}

gene_prop_plot <- function(sce, gene, group_by, cols=NULL, threshold=0) {
	gene.d <- data.frame(
		group = colData(sce)[[group_by]],
		logcount = as.numeric(logcounts(sce[gene, ]))
	);

	if (is.null(cols)) {
		cols <- pal_npg()(length(levels(gene.d$group)));
	}

	p <- with(gene.d, fisher.test(table(group, logcount > threshold))$p.value);

	prop <- function(y) {
		r <- prop.test(sum(y > threshold, na.rm=TRUE), sum(!is.na(y)));
		data.frame(
			y = r$estimate * 100,
			ymin = max(0, r$conf.int[1]) * 100,
			ymax = min(1, r$conf.int[2]) * 100
		)
	}

	summary.d <- group_by(gene.d, group) %>% 
		summarize(prop(logcount));

	ggplot(summary.d, aes(x=group, y=y, ymin=ymin, ymax=ymax, fill=group)) +
		theme_classic() +
		geom_col() + 
		geom_errorbar(width=0.2) +
		scale_fill_manual(values=cols, guide="none") +
		xlab("") + ylab("% expressed") +
		ggtitle(sprintf("*p* = %s", format(p, digits=2))) +
		theme(
			plot.title = element_markdown(size=10),
			axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
		)
}

