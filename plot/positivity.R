library(io)
library(ggplot2)
library(ggsci)
library(binom)

d <- qread("positivity.tsv")
ci <- binom.confint(d$positive, d$total, method="agresti-coull", conf.level=0.8);
ci$lower <- pmax(ci$lower, 0);
ci$upper <- pmin(ci$upper, 1);
d <- cbind(d, ci[, c("mean", "lower", "upper")]);
d$group <- factor(d$group, levels=unique(d$group));

qdraw(
	ggplot(d, aes(x=group, y=mean, ymin=lower, ymax=upper, fill=antigen)) +
		theme_classic() +
		geom_col() + geom_errorbar(width=0.3) +
		facet_grid(tcell ~ antigen) +
		scale_fill_manual(values=c("IE-1"="#BC3C29", "pp65"="#0072B5"), guide=FALSE) +
		xlab("") + ylab("proportion of responders") +
		theme(
			strip.background = element_blank(),
			axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
			strip.text.y = element_text(angle=0)
		) +
		ylim(0, 1)
	,
	width = 3, height = 4,
	file = "v160_prop-responder.pdf"
)
