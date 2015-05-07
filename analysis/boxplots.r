#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))


theme_set(theme_bw(16) + theme(strip.background = element_blank()))

agg <- read.csv("../sims/tmp.csv", as.is = TRUE, comment='#')


ratio.by.replicate.bxplot <- function(agg) {
    df <- agg %>%
        extract(root, c('replicate', 'tevent', 'donor', 'recipient'), 'replicate_(\\d+)/(.*)/(.*)/(.*)', remove=TRUE, convert=TRUE) %>%
        mutate(score_ratio = p_score / b_score)
    

    ggplot(df, aes(factor(replicate), score_ratio)) +
        geom_boxplot() +
        scale_x_discrete(labels=sprintf('replicate %s', levels(factor(df$replicate)))) +
        xlab(NULL) +
        ylab("Ratio N-W score: PRANK / Beast") +
        ggtitle("PRANK/Beast ratio of founder inference score\nacross replicated simulations") +
        theme(plot.title = element_text(lineheight=.8, face="bold"))
}



ratio.by.xmit.bxplot <- function(agg) {
    # boxplot of prank/beast ratio at different transmission points,
    # averaged across replicates

    df <- agg %>%
        extract(root, c('replicate', 'tevent', 'donor', 'recipient'), 'replicate_(\\d+)/(.*)/(.*)/(.*)', remove=TRUE, convert=TRUE) %>%
        group_by(tevent, donor, recipient) %>%
        summarise_each(funs(mean)) %>%
        mutate(score_ratio = p_score / b_score)


    ggplot(df, aes(factor(tevent), score_ratio)) +
        geom_boxplot() +
        scale_x_discrete(labels=sprintf("@%s generations", levels(factor(df$tevent)))) +
        xlab("Transmission event") +
        ylab("Ratio N-W score: PRANK / Beast") +
        ggtitle("PRANK/Beast score ratio\naveraged over replicates")
}

consensus_ratio.by.xmit.bxplot <- function(agg) {
    # boxplot of prank/beast ratio at different transmission points,
    # averaged across replicates

    df <- agg %>%
        extract(root, c('replicate', 'tevent', 'donor', 'recipient'), 'replicate_(\\d+)/(.*)/(.*)/(.*)', remove=TRUE, convert=TRUE) %>%
        group_by(tevent, donor, recipient) %>%
        summarise_each(funs(mean)) %>%
        mutate(score_ratio = p_score / c_score)


    ggplot(df, aes(factor(tevent), score_ratio)) +
        geom_boxplot() +
        scale_x_discrete(labels=sprintf("@%s generations", levels(factor(df$tevent)))) +
        xlab("Transmission event") +
        ylab("Ratio N-W score: PRANK / Consensus") +
        ggtitle("PRANK/Consensus score ratio\naveraged over replicates")
}

scores.by.method.density <- function(agg) {
    # plot density of scores
    agg %>%
        extract(root, c('replicate', 'tevent', 'donor', 'recipient'), 'replicate_(\\d+)/(.*)/(.*)/(.*)', remove=TRUE, convert=TRUE) %>%
        gather(measure, score, c(p_score, b_score, c_score)) %>%
        ggplot(aes(x=score)) +
        ggtitle("Density of scores by inference method") +
        geom_density(aes(fill=factor(measure)), size=1, alpha=0.3) +
        scale_fill_discrete(name="Inference\nmethod", labels=c("PRANK","Beast", "Consensus"))
}

gaps.by.method.density <- function(agg) {
    # plot density of scores
    agg %>%
        extract(root, c('replicate', 'tevent', 'donor', 'recipient'), 'replicate_(\\d+)/(.*)/(.*)/(.*)', remove=TRUE, convert=TRUE) %>%
        # REMIND - swap these when the distance calculation is fixed
        gather(measure, score, c(p_gaps, b_gaps, c_gaps)) %>%
        ggplot(aes(x=score)) +
        ggtitle("Proportion gap sites\nbetween inferred and actual founder") +
        geom_density(aes(fill=factor(measure)), size=1, alpha=0.3) +
        xlab("Proportion of gapped sites") +
        scale_fill_discrete(name="Inference\nmethod", breaks=c('p_gaps', 'b_gaps', 'c_gaps'),labels=c("PRANK","Beast", "Consensus"))
}

match.by.method.density <- function(agg) {
    agg %>%
        extract(root, c('replicate', 'tevent', 'donor', 'recipient'), 'replicate_(\\d+)/(.*)/(.*)/(.*)', remove=TRUE, convert=TRUE) %>%
        # REMIND - swap these when the distance calculation is fixed
        gather(measure, score, c(p_identity, b_identity, c_identity)) %>%
        ggplot(aes(x=score)) +
        ggtitle("Proportion of matching sites\nbetween inferred and actual founder") +
        geom_density(aes(fill=factor(measure)), size=1, alpha=0.3) +
        xlab("Proportion of identical sites") +
        scale_fill_discrete(name="Inference\nmethod", breaks=c('p_identity', 'b_identity', 'c_identity'),labels=c("PRANK","Beast", "Consensus"))
}

# create all plot objects.
# we will arrange print them afterward.
plot.ratio.by.replicate <- ratio.by.replicate.bxplot(agg)
plot.ratio.by.xmit <- ratio.by.xmit.bxplot(agg)
plot.scores.by.method <- scores.by.method.density(agg)
plot.gaps.by.method <- gaps.by.method.density(agg)
plot.match.by.method <- match.by.method.density(agg)
plot.consensus.ratio.by.xmit <- consensus_ratio.by.xmit.bxplot(agg)




# find all plot objects and print them!
# http://stackoverflow.com/a/20502085/1135316
plots <- ls(pattern=glob2rx('plot.*'))
l = mget(plots)

# save one plot per file, in different formats
bquiet = mapply(ggsave, file=paste0(names(l), ".pdf"), plot=l, width = 7, height = 5)
bquiet = mapply(ggsave, file=paste0(names(l), ".svg"), plot=l, width = 7, height = 5)
bquiet = mapply(ggsave, file=paste0(names(l), ".png"), plot=l, width = 7, height = 5)

