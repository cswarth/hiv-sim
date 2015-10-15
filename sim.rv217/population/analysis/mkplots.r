#!/usr/bin/env Rscript 
# expect to be in ~/src/matsen/hiv-sim/sim.rv217/population/

source('../../R/utils.r')

load(file="distances.Rda")

make_plot <- function(tmp, allow.indels) {
    .facet_labeller <- function(labels, sep = ", ") {
        var2str <- function(variable, value) {
            switch(variable,
                   indel=ifelse(value, "With indels", "No indels"),
                   fitness=ifelse(value =='noselection', 'no fitness', paste(value, 'selection')),
                   as.character(value))
        }

        mapply(var2str, names(labels), labels) %>% apply(1, function(s) paste0(s, collapse=sep)) %>% list()
    }


    tmp %>% filter(indel==allow.indels) %>% 
        gather(method, dist, starts_with("dist.")) %>%
        mutate(method=sub('dist.','',method)) %>%
        mutate(method=factor(method)) %>%
        mutate(nseq=factor(nseq)) %>%
        filter(method %in% c('unguided')) %>%
        #filter(fitness=='noselection') %>%

        ggplot(aes(x=wpi, y=dist, fill=population)) +
        facet_wrap(indel ~ fitness, labeller=.facet_labeller, ncol=1) +
        #facet_grid(indel ~ fitness, labeller=facet_labeller) +
        theme(axis.text.x = element_text(size = 8, colour = "red", angle = 45)) +
        geom_boxplot() +
        xlab("Generations") +
        ylab("N-W Pairwise Distance") +
        ggtitle("Unguided Prank") +
        scale_fill_discrete(name="Grouping",
                            breaks=c('10', '20', '40', 
                                '1000', '3000', '5000', '7000', '10000',
                                "beast", "consensus", 'pcodon', 'pdna', 'unguided', 
                                'relaxed', 'strict'),
                            labels=c('20 sequences', '40 sequences', '80 sequences', 
                                'population 1000', 'population 3000', 'population 5000', 'population 7000', 'population 10000',
                                "Beast", "Consensus", 'Prank codon', 'Prank dna', 'Prank unguided', 
                                'Beast relaxed', 'Beast strict')) 


}

save.plots <- function(g, show.indels, basename='population', width=8, height=11) {
    indelstr <- ifelse(show.indels,'indels','noindels')
    png(file=sprintf('%s_%s.png', basename, indelstr), width=width, height=height, units = "in", res=300)
    print(g)
    invisible(dev.off())

    pdf(file=sprintf('%s_%s.pdf', basename, indelstr), width=width, height=height)
    print(g)
    invisible(dev.off())
}


save.plots(make_plot(tmp, FALSE), FALSE)
save.plots(make_plot(tmp, TRUE), TRUE)
