#!/usr/bin/env Rscript

# plot distances form founder infered by prank to actual founder.
#
# Works of files created by a simulation process to mimic some of
# those that we obtain from the RV217 project.  The RV217 sample ahave
# around 10-12 sequences per sample, all taken from the c2v3c3 region
# gp-120 from the HIV-1 genome.

setwd("~/src/matsen/hiv-sim/sim.rv217/analysis")

source("utils.r")


png(file="lines_by_method.png", width=18, height=12, units = "in", res=300)
lines_by_method(tmp)
dev.off()

founder <- founder.sequence()

# recursively list files matching 'sample.fa'
files <- list.files(path='../build', pattern="sample\\.fa$", recursive=TRUE)

# The files look like,
# 	"1.0E-5/frequency_noindel/replicate_0/5000/sample.fa"
# break down the path components into seperate columns so the data
# can be filtered appropriately.

df <- as_data_frame(list(dir=dirname(files))) %>% 
    mutate(rate=factor(as.numeric(str_match(dir, "^([^/]+)/")[,2]))) %>%
    mutate(rate=reorder(rate, as.numeric(as.character(rate)), FUN=max)) %>%
    mutate(fitness=str_match(dir, "/([^/]+)_[^/]+/replicate")[,2]) %>%
    mutate(fitness=ifelse(is.na(fitness), 'none', fitness)) %>%
    mutate(fitness=factor(fitness)) %>%
    mutate(indel=is.na(str_extract(dir, "noindel"))) %>%
    mutate(replicate=factor(str_match(dir, "replicate_([0-9]*)")[,2])) %>%
    mutate(replicate=reorder(replicate, as.numeric(as.character(replicate)), FUN=max)) %>%
    mutate(wpi=factor(str_match(dir, "/([0-9]+)$")[,2])) %>%
    mutate(wpi=reorder(wpi, as.integer(as.character(wpi)), FUN=max)) %>%
    mutate(dir=file.path('..', 'build', dir))

tmp <- df %>% rowwise() %>% do({data.frame(.,distances.pwscore(.$dir), stringsAsFactors=FALSE)}) %>% ungroup()
save(tmp,file="distances.Rda")
load(file="distances.Rda")

png(file="distances.png", width=18, height=12, units = "in", res=300)
boxplot(tmp)
dev.off()

    
cat(sprintf("%.2f%% both PRANK models have the same score\n", nrow(tmp %>% filter(dist.pcodon == dist.pdna))/nrow(tmp) * 100))
cat(sprintf("%.2f%% consensus does better than both prank models\n", nrow(tmp %>% filter(dist.consensus > dist.pdna & dist.consensus > dist.pcodon ))/nrow(tmp) * 100))

# which sample has the greatest difference between consensus score and prank dna score?
# this would be the sample that does best on consensus compared to prank
tmp %>% ungroup() %>%
    filter(indel & as.numeric(as.character(wpi)) < 5000) %>%
    mutate(difference=dist.consensus - dist.pdna) %>%
    arrange(desc(difference)) %>%
    head(1) %>% write.alignment()

# which sample has the smallest difference between consensus score and prank dna score?
# this would be the sample that does better with prank dna compared to consesnsus.
tmp %>% ungroup() %>%
    filter(indel & as.numeric(as.character(wpi)) < 5000) %>%
    mutate(difference=dist.consensus - dist.pdna) %>%
    arrange(desc(difference)) %>%
    tail(1)  %>% write.alignment()

# which sample has the worst beast performance?
tmp %>% ungroup() %>%
    filter(!indel & as.numeric(as.character(wpi)) < 5000) %>%
    arrange(dist.beast) %>%
    head(1)  %>% write.alignment()




aln <- alignment(dir.biggest)
print(aln)

# who has the highest disparity between consensus score and prank_dna?

