#!/usr/bin/env Rscript 
# expect to be in ~/src/matsen/hiv-sim/sim.rv217/population/

source('../../R/utils.r')

founder <- founder.sequence(path=file.path('..', '..', '..', 'templates', 'HIV1C2C3.fasta'))



# recursively list files matching 'sample.fa'
files <- list.files(path='../build', pattern="sample\\.fa$", recursive=TRUE)

# The file names look like,
# 	"1.0E-5/frequency_noindel/replicate_0/5000/sample.fa"
# break down the path components into seperate columns so the data
# can be filtered appropriately.

df <- as_data_frame(list(dir=dirname(files))) %>%
    separate(dir, c('rate', 'fitness', 'indel', 'population', 'replicate', 'wpi', 'nseq'), "/", remove=FALSE) %>%
    mutate(rate=reorder(rate, as.numeric(as.character(rate)), FUN=max)) %>%
    mutate(fitness=factor(fitness)) %>%
    mutate(indel=is.na(str_extract(indel, "noindel"))) %>%
    mutate(replicate=factor(str_match(replicate, "rep=([0-9]*)")[,2])) %>%
    mutate(wpi=factor(str_match(wpi, "gen=([0-9]*)")[,2])) %>%
    mutate(population=factor(str_match(population, "pop=([0-9]*)")[,2])) %>%
    mutate(nseq=as.integer(str_match(nseq, "N=([0-9]*)")[,2])) %>%

    mutate(replicate=reorder(replicate, as.numeric(as.character(replicate)), FUN=max)) %>%
    mutate(wpi=reorder(wpi, as.integer(as.character(wpi)), FUN=max)) %>%
    mutate(population=reorder(population, as.integer(as.character(population)), FUN=max)) %>%
    mutate(dir=file.path('..', 'build', dir))

get_sequences <- function(dir) {
    c(strict=beast.sequence(file.path(dir,"strict")),
      unguided=prank.sequence(dir, 'unguided'),
      consensus=consensus.sequence(dir))
}


cat(sprintf('Number of samples to process: %d\n', nrow(df)))
tmp <- df %>% rowwise() %>% do({data.frame(.,distances.pwscore(.$dir), stringsAsFactors=FALSE)}) %>% ungroup()

save(tmp,file="distances.Rda")

