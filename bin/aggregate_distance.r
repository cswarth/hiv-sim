#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

main <- function(args) {
    option_list <- list(
        ## make_option("--normalization", default="area",
        ##             help = "Function to normalize data, \"quantile\" or \"area\" [default \"%default\"]"),
        ## make_option("--genefile", default="",
        ##             help="file containing name genes to monitor. [default \"%default\"]"),
        ## make_option("--sample", default="0",
        ##             help="number of genes to sample from genefile. [default \"%default\"]"),
        ## make_option("--genes", default="",
        ##             help="gene or genes to monitor. [default \"%default\"]"),
        ## make_option("--seed", default="1234",
        ##             help="random seed. [default \"%default\"]"),
        ## make_option("--restart", default="1",
        ##             help="restart batch. [default \"%default\"]"),
        ## make_option("--sqlite", default="ProfileCache.sqlite",
        ##             help="sqlite database. [default \"%default\"]")
        )

    # get command line options, if help option encountered print help and exit,
    # otherwise if options not found on command line then set defaults,
    parser <- OptionParser(usage = "%prog [options] files [...]", option_list=option_list)

    arguments <- parse_args(parser, args=args, positional_arguments = c(1, Inf))
    opt <- arguments$options
    files <- arguments$args

    df <- t(sapply(files, read.delim, comment.char = "#"))
    names <- do.call(rbind, strsplit(rownames(df), "/"))

    # Add portions of the file path as columns in the table
    df <- cbind(left=as.numeric(names[,2]), right=as.numeric(names[,3]),  clock=as.character(names[,4]), df )

    cat("# use read.csv(file, comment='#') to read this file into R\n")
    cat("#\n")
    cat(paste0("# date: ", Sys.time(), "\n"))
    cat(paste0("# working directory: ", getwd(), "\n"))
    cat(paste0("# command line: ", paste(commandArgs(trailingOnly = TRUE), collapse=" "), "\n"))

    write.csv(df, row.names = FALSE)

}

main(commandArgs(trailingOnly = TRUE))
