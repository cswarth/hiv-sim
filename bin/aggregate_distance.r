#!/usr/bin/env Rscript

# Aggregate individual distance files into one distance file.
#
# This script expects to be invoked with a bunch of paths to individual
# CSV files.
#
# e.g. aggregate_distance.r runs/300/100/500/relaxed/distance.csv runs/500/300/1500/relaxed/distance.csv
#
# except that one should expect dozens of filenames on the
# commandline.
#
# The directory paths to the distance.csv files are significant!
# They must have at least 4 components in the path and each component
# has a particular interpretation.  For example, the first numeric
# component is the time since infection for patient-1, the second is
# time of transmission, and the third is time since infection for
# patient-2


suppressPackageStartupMessages(library("optparse"))

main <- function(args) {
    option_list <- list()

    # get command line options, if help option encountered print help and exit,
    # otherwise if options not found on command line then set defaults,
    parser <- OptionParser(usage = "%prog [options] files [...]", option_list=option_list)

    arguments <- parse_args(parser, args=args, positional_arguments = c(1, Inf))
    opt <- arguments$options
    files <- arguments$args

    ll <- lapply(files, read.delim, comment.char = "#")
    df <- do.call(rbind, ll)

    # Add portions of the file path as columns in the table.
    #
    # The last four path components provide information about the
    # patient sampling dates and which molecular clock model was used
    # by BEAST.  We have to count components from the right because an
    # unknown number of path elements may precede the filename.
    #
    m <- do.call(rbind, lapply(strsplit(files, "/"), function(s) head(tail(s, 5), 4)))
    m <- data.frame(m, stringsAsFactors=F)

    # Tinf.p1 - Time since infection, patient 1  (relative to origin)
    # Ttrans - Transmission time (relative to origin)
    # Tinf.p2 - Time since infection, patien2 (relative to transmission)
    names(m) <- c("Tinf.p1", "Ttrans", "Tinf.p2", "clock")
    m <- within(m, {
        Tinf.p1 <- as.numeric(Tinf.p1)
        Ttrans <- as.numeric(Ttrans)
        Tinf.p2 <- as.numeric(Tinf.p2)
    })
    

    df <- cbind(m, df)

    cat("# use read.csv(file, comment='#') to read this file into R\n")
    cat("#\n")
    cat(paste0("# date: ", Sys.time(), "\n"))
    cat(paste0("# working directory: ", getwd(), "\n"))
    cat(paste0("# command line: ", paste(commandArgs(trailingOnly = FALSE), collapse=" "), "\n"))

    write.csv(df, row.names = FALSE)

}

main(commandArgs(trailingOnly = TRUE))
