#!/usr/bin/env Rscript
library("Biostrings")
library(tidyverse)
library(argparser)

# example command line input to run script:
# Rscript kmer_argparser.R 'takehome/challenge1/experiment1.fasta' 3
# Make sure to change the permissions on the script so that it's excutable:
# chmod +x kmer_argparser.R

p <- arg_parser("input a fasta file and a kmer length")
p <- add_argument(p, "input_file", help = "fasta input file")
p <- add_argument(p, "kmer_length", help = "select the length of kmer")

p <- add_argument(p, "--output_file", help = "returns a tsv file with kmer counts sorted by abundance")
argv <- parse_args(p)

seq1 = readDNAStringSet(argv$input_file) %>% as.list()
#ex1 = readDNAStringSet('takehome/challenge1/experiment1.fasta') %>% as.list() %>% tibble()
seq1

convert_to_df <- function(s) {
sequence = paste(s)
df <- tibble(sequence)
}

fasta_df <- map_df(seq1, convert_to_df)

# kmer_len = function(string, length) {
# l <- str_length(string)
# str_sub(string, start = seq(1L, l - length + 1, length), end = seq(length , l, length))
# }

kmer_len <- function(s,width) { # does not return list
  L <- str_length(s)
  str_sub(s, start=seq(1L,L-width+1,width), end=seq(width,L,width))
}

paste0("The kmer's length is ", argv$kmer_length)

fasta_df$sequence %>%
  kmer_len(width = as.numeric(argv$kmer_length)) %>%
  tibble() %>%
  rename(kmer = ".") %>%
  count(kmer, sort = T)

write_tsv(argv$output_file)

