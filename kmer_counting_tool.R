#!/usr/bin/env Rscript
suppressMessages(library("Biostrings"))
suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
options(warn=-1)

# example command line input to run script with kmer length=4:
# Rscript kmer_counting_tool.R 'takehome/challenge1/experiment1.fasta' 4 --output_file 'output-kmer-4.tsv'
# Ex2: Rscript kmer_counting_tool.R 'nonstandard_nucs.fasta' 4 --output_file 'output-kmer-4-NS.tsv'
# Change the permissions for this script locally to make it executable:
# chmod +x kmer_counting_tool.R

#define arguments for argparser
p <- arg_parser("Input: a fasta file and kmer length. Output: A tab separated file with counts per kmer")
p <- add_argument(p, "input_file", help = "fasta input file")
p <- add_argument(p, "kmer_length", help = "Choose kmer length")

p <- add_argument(p, "--output_file", help = "returns a tsv file with kmer counts sorted by abundance")

argv <- parse_args(p)

#use Biostrings to import fasta file
seq1 = readDNAStringSet(argv$input_file) %>% as.list()
seq1

#Function:
#convert_to_tibble()
convert_to_tibble <- function(s) {
  sequence = paste(s)
  df <- tibble(sequence)
}

fasta_df <- map_df(seq1, convert_to_tibble)

paste0("The input sequence is ", nchar(fasta_df$sequence), " bp's in length")

#Function:
#kmer_length()
kmer_len <- function(s, length) {
  L <- str_length(s)
  str_sub(s, start=seq(1L,L-length+1,length), end=seq(length,L,length))
}

paste0("The kmer length is ", argv$kmer_length)

#Message for out-of-bound kmers (longer or shorter than the input sequence)
if (as.integer(argv$kmer_length) <= 0) {
  print("The kmer length must be a postive integer greater than 0. Choose a longer kmer.")
} else
  if (as.integer(argv$kmer_length) > nchar(fasta_df$sequence)) {#the length of the FASTA record
    print("The kmer length is greater than the length of the input sequence. Choose a shorter kmer.")
  } else {
    print("Good choice! The kmer length is within the sequence range.")
  }

#non_stand_nucs
non_stand_nucs <- paste0("[", paste(letters[-c(1,3,7,20)], collapse = ""), paste(toupper(letters[-c(1,3,7,20)]), collapse = ""), "]") #bdefhijklmnopqrsuvwxyzBDEFHIJKLMNOPQRSUVWXYZ
paste0("This script checks for the following non-standard nucleotides ", non_stand_nucs)

#Group and count each kmer, sort by abundance
a <- fasta_df$sequence %>%
  kmer_len(length = as.numeric(argv$kmer_length)) %>%
  tibble() %>%
  rename(kmer = ".") %>%
  count(kmer, name = 'length', sort = T) %>%
  mutate(standard_nucs = str_detect(kmer, non_stand_nucs, negate = T))

a %>%
  select(kmer, length)

#Check for nonstandard nucs (aka, not 'actgACTG')
#a$standard_nucs[30] = FALSE
if (any(FALSE) %in% a$standard_nucs) {
  print("Warning: Your fasta sequence includes nonstandard nucleic acids")
} else {
  print("Great news: Your fasta sequence contains only A, C, G, and T")
}

write_delim(a, argv$output_file, delim = "\t")
