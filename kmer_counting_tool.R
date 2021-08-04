#!/usr/bin/env Rscript
suppressMessages(library("Biostrings"))
suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
options(warn=-1)

# Make sure to change the permissions for this script in the command line to make it executable:
# chmod +x kmer_counting_tool.R
# example command line input to run script with k-mer length=4
# Ex1:
# Rscript kmer_counting_tool.R 'takehome/challenge1/experiment1.fasta' 4 --output_file 'output-files/exp1-kmer-4.tsv'
# Ex2:
# Rscript kmer_counting_tool.R 'takehome/challenge1/nonstandard_nucs.fasta' 4 --output_file 'output-files/exp1-kmer-4-nonstandard-nucs.tsv'


#define arguments for argparser
p <- arg_parser("Input: a FASTA file and k-mer length. Output: A tab separated file with counts per k-mer")
p <- add_argument(p, "input_file", help = "FASTA input file")
p <- add_argument(p, "kmer_length", help = "Choose k-mer length")

p <- add_argument(p, "--output_file", help = "returns a tsv file with k-mer counts arranged in descending order")

argv <- parse_args(p)

#use Biostrings to import FASTA file
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

paste0("The k-mer length is ", argv$kmer_length)

#Message for out-of-bound k-mers (longer or shorter than the input sequence)
if (as.integer(argv$kmer_length) <= 0) {
  print("The k-mer length must be a postive integer greater than 0. Choose a longer k-mer.")
} else
  if (as.integer(argv$kmer_length) > nchar(fasta_df$sequence)) {#the length of the FASTA record
    print("The k-mer length is greater than the length of the input sequence. Choose a shorter k-mer.")
  } else {
    print("Good choice! The k-mer length is within the sequence range.")
  }

#Nonstandard nucleotide list
non_stand_nucs <- paste0("[", paste(letters[-c(1,3,7,20)], collapse = ""), paste(toupper(letters[-c(1,3,7,20)]), collapse = ""), "]") #bdefhijklmnopqrsuvwxyzBDEFHIJKLMNOPQRSUVWXYZ
paste0("This script checks for the following non-standard nucleotides ", non_stand_nucs)

#Group and count each k-mer, sort by abundance
a <- fasta_df$sequence %>%
  kmer_len(length = as.numeric(argv$kmer_length)) %>%
  tibble() %>%
  rename(kmer = ".") %>%
  count(kmer, name = 'length', sort = T) %>%
  mutate(standard_nucs = str_detect(kmer, non_stand_nucs, negate = T))

a %>%
  select(kmer, length)

#Check for nonstandard nucleotides (aka, not 'actgACTG')
if (any(FALSE) %in% a$standard_nucs) {
  print("Warning: Your FASTA sequence includes nonstandard nucleic acids")
} else {
  print("Great news: Your FASTA sequence contains only A, C, G, and T")
}

write_delim(a, argv$output_file, delim = "\t")
