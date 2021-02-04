## Generate freebayes params ##
# make bed for parallelising freebayes

#' This script reads a .fai index and a set of chroms and outputs a series of BED files
#' which split the genome into separate regions.

library(dplyr)
library(data.table)
library(glue)

# read inputs
chroms = snakemake@params[['chroms']]
chunks = snakemake@params[['chunks']]
fai = fread(snakemake@input[['index']])

# select chroms we want, and start, end columns
fai = fai[fai$V1 %in% chroms, c(1,2)]

# for each chromsome
for (chrom in chroms){
   #subset index to desired chrom
   f = fai[fai$V1 == chrom]
   #get sequence of n chunks from 0 to length of chrom
   bedseq = round(seq(0, f$V2, length.out = chunks))
   
   #for each chunk
   for (i in 1:(chunks-1)){
      #write bed file, one for each chunk/interval, which will be passed as input to freebayes
      row = c(chrom, bedseq[i], bedseq[i+1])
      data.frame(row) %>% t() %>% fwrite(., glue("resources/regions/genome.{chrom}.region.{i}.bed"), sep="\t", col.names = FALSE)
   }
}

sessionInfo()