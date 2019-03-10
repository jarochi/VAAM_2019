library(tidyverse)
library(biogram)

csgc <- biogram::read_fasta("data/CsgC.fasta")

ecoli_names <- grep("Escherichia", names(csgc))


for (i in ecoli_names) {
  print(csgc[i])
  CsgC_ecoli_only <- rbind(csgc[i])
}


CsgC_ecoli_only <- lapply(csgc, rbind(csgc[i]))


CsgC_ecoli_only <- lapply(ecoli_names, function(i){
  csgc[1]
})

biogram::write_fasta(CsgC_ecoli_only, "CsgC_ecoli_only.fasta")
biogram::write_encoding(CsgC_ecoli_only, "CsgC_ecoli_only.fasta")    # w opisie że czyta | Żle zapisuje obiekt po read_fasta !!
biogram::