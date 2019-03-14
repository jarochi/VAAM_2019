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
  
  
  write_fasta <- function(seq, file, nchar = 80) {
    char_vec <- unlist(lapply(1L:length(seq), function(ith_id) 
      c(paste0(">", names(seq[ith_id])), 
        lapply(split(seq[[ith_id]], floor(seq_along(seq[[ith_id]])/nchar)), paste0, collapse = ""))))
    writeLines(text = char_vec, con = file)
  }



zz <- csgc %>% unlist




char_vec <- unlist(lapply(1L:length(csgc), function(ith_id){
  c(paste0(">", names(csgc[ith_id])), 
    lapply(split(csgc[[ith_id]], floor(seq_along(csgc[[ith_id]])/nchar)), paste0, collapse = ""))
}))
