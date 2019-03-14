library(dplyr)
library(ggplot2)
library(ggrepel)

if(Sys.info()[["nodename"]] %in% c("amyloid", "lori")) {
  seq_path_CsgA <- "/home/michal/Dropbox/dropbox-amylogram/PSI-blast/CsgA_muscle.fas"
  seq_path_CsgB <- "/home/michal/Dropbox/dropbox-amylogram/PSI-blast/CsgB_muscle.fas"
}else{
  seq_path_CsgA <- "data/CsgA-only-coli_mod.fas"
  seq_path_CsgB <- "data/CsgB_muscle.fas"
  # seq_path_CsgA <- "data/CsgA-only-coli_bez_ostatnich.fas"
}

aas <- c(seqinr::a()[-1], "-")

CsgA_regions <- list(R1 = 43L:65, 
                     R2 = 66L:87,
                     R3 = 88L:110,
                     R4 = 111L:132,
                     R5 = 133L:151)

seq_path <- seq_path_CsgA

all_lines <- readLines(seq_path)

prot_id <- cumsum(grepl("^>", all_lines))

all_prots <- split(all_lines, prot_id)

aln_dat <- lapply(all_prots, function(ith_prot) {
  strsplit(paste0(ith_prot[-1], collapse = ""), "")[[1]]
}) %>% 
  do.call(rbind, .)

real_positions <- cumsum(aln_dat[1, ] != "-")

only_regions <- lapply(CsgA_regions, function(ith_region){
  aln_dat[, real_positions %in% ith_region]})

res_in_aln <- list(R1 = c(S = 1, Q = 7, G = 9, G = 11, N = 12, A = 14, Q = 18),
                   R2 = c(S = 1, Q = 7, G = 9, G = 11, N = 12, A = 14, Q = 18),
                   R3 = c(S = 1, Q = 7, G = 9, G = 11, N = 12, A = 14, Q = 18),
                   R4 = c(S = 1, Q = 7, G = 9, G = 11, N = 12, A = 14, Q = 18),
                   R5 = c(S = 1, Q = 7, G = 9, G = 11, N = 12, A = 14, Q = 18))

CsgA_pattern <- lapply(1L:length(res_in_aln), function(ith_region) {
  only_regions[[ith_region]][, unname(res_in_aln[[ith_region]])]  
})


CsgA_motifs_per_species <- lapply(1L:nrow(CsgA_pattern[[1]]), function(ith_specimen) 
  sapply(1L:length(CsgA_pattern), function(ith_region) {
    CsgA_pattern[[ith_region]][ith_specimen, ]
  })
) 

mutants <- sapply(CsgA_motifs_per_species, function(ith_specimen)
  apply(ith_specimen, 1, function(i) length(unique(i))) %>% 
    sum
) != 7

mutation_df <- lapply(CsgA_motifs_per_species[mutants], function(single_mutant) {
  mutation_count <- apply(single_mutant, 1, function(i) length(unique(i)))

  data.frame(nonmutated = sum(mutation_count == 1), 
             frac_nonmutated = mean(mutation_count == 1),
             max_mutations = max(mutation_count),
             total_mutations = sum(mutation_count[mutation_count != 1]),
             pos_max_mutations = which.max(mutation_count))
}) %>% 
  bind_rows()

save(mutation_df, file = "mutation_df.RData")

group_by(mutation_df, pos_max_mutations, nonmutated, max_mutations) %>% 
  summarise(count = length(pos_max_mutations)) %>%  
  ungroup %>% 
  mutate(pos_max_mutations = factor(pos_max_mutations, 
                                    labels = paste0(paste(c("S", "Q", "G", "G", "N", "A", "Q"), res_in_aln[[1]], sep = " ("), ")")),
         nonmutated = paste0(nonmutated, " nonmutated resiues in the motif"),
         max_mutations = max_mutations - 1) %>% 
  ggplot(aes(x = factor(pos_max_mutations), y = count, color = factor(max_mutations), label = count)) +
  geom_point(size = 4) +
  geom_text_repel(color = "black", force = 10) +
  scale_color_discrete("Number of other amino acids\non the given position") + 
  guides(colour = guide_legend(nrow = 1)) +
  scale_x_discrete("Position in the motif with the largest number of mutations") +
  scale_y_continuous("Number of strains") +
  facet_wrap(~ nonmutated) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))

sapply(CsgA_motifs_per_species[mutants][mutation_df[["nonmutated"]] == 5], function(single_mutant) {
  paste0(which(apply(single_mutant, 1, function(i) length(unique(i))) != 1), collapse = "")
}) %>% table
