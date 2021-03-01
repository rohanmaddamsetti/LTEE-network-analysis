## make-MAE-KO-mutation-table.R by Rohan Maddamsetti.
## genes that are inactivated are now annotated by breseq 0.35.5.

library(tidyverse)

filter.for.KO.muts <- function(df) {
    df %>% rename(strain = title) %>% select(strain, mutation_category, genes_inactivated, locus_tags_inactivated) %>% filter(str_length(locus_tags_inactivated) > 0)
}

JEB807.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB807.tsv", sep = "\t") %>% filter.for.KO.muts
JEB808.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB808.tsv", sep = "\t") %>% filter.for.KO.muts
JEB809.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB809.tsv", sep = "\t") %>% filter.for.KO.muts
JEB810.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB810.tsv", sep = "\t") %>% filter.for.KO.muts
JEB811.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB811.tsv", sep = "\t") %>% filter.for.KO.muts
JEB812.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB812.tsv", sep = "\t") %>% filter.for.KO.muts
JEB813.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB813.tsv", sep = "\t") %>% filter.for.KO.muts
JEB814.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB814.tsv", sep = "\t") %>% filter.for.KO.muts
JEB815.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB815.tsv", sep = "\t") %>% filter.for.KO.muts
JEB816.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB816.tsv", sep = "\t") %>% filter.for.KO.muts
JEB817.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB817.tsv", sep = "\t") %>% filter.for.KO.muts
JEB818.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB818.tsv", sep = "\t") %>% filter.for.KO.muts
JEB819.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB819.tsv", sep = "\t") %>% filter.for.KO.muts
JEB820.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB820.tsv", sep = "\t") %>% filter.for.KO.muts
JEB821.df <- read.csv("../results/resilience/annotated-MAE-clone-curated/JEB821.tsv", sep = "\t") %>% filter.for.KO.muts

MAE.KO.mutations.df <- rbind(JEB807.df, JEB808.df, JEB809.df, JEB810.df, JEB811.df,
                             JEB812.df, JEB813.df, JEB814.df, JEB815.df, JEB816.df,
                             JEB817.df, JEB818.df, JEB819.df, JEB820.df, JEB821.df)
write.csv(MAE.KO.mutations.df,
          file="../results/resilience/annotated-MAE-clone-curated/all-MAE-KO-mutations.csv")
                             
