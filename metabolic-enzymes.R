## metabolic-enzymes.R by Rohan Maddamsetti
## In short, there looks like there is idiosyncratic purifying selection on
## core metabolic networks (superessential reactions and specialist enzymes)
## in the LTEE.

source("metagenomics-library.R")
library(UpSetR)
####################
## (METAGENOMIC) DATA PREPROCESSING

## get the lengths of all genes in REL606.
## This excludes genes in repetitive regions of the genome.
## See Section 4.3.1
## "Removing mutations in repetitive regions of the genome"
## in Ben Good's LTEE metagenomics paper for more details.
## This filtering is done in my python script printEcoliIDs.py.
## Do by running:
## python printEcoliIDs.py -i ../data/REL606.8-RM.gbk > ../results/REL606_IDs.csv.

REL606.genes <- read.csv('../results/REL606_IDs.csv',as.is=TRUE) %>%
    mutate(gene_length=strtoi(gene_length))
## It turns out that some gene names map to multiple genes.
duplicate.genes <- REL606.genes %>%
    group_by(Gene) %>%
    summarize(copies=length(unique(gene_length))) %>%
    filter(copies>1)
## There are 3 such cases: alr, bioD, maf (3 each for 6 loci total).
## filter them from this analysis.
REL606.genes <- REL606.genes %>%
    filter(!(Gene %in% duplicate.genes$Gene))

## Order nonmutator pops, then hypermutator pops by converting Population to
## factor type and setting the levels.
nonmutator.pops <- c("Ara-5", "Ara-6", "Ara+1", "Ara+2", "Ara+4", "Ara+5")
hypermutator.pops <- c("Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara+3", "Ara+6")
LTEE.pop.vec <- c(nonmutator.pops, hypermutator.pops)

## import genes in the table created by:
## conda activate ltee-metagenomix
## cd LTEE-metagenomic-repo
## python rohan-write-csv.py > ../results/LTEE-metagenome-mutations.csv
gene.mutation.data <- read.csv(
    '../results/LTEE-metagenome-mutations.csv',
    header=TRUE,as.is=TRUE) %>%
    mutate(Generation=t0/10000) %>%
    ## This for changing the ordering of populations in plots.
    mutate(Population = factor(
               Population,
               levels = LTEE.pop.vec)) %>%
    inner_join(REL606.genes) %>%
    filter(Gene!='intergenic')


## To simplify figure presentation, put the hypermutator results
## in the main text, and the nonmutator results in the supplement.
hypermutator.data <- gene.mutation.data %>%
    filter(Population %in% hypermutator.pops) %>%
    ## This for changing the ordering of populations in plots.
    mutate(Population = factor(Population,
                             levels = hypermutator.pops))

nonmutator.data <- gene.mutation.data %>%
    filter(Population %in% nonmutator.pops) %>%
    ## This for changing the ordering of populations in plots.
    mutate(Population=factor(Population,
                             levels = nonmutator.pops))

########################################################
## METABOLIC ENZYME DATASETS
########################################################

## BiGG Models E. coli core.
BiGG.core <- read.csv("../results/metabolic-enzymes/BiGG-Model-Ecoli_core.csv") %>%
    mutate(blattner = BiGG.ID) %>%
    left_join(REL606.genes)

## superessential metabolic reactions reported by Barve and Wagner (2012).
superessential.rxns.df <- read.csv("../results/metabolic-enzymes/Barve2012-S6-superessential.csv") %>%
    inner_join(REL606.genes)

## Look at specialist and generalist enzymes in Nam et al. (2012):
## Network context and selection in the evolution of enzyme specificity.

## 1157 genes.
Nam.df <- read.csv("../results/metabolic-enzymes/Nam2012_Database_S1.csv") %>%
    left_join(REL606.genes) %>% filter(!is.na(Gene))

specialist.enzymes <- Nam.df %>% filter(Class=="Spec.")
generalist.enzymes <- Nam.df %>% filter(Class=="Gen.")

########################################################
## KNOCKOUTS OF METABOLIC ENZYMES IN 50K GENOMES
########################################################

## Question: have any of the BiGG core, superessential,
## or specialist/generalist enzymes been knocked out in any of the
## 50K LTEE A clones?

make.list.of.strain.to.KOed.genes <- function(LTEE.KO.data) {
## make a list of strain to vectors of KO'ed genes.

    strip.and.split <- function(x) {
        ## remove [square brackets] and split on commas.
        str_replace_all(x, c("\\[" = "",  "\\]" = "")) %>%
            str_split(",")        
    }
    
    ## This is a helper function that does the string manipulation
    ## to extract and concatenate genes from rows of the dataframe
    ## to return a vector of genes.
    KO.df.to.KO.vec <- function(strain.df) {
        parsed.gene.list <- sapply(strain.df$gene_list, strip.and.split)
        reduce(parsed.gene.list, .f = c ) %>% unique()
    }
    
    LTEE.KO.data %>%
        ## for easy plotting, split on population.
        ## CRITICAL ASSUMPTION: one strain per population.
        split(.$population) %>%
        map(.f = KO.df.to.KO.vec)
}


LTEE.50K.A.clone.KO.muts <- read.csv(
    "../data/LTEE-264-genomes-SNP-nonsense-small-indel-MOB-large-deletions.csv") %>%
    select(-X) %>% ## drop this indexing column.
    ## filter out intergenic mutations.
    filter(!str_detect(gene_position, "intergenic")) %>%
    ## 50,000 generations only
    filter(time == 50000) %>%
    ## and A clone only.
    filter(clone == 'A')

LTEE.50K.A.clone.KO.metadata <- LTEE.50K.A.clone.KO.muts %>%
    select(population,time,strain,clone,mutator_status) %>%
    distinct() %>%
    mutate(Generation=time/10000) %>%
    ## This for changing the ordering of populations in plots.
    mutate(population=factor(population,levels = LTEE.pop.vec))


## make a data structure for genes affected by KO mutations in 50K A clones.
KOed.genes.in.LTEE.50K.A.clones <- LTEE.50K.A.clone.KO.muts %>%
    select(population, strain, clone, gene_list) %>%
    make.list.of.strain.to.KOed.genes()

## turn this data structure into a data.frame.
KOed.50K.A.clone.genes.df <- KOed.genes.in.LTEE.50K.A.clones %>%
    enframe() %>%
    unnest_longer(value) %>%
    rename(Population = name) %>%
    rename(Gene = value)
    

## KO'ed BiGG core genes
KOed.50K.BiGG.core.genes <- KOed.50K.A.clone.genes.df %>%
    filter(Gene %in% BiGG.core$Gene) %>%
    mutate(MetabolicClass = "BiGG_core")

## KO'ed superessential genes
KOed.50K.superessential.genes <- KOed.50K.A.clone.genes.df %>%
    filter(Gene %in% superessential.rxns.df$Gene) %>%
    mutate(MetabolicClass = "Superessential")

## KO'ed specialist and generalist enzymes in Nam et al. (2012):
## Network context and selection in the evolution of enzyme specificity.

KOed.50K.specialist.genes <- KOed.50K.A.clone.genes.df %>%
    filter(Gene %in% specialist.enzymes$Gene) %>%
    mutate(MetabolicClass = "Specialist")

KOed.50K.generalist.genes <- KOed.50K.A.clone.genes.df %>%
    filter(Gene %in% generalist.enzymes$Gene) %>%
    mutate(MetabolicClass = "Generalist")

## combine these tables, add metadata, and write to file, so that these
## gene knockouts can be used for dynamic FBA analysis with COMETS.
KOed.50K.metabolic.enzymes <- rbind(KOed.50K.BiGG.core.genes,
                                KOed.50K.superessential.genes,
                                KOed.50K.specialist.genes,
                                KOed.50K.generalist.genes) %>%
    inner_join(REL606.genes)

write.csv(KOed.50K.metabolic.enzymes,
          "../results/metabolic-enzymes/KOed-metabolic-enzymes-in-LTEE-50K-A-clones.csv",
          row.names = FALSE)

########################################################
## METABOLIC ENZYME STIMS ANALYSIS.
########################################################

## Run STIMS on BiGG Models E. coli core.
## Hypothesis: E. coli core metabolism
## is evolving under purifying selection in the LTEE.

BiGG.core <- read.csv("../results/metabolic-enzymes/BiGG-Model-Ecoli_core.csv") %>%
    mutate(blattner = BiGG.ID) %>%
    left_join(REL606.genes)

## plot just the hypermutator populations.
BiGG.core.hypermut.data <- hypermutator.data %>%
    filter(Gene %in% BiGG.core$Gene)

c.BiGG.core.hypermut <- calc.cumulative.muts(
    BiGG.core.hypermut.data,
    BiGG.core,
    manual.pop.levels.vec = hypermutator.pops)

Fig1A <- plot.base.layer(
    hypermutator.data,
    REL606.genes,
    subset.size = length(unique(BiGG.core$Gene)),
    manual.pop.levels.vec = hypermutator.pops,
    plot.rows = 2) %>%
    add.cumulative.mut.layer(c.BiGG.core.hypermut, my.color="black") +
    ggtitle("BiGG core metabolic enzymes")

## plot just the nonmutator populations.
BiGG.core.nonmut.data <- nonmutator.data %>%
    filter(Gene %in% BiGG.core$Gene)

c.BiGG.core.nonmut <- calc.cumulative.muts(
    BiGG.core.nonmut.data,
    BiGG.core,
    manual.pop.levels.vec = nonmutator.pops)

S1FigA <- plot.base.layer(
    nonmutator.data,
    REL606.genes,
    subset.size = length(unique(BiGG.core$Gene)),
    manual.pop.levels.vec = nonmutator.pops,
    plot.rows = 2) %>%
    add.cumulative.mut.layer(c.BiGG.core.nonmut, my.color="black") +
    ggtitle("BiGG core metabolic enzymes")
    
## calculate formal p-values.
BiGG.core.pvals <- calc.traj.pvals(gene.mutation.data, REL606.genes, unique(BiGG.core$Gene))

## results:
##  > BiGG.core.pvals
##  # A tibble: 12 x 3
##    Population count p.val
##    <fct>      <int> <dbl>
##  1 Ara-5       3118 0.312
##  2 Ara-6       2079 0.208
##  3 Ara+1       7836 0.784
##  4 Ara+2       4072 0.407
##  5 Ara+4       6300 0.63 
##  6 Ara+5       2438 0.244
##  7 Ara-1       9502 0.950
##  8 Ara-2       7084 0.708
##  9 Ara-3       1340 0.134
## 10 Ara-4       9342 0.934
## 11 Ara+3       9647 0.965
## 12 Ara+6       9967 0.997

#####################################################################################
## examine superessential metabolic reactions reported by Barve and Wagner (2012).

## plot just the hypermutator populations.
superessential.hypermut.data <- hypermutator.data %>%
    filter(Gene %in% superessential.rxns.df$Gene)

c.superessential.hypermut <- calc.cumulative.muts(
    superessential.hypermut.data,
    superessential.rxns.df,
    manual.pop.levels.vec = hypermutator.pops)

Fig1B <- plot.base.layer(
    hypermutator.data,
    REL606.genes,
    subset.size = length(unique(superessential.rxns.df$Gene)),
    manual.pop.levels.vec = hypermutator.pops,
    plot.rows = 2,
    my.color = "pink") %>%
    add.cumulative.mut.layer(c.superessential.hypermut, my.color="red") +
    ggtitle("Superessential metabolic enzymes")

ggsave("../results/metabolic-enzymes/Fig1B.pdf",Fig1B)

## plot just the nonmutator populations.
superessential.nonmut.data <- nonmutator.data %>%
    filter(Gene %in% superessential.rxns.df$Gene)

c.superessential.nonmut <- calc.cumulative.muts(
    superessential.nonmut.data,
    superessential.rxns.df,
    manual.pop.levels.vec = nonmutator.pops)

S1FigB <- plot.base.layer(
    nonmutator.data,
    REL606.genes,
    subset.size = length(unique(superessential.rxns.df$Gene)),
    manual.pop.levels.vec = nonmutator.pops,
    plot.rows = 2,
    my.color = "pink") %>%
    add.cumulative.mut.layer(c.superessential.nonmut, my.color="red") +
    ggtitle("Superessential metabolic enzymes")

## calculate formal p-values.
superessential.pvals <- calc.traj.pvals(gene.mutation.data, REL606.genes, unique(superessential.rxns.df$Gene))

## results:
## > superessential.pvals
## # A tibble: 12 x 3
##    Population count p.val
##    <fct>      <int> <dbl>
##  1 Ara-5       6279 0.628
##  2 Ara-6       9133 0.913
##  3 Ara+1       8249 0.825
##  4 Ara+2       6725 0.672
##  5 Ara+4       8841 0.884
##  6 Ara+5       1976 0.198
##  7 Ara-1       9999 1.00 
##  8 Ara-2       6319 0.632
##  9 Ara-3       5968 0.597
## 10 Ara-4       8560 0.856
## 11 Ara+3       5476 0.548
## 12 Ara+6       9902 0.990

################################################################
## Figure 1 combines the BiGG core and superessential gene results
## for hypermutators.
Fig1 <- plot_grid(Fig1A, Fig1B, labels=c('A','B'),nrow=2)
save_plot("../results/metabolic-enzymes/Fig1.pdf",Fig1, base_height=7,base_asp=1)

## Supplementary Figure S1 combines the BiGG and superessential gene
## results for nonmutators.
S1Fig <- plot_grid(S1FigA, S1FigB, labels=c('A','B'),nrow=2)
save_plot("../results/metabolic-enzymes/S1Fig.pdf",S1Fig, base_height=7,base_asp=1)

################################################################
## Look at specialist and generalist enzymes in Nam et al. (2012):
## Network context and selection in the evolution of enzyme specificity.

## plot just the hypermutator populations.
specialist.hypermut.data <- hypermutator.data %>%
    filter(Gene %in% specialist.enzymes$Gene)

c.specialists.hypermut <- calc.cumulative.muts(
    specialist.hypermut.data,
    specialist.enzymes,
    manual.pop.levels.vec = hypermutator.pops)

generalist.hypermut.data <- hypermutator.data %>%
    filter(Gene %in% generalist.enzymes$Gene)

c.generalists.hypermut <- calc.cumulative.muts(
    generalist.hypermut.data,
    generalist.enzymes,
    manual.pop.levels.vec = hypermutator.pops)

Fig2A <- plot.base.layer(
    hypermutator.data,
    REL606.genes,
    subset.size=length(unique(specialist.enzymes$Gene)),
    manual.pop.levels.vec = hypermutator.pops,
    plot.rows = 2,
    my.color = "plum1") %>%
    add.cumulative.mut.layer(c.specialists.hypermut, my.color="darkorchid4") +
    ggtitle("Specialist enzymes")

Fig2B <- plot.base.layer(
    hypermutator.data,
    REL606.genes,
    subset.size=length(unique(generalist.enzymes$Gene)),
    manual.pop.levels.vec = hypermutator.pops,
    plot.rows = 2,
    my.color = "darkolivegreen1") %>%
    add.cumulative.mut.layer(c.generalists.hypermut, my.color="springgreen4") +
    ggtitle("Generalist enzymes")

Fig2 <- plot_grid(Fig2A, Fig2B, labels=c('A','B'),nrow=2)
save_plot("../results/metabolic-enzymes/Fig2.pdf",Fig2, base_height=7,base_asp=1)

## plot just the nonmutator populations.
specialist.nonmut.data <- nonmutator.data %>%
    filter(Gene %in% specialist.enzymes$Gene)

c.specialists.nonmut <- calc.cumulative.muts(
    specialist.nonmut.data,
    specialist.enzymes,
    manual.pop.levels.vec = nonmutator.pops)

generalist.nonmut.data <- nonmutator.data %>%
    filter(Gene %in% generalist.enzymes$Gene)

c.generalists.nonmut <- calc.cumulative.muts(
    generalist.nonmut.data,
    generalist.enzymes,
    manual.pop.levels.vec = nonmutator.pops)

S2FigA <- plot.base.layer(
    nonmutator.data,
    REL606.genes,
    subset.size=length(unique(specialist.enzymes$Gene)),
    manual.pop.levels.vec = nonmutator.pops,
    plot.rows = 2,
    my.color = "plum1") %>%
    add.cumulative.mut.layer(c.specialists.nonmut, my.color="darkorchid4") +
    ggtitle("Specialist enzymes")

S2FigB <- plot.base.layer(
    nonmutator.data,
    REL606.genes,
    subset.size=length(unique(generalist.enzymes$Gene)),
    manual.pop.levels.vec = nonmutator.pops,
    plot.rows = 2,
    my.color = "darkolivegreen1") %>%
    add.cumulative.mut.layer(c.generalists.nonmut, my.color="springgreen4") +
    ggtitle("Generalist enzymes")

S2Fig <- plot_grid(S2FigA, S2FigB, labels=c('A','B'),nrow=2)
save_plot("../results/metabolic-enzymes/S2Fig.pdf",S2Fig, base_height=7,base_asp=1)

## calculate formal p-values.
specialist.pvals <- calc.traj.pvals(gene.mutation.data, REL606.genes, unique(specialist.enzymes$Gene))
## results:
## > specialist.pvals
## # A tibble: 12 x 3
##    Population count p.val
##    <fct>      <int> <dbl>
##  1 Ara-5       2349 0.235
##  2 Ara-6       8521 0.852
##  3 Ara+1       9900 0.99 
##  4 Ara+2       9649 0.965
##  5 Ara+4       8247 0.825
##  6 Ara+5       9379 0.938
##  7 Ara-1       5472 0.547
##  8 Ara-2       4766 0.477
##  9 Ara-3       9266 0.927
## 10 Ara-4       9506 0.951
## 11 Ara+3       7796 0.780
## 12 Ara+6      10000 1    

generalist.pvals <- calc.traj.pvals(gene.mutation.data, REL606.genes, unique(generalist.enzymes$Gene))
## results:
## > generalist.pvals
## # A tibble: 12 x 3
##    Population count  p.val
##    <fct>      <int>  <dbl>
##  1 Ara-5       7166 0.717 
##  2 Ara-6       4566 0.457 
##  3 Ara+1       6108 0.611 
##  4 Ara+2       4316 0.432 
##  5 Ara+4       4127 0.413 
##  6 Ara+5        277 0.0277
##  7 Ara-1       9993 0.999 
##  8 Ara-2       7243 0.724 
##  9 Ara-3       3413 0.341 
## 10 Ara-4       6764 0.676 
## 11 Ara+3       7066 0.707 
## 12 Ara+6       8988 0.899 

################
## Make an UpSet plot to examine set overlap.

## Cite: Jake R Conway, Alexander Lex, Nils Gehlenborg UpSetR:
## An R Package for the Visualization of Intersecting Sets and their
## Properties doi: https://doi.org/10.1093/bioinformatics/btx364

## The as.numeric() calls turn TRUE/FALSE to 1/0.
REL606.UpSet.data <- REL606.genes %>%
    mutate(`BiGG Core` = Gene %in% BiGG.core$Gene) %>%
    mutate(`BiGG Core` = as.numeric(`BiGG Core`)) %>%
    mutate(Superessential = Gene %in% superessential.rxns.df$Gene) %>%
    mutate(Superessential = as.numeric(Superessential)) %>%
    mutate(Specialist = Gene %in% specialist.enzymes$Gene) %>%
    mutate(Specialist = as.numeric(Specialist)) %>%
    mutate(Generalist = Gene %in% generalist.enzymes$Gene) %>%
    mutate(Generalist = as.numeric(Generalist))

pdf("../results/metabolic-enzymes/Fig3.pdf",onefile=FALSE)
upset(REL606.UpSet.data, sets = c("BiGG Core", "Superessential", "Specialist", "Generalist"), mb.ratio = c(0.7, 0.3), order.by = "freq", text.scale=2)
dev.off()

write.csv(REL606.UpSet.data,file="../results/metabolic-enzymes/TableS1.csv")
