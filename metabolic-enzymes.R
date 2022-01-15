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

######################################################################
## Figure 1: Make an UpSet plot to examine set overlap.

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

pdf("../results/metabolic-enzymes/Fig1.pdf",onefile=FALSE)
upset(REL606.UpSet.data,
      sets = c("BiGG Core", "Superessential",
               "Specialist", "Generalist"),
      mb.ratio = c(0.7, 0.3), order.by = "freq", text.scale=2)
dev.off()

write.csv(REL606.UpSet.data,file="../results/metabolic-enzymes/TableS1.csv")

########################################################
## KNOCKOUTS OF METABOLIC ENZYMES IN 50K LTEE GENOMES
########################################################

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

## We need this to remove "KO'ed genes" that were deleted after an amplification,
## preserving overall copy number.
AMPed.50K.A.clone.genes.df <- read.csv(
    "../data/LTEE-264-genomes-large-amplifications.csv") %>%
    select(-X) %>% ## drop this indexing column.
    ## 50,000 generations only
    filter(time == 50000) %>%
    ## and A clone only.
    filter(clone == 'A') %>%
    ## make a dataframe for genes affected by AMP mutations in 50K A clones.
    select(population, strain, clone, gene_list) %>%
    ## make a data structure for genes affected by AMP mutations in 50K A clones.
    ## This function should still work even though its a misnomer here.
    make.list.of.strain.to.KOed.genes() %>%
    ## turn this data structure into a data.frame.
    enframe() %>%
    unnest_longer(value) %>%
    rename(Population = name) %>%
    rename(Gene = value) %>%
    ## add metadata.
    inner_join(REL606.genes)

## get candidate KO mutations, process, and filter for those in amplifications.
KOed.50K.A.clone.genes.df <- read.csv(
    "../data/LTEE-264-genomes-SNP-nonsense-small-indel-MOB-large-deletions.csv") %>%
    select(-X) %>% ## drop this indexing column.
    ## filter out intergenic mutations.
    filter(!str_detect(gene_position, "intergenic")) %>%
    ## 50,000 generations only
    filter(time == 50000) %>%
    ## and A clone only.
    filter(clone == 'A') %>%
## make a dataframe for genes affected by KO mutations in 50K A clones.
    select(population, strain, clone, gene_list) %>%
    ## make a data structure for genes affected by KO mutations in 50K A clones.
    make.list.of.strain.to.KOed.genes() %>%
    ## turn this data structure into a data.frame.
    enframe() %>%
    unnest_longer(value) %>%
    rename(Population = name) %>%
    rename(Gene = value) %>%
    ## add metadata.
    inner_join(REL606.genes) %>%
    ## ABSOLUTELY CRITICAL: remove all "KO" mutations that were previously duplicated.
    anti_join(AMPed.50K.A.clone.genes.df)


## Use the Favate et al. (2021) Riboseq/RNAseq data to cross-check.
## Import separately for cross-checking in the end.
favate.data <- read.csv("../data/Favate2021_table_s1_read_counts.csv") %>%
    ## ignore REL606 and REL607 data.
    filter(!(line %in% c("REL606", "REL607"))) %>%
    group_by(line, target_id, eff_length, length) %>%
    summarize(total_est_counts = sum(est_counts)) %>%
    rename(Population = line) %>%
    rename(locus_tag = target_id)

## IMPORTANT: Ara+6 is not represented in these data, due to contamination.
favate.no.expression.df <- favate.data %>%
    ## require no transcription or translation in any replicate.
    filter(total_est_counts == 0) %>%
    as_tibble() %>%
    ## I'm not sure what the ERC_00XXX IDs are. Ignore those and the tRNAs for now.
    filter(str_detect(locus_tag,"^ECB")) %>%
    inner_join(REL606.genes) %>%
    ## remove unnecessary columns before joining to the KO dataset from genomics.
    select(-eff_length, -length, -total_est_counts)

## combine the genomics, RNAseq, and Riboseq data to generate the table of genes
## to remove from the 50K metabolic networks.
## There must be no expression, and use the KOs to get Ara+6.
inactive.50K.A.clone.genes.df <- left_join(
    favate.no.expression.df, KOed.50K.A.clone.genes.df)

## write to file, so that these gene knockouts can be used for FBA.
write.csv(inactive.50K.A.clone.genes.df,
          "../results/metabolic-enzymes/inactive-genes-in-LTEE-50K-A-clones.csv",
          row.names = FALSE)

## This is quite interesting: several genes with premature stops are still
## being expressed, some at quite high levels. This could represent the
## evolution of genes by removing extraneous C-terminal domains.
## worth exploring in future work.
expressed.putative.KOs <- anti_join(
    KOed.50K.A.clone.genes.df, favate.no.expression.df) %>%
    filter(Population != "Ara+6") %>%
    left_join(favate.data) %>%
    select(-eff_length, -length) %>%
    ## reorder the columns
    select(Population, Gene, total_est_counts, product, locus_tag,
           blattner, gene_length, start, end, strand)


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

Fig2A <- plot.base.layer(
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

Fig2B <- plot.base.layer(
    hypermutator.data,
    REL606.genes,
    subset.size = length(unique(superessential.rxns.df$Gene)),
    manual.pop.levels.vec = hypermutator.pops,
    plot.rows = 2,
    my.color = "pink") %>%
    add.cumulative.mut.layer(c.superessential.hypermut, my.color="red") +
    ggtitle("Superessential metabolic enzymes")

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
## Figure 2 combines the BiGG core and superessential gene results
## for hypermutators.
Fig2 <- plot_grid(Fig2A, Fig2B, labels=c('A','B'),nrow=2)
save_plot("../results/metabolic-enzymes/Fig1.pdf",Fig2, base_height=7,base_asp=1)

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

Fig3A <- plot.base.layer(
    hypermutator.data,
    REL606.genes,
    subset.size=length(unique(specialist.enzymes$Gene)),
    manual.pop.levels.vec = hypermutator.pops,
    plot.rows = 2,
    my.color = "plum1") %>%
    add.cumulative.mut.layer(c.specialists.hypermut, my.color="darkorchid4") +
    ggtitle("Specialist enzymes")

Fig3B <- plot.base.layer(
    hypermutator.data,
    REL606.genes,
    subset.size=length(unique(generalist.enzymes$Gene)),
    manual.pop.levels.vec = hypermutator.pops,
    plot.rows = 2,
    my.color = "darkolivegreen1") %>%
    add.cumulative.mut.layer(c.generalists.hypermut, my.color="springgreen4") +
    ggtitle("Generalist enzymes")

Fig3 <- plot_grid(Fig3A, Fig3B, labels=c('A','B'),nrow=2)
save_plot("../results/metabolic-enzymes/Fig3.pdf",Fig3, base_height=7,base_asp=1)

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

###############################################################################
## Figure 4. Genome Jenga analysis. Panels A and B are illustrations
## of the Jenga Hypothesis and the Genome Jenga algorithm, respectively.

minimal.genome.data <- read.csv(
    "../results/metabolic-enzymes/jenga_minimal_genomes.csv")
essential.genes.data <- read.csv(
    "../results/metabolic-enzymes/jenga_essential_genes.csv")
minimal.rxns.data <- read.csv(
    "../results/metabolic-enzymes/jenga_minimal_reactions.csv")


##### Figure 4, panels C, D, E.
minimal.genome.count.df <- minimal.genome.data %>%
    group_by(locus_tag) %>%
    summarize(Count = n()) %>%
    arrange(desc(Count)) %>%
    mutate(locus_tag = factor(locus_tag,levels=.$locus_tag)) %>%
    mutate(rank = row_number(locus_tag))

essential.gene.count.df <- essential.genes.data %>%
    group_by(locus_tag) %>%
    summarize(Count = n()) %>%
    arrange(desc(Count)) %>%
    mutate(locus_tag = factor(locus_tag,levels=.$locus_tag)) %>%
    mutate(rank = row_number(locus_tag))

minimal.rxn.count.df <- minimal.rxns.data %>%
    group_by(Reaction) %>%
    summarize(Count = n()) %>%
    arrange(desc(Count)) %>%
    mutate(Reaction = factor(Reaction,levels=.$Reaction)) %>%
    mutate(rank = row_number(Reaction))

Fig4C <- minimal.genome.count.df %>%
    ggplot(aes(x = rank, y = Count)) +
    geom_point(size=0.5) +
    xlab("Genes") +
    ylab("Number of genomes") +
    theme_classic()

Fig4D <- essential.gene.count.df %>% 
    ggplot(aes(x = rank, y = Count)) +
    geom_point(size=0.5) +
    geom_vline( xintercept = 202, color="red", linetype="dashed") +
    xlim(0,405) +
    xlab("Essential genes") +
    ylab("Number of genomes") +
    theme_classic()

Fig4E <- minimal.rxn.count.df %>%
    ggplot(aes(x = rank, y = Count)) +
    geom_point(size=0.5) +
    xlab("Metabolic reactions") +
    ylab("Number of genomes") +
    theme_classic()

Fig4CDE <- plot_grid(Fig4C, Fig4D, Fig4E, labels=c('C','D','E'),nrow=1)

##### Figure 4, panels F, G, H.

minimal.genome.sizes.df <- minimal.genome.data %>%
    group_by(Replicate) %>%
    summarize(NumGenes = n()) %>%
    arrange(desc(NumGenes))

essentialome.sizes.df <- essential.genes.data %>%
    group_by(Replicate) %>%
    summarize(NumEssentialGenes = n()) %>%
    arrange(desc(NumEssentialGenes))

minimal.rxn.network.sizes.df <- minimal.rxns.data %>%
    group_by(Replicate) %>%
    summarize(NumRxns = n()) %>%
    arrange(desc(NumRxns))

## Fig 4F. Minimal genomes have variable numbers of genes.
Fig4F <- minimal.genome.sizes.df %>%
    ggplot(aes(x = NumGenes)) +
    geom_histogram(binwidth=1) +
    theme_classic() +
    ylab("Count") +
    xlab("Genes per genome")

## Fig 4G. Evolved genomes are more fragile, and the number of essential
## genes varies across the minimal genomes.
## IMPORTANT RESULT! The predicted essentialomes are
## MUCH larger than the essentialome of REL606: 202 genes.
Fig4G <- essentialome.sizes.df %>%
    ggplot(aes(x = NumEssentialGenes)) +
    geom_histogram(binwidth=1) +
    theme_classic() +
    ylab("Count") +
    xlab("Essential genes per genome") +
    geom_vline( xintercept=202, color="red", linetype="dashed") +
    geom_label(x = 242, y = 50, size = 3,
               label="202 essential genes\nin the ancestral\nmetabolic network",
               label.size = 0, color = "red")

## Fig 4H. Idiosyncratic variation in reaction network size.
## This distribution has multiple modes.
Fig4H <- minimal.rxn.network.sizes.df %>%
    ggplot(aes(x = NumRxns)) +
    geom_histogram(binwidth=1) +
    theme_classic() +
    ylab("Count") +
    xlab("Reactions per network")

Fig4FGH <- plot_grid(Fig4F, Fig4G, Fig4H, labels=c('F','G','H'),nrow=1)

Fig4CDEFGH <- plot_grid(Fig4CDE, Fig4FGH, nrow = 2)
ggsave("../results/metabolic-enzymes/Fig4CDEFGH.pdf", height = 5, width = 7.5)
###############################################################################
## Make STIMS figures for core and essential genes in the 1000 minimal genomes.

## core genome of the minimal genomes.
minimal.core <- minimal.genome.count.df %>%
    filter(Count == 1000) %>%
    inner_join(REL606.genes)
## Write the minimal genome core genes to file.
jenga.genome.core.csv <- "../results/metabolic-enzymes/jenga-genome-core.csv"
write.csv(minimal.core, file = jenga.genome.core.csv)

## core essential genes.
minimal.essential <- essential.gene.count.df %>%
    filter(Count == 1000) %>%
    inner_join(REL606.genes)
## Write the essential core genes to file.
jenga.essential.core.csv = "../results/metabolic-enzymes/jenga-essential-core.csv"
write.csv(minimal.essential, file = jenga.essential.core.csv)


## plot just the hypermutator populations.
minimal.core.hypermut.data <- hypermutator.data %>%
    filter(Gene %in% minimal.core$Gene)

c.minimal.core.hypermut <- calc.cumulative.muts(
    minimal.core.hypermut.data,
    minimal.core,
    manual.pop.levels.vec = hypermutator.pops)

minimal.essential.hypermut.data <- hypermutator.data %>%
    filter(Gene %in% minimal.essential$Gene)

c.minimal.essential.hypermut <- calc.cumulative.muts(
    minimal.essential.hypermut.data,
    minimal.essential,
    manual.pop.levels.vec = hypermutator.pops)

Fig5A <- plot.base.layer(
    hypermutator.data,
    REL606.genes,
    subset.size=length(unique(minimal.core$Gene)),
    manual.pop.levels.vec = hypermutator.pops,
    plot.rows = 2,
    my.color = "plum1") %>%
    add.cumulative.mut.layer(c.minimal.core.hypermut, my.color="darkorchid4") +
    ggtitle("Core genes found in all 1000 minimal genomes")

Fig5B <- plot.base.layer(
    hypermutator.data,
    REL606.genes,
    subset.size=length(unique(minimal.essential$Gene)),
    manual.pop.levels.vec = hypermutator.pops,
    plot.rows = 2,
    my.color = "darkolivegreen1") %>%
    add.cumulative.mut.layer(c.minimal.essential.hypermut, my.color="springgreen4") +
    ggtitle("Essential genes found in all 1000 minimal genomes")

Fig5 <- plot_grid(Fig5A, Fig5B, labels=c('A','B'),nrow=2)
save_plot("../results/metabolic-enzymes/Fig5.pdf",Fig5, base_height=7,base_asp=1)

## plot just the nonmutator populations.
minimal.core.nonmut.data <- nonmutator.data %>%
    filter(Gene %in% minimal.core$Gene)

c.minimal.core.nonmut <- calc.cumulative.muts(
    minimal.core.nonmut.data,
    minimal.core,
    manual.pop.levels.vec = nonmutator.pops)

minimal.essential.nonmut.data <- nonmutator.data %>%
    filter(Gene %in% minimal.essential$Gene)

c.minimal.essential.nonmut <- calc.cumulative.muts(
    minimal.essential.nonmut.data,
    minimal.essential,
    manual.pop.levels.vec = nonmutator.pops)

S4FigA <- plot.base.layer(
    nonmutator.data,
    REL606.genes,
    subset.size=length(unique(minimal.core$Gene)),
    manual.pop.levels.vec = nonmutator.pops,
    plot.rows = 2,
    my.color = "plum1") %>%
    add.cumulative.mut.layer(c.minimal.core.nonmut, my.color="darkorchid4") +
    ggtitle("Core genes found in all 1000 minimal genomes")

S4FigB <- plot.base.layer(
    nonmutator.data,
    REL606.genes,
    subset.size=length(unique(minimal.essential$Gene)),
    manual.pop.levels.vec = nonmutator.pops,
    plot.rows = 2,
    my.color = "darkolivegreen1") %>%
    add.cumulative.mut.layer(c.minimal.essential.nonmut, my.color="springgreen4") +
    ggtitle("Essential genes found in all 1000 minimal genomes")

S4Fig <- plot_grid(S4FigA, S4FigB, labels=c('A','B'),nrow=2)
save_plot("../results/metabolic-enzymes/S4Fig.pdf",S4Fig, base_height=7,base_asp=1)


## Make STIMS figures for essential genes for glucose and citrate growth
## in the REL606 FBA metabolic model.

FBA.glucose.essential <- read.csv(
    "../results/metabolic-enzymes/glucose_FBA_essential.csv")

FBA.citrate.essential <- read.csv(
    "../results/metabolic-enzymes/citrate_FBA_essential.csv")

## plot just the hypermutator populations.
glucose.essential.hypermut.data <- hypermutator.data %>%
    filter(Gene %in% FBA.glucose.essential$Gene)

c.glucose.essential.hypermut <- calc.cumulative.muts(
    glucose.essential.hypermut.data,
    FBA.glucose.essential,
    manual.pop.levels.vec = hypermutator.pops)

citrate.essential.hypermut.data <- hypermutator.data %>%
    filter(Gene %in% FBA.citrate.essential$Gene)

c.citrate.essential.hypermut <- calc.cumulative.muts(
    citrate.essential.hypermut.data,
    FBA.citrate.essential,
    manual.pop.levels.vec = hypermutator.pops)

Fig6A <- plot.base.layer(
    hypermutator.data,
    REL606.genes,
    subset.size=length(unique(FBA.glucose.essential$Gene)),
    manual.pop.levels.vec = hypermutator.pops,
    plot.rows = 2,
    my.color = "plum1") %>%
    add.cumulative.mut.layer(c.glucose.essential.hypermut, my.color="darkorchid4") +
    ggtitle("Genes essential for growth on glucose in the REL606 metabolic model")

Fig6B <- plot.base.layer(
    hypermutator.data,
    REL606.genes,
    subset.size=length(unique(FBA.citrate.essential$Gene)),
    manual.pop.levels.vec = hypermutator.pops,
    plot.rows = 2,
    my.color = "darkolivegreen1") %>%
    add.cumulative.mut.layer(c.citrate.essential.hypermut, my.color="springgreen4") +
    ggtitle("Genes essential for growth on citrate in the REL606 metabolic model")

Fig6 <- plot_grid(Fig6A, Fig6B, labels=c('A','B'),nrow=2)
save_plot("../results/metabolic-enzymes/Fig6.pdf",Fig6, base_height=7,base_asp=1)

## plot just the nonmutator populations.
glucose.essential.nonmut.data <- nonmutator.data %>%
    filter(Gene %in% FBA.glucose.essential$Gene)

c.glucose.essential.nonmut <- calc.cumulative.muts(
    glucose.essential.nonmut.data,
    FBA.glucose.essential,
    manual.pop.levels.vec = nonmutator.pops)

citrate.essential.nonmut.data <- nonmutator.data %>%
    filter(Gene %in% FBA.citrate.essential$Gene)

c.citrate.essential.nonmut <- calc.cumulative.muts(
    citrate.essential.nonmut.data,
    FBA.citrate.essential,
    manual.pop.levels.vec = nonmutator.pops)

S5FigA <- plot.base.layer(
    nonmutator.data,
    REL606.genes,
    subset.size=length(unique(FBA.glucose.essential$Gene)),
    manual.pop.levels.vec = nonmutator.pops,
    plot.rows = 2,
    my.color = "plum1") %>%
    add.cumulative.mut.layer(c.glucose.essential.nonmut, my.color="darkorchid4") +
    ggtitle("Genes essential for growth on glucose in the REL606 metabolic model")

S5FigB <- plot.base.layer(
    nonmutator.data,
    REL606.genes,
    subset.size=length(unique(FBA.citrate.essential$Gene)),
    manual.pop.levels.vec = nonmutator.pops,
    plot.rows = 2,
    my.color = "darkolivegreen1") %>%
    add.cumulative.mut.layer(c.citrate.essential.nonmut, my.color="springgreen4") +
    ggtitle("Genes essential for growth on citrate in the REL606 metabolic model")

S5Fig <- plot_grid(S5FigA, S5FigB, labels=c('A','B'),nrow=2)
save_plot("../results/metabolic-enzymes/S5Fig.pdf",S5Fig, base_height=7,base_asp=1)
