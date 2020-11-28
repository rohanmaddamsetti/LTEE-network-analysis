## metabolic-enzymes.R by Rohan Maddamsetti

source("metagenomics-library.R")

####################
## DATA PREPROCESSING

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

## import genes in the table created by:
## conda activate ltee-metagenomix
## cd LTEE-metagenomic-repo
## python rohan-write-csv.py > ../results/LTEE-metagenome-mutations.csv
gene.mutation.data <- read.csv(
    '../results/LTEE-metagenome-mutations.csv',
    header=TRUE,as.is=TRUE) %>%
    mutate(Generation=t0/10000) %>%
    ## This for changing the ordering of populations in plots.
    mutate(Population=factor(Population,
                             levels=c(nonmutator.pops,hypermutator.pops))) %>%
    inner_join(REL606.genes) %>%
    filter(Gene!='intergenic')

KO.gene.mutation.data <- gene.mutation.data %>%
    filter(Annotation %in% c("sv","indel","nonsense"))

#######################################################################
## METABOLIC ENZYME ANALYSIS.
########################################################

## Some analyses in this vein have already been published in:
## Metabolic Determinants of Enzyme Evolution in a Genome-Scale Bacterial Metabolic Network
## by Jose Aguilar-Rodriguez and Andreas Wagner.

## Only Ara-1 and Ara+6 show evidence of purifying selection on
## superessential metabolic enzymes. Why only these two? No idea why.

## In short, there looks like there is idiosyncratic purifying selection on
## core metabolic networks (superessential reactions and specialist enzymes)
## in the LTEE. Results depend on the population, and are not consistent across
## populations.

## Report these findings in the main text, but keep it brief, and put all STIMS
## figures into the Supplement.

## examine superessential metabolic reactions reported by Barve and Wagner (2012).
superessential.rxns.df <- read.csv("../results/metabolic-enzymes/Barve2012-S6-superessential.csv")

superessential.mut.data <- gene.mutation.data %>%
    filter(Gene %in% superessential.rxns.df$Gene)

c.superessential <- calc.cumulative.muts(superessential.mut.data)

superessential.base.layer <- plot.base.layer(
    gene.mutation.data,
    subset.size=length(unique(superessential.rxns.df$Gene)))

## plot of superessential metabolic enzymes analysis
superessential.fig <- superessential.base.layer %>% 
    add.cumulative.mut.layer(c.superessential, my.color="black")
ggsave("../results/metabolic-enzymes/superessential.pdf", superessential.fig)

## calculate formal p-values.
superessential.pvals <- calc.traj.pvals(gene.mutation.data, unique(superessential.rxns.df$Gene))
## results:
## A tibble: 12 x 3
##   Population count p.val
##   <fct>      <int> <dbl>
## 1 Ara-5       6363 0.636
## 2 Ara-6       9150 0.915
## 3 Ara+1       8365 0.836
## 4 Ara+2       6800 0.68 
## 5 Ara+4       8913 0.891
## 6 Ara+5       1892 0.189
## 7 Ara-1       9999 1.00 
## 8 Ara-2       6484 0.648
## 9 Ara-3       6178 0.618
##10 Ara-4       8766 0.877
##11 Ara+3       5669 0.567
##12 Ara+6       9933 0.993

################################################################
## look at specialist and generalist enzymes in Nam et al. (2012):
## Network context and selection in the evolution of enzyme specificity.

## 1157 genes.
Nam.df <- read.csv("../results/metabolic-enzymes/Nam2012_Database_S1.csv") %>%
    left_join(REL606.genes) %>% filter(!is.na(Gene))

specialist.enzymes <- Nam.df %>% filter(Class=="Spec.")
generalist.enzymes <- Nam.df %>% filter(Class=="Gen.")

specialist.mut.data <- gene.mutation.data %>%
    filter(Gene %in% specialist.enzymes$Gene)
c.specialists <- calc.cumulative.muts(specialist.mut.data)

generalist.mut.data <- gene.mutation.data %>%
    filter(Gene %in% generalist.enzymes$Gene)
c.generalists <- calc.cumulative.muts(generalist.mut.data)

specialist.base.layer <- plot.base.layer(
    gene.mutation.data,
    subset.size=length(unique(specialist.enzymes$Gene)))

generalist.base.layer <- plot.base.layer(
    gene.mutation.data,
    subset.size=length(unique(generalist.enzymes$Gene)))

## specialist fig.
specialist.fig <- specialist.base.layer %>% ## null for specialists
    add.cumulative.mut.layer(c.specialists, my.color="black")
ggsave("../results/metabolic-enzymes/specialist.pdf", specialist.fig)

generalist.fig <- generalist.base.layer %>% ## null for generalists
    add.cumulative.mut.layer(c.generalists, my.color="black")
ggsave("../results/metabolic-enzymes/generalist.pdf", generalist.fig)

#####################################################################################
## Run STIMS on BiGG Models E. coli core. Hypothesis: E. coli core metabolism
## is evolving under purifying selection in the LTEE.

BiGG.core <- read.csv("../results/metabolic-enzymes/BiGG-Model-Ecoli_core.csv") %>%
    mutate(blattner = BiGG.ID) %>%
    left_join(REL606.genes)

BiGG.core.mut.data <- gene.mutation.data %>%
    filter(Gene %in% BiGG.core$Gene)
c.BiGG.core <- calc.cumulative.muts(BiGG.core.mut.data)

BiGG.core.base.layer <- plot.base.layer(
    gene.mutation.data,
    subset.size=length(unique(BiGG.core$Gene)))

## BiGG core fig.
BiGG.core.fig <- BiGG.core.base.layer %>% ## null for BiGG core
    add.cumulative.mut.layer(c.BiGG.core, my.color="black")
ggsave("../results/metabolic-enzymes/BiGG-core.pdf", BiGG.core.fig)

## can I come up with a good explanation for the idiosyncratic patterns
## of purifying selection on these sets of metabolic enzymes, across LTEE
## populations? Anything that's testable? Check out and cite Tim Cooper's
## recent paper on pykF and epistasis.

## One guess: pykF knockout is why superessential genes are under strong
## purifying selection in Ara-1, but not the other pops.
