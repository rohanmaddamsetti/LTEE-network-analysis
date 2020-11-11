## resilience-analysis.R by Rohan Maddamsetti.

## CRITICAL TODO: use the REL606 resilience value for the 14 genomes without
## any KO mutations, then re-calculate statistics.

## CRITICAL TODO: for a really fair comparison, I should fix the y-intercept for the linear regressions to
## the starting resilience value for REL606.

## CRITICAL TODO: As an additional control analysis,
## bootstrap resilience using KO mutations rather than genes, in some way.
## The idea here is that genes that are over-represented are probably under
## positive selection. So choose those more often. Then, is the
## realized resilience trajectory more shallow than this bootstrapped
## trajectory?

library(tidyverse)

################################################################
## PPI NETWORK RESILIENCE RESULTS.

## Order nonmutator pops, then hypermutator pops by converting Population to
## factor type and setting the levels.
nonmutator.pops <- c("Ara-5", "Ara-6", "Ara+1", "Ara+2", "Ara+4", "Ara+5")
hypermutator.pops <- c("Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara+3", "Ara+6")


## -- add LTEE strains with no KO to table, with ancestral resilience.

LTEE.genomes.KO.metadata <- read.csv("../data/LTEE-264-genomes-SNP-nonsense-small-indel-MOB-large-deletions.csv") %>%
    select(population,time,strain,clone,mutator_status) %>%
    distinct() %>%
    mutate(Generation=time/10000) %>%
    ## This for changing the ordering of populations in plots.
    mutate(population=factor(population,levels=c(nonmutator.pops,hypermutator.pops)))


zitnik.network.resilience.df <- read.csv("../results/resilience/Zitnik_PPI_LTEE_genome_resilience.csv") %>%
    full_join(LTEE.genomes.KO.metadata) %>%
    mutate(log.resilience=log(resilience))

zitnik.regression <- lm(data=zitnik.network.resilience.df,log.resilience~time)
cor.test(zitnik.network.resilience.df$time,zitnik.network.resilience.df$log.resilience)
cor.test(zitnik.network.resilience.df$time,zitnik.network.resilience.df$resilience)

zitnik.resilience.plot <- ggplot(zitnik.network.resilience.df,
                                 aes(x=Generation,y=resilience,color=population)) +
    facet_wrap(.~population) +
    geom_smooth(method="lm") +
    theme_classic() + geom_point()
ggsave("../results/resilience/zitnik-resilience.pdf", zitnik.resilience.plot)

cong.network.resilience.df <- read.csv("../results/resilience/Cong_PPI_LTEE_genome_resilience.csv") %>%
    full_join(LTEE.genomes.KO.metadata) %>%
    mutate(log.resilience=log(resilience))

cong.regression <- lm(data=cong.network.resilience.df,log.resilience~time)
cor.test(cong.network.resilience.df$time,cong.network.resilience.df$log.resilience)
cor.test(cong.network.resilience.df$time,cong.network.resilience.df$resilience)

cong.resilience.plot <- ggplot(cong.network.resilience.df,
                                 aes(x=Generation,y=resilience,color=population)) +
    facet_wrap(.~population) +
    geom_smooth(method="lm") +
    theme_classic() + geom_point()
ggsave("../results/resilience/cong-resilience.pdf", cong.resilience.plot)

#######################################################################
## Now examine the preliminary randomized results.

################### within populations.

## for Zitnik PPI network.
zitnik.randomized.within <- read.csv("../results/resilience/Zitnik_PPI_within_pops_randomized_resilience.csv") %>%
    full_join(LTEE.genomes.KO.metadata) %>%
    mutate(log.resilience=log(randomized_resilience))

zitnik.randomized.within.lm <- lm(data=zitnik.randomized.within,log.resilience~time)
cor.test(zitnik.randomized.within$time,zitnik.randomized.within$log.resilience)
cor.test(zitnik.randomized.within$time,zitnik.randomized.within$randomized_resilience)

zitnik.randomized.within.plot <- ggplot(zitnik.randomized.within,
                                 aes(x=Generation,y=randomized_resilience,color=population)) +
    facet_wrap(.~population) +
    geom_smooth(method="lm") +
    theme_classic() + geom_point()
ggsave("../results/resilience/zitnik-randomized-within.pdf", zitnik.randomized.within.plot)

## for Cong PPI network.
cong.randomized.within <- read.csv("../results/resilience/Cong_PPI_within_pops_randomized_resilience.csv") %>%
    full_join(LTEE.genomes.KO.metadata) %>%
    mutate(log.resilience=log(randomized_resilience))

cong.randomized.within.lm <- lm(data=cong.randomized.within,log.resilience~time)
cor.test(cong.randomized.within$time,cong.randomized.within$log.resilience)
cor.test(cong.randomized.within$time,cong.randomized.within$randomized_resilience)

cong.randomized.within.plot <- ggplot(cong.randomized.within,
                                 aes(x=Generation,y=randomized_resilience,color=population)) +
    facet_wrap(.~population) +
    geom_smooth(method="lm") +
    theme_classic() + geom_point()
ggsave("../results/resilience/cong-randomized-within.pdf", cong.randomized.within.plot)


###################### across populations.
# for Zitnik PPI network.
zitnik.randomized.across <- read.csv("../results/resilience/Zitnik_PPI_across_pops_randomized_resilience.csv") %>%
    full_join(LTEE.genomes.KO.metadata) %>%
    mutate(log.resilience=log(randomized_resilience))

zitnik.randomized.across.lm <- lm(data=zitnik.randomized.across,log.resilience~time)
cor.test(zitnik.randomized.across$time,zitnik.randomized.across$log.resilience)
cor.test(zitnik.randomized.across$time,zitnik.randomized.across$randomized_resilience)

zitnik.randomized.across.plot <- ggplot(zitnik.randomized.across,
                                 aes(x=Generation,y=randomized_resilience,color=population)) +
    facet_wrap(.~population) +
    geom_smooth(method="lm") +
    theme_classic() + geom_point()
ggsave("../results/resilience/zitnik-randomized-across.pdf", zitnik.randomized.across.plot)

## for Cong PPI network.
cong.randomized.across <- read.csv("../results/resilience/Cong_PPI_across_pops_randomized_resilience.csv") %>%
    full_join(LTEE.genomes.KO.metadata) %>%
    mutate(log.resilience=log(randomized_resilience))

cong.randomized.across.lm <- lm(data=cong.randomized.across,log.resilience~time)
cor.test(cong.randomized.across$time,cong.randomized.across$log.resilience)
cor.test(cong.randomized.across$time,cong.randomized.across$randomized_resilience)

cong.randomized.across.plot <- ggplot(cong.randomized.across,
                                 aes(x=Generation,y=randomized_resilience,color=population)) +
    facet_wrap(.~population) +
    geom_smooth(method="lm") +
    theme_classic() + geom_point()
ggsave("../results/resilience/cong-randomized-across.pdf", cong.randomized.across.plot)

############################# all genes allowed.

# for Zitnik PPI network.
zitnik.randomized.all <- read.csv("../results/resilience/Zitnik_PPI_all_genes_randomized_resilience.csv") %>%
    full_join(LTEE.genomes.KO.metadata) %>%
    mutate(log.resilience=log(randomized_resilience))

zitnik.randomized.all.lm <- lm(data=zitnik.randomized.all,log.resilience~time)
cor.test(zitnik.randomized.all$time,zitnik.randomized.all$log.resilience)
cor.test(zitnik.randomized.all$time,zitnik.randomized.all$randomized_resilience)

zitnik.randomized.all.plot <- ggplot(zitnik.randomized.all,
                                 aes(x=Generation,y=randomized_resilience,color=population)) +
    facet_wrap(.~population) +
    geom_smooth(method="lm") +
    theme_classic() + geom_point()
ggsave("../results/resilience/zitnik-randomized-all.pdf", zitnik.randomized.all.plot)

## for Cong PPI network.
cong.randomized.all <- read.csv("../results/resilience/Cong_PPI_all_genes_randomized_resilience.csv") %>%
    full_join(LTEE.genomes.KO.metadata) %>%
    mutate(log.resilience=log(randomized_resilience))

cong.randomized.all.lm <- lm(data=cong.randomized.all,log.resilience~time)
cor.test(cong.randomized.all$time,cong.randomized.all$log.resilience)
cor.test(cong.randomized.all$time,cong.randomized.all$randomized_resilience)

cong.randomized.all.plot <- ggplot(cong.randomized.all,
                                 aes(x=Generation,y=randomized_resilience,color=population)) +
    facet_wrap(.~population) +
    geom_smooth(method="lm") +
    theme_classic() + geom_point()
ggsave("../results/resilience/cong-randomized-all.pdf", cong.randomized.all.plot)

## although more simulation runs are needed, these results suggest purifying selection on
## which genes are affected by KO mutations in each population,
## as resilience falls more slowly in the real data, than in the simulated data.

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
superessential.rxns.df <- read.csv("../results/thermostability/Barve2012-S6-superessential.csv")

superessential.mut.data <- gene.mutation.data %>%
    filter(Gene %in% superessential.rxns.df$Gene)

c.superessential <- calc.cumulative.muts(superessential.mut.data)

superessential.base.layer <- plot.base.layer(
    gene.mutation.data,
    subset.size=length(unique(superessential.rxns.df$Gene)))

## plot of superessential metabolic enzymes analysis
superessential.fig <- superessential.base.layer %>% 
    add.cumulative.mut.layer(c.superessential, my.color="black")
ggsave("../results/thermostability/figures/superessential.pdf", superessential.fig)

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
Nam.df <- read.csv("../results/thermostability/Nam2012_Database_S1.csv") %>%
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
ggsave("../results/thermostability/figures/specialist.pdf", specialist.fig)

generalist.fig <- generalist.base.layer %>% ## null for generalists
    add.cumulative.mut.layer(c.generalists, my.color="black")
ggsave("../results/thermostability/figures/generalist.pdf", generalist.fig)

######### TODO: Do network resilience analysis on the E. coli metabolic network.
## can I come up with a good explanation for the idiosyncratic patterns
## of purifying selection on these sets of metabolic enzymes, across LTEE
## populations?
