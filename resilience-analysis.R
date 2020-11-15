## resilience-analysis.R by Rohan Maddamsetti.

## CRITICAL TODO: for a really fair comparison,
## I should fix the y-intercept for the linear regressions to
## the starting resilience value for REL606.

## CRITICAL TODO: As an additional control analysis,
## bootstrap resilience using KO mutations rather than genes, in some way.
## The idea here is that genes that are over-represented are probably under
## positive selection. So choose those more often. Then, is the
## realized resilience trajectory more shallow than this bootstrapped
## trajectory?

library(tidyverse)
library(cowplot)

################################################################
## PPI NETWORK RESILIENCE RESULTS.

set.no.KO.strains.resilience.to.REL606 <- function(df) {
    ## strains with no KO mutations get the ancestral resilience value.
    ancestral.resilience <- unique(filter(df, strain == 'REL606')$resilience)
    new.df <- df %>% mutate(resilience = ifelse(is.na(resilience),
                                      ancestral.resilience,
                                      resilience))
    return(new.df)
}

read.network.resilience.df <- function(path, LTEE.genome.metadata) {
    ## This function makes sure that metadata and resilience values are set
    ## consistently, by wrapping the stuff that needs to happen in the right order
    ## of function calls.
    resilience.df <- read.csv(path) %>%
        select(-X) %>% ## remove pandas crud
        full_join(LTEE.genome.metadata) %>%
        set.no.KO.strains.resilience.to.REL606() %>%
        mutate(log.resilience=log(resilience))
    
    return(resilience.df)
}

make.resilience.plot <- function(resilience.df) {
    p <- ggplot(resilience.df,
                aes(x=Generation,
                    y=resilience,
                    color=population)) +
        facet_wrap(.~population,nrow=4) +
        geom_smooth(method="lm") +
        theme_classic() + geom_point() + guides(color=FALSE) +
        ylab("Network resilience")
    return(p)
}

make.Cong.resilience.plot <- function(resilience.df) {
    p <- make.resilience.plot(resilience.df) +
        ## IMPORTANT: axes need to be consistent when comparing plots.
        ylim(0.15,0.18)
    return(p)
}

make.Zitnik.resilience.plot <- function(resilience.df) {
    p <- make.resilience.plot(resilience.df) +
        ## IMPORTANT: axes need to be consistent when comparing plots.
        ylim(0.40,0.414)
    return(p)

}

make.big.resilience.plot <- function(big.resilience.df, plot.legend=FALSE) {
    ## plot data and randomized data on the same figure.
    p <- ggplot(big.resilience.df,
                aes(x = Generation,
                    y = resilience,
                    color = timeseries)) +
        facet_wrap(.~population,nrow=4) +
        theme_classic() +
        geom_point(alpha = 0.2, size = 0.5) +
        geom_smooth(alpha = 0.2, method = "lm") +
        ylab("Network resilience")
    if (!plot.legend) {
        p <- p + guides(color=FALSE)
    }
    return(p)
}

make.big.Cong.resilience.plot <- function(big.resilience.df) {
    p <- make.big.resilience.plot(big.resilience.df) +
        ## IMPORTANT: axes need to be consistent when comparing plots.
        ylim(0.15,0.18) +
        ggtitle("Cong PPI dataset")
    return(p)
}

make.big.Zitnik.resilience.plot <- function(big.resilience.df) {
    p <- make.big.resilience.plot(big.resilience.df) +
        ## IMPORTANT: axes need to be consistent when comparing plots.
        ylim(0.40,0.414) +
        ggtitle("Zitnik PPI dataset")
    return(p)

}

calc.regression.and.correlations <- function(resilience.df) {
    ## TODO: fix the y-intercept at the ancestral REL606 resilience.

    ## these two correlations should be the same, within numerical error.
    print(cor.test(resilience.df$time, resilience.df$log.resilience))
    print(cor.test(resilience.df$time, resilience.df$resilience))
    my.regression <- lm(data=resilience.df, log.resilience~time)
    return(my.regression)
}

## Order nonmutator pops, then hypermutator pops by converting Population to
## factor type and setting the levels.
nonmutator.pops <- c("Ara-5", "Ara-6", "Ara+1", "Ara+2", "Ara+4", "Ara+5")
hypermutator.pops <- c("Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara+3", "Ara+6")

## we will add REL606 as ancestor to all 12 populations.
REL606.metadata <- data.frame(strain = rep('REL606', 12),
                              population = c(nonmutator.pops, hypermutator.pops),
                              time = rep(0, 12),
                              clone = rep('A', 12),
                              mutator_status = rep('non-mutator', 12))

LTEE.genomes.KO.metadata <- read.csv(
    "../data/LTEE-264-genomes-SNP-nonsense-small-indel-MOB-large-deletions.csv") %>%
    select(population,time,strain,clone,mutator_status) %>%
    distinct() %>%
    ## let's add REL606.
    full_join(REL606.metadata) %>%
    mutate(Generation=time/10000) %>%
    ## This for changing the ordering of populations in plots.
    mutate(population=factor(population,levels=c(nonmutator.pops,hypermutator.pops)))

#######################################################################
## PPI network resilience results.

## Now, Zitnik PPI resilience analysis.
zitnik.network.resilience.df <- read.network.resilience.df(
    "../results/resilience/Zitnik_PPI_LTEE_genome_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "actual_data")

zitnik.randomized.within <- read.network.resilience.df(
    "../results/resilience/Zitnik_PPI_within_pops_randomized_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "randomized_within")

zitnik.randomized.across <- read.network.resilience.df(
    "../results/resilience/Zitnik_PPI_across_pops_randomized_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "randomized_across")

zitnik.randomized.all <- read.network.resilience.df(
    "../results/resilience/Zitnik_PPI_all_genes_randomized_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "randomized_all")

## make a big data frame to plot all the data in the same figure.
big.zitnik.df <- zitnik.network.resilience.df %>%
    full_join(zitnik.randomized.within) %>%
    full_join(zitnik.randomized.across) %>%
    full_join(zitnik.randomized.all)
## make the big plot for Zitnik data and randomized data.
big.zitnik.plot <- make.big.Zitnik.resilience.plot(big.zitnik.df)

zitnik.regression <- calc.regression.and.correlations(zitnik.network.resilience.df)
zitnik.randomized.within.lm <- calc.regression.and.correlations(zitnik.randomized.within)
zitnik.randomized.across.lm <- calc.regression.and.correlations(zitnik.randomized.across)
zitnik.randomized.all.lm <- calc.regression.and.correlations(zitnik.randomized.all)

zitnik.resilience.plot <- make.Zitnik.resilience.plot(zitnik.network.resilience.df)
zitnik.randomized.within.plot <- make.Zitnik.resilience.plot(zitnik.randomized.within)
zitnik.randomized.across.plot <- make.Zitnik.resilience.plot(zitnik.randomized.across)
zitnik.randomized.all.plot <- make.Zitnik.resilience.plot(zitnik.randomized.all)

ggsave("../results/resilience/zitnik-resilience.pdf", zitnik.resilience.plot)
ggsave("../results/resilience/zitnik-randomized-within.pdf", zitnik.randomized.within.plot)
ggsave("../results/resilience/zitnik-randomized-across.pdf", zitnik.randomized.across.plot)
ggsave("../results/resilience/zitnik-randomized-all.pdf", zitnik.randomized.all.plot)

## Now, Cong PPI resilience analysis.
cong.network.resilience.df <- read.network.resilience.df(
    "../results/resilience/Cong_PPI_LTEE_genome_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "actual_data")

cong.randomized.within <- read.network.resilience.df(
    "../results/resilience/Cong_PPI_within_pops_randomized_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "randomized_within")

cong.randomized.across <- read.network.resilience.df(
    "../results/resilience/Cong_PPI_across_pops_randomized_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "randomized_across")

cong.randomized.all <- read.network.resilience.df(
    "../results/resilience/Cong_PPI_all_genes_randomized_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "randomized_all")

big.cong.df <- cong.network.resilience.df %>%
    full_join(cong.randomized.within) %>%
    full_join(cong.randomized.across) %>%
    full_join(cong.randomized.all)
## make the big plot for Cong data and randomized data.
big.cong.plot <- make.big.Cong.resilience.plot(big.cong.df)
Fig1.legend <- get_legend(make.big.resilience.plot(big.cong.df, plot.legend=TRUE))

Fig1 <- plot_grid(big.zitnik.plot,big.cong.plot,Fig1.legend,
                  labels=c('A','B',NULL),nrow=1, rel_widths=c(1,1,0.4))
ggsave("../results/resilience/Fig1.pdf",height=6, width=12)

cong.regression <- calc.regression.and.correlations(cong.network.resilience.df)
cong.randomized.within.lm <- calc.regression.and.correlations(cong.randomized.within)
cong.randomized.across.lm <- calc.regression.and.correlations(cong.randomized.across)
cong.randomized.all.lm <- calc.regression.and.correlations(cong.randomized.all)

cong.resilience.plot <- make.Cong.resilience.plot(cong.network.resilience.df)
cong.randomized.within.plot <- make.Cong.resilience.plot(cong.randomized.within)
cong.randomized.across.plot <- make.Cong.resilience.plot(cong.randomized.across)
cong.randomized.all.plot <- make.Cong.resilience.plot(cong.randomized.all)

ggsave("../results/resilience/cong-resilience.pdf", cong.resilience.plot)
ggsave("../results/resilience/cong-randomized-within.pdf", cong.randomized.within.plot)
ggsave("../results/resilience/cong-randomized-across.pdf", cong.randomized.across.plot)
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
