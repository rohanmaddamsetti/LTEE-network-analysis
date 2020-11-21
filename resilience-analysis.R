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
        ylab("Network resilience") +
        theme(legend.position="bottom")
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

calc.grand.correlation <- function(resilience.df) {
    print(cor.test(resilience.df$time, resilience.df$resilience))
}

calc.resilience.regression.df <- function(resilience.df) {
    ## returns a dataframe with the linear regression coefficients
    ## per population as a column.

    calc.resilience.slope.per.pop.df <- function(pop.df) {
        ## get the ancestral resilience.
        REL606.df <- pop.df %>%
            filter(Generation == 0)
        ancestral.resilience <- unique(REL606.df$resilience)
        ## fix the y-intercept at the ancestral REL606 resilience.
        my.regression <- lm(data=pop.df, I(resilience - ancestral.resilience) ~ 0 + Generation)
        resilience.coefficient <- my.regression$coefficients[[1]]
        pop.df.with.slope <- pop.df %>%
            select(population, timeseries) %>%
            mutate(resilience.slope = resilience.coefficient) %>%
            distinct()
        return(pop.df.with.slope)
    }

    summary.df <- resilience.df %>%
        split(.$population) %>%
        map_dfr(.f = calc.resilience.slope.per.pop.df)
    return(summary.df)
}

regression.slope.test <- function(full.slope.df, null.comp.string) {

    if (null.comp.string == "randomized_within") {
        null.df <- full.slope.df %>%
            filter(timeseries == "randomized_within")
    } else if (null.comp.string == "randomized_across") {
        null.df <- full.slope.df %>%
            filter(timeseries == "randomized_across")
    } else if (null.comp.string == "randomized_all") {
        null.df <- full.slope.df %>%
            filter(timeseries == "randomized_all")
    } else {
        print("ERROR: null.comp.string is not recognized.")
        return(NULL)
    }

    data.df <- full.slope.df %>%
        filter(timeseries == "actual_data")

    ## a critical assumption of the statistical test is that entries of the pair of
    ## vectors being compared correspond to matching populations.
    ## this assertion checks this critical invariant.
    stopifnot(data.df$population == null.df$population)

    ## print out the data going into the test for error checking.
    difference.vec <- data.df$resilience.slope - null.df$resilience.slope
    print("DATA VECTOR:")
    print(data.df$resilience.slope)
    print("NULL VECTOR:")
    print(null.df$resilience.slope)
    print("DATA - NULL:")
    print(difference.vec)

    ## calculate Wilcoxon signed-rank test to do a nonparametric paired sample test.
    wilcox.test(data.df$resilience.slope,
                null.df$resilience.slope,
                paired = TRUE, alternative = "greater")
}

##########################################################################
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

zitnik.slope.df <- calc.resilience.regression.df(zitnik.network.resilience.df)
zitnik.randomized.within.slope.df <- calc.resilience.regression.df(zitnik.randomized.within)
zitnik.randomized.across.slope.df <- calc.resilience.regression.df(zitnik.randomized.across)
zitnik.randomized.all.slope.df <- calc.resilience.regression.df(zitnik.randomized.all)

full.zitnik.slope.df <- zitnik.slope.df %>%
    full_join(zitnik.randomized.within.slope.df) %>%
    full_join(zitnik.randomized.across.slope.df) %>%
    full_join(zitnik.randomized.all.slope.df)

regression.slope.test(full.zitnik.slope.df, "randomized_within")
regression.slope.test(full.zitnik.slope.df, "randomized_across")
regression.slope.test(full.zitnik.slope.df, "randomized_all")

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

Fig1 <- plot_grid(plot_grid(big.zitnik.plot,big.cong.plot,
                            labels=c('A','B',NULL),nrow=1, rel_widths=c(1,1,0.4)),
                  Fig1.legend,ncol=1, rel_heights=c(1,0.05))
ggsave("../results/resilience/Fig1.pdf",height=7,width=6)

cong.slope.df <- calc.resilience.regression.df(cong.network.resilience.df)
cong.randomized.within.slope.df <- calc.resilience.regression.df(cong.randomized.within)
cong.randomized.across.slope.df <- calc.resilience.regression.df(cong.randomized.across)
cong.randomized.all.slope.df <- calc.resilience.regression.df(cong.randomized.all)

full.cong.slope.df <- cong.slope.df %>%
    full_join(cong.randomized.within.slope.df) %>%
    full_join(cong.randomized.across.slope.df) %>%
    full_join(cong.randomized.all.slope.df)

regression.slope.test(full.cong.slope.df, "randomized_within")
regression.slope.test(full.cong.slope.df, "randomized_across")
regression.slope.test(full.cong.slope.df, "randomized_all")


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
