## resilience-analysis.R by Rohan Maddamsetti.

library(tidyverse)
library(cowplot)

################################################################
## FUNCTIONS.

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
## FUNCTIONS for big resilience data analysis.

add_timeseries_column <- function(the.big.df) {

    ## there are six possible values for the timeseries column,
    ## based on the value in the run_type column.
    ## actual_data, randomized_within, weighted_randomized_within,
    ## randomized_across, weighted_randomized_across, randomized_all.
    the.big.df %>%
        mutate(
            timeseries = case_when(
                run_type == "LTEE_genome_resilience" ~ "actual_data",
                run_type == "across_pops_randomized_resilience"
                ~ "randomized_across",
                run_type == "across_pops_weighted_randomized_resilience"
                ~ "weighted_randomized_across",
                run_type == "all_genes_randomized_resilience" ~ "randomized_all"
            )
        )
}

##########################################################################
## IMPORT METADATA.

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
## PPI network resilience results, using Duke Compute Cluster runs.

## IMPORTANT BUG TO FIX: REL958B and REL1066B have NA values in their rows. Why?
filter(big.resilience.df,is.na(replicate))

big.resilience.df <- read.csv("../results/resilience/collated-resilience-runs.csv") %>%
    full_join(LTEE.genomes.KO.metadata) %>%
    set.no.KO.strains.resilience.to.REL606() %>%
    mutate(log.resilience=log(resilience)) %>%
    add_timeseries_column() %>%
    ## HACK TO REMOVE THE TWO PROBLEMATIC CASES FOR NOW (BEFORE DEBUGGING):
    filter(!is.na(replicate))

big.cong.resilience.df <- big.resilience.df %>%
    filter(dataset=="Cong")

big.zitnik.resilience.df <- big.resilience.df %>%
    filter(dataset=="Zitnik")


big.cong.resilience.plot <- make.big.Cong.resilience.plot(big.cong.resilience.df)
big.zitnik.resilience.plot <- make.big.Zitnik.resilience.plot(big.zitnik.resilience.df)

Fig1.legend <- get_legend(make.big.resilience.plot(big.cong.resilience.df,
                                                   plot.legend=TRUE))

Fig1 <- plot_grid(plot_grid(big.zitnik.resilience.plot,big.cong.resilience.plot,
                            labels=c('A','B',NULL),nrow=1, rel_widths=c(1,1,0.4)),
                  Fig1.legend,ncol=1, rel_heights=c(1,0.05))
ggsave("../results/resilience/figures/Fig1.pdf",height=7,width=6)


## Calculate statistics on differences between slopes.
big.cong.slope.df <- big.cong.resilience.df %>%
    split(.$timeseries) %>%
    map_dfr(.f = calc.resilience.regression.df)

big.zitnik.slope.df <- big.zitnik.resilience.df %>%
    split(.$timeseries) %>%
    map_dfr(.f = calc.resilience.regression.df)

regression.slope.test(big.zitnik.slope.df, "randomized_across")
regression.slope.test(big.zitnik.slope.df, "randomized_all")

regression.slope.test(big.cong.slope.df, "randomized_across")
regression.slope.test(big.cong.slope.df, "randomized_all")

#######################################################################
## ORIGINAL PPI network resilience results.

## Now, Zitnik PPI resilience analysis.
zitnik.network.resilience.df <- read.network.resilience.df(
    "../results/resilience/original-run-results/Zitnik_PPI_LTEE_genome_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "actual_data")

zitnik.randomized.within <- read.network.resilience.df(
    "../results/resilience/original-run-results/Zitnik_PPI_within_pops_randomized_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "randomized_within")

zitnik.randomized.across <- read.network.resilience.df(
    "../results/resilience/original-run-results/Zitnik_PPI_across_pops_randomized_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "randomized_across")

zitnik.randomized.all <- read.network.resilience.df(
    "../results/resilience/original-run-results/Zitnik_PPI_all_genes_randomized_resilience.csv",
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
    "../results/resilience/original-run-results/Cong_PPI_LTEE_genome_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "actual_data")

cong.randomized.within <- read.network.resilience.df(
    "../results/resilience/original-run-results/Cong_PPI_within_pops_randomized_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "randomized_within")

cong.randomized.across <- read.network.resilience.df(
    "../results/resilience/original-run-results/Cong_PPI_across_pops_randomized_resilience.csv",
    LTEE.genomes.KO.metadata) %>%
    mutate(timeseries = "randomized_across")

cong.randomized.all <- read.network.resilience.df(
    "../results/resilience/original-run-results/Cong_PPI_all_genes_randomized_resilience.csv",
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
ggsave("../results/resilience/figures/oldFig1.pdf",height=7,width=6)

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

