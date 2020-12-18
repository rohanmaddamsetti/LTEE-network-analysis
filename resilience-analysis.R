## resilience-analysis.R by Rohan Maddamsetti.

library(tidyverse)
library(cowplot)
library(circlize)

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

    if (null.comp.string == "randomized_across") {
        null.df <- full.slope.df %>%
            filter(timeseries == "randomized_across")
    } else if (null.comp.string == "weighted_randomized_across") {
        null.df <- full.slope.df %>%
            filter(timeseries == "weighted_randomized_across")
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
regression.slope.test(big.zitnik.slope.df, "weighted_randomized_across")
regression.slope.test(big.zitnik.slope.df, "randomized_all")

regression.slope.test(big.cong.slope.df, "randomized_across")
regression.slope.test(big.cong.slope.df, "weighted_randomized_across")
regression.slope.test(big.cong.slope.df, "randomized_all")

########################################################################

## Make a circos plot for the zitnik and cong networks, as affected by KO
## mutations in the 50K A clones.

LTEE.50K.A.clone.metadata <- LTEE.genomes.KO.metadata %>%
    filter(time == 50000) %>%
    filter(clone == 'A')


## First, make a circos plot for each of the PPI networks.

## Then make a facet style plot, showing the edges that are deleted due to KO
## mutations in each of the 12 50K A clones.


REL606.genes <- read.csv("../results/REL606_IDs.csv")

## get edges in Zitnik PPI network.
raw.zitnik.edges <- read.csv("../results/thermostability/Ecoli-Zitnik-data/Ecoli-treeoflife.interactomes/511145.txt", sep=" ", header=FALSE)

## get edges in Cong PPI network.
raw.cong.edges <- read.csv("../results/thermostability/Cong-good-interaction-set.tsv",sep="\t",header=FALSE) %>%
    mutate(edge.number = row_number())

cong.node1 <- raw.cong.edges %>% 
    mutate(Gene = V1) %>%
    select(Gene, edge.number) %>%
    left_join(REL606.genes) %>%
    ## remove edges with NA values.
    filter(complete.cases(.))

cong.node2 <- raw.cong.edges %>% 
    mutate(Gene = V2) %>%
    select(Gene, edge.number) %>%
    left_join(REL606.genes) %>%
    ## remove edges with NA values.
    filter(complete.cases(.))

## filter for edges between nodes that have not been removed.
cong.node1 <- cong.node1 %>%
    filter(edge.number %in% cong.node2$edge.number) %>%
    mutate(chr = "REL606") %>%
    relocate(chr) %>%
    relocate(Gene, .after = last_col())

cong.node2 <- cong.node2 %>%
    filter(edge.number %in% cong.node1$edge.number) %>%
    mutate(chr = "REL606") %>%
    relocate(chr) %>%
    relocate(Gene, .after = last_col())


## just use the essential columns for plotting.
cong.node1.bed <- cong.node1 %>%
    select(chr, start, end, Gene)
cong.node2.bed <- cong.node2 %>%
    select(chr, start, end, Gene)

## one row df to initialize plot with one sector.
REL606.sector.df <- data.frame(chr="REL606",start=1, end=4629812)
circos.genomicInitialize(REL606.sector.df)


## now plot the cong PPI network.
circos.genomicLink(cong.node1.bed, cong.node2.bed,
                   col = rep("light gray", nrow(cong.node1.bed)),
                   border = NA)


## now make a circos plot for each of the 12 50K A clones,
## showing the edges that are deleted due to KO mutations.



## Then make a facet style plot, showing the edges that are deleted due to KO
## mutations in each of the 12 50K A clones.
