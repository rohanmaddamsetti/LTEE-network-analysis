## resilience-analysis.R by Rohan Maddamsetti.

library(tidyverse)
library(cowplot)
###########################################################################

REL606.CHR.LENGTH <- 4629812 ## constant used in a few places in the code.

## FUNCTIONS.

set.no.KO.strains.resilience.to.REL606 <- function(df) {
    
    ## get the average ancestral resilience.
    REL606.df <- df %>%
        filter(Generation == 0)
    ancestral.resilience <- mean(unique(REL606.df$resilience))

    ## strains with no KO mutations get the ancestral resilience value.
    new.df <- df %>% mutate(resilience = ifelse(is.na(resilience),
                                      ancestral.resilience,
                                      resilience))
    return(new.df)
}


add_Treatment_column <- function(the.big.df) {

    ## there are six possible values for the Treatment column,
    ## based on the value in the run_type column.
    ## actual data, randomized_within, weighted_randomized_within,
    ## randomized_across, randomized over LTEE, randomized over genome.
    the.big.df %>%
        mutate(
            Treatment = case_when(
                run_type == "LTEE_genome_resilience" ~ "actual data",
                run_type == "across_pops_randomized_resilience"
                ~ "randomized_across",
                run_type == "across_pops_weighted_randomized_resilience"
                ~ "randomized over LTEE",
                run_type == "all_genes_randomized_resilience" ~ "randomized over genome"
            )
        )
}


make.big.resilience.plot <- function(big.resilience.df, plot.legend=FALSE) {

    ## for plotting linear regressions.
    lm.df <- calc.resilience.regression.df(big.resilience.df)

    ## plot data and randomized data on the same figure.
    p <- ggplot(big.resilience.df,
                aes(x = Generation,
                    y = resilience,
                    color = Treatment)) +
        facet_wrap(.~population,nrow=4) +
        theme_classic() +
        ylab("PPI network resilience") +
        xlab("Time (x 10,000 generations)") +
        theme(legend.position="bottom") +
        geom_point(alpha = 0.2, size = 0.5) +
        geom_abline(data = lm.df,
                    aes(intercept = resilience.yintercept,
                        slope = resilience.slope,
                        color = Treatment)) +
        scale_color_manual(values = c("#D55E00", "#F0E442", "#0072B2"))
    
    if (!plot.legend) {
        p <- p + guides(color=FALSE)
    }
    return(p)
}


make.big.Cong.resilience.plot <- function(big.resilience.df) {
    p <- make.big.resilience.plot(big.resilience.df) +
        ## IMPORTANT: axes need to be consistent when comparing plots.
        ylim(0.14, 0.18) +
        ggtitle("Cong PPI dataset")
    return(p)
}


make.big.Zitnik.resilience.plot <- function(big.resilience.df) {
    p <- make.big.resilience.plot(big.resilience.df) +
        ## IMPORTANT: axes need to be consistent when comparing plots.
        ylim(0.394, 0.414) +
        ggtitle("Zitnik PPI dataset")
    return(p)

}


calc.resilience.regression.df <- function(resilience.df) {
    ## returns a dataframe with the linear regression coefficients
    ## per population as a column.

    ## get the average ancestral resilience.
    REL606.df <- resilience.df %>%
        filter(Generation == 0)
    ancestral.resilience <- mean(unique(REL606.df$resilience))
    
    calc.resilience.slope.per.pop.df <- function(pop.df) {
        ## fix the y-intercept at the ancestral REL606 resilience.
        my.regression <- lm(data=pop.df, I(resilience - ancestral.resilience) ~ 0 + Generation)
        resilience.coefficient <- my.regression$coefficients[[1]]
        pop.df.with.slope <- pop.df %>%
            select(population, Treatment) %>%
            mutate(resilience.yintercept = ancestral.resilience) %>%
            mutate(resilience.slope = resilience.coefficient) %>%
            distinct()
        return(pop.df.with.slope)
    }

    summary.df <- resilience.df %>%
        split(list(.$population,.$Treatment)) %>%
        map_dfr(.f = calc.resilience.slope.per.pop.df)
    return(summary.df)
}


regression.slope.test <- function(full.slope.df, null.comp.string) {

    if (null.comp.string == "randomized_across") {
        null.df <- full.slope.df %>%
            filter(Treatment == "randomized_across")
    } else if (null.comp.string == "randomized over LTEE") {
        null.df <- full.slope.df %>%
            filter(Treatment == "randomized over LTEE")
    } else if (null.comp.string == "randomized over genome") {
        null.df <- full.slope.df %>%
            filter(Treatment == "randomized over genome")
    } else {
        print("ERROR: null.comp.string is not recognized.")
        return(NULL)
    }

    data.df <- full.slope.df %>%
        filter(Treatment == "actual data")

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


calc.REL606.PPI.gene.distances <- function(node1.df, node2.df) {
    ## Calculate the spatial distribution of edges to genomic distances
    ## in the REL606 chromosome.
    ## This is for comparing the spatial distribution of KO'ed edges to
    ## genomic distances, with and without large deletions.

    ## helper function that calculates the distance between two genes on
    ## the REL606 chromosome.
    calc.PPI.gene.distance <- function(gene1.start, gene2.start) {
        
        if (gene1.start < gene2.start) {
            x1 <- gene1.start
            x2 <- gene2.start
        } else {
            x1 <- gene2.start
            x2 <- gene1.start
        }

        d1 <- x2 - x1
        d2 <- REL606.CHR.LENGTH - x2 + x1
        return(min(d1, d2))
    }
    
    PPI.dist.vec <- map2_dbl(.x = node1.df$start, .y = node2.df$start,
                         .f = calc.PPI.gene.distance)
    PPI.dist.df <- data.frame(PPI.dist=PPI.dist.vec)
   return(PPI.dist.df) 
}


get.KOed.edges <- function(KO.gene.vec, node1.df, node2.df) {
    ## returns a vector for edges that are knocked out,
    ## based on whether the nodes for the link are in KO.gene.vec.
    
    node1.vec <- node1.df$Gene
    node2.vec <- node2.df$Gene
    
    boolvec1 <- sapply(node1.vec, function(x) ifelse(x %in% KO.gene.vec,TRUE,FALSE))
    boolvec2 <- sapply(node2.vec, function(x) ifelse(x %in% KO.gene.vec,TRUE,FALSE))
    or.boolvec <- boolvec1 | boolvec2
    return(or.boolvec)
}


pop.to.KO.edge.distances <- function(cur.pop, clone.to.KOed.genes.list, node1.df, node2.df) {

    KO.gene.vec <- clone.to.KOed.genes.list[[cur.pop]]
    KO.vec <- get.KOed.edges(KO.gene.vec, node1.df, node2.df)

    KOed.node1.df <- node1.df %>%
        mutate(KOed.edge = KO.vec) %>%
        filter(KOed.edge == TRUE)
    KOed.node2.df <- node2.df %>%
        mutate(KOed.edge = KO.vec) %>%
        filter(KOed.edge == TRUE)
    
    KOed.edge.distance.df <- calc.REL606.PPI.gene.distances(KOed.node1.df, KOed.node2.df)
    return(KOed.edge.distance.df)
}


pop.to.KO.degrees <- function(cur.pop, clone.to.KOed.genes.list, degree.df) {
    ## calculate degree distribution for KO'ed genes of various stripes.
    KO.gene.vec <- clone.to.KOed.genes.list[[cur.pop]]
    degree.df %>%
        filter(Gene %in% KO.gene.vec)
}


##########################################################################
## IMPORT METADATA.

## Order nonmutator pops, then hypermutator pops by converting Population to
## factor type and setting the levels.
nonmutator.pops <- c("Ara-5", "Ara-6", "Ara+1", "Ara+2", "Ara+4", "Ara+5")
hypermutator.pops <- c("Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara+3", "Ara+6")

LTEE.pop.vec <- c(nonmutator.pops, hypermutator.pops)

## we will add REL606 as ancestor to all 12 populations.
REL606.metadata <- data.frame(strain = rep('REL606', 12),
                              population = LTEE.pop.vec,
                              time = rep(0, 12),
                              clone = rep('A', 12),
                              mutator_status = rep('non-mutator', 12))

LTEE.genomes.KO.muts <- read.csv(
    "../data/LTEE-264-genomes-SNP-nonsense-small-indel-MOB-large-deletions.csv") %>%
    ## filter out intergenic mutations.
    filter(!str_detect(gene_position, "intergenic"))

LTEE.genomes.KO.metadata <- LTEE.genomes.KO.muts %>%
    select(population,time,strain,clone,mutator_status) %>%
    distinct() %>%
    ## let's add REL606.
    full_join(REL606.metadata) %>%
    mutate(Generation=time/10000) %>%
    ## This for changing the ordering of populations in plots.
    mutate(population=factor(population,levels = LTEE.pop.vec))

## import metadata for genes in the REL606 genomes.
REL606.genes <- read.csv("../results/REL606_IDs.csv")

## make a data structure for genes affected by KO mutations in 50K genomes.
LTEE.50K.KO.data <- LTEE.genomes.KO.muts %>%
    ## remove intergenic mutations.
    filter(!str_detect(gene_position, "intergenic")) %>%
    filter(time == 50000) %>%
    select(population, strain, clone, gene_list)

LTEE.50K.A.clone.KO.data <- LTEE.50K.KO.data %>%
    filter(clone == 'A')

KOed.genes.in.LTEE.50K.A.clones <- make.list.of.strain.to.KOed.genes(
    LTEE.50K.A.clone.KO.data)


#######################################################################
## PPI network resilience results, using Duke Compute Cluster runs.

big.resilience.df <- read.csv("../results/resilience/collated-resilience-runs.csv") %>%
    full_join(LTEE.genomes.KO.metadata) %>%
    set.no.KO.strains.resilience.to.REL606() %>%
    mutate(log.resilience=log(resilience)) %>%
    ## remove the runs for the set of genes KO'ed across the LTEE populations.
    filter(run_type != "across_pops_randomized_resilience") %>%
    add_Treatment_column()

big.cong.resilience.df <- big.resilience.df %>%
    filter(dataset=="Cong")

big.zitnik.resilience.df <- big.resilience.df %>%
    filter(dataset=="Zitnik")

big.cong.resilience.plot <- make.big.Cong.resilience.plot(big.cong.resilience.df)
big.zitnik.resilience.plot <- make.big.Zitnik.resilience.plot(big.zitnik.resilience.df)

Fig2.legend <- get_legend(make.big.resilience.plot(big.cong.resilience.df,
                                                   plot.legend=TRUE))

Fig2 <- plot_grid(plot_grid(big.zitnik.resilience.plot,big.cong.resilience.plot,
                            labels=c('A','B',NULL),nrow=1, rel_widths=c(1,1,0.4)),
                  Fig2.legend,ncol=1, rel_heights=c(1,0.05))
ggsave("../results/resilience/figures/Fig2.pdf",height=7,width=6)


## calculate the slopes of the resilience regressions.
big.cong.slope.df <- calc.resilience.regression.df(big.cong.resilience.df)
big.zitnik.slope.df <- calc.resilience.regression.df(big.zitnik.resilience.df)

cong.slopes.from.data <- big.cong.slope.df %>%
    filter(Treatment == "actual data")

zitnik.slopes.from.data <- big.zitnik.slope.df %>%
    filter(Treatment == "actual data")

## Calculate statistics on differences between slopes.
regression.slope.test(big.zitnik.slope.df, "randomized over genome")
regression.slope.test(big.cong.slope.df, "randomized over genome")

regression.slope.test(big.zitnik.slope.df, "randomized over LTEE")
regression.slope.test(big.cong.slope.df, "randomized over LTEE")


########################################################################
## Analysis to see whether large deletions systematically bias
## the result reported in Figure 1.
########################################################################

## get edges in Zitnik PPI network.
raw.zitnik.edges <- read.csv("../results/thermostability/Ecoli-Zitnik-data/Ecoli-treeoflife.interactomes/511145.txt", sep=" ", header=FALSE) %>%
    mutate(edge.number = row_number())

zitnik.node1 <- raw.zitnik.edges %>% 
    mutate(blattner = V1) %>%
    select(blattner, edge.number) %>%
    left_join(REL606.genes) %>%
    ## remove edges with NA values.
    filter(complete.cases(.))

zitnik.node2 <- raw.zitnik.edges %>% 
    mutate(blattner = V2) %>%
    select(blattner, edge.number) %>%
    left_join(REL606.genes) %>%
    ## remove edges with NA values.
    filter(complete.cases(.))

## filter for edges between nodes that have not been removed.
zitnik.node1 <- zitnik.node1 %>%
    filter(edge.number %in% zitnik.node2$edge.number) %>%
    select(Gene, start, edge.number)

zitnik.node2 <- zitnik.node2 %>%
    filter(edge.number %in% zitnik.node1$edge.number) %>%
    select(Gene, start, edge.number)

#############################

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
    select(Gene, start, edge.number)

cong.node2 <- cong.node2 %>%
    filter(edge.number %in% cong.node1$edge.number) %>%
    select(Gene, start, edge.number)

#############################

## Dataframes of the distance between interacting genes.
cong.dist.df <- calc.REL606.PPI.gene.distances(cong.node1, cong.node2)
zitnik.dist.df <- calc.REL606.PPI.gene.distances(zitnik.node1, zitnik.node2)

## Calculate the genomic distance for KO'ed interactions across all 12 LTEE pops.

## use partial function application so that cur.pop is the only free parameter.
pop.to.clone.A.cong.PPI.gene.dists <- partial(
    .f = pop.to.KO.edge.distances,
    clone.to.KOed.genes.list = KOed.genes.in.LTEE.50K.A.clones,
    node1.df = cong.node1,
    node2.df = cong.node2)

pop.to.clone.A.zitnik.PPI.gene.dists <- partial(
    .f = pop.to.KO.edge.distances,
    clone.to.KOed.genes.list = KOed.genes.in.LTEE.50K.A.clones,
    node1.df = zitnik.node1,
    node2.df = zitnik.node2)

## IMPORTANT NOTE: These distributions involve
## parallel evolution. So edges that are knocked out in multiple populations
## are counted multiple times.

KOed.50K.clone.A.cong.PPI.dists <- map_dfr(
    .x = LTEE.pop.vec,
    .f = pop.to.clone.A.cong.PPI.gene.dists)

KOed.50K.clone.A.zitnik.PPI.dists <- map_dfr(
    .x = LTEE.pop.vec,
    .f = pop.to.clone.A.zitnik.PPI.gene.dists)

#######################################################
## redo, now omitting large deletions from the picture.

NoDel.LTEE.50K.A.clone.KO.data <- LTEE.genomes.KO.muts %>%
    ## remove intergenic mutations.
    filter(!str_detect(gene_position, "intergenic")) %>%
    filter(time == 50000) %>%
    select(population, strain, clone, gene_list) %>%
    filter(clone == 'A') %>%
    ## filter out multigene deletions.
    filter(!str_detect(gene_list, ','))

OnlyDel.LTEE.50K.A.clone.KO.data <- LTEE.genomes.KO.muts %>%
    ## remove intergenic mutations.
    filter(!str_detect(gene_position, "intergenic")) %>%
    filter(time == 50000) %>%
    select(population, strain, clone, gene_list) %>%
    filter(clone == 'A') %>%
    ## only multigene deletions.
    filter(str_detect(gene_list, ','))


NoDel.KOed.genes.in.LTEE.50K.A.clones <- make.list.of.strain.to.KOed.genes(
    NoDel.LTEE.50K.A.clone.KO.data)


OnlyDel.KOed.genes.in.LTEE.50K.A.clones <- make.list.of.strain.to.KOed.genes(
    OnlyDel.LTEE.50K.A.clone.KO.data)


## use partial function application so that cur.pop is the only free parameter.
NoDel.pop.to.clone.A.cong.PPI.gene.dists <- partial(.f = pop.to.KO.edge.distances,
                         clone.to.KOed.genes.list = NoDel.KOed.genes.in.LTEE.50K.A.clones,
                         node1.df = cong.node1,
                         node2.df = cong.node2)

OnlyDel.pop.to.clone.A.cong.PPI.gene.dists <- partial(.f = pop.to.KO.edge.distances,
                         clone.to.KOed.genes.list = OnlyDel.KOed.genes.in.LTEE.50K.A.clones,
                         node1.df = cong.node1,
                         node2.df = cong.node2)

NoDel.pop.to.clone.A.zitnik.PPI.gene.dists <- partial(.f = pop.to.KO.edge.distances,
                         clone.to.KOed.genes.list = NoDel.KOed.genes.in.LTEE.50K.A.clones,
                         node1.df = zitnik.node1,
                         node2.df = zitnik.node2)

OnlyDel.pop.to.clone.A.zitnik.PPI.gene.dists <- partial(.f = pop.to.KO.edge.distances,
                         clone.to.KOed.genes.list = OnlyDel.KOed.genes.in.LTEE.50K.A.clones,
                         node1.df = zitnik.node1,
                         node2.df = zitnik.node2)


NoDel.KOed.50K.clone.A.cong.PPI.dists <- map_dfr(
    .x = LTEE.pop.vec,
    .f = NoDel.pop.to.clone.A.cong.PPI.gene.dists)

OnlyDel.KOed.50K.clone.A.cong.PPI.dists <- map_dfr(
    .x = LTEE.pop.vec,
    .f = OnlyDel.pop.to.clone.A.cong.PPI.gene.dists)


NoDel.KOed.50K.clone.A.zitnik.PPI.dists <- map_dfr(
    .x = LTEE.pop.vec,
    .f = NoDel.pop.to.clone.A.zitnik.PPI.gene.dists)

OnlyDel.KOed.50K.clone.A.zitnik.PPI.dists <- map_dfr(
    .x = LTEE.pop.vec,
    .f = OnlyDel.pop.to.clone.A.zitnik.PPI.gene.dists)


## surprising: the patterns are actually in the opposite direction than I expected!
wilcox.test(NoDel.KOed.50K.clone.A.zitnik.PPI.dists$PPI.dist,
            OnlyDel.KOed.50K.clone.A.zitnik.PPI.dists$PPI.dist)

mean(OnlyDel.KOed.50K.clone.A.zitnik.PPI.dists$PPI.dist)
mean(NoDel.KOed.50K.clone.A.zitnik.PPI.dists$PPI.dist)

wilcox.test(NoDel.KOed.50K.clone.A.cong.PPI.dists$PPI.dist,
            OnlyDel.KOed.50K.clone.A.cong.PPI.dists$PPI.dist)

mean(OnlyDel.KOed.50K.clone.A.cong.PPI.dists$PPI.dist)
mean(NoDel.KOed.50K.clone.A.cong.PPI.dists$PPI.dist)


########################
## Let's examine the degree distribution of KO'ed genes,
## with and without large deletions.

zitnik.degree.df <- read.csv("../results/thermostability/Zitnik_network_statistics.csv")
cong.degree.df <- read.csv("../results/thermostability/Cong_network_statistics.csv")


OnlyDel.pop.to.clone.A.cong.PPI.degrees <- partial(.f = pop.to.KO.degrees,
                         clone.to.KOed.genes.list = OnlyDel.KOed.genes.in.LTEE.50K.A.clones,
                         degree.df = cong.degree.df)

NoDel.pop.to.clone.A.cong.PPI.degrees <- partial(.f = pop.to.KO.degrees,
                         clone.to.KOed.genes.list = NoDel.KOed.genes.in.LTEE.50K.A.clones,
                         degree.df = cong.degree.df)


OnlyDel.pop.to.clone.A.zitnik.PPI.degrees <- partial(.f = pop.to.KO.degrees,
                         clone.to.KOed.genes.list = OnlyDel.KOed.genes.in.LTEE.50K.A.clones,
                         degree.df = zitnik.degree.df)

NoDel.pop.to.clone.A.zitnik.PPI.degrees <- partial(.f = pop.to.KO.degrees,
                         clone.to.KOed.genes.list = NoDel.KOed.genes.in.LTEE.50K.A.clones,
                         degree.df = zitnik.degree.df)


OnlyDel.KOed.50K.clone.A.cong.PPI.degrees <- map_dfr(
    .x = LTEE.pop.vec,
    .f = OnlyDel.pop.to.clone.A.cong.PPI.degrees) 


NoDel.KOed.50K.clone.A.cong.PPI.degrees <- map_dfr(
    .x = LTEE.pop.vec,
    .f = NoDel.pop.to.clone.A.cong.PPI.degrees)


OnlyDel.KOed.50K.clone.A.zitnik.PPI.degrees <- map_dfr(
    .x = LTEE.pop.vec,
    .f = OnlyDel.pop.to.clone.A.zitnik.PPI.degrees)

NoDel.KOed.50K.clone.A.zitnik.PPI.degrees <- map_dfr(
    .x = LTEE.pop.vec,
    .f = NoDel.pop.to.clone.A.zitnik.PPI.degrees)

## no significant different in degree between KO'ed genes within and outside of large
## deletions in the LTEE.
wilcox.test(OnlyDel.KOed.50K.clone.A.zitnik.PPI.degrees$Degree,
            NoDel.KOed.50K.clone.A.zitnik.PPI.degrees$Degree)

mean(OnlyDel.KOed.50K.clone.A.zitnik.PPI.degrees$Degree)
mean(NoDel.KOed.50K.clone.A.zitnik.PPI.degrees$Degree)


wilcox.test(OnlyDel.KOed.50K.clone.A.cong.PPI.degrees$Degree,
            NoDel.KOed.50K.clone.A.cong.PPI.degrees$Degree)

mean(OnlyDel.KOed.50K.clone.A.cong.PPI.degrees$Degree)
mean(NoDel.KOed.50K.clone.A.cong.PPI.degrees$Degree)

########################################################
## Compare population fitness (Wiser et al. 2013)
## to the resilience of genomes isolated at those time points.
## Overall, I find this result inconclusive.
## measuring the actual fitnesses of the relevant clones
## would give a better picture.


wiser.data <- read.csv("../data/Concatenated.LTEE.data.all.csv") %>%
    ## get time units to match the resilience dataframes.
    mutate(Generation = Generation/10000) %>%
    ## make population column name match the resilience dataframes.
    mutate(population = Population) %>%
    select(Generation, Red.Pop, White.Pop, population, Rep, Fitness, Complete) %>%
    group_by(Generation, Red.Pop, White.Pop, population) %>%
    ## calculate the average fitness over replicates
    summarize(mean.Fitness = mean(Fitness)) %>%
    ungroup() %>%
    select(Generation, population, mean.Fitness)


cong.resilience.with.fitness.df <- big.cong.resilience.df %>%
    filter(run_type == "LTEE_genome_resilience") %>%
    filter(Treatment == "actual data") %>%
    select(dataset, replicate, strain, resilience,
           population, Generation) %>%
    group_by(dataset, strain, population, Generation) %>%
    summarize(mean.resilience = mean(resilience)) %>%
    ungroup() %>%
    left_join(wiser.data) %>%
    filter(complete.cases(.))

zitnik.resilience.with.fitness.df <- big.zitnik.resilience.df %>%
    filter(run_type == "LTEE_genome_resilience") %>%
    filter(Treatment == "actual data") %>%
    select(dataset, replicate, strain, resilience,
           population, Generation) %>%
    group_by(dataset, strain, population, Generation) %>%
    summarize(mean.resilience = mean(resilience)) %>%
    ungroup() %>%
    left_join(wiser.data) %>%
    filter(complete.cases(.))


cong.fitness.plot.with.legend <- ggplot(cong.resilience.with.fitness.df,
                            aes(x=mean.Fitness,y=mean.resilience,color=population)) +
    geom_point() +
    theme_classic() +
    xlab("Mean population fitness") +
    ylab("PPI network resilience") +
    theme(legend.position="bottom")

cong.fitness.plot <- cong.fitness.plot.with.legend + guides(color=FALSE)
fitness.resilience.legend <- get_legend(cong.fitness.plot.with.legend)

zitnik.fitness.plot <- ggplot(zitnik.resilience.with.fitness.df,
                            aes(x=mean.Fitness,y=mean.resilience,color=population)) +
    geom_point() +
    theme_classic() +
    xlab("Mean population fitness") +
    ylab("PPI network resilience") +
    guides(color=FALSE)

cong.fitness.curve <- ggplot(cong.resilience.with.fitness.df,
                             aes(x = Generation,
                                 y = mean.Fitness,
                                 color = mean.resilience)) +
    geom_point() +
    theme_classic() +
    xlab("Time (x 10,000 generations)") +
    ylab("Mean population fitness") +
    scale_color_viridis_c(name = "resilience", option = "plasma", direction = -1) +
    ggtitle("Cong PPI dataset")

zitnik.fitness.curve <- ggplot(zitnik.resilience.with.fitness.df,
                               aes(x = Generation,
                                   y = mean.Fitness,
                                   color = mean.resilience)) +
    geom_point() +
    theme_classic() +
    xlab("Time (x 10,000 generations)") +
    ylab("Mean population fitness") +
    scale_color_viridis_c(name = "resilience", option = "plasma", direction = -1) +
    ggtitle("Zitnik PPI dataset")

Fig3 <- plot_grid(
    plot_grid(
        zitnik.fitness.curve,
        cong.fitness.curve,
        zitnik.fitness.plot,
        cong.fitness.plot,
        labels=c('A','B','',''),
        nrow=2),
    fitness.resilience.legend,
    ncol=1,
    rel_heights = c(1,0.15))
ggsave("../results/resilience/figures/Fig3.pdf", Fig3, height=6, width=7)

## calculate correlation coefficients and p-values.
cor.test(zitnik.resilience.with.fitness.df$mean.resilience,
         zitnik.resilience.with.fitness.df$mean.Fitness, method = "spearman")

cor.test(cong.resilience.with.fitness.df$mean.resilience,
         cong.resilience.with.fitness.df$mean.Fitness, method = "spearman")

#########################################################################################

## Analysis of essential genes and PPI network resilience.

## 1) take essential genes in REL606 from Couce paper,
## and count the number of knockout mutations in those genes
## in the 50K clone A genomes. do a one-sided binomial test
## for purifying selection.

## 2) compare the PPI degree for essential genes,
## compared to the remainder of the genes in the genome.

## 3) calculate how resilience changes for all single gene
## knockouts in the REL606 genome. make a histogram of resilience
## for the genes, and color essential genes in red, and do a wilcox
## test on a difference between the two distributions.

#########################################################################################

## data structures for KO mutations in the LTEE genomes.
## LTEE.genomes.KO.muts
## LTEE.genomes.KO.metadata


## Get essential and near-essential genes reported in
## Supplementary Table 1 of Couce et al. 2017.
## I manually fixed the names of a couple genes in this dataset.
## The original names are in the "Name" column, and updated names
## are in the "Gene" column.
essential.genes <- read.csv("../data/Couce2017-LTEE-essential.csv") %>%
    inner_join(REL606.genes) %>% filter(!(is.na(locus_tag)))

nonessential.genes <- REL606.genes %>%
    filter(!(Gene %in% essential.genes$Gene))

essential.gene.regex <- reduce(essential.genes$Gene, .f = partial(paste, sep = '|'))

## 44 KO mutations affected essential genes in 50K clone A genomes.
essential.50K.A.clone.KO.data <- LTEE.50K.A.clone.KO.data %>%
    filter(str_detect(gene_list, essential.gene.regex))

## 897 KO mutations did not affect any essential genes in 50K clone A genomes.
nonessential.50K.A.clone.KO.data <- LTEE.50K.A.clone.KO.data %>%
    filter(!str_detect(gene_list, essential.gene.regex))

essential.length <- sum(essential.genes$gene_length) ## 499180 bp
nonessential.length <- sum(nonessential.genes$gene_length) ## 3462963 bp.
essential.target.prob <- essential.length/(essential.length + nonessential.length)

## 1) essential genes are rarely disrupted in the 50K genomes (p = 1.27e-16).
essential.binom.test <- binom.test(x = nrow(essential.50K.A.clone.KO.data),
                                   n = nrow(essential.50K.A.clone.KO.data) +
                                       nrow(nonessential.50K.A.clone.KO.data),
                                   p = essential.target.prob,
                                   alternative="less")

## Nevertheless, a significant proportion of genes under positive selection
## in the LTEE are essential genes, as reported in Maddamsetti et al. (2017).

## get mutation parallelism in the LTEE genomes published in Tenaillon et al. (2016).
## for comparison to essential genes.
## order by positive G scores.
nonmut.genomics <- read.csv('../data/tenaillon2016-nonmutator-parallelism.csv') %>%
    ## make sure these genes pass the filters on REL606.genes.
    filter(Gene.name %in% REL606.genes$Gene) %>%
    mutate(isEssential = Gene.name %in% essential.genes$Gene) %>%
    arrange(desc(G.score))

## genes with negative G-scores are the most interesting--
## could indicate purifying selection.
hypermut.genomics <- read.csv('../data/tenaillon2016-mutator-parallelism.csv') %>%
    ## make sure these genes pass the filters on REL606.genes.
    filter(Gene.name %in% REL606.genes$Gene) %>%
    mutate(isEssential = Gene.name %in% essential.genes$Gene) %>%
    arrange(G.score)

## These are the 54 genes with highest G-score and 2 or more observed nonsynonymous
## mutations
top.nonmut.genomics <- slice_max(nonmut.genomics, n = 54, order_by = G.score)
## no criterion here, so just take 50.
top.hypermut.genomics <- slice_max(hypermut.genomics, n = 50, order_by = G.score)

## 22 out of 54 top non-mut genes are essential.
nonmut.top.hit.essential <- essential.genes %>%
    filter(Gene %in% top.nonmut.genomics$Gene.name)
## calculate the probability of this outcome by chance (hypergeometric distribution),
## using a 2x2 contingency table and Fisher's exact test.
## 541 essential and nearly-essential genes, out of 4112 in REL606.genes.
## the p-value (tail probability).
fisher.test(matrix(c(22,541-22,54-22,4112-541-54+22),2))


## what about the hypermutators? 3 out of 50 top hypermut genes.
hypermut.top.hit.essential <- essential.genes %>%
    filter(Gene %in% top.hypermut.genomics$Gene.name)

## what about the bottom 50 genes?
bottom.nonmut.genomics <- slice_min(nonmut.genomics, n = 50, order_by = G.score)
bottom.hypermut.genomics <- slice_min(hypermut.genomics, n = 50, order_by = G.score)

nonmut.bottom.hit.essential <- essential.genes %>%
    filter(Gene %in% bottom.nonmut.genomics$Gene.name)
## what about the hypermutators? 3 out of 50 top hypermut genes.
hypermut.bottom.hit.essential <- essential.genes %>%
    filter(Gene %in% bottom.hypermut.genomics$Gene.name)

hypermut.negative.G <- hypermut.genomics %>% filter(G.score < 0)


## make a plot or report the G-score distribution for essential genes?
## interesting possibility-- look at genes with negative G scores in the
## hypermutators.

nonmut.G.essentiality.plot <- ggplot(nonmut.genomics,
                                     aes(x = isEssential,
                                         y = G.score)) +
    geom_violin() + theme_classic()

hypermut.G.essentiality.plot <- ggplot(hypermut.genomics,
                                     aes(x = isEssential,
                                         y = G.score)) +
    geom_violin() + theme_classic()

hypermut.negativeG.essentiality.plot <- ggplot(hypermut.negative.G,
                                     aes(x = isEssential,
                                         y = G.score)) +
    geom_boxplot() + theme_classic()


## 2) compare the PPI degree for essential genes,
## compared to the remainder of the genes in the genome.
essentiality.zitnik.degree.df <- zitnik.degree.df %>%
    mutate(isEssential = Gene %in% essential.genes$Gene)
essentiality.cong.degree.df <- cong.degree.df %>%
    mutate(isEssential = Gene %in% essential.genes$Gene)

zitnik.essentiality.plot <- ggplot(essentiality.zitnik.degree.df,
                                   aes(x = isEssential, y = Degree)) +
    geom_boxplot() + theme_classic()

cong.essentiality.plot <- ggplot(essentiality.cong.degree.df,
                                   aes(x = isEssential, y = Degree)) +
    geom_boxplot() + theme_classic()

wilcox.test(filter(essentiality.zitnik.degree.df,isEssential==TRUE)$Degree,
            filter(essentiality.zitnik.degree.df,isEssential==FALSE)$Degree,
            alternative="greater")$p.value

wilcox.test(filter(essentiality.cong.degree.df,isEssential==TRUE)$Degree,
            filter(essentiality.cong.degree.df,isEssential==FALSE)$Degree,
            alternative="greater")$p.value

## 3) calculate how resilience changes for all single gene
## knockouts in the REL606 genome. make a histogram of resilience
## for the genes, and do a wilcoxon test on a difference between the two distributions.

zitnik.single.KO.resilience.data <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_single_KO_resilience_Rep999.csv")

cong.single.KO.resilience.data <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_single_KO_resilience_Rep999.csv")


single.KO.analysis <- function(zitnik.single.KO.data, cong.single.KO.data,
                               base.strain, outf) {
    
    zitnik.base.resilience <- zitnik.single.KO.data %>% filter(strain == base.strain)
    cong.base.resilience <- cong.single.KO.data %>% filter(strain == base.strain)

    zitnik.single.KO.resilience.df <- zitnik.single.KO.data %>%
        filter(strain != base.strain) %>%
        mutate(Gene = str_split_fixed(strain, "_", n = 2)[,1]) %>%
        mutate(isEssentialKO = Gene %in% essential.genes$Gene) %>%
        mutate(Essentiality = ifelse(isEssentialKO, "Essential gene disruption", "Non-essential gene disruption")) %>%
        mutate(Dataset = "Zitnik")
    
    cong.single.KO.resilience.df <- cong.single.KO.data %>%
        filter(strain != base.strain) %>%
        mutate(Gene = str_split_fixed(strain, "_", n = 2)[,1]) %>%
        mutate(isEssentialKO = Gene %in% essential.genes$Gene) %>%
        mutate(Essentiality = ifelse(isEssentialKO, "Essential gene disruption", "Non-essential gene disruption")) %>%
        mutate(Dataset = "Cong")
    
    title <- ggdraw() + 
        draw_label(
            paste0("Single-gene disruptions of the ", base.strain,
                     " PPI network"),
            fontface = 'bold',
            x = 0,
            hjust = 0
        ) +
        theme(
            ## add margin on the left of the drawing canvas,
            ## so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
        )
    
    panelA <- ggplot(zitnik.single.KO.resilience.df,
                     aes(x = Essentiality, y = resilience)) +
        theme_classic() +
        xlab("") + ## remove x-axis label.
        ylab("PPI network resilience") +
        geom_jitter(alpha = 0.2, height = 0) +
        geom_hline(yintercept = zitnik.base.resilience$resilience,
                   color = 'red', linetype = "dashed") +
        ggtitle("Zitnik PPI dataset")
    
    panelB <- ggplot(cong.single.KO.resilience.df,
                     aes(x = Essentiality,
                         y = resilience)) +
        xlab("") + ## remove x-axis label.
        ylab("PPI network resilience") +
        theme_classic() +
        geom_jitter(alpha = 0.2, height = 0) +
        geom_hline(yintercept = cong.base.resilience$resilience,
                   color = 'red', linetype = "dashed") +
        ggtitle("Cong PPI dataset")
    
    panels <- plot_grid(panelA, panelB, ncol = 1)
    
    Fig <- plot_grid(
        title, panels,
        ncol = 1,
        ## rel_heights values control vertical title margins
        rel_heights = c(0.1, 1)
    )
    
    ggsave(outf, Fig, height = 6, width = 5)
    

    ## For REL606:
    ## single-gene knockouts are much more likely to INCREASE network resilience than
    ## decrease network resilience. In addition, a single knockout of an essential genes
    ## has on average a stronger positive effect on network resilience than knocking out
    ## a non-essential gene. This pattern holds for both datasets.

    ## The patterns are idiosyncratic and not always consistent across datasets for the
    ## 50,000 generations LTEE clones.
    print("mean of Zitnik single essential KO resilience")
    print(mean(filter(zitnik.single.KO.resilience.df, isEssentialKO == TRUE)$resilience))
    print("mean of Zitnik single non-essential KO resilience")
    print(mean(filter(zitnik.single.KO.resilience.df, isEssentialKO == FALSE)$resilience))
    print(wilcox.test(x = filter(zitnik.single.KO.resilience.df, isEssentialKO == TRUE)$resilience, y = filter(zitnik.single.KO.resilience.df, isEssentialKO == FALSE)$resilience))

    print("mean of Cong single essential KO resilience")
    print(mean(filter(cong.single.KO.resilience.df, isEssentialKO == TRUE)$resilience))
    print("mean of Cong single non-essential KO resilience")
    print(mean(filter(cong.single.KO.resilience.df, isEssentialKO == FALSE)$resilience))
    print(
        wilcox.test(x = filter(cong.single.KO.resilience.df, isEssentialKO == TRUE)$resilience, y = filter(cong.single.KO.resilience.df, isEssentialKO == FALSE)$resilience))
    
}

## run single KO analysis of REL606 PPI network.
single.KO.analysis(zitnik.single.KO.resilience.data,
                   cong.single.KO.resilience.data,
                   "REL606",
                   "../results/resilience/figures/Fig4.pdf")

##################################################
## Supplementary Figure S1: repeat this analysis, on the 50K clone A genomes.
## do we see the same trend, or is it reversed (i.e. essential genes tend to reduce
## resilience when removed).
##################################################

## Ara+1 REL11392 PPI network.
zitnik.single.KO.Ara.plus1.data <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_single_KO_resilience_Ara+1_Rep999.csv") %>%
    mutate(population = "Ara+1") %>%
    mutate(starting_clone = "REL11392")

cong.single.KO.Ara.plus1.data <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_single_KO_resilience_Ara+1_Rep999.csv") %>%
    mutate(population = "Ara+1") %>%
    mutate(starting_clone = "REL11392")

## Ara+2 REL11342 PPI network.
zitnik.single.KO.Ara.plus2.data <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_single_KO_resilience_Ara+2_Rep999.csv") %>%
    mutate(population = "Ara+2") %>%
    mutate(starting_clone = "REL11342")

cong.single.KO.Ara.plus2.data <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_single_KO_resilience_Ara+2_Rep999.csv") %>%
    mutate(population = "Ara+2") %>%
    mutate(starting_clone = "REL11342")

## Ara+3 REL11345 PPI network.
zitnik.single.KO.Ara.plus3.data <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_single_KO_resilience_Ara+3_Rep999.csv") %>%
    mutate(population = "Ara+3") %>%
    mutate(starting_clone = "REL11345")

cong.single.KO.Ara.plus3.data <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_single_KO_resilience_Ara+3_Rep999.csv") %>%
    mutate(population = "Ara+3") %>%
    mutate(starting_clone = "REL11345")

## Ara+4 REL11348 PPI network.
zitnik.single.KO.Ara.plus4.data <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_single_KO_resilience_Ara+4_Rep999.csv") %>%
    mutate(population = "Ara+4") %>%
    mutate(starting_clone = "REL11348")

cong.single.KO.Ara.plus4.data <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_single_KO_resilience_Ara+4_Rep999.csv") %>%
    mutate(population = "Ara+4") %>%
    mutate(starting_clone = "REL11348")

## Ara+5 REL11367 PPI network.
zitnik.single.KO.Ara.plus5.data <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_single_KO_resilience_Ara+5_Rep999.csv") %>%
    mutate(population = "Ara+5") %>%
    mutate(starting_clone = "REL11367")

cong.single.KO.Ara.plus5.data <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_single_KO_resilience_Ara+5_Rep999.csv") %>%
    mutate(population = "Ara+5") %>%
    mutate(starting_clone = "REL11367")

## Ara+6 REL11370 PPI network.
zitnik.single.KO.Ara.plus6.data <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_single_KO_resilience_Ara+6_Rep999.csv") %>%
    mutate(population = "Ara+6") %>%
    mutate(starting_clone = "REL11370")

cong.single.KO.Ara.plus6.data <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_single_KO_resilience_Ara+6_Rep999.csv") %>%
    mutate(population = "Ara+6") %>%
    mutate(starting_clone = "REL11370")

## Ara-1 REL11330 PPI network.
zitnik.single.KO.Ara.minus1.data <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_single_KO_resilience_Ara-1_Rep999.csv") %>%
    mutate(population = "Ara-1") %>%
    mutate(starting_clone = "REL11330")

cong.single.KO.Ara.minus1.data <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_single_KO_resilience_Ara-1_Rep999.csv") %>%
    mutate(population = "Ara-1") %>%
    mutate(starting_clone = "REL11330")

## Ara-2 REL11333 PPI network.
zitnik.single.KO.Ara.minus2.data <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_single_KO_resilience_Ara-2_Rep999.csv") %>%
    mutate(population = "Ara-2") %>%
    mutate(starting_clone = "REL11333")

cong.single.KO.Ara.minus2.data <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_single_KO_resilience_Ara-2_Rep999.csv") %>%
    mutate(population = "Ara-2") %>%
    mutate(starting_clone = "REL11333")

## Ara-3 REL11364 PPI network.
zitnik.single.KO.Ara.minus3.data <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_single_KO_resilience_Ara-3_Rep999.csv") %>%
    mutate(population = "Ara-3") %>%
    mutate(starting_clone = "REL11364")

cong.single.KO.Ara.minus3.data <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_single_KO_resilience_Ara-3_Rep999.csv") %>%
    mutate(population = "Ara-3") %>%
    mutate(starting_clone = "REL11364")

## Ara-4 REL11336 PPI network.
zitnik.single.KO.Ara.minus4.data <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_single_KO_resilience_Ara-4_Rep999.csv") %>%
    mutate(population = "Ara-4") %>%
    mutate(starting_clone = "REL11336")

cong.single.KO.Ara.minus4.data <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_single_KO_resilience_Ara-4_Rep999.csv") %>%
    mutate(population = "Ara-4") %>%
    mutate(starting_clone = "REL11336")

## Ara-5 REL11339 PPI network.
zitnik.single.KO.Ara.minus5.data <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_single_KO_resilience_Ara-5_Rep999.csv") %>%
    mutate(population = "Ara-5") %>%
    mutate(starting_clone = "REL11339")

cong.single.KO.Ara.minus5.data <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_single_KO_resilience_Ara-5_Rep999.csv") %>%
    mutate(population = "Ara-5") %>%
    mutate(starting_clone = "REL11339")

## Ara-6 REL11389 PPI network.
zitnik.single.KO.Ara.minus6.data <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_single_KO_resilience_Ara-6_Rep999.csv") %>%
    mutate(population = "Ara-6") %>%
    mutate(starting_clone = "REL11389")

cong.single.KO.Ara.minus6.data <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_single_KO_resilience_Ara-6_Rep999.csv") %>%
    mutate(population = "Ara-6") %>%
    mutate(starting_clone = "REL11389")

zitnik.single.KO.50K.data <- rbind(
    zitnik.single.KO.Ara.plus1.data,
    zitnik.single.KO.Ara.plus2.data,
    zitnik.single.KO.Ara.plus3.data,
    zitnik.single.KO.Ara.plus4.data,
    zitnik.single.KO.Ara.plus5.data,
    zitnik.single.KO.Ara.plus6.data,
    zitnik.single.KO.Ara.minus1.data,
    zitnik.single.KO.Ara.minus2.data,
    zitnik.single.KO.Ara.minus3.data,
    zitnik.single.KO.Ara.minus4.data,
    zitnik.single.KO.Ara.minus5.data,
    zitnik.single.KO.Ara.minus6.data
)

cong.single.KO.50K.data <- rbind(
    cong.single.KO.Ara.plus1.data,
    cong.single.KO.Ara.plus2.data,
    cong.single.KO.Ara.plus3.data,
    cong.single.KO.Ara.plus4.data,
    cong.single.KO.Ara.plus5.data,
    cong.single.KO.Ara.plus6.data,
    cong.single.KO.Ara.minus1.data,
    cong.single.KO.Ara.minus2.data,
    cong.single.KO.Ara.minus3.data,
    cong.single.KO.Ara.minus4.data,
    cong.single.KO.Ara.minus5.data,
    cong.single.KO.Ara.minus6.data
)


## get the initial resilience values for each 50K clone.
zitnik.base.resilience.df <- zitnik.single.KO.50K.data %>%
    filter(!str_detect(strain, "_knockout")) %>%
    ## This for changing the ordering of populations in plots.
    mutate(population=factor(population,levels = LTEE.pop.vec))

## get the initial resilience values for each 50K clone.
cong.base.resilience.df <- cong.single.KO.50K.data %>%
    filter(!str_detect(strain, "_knockout")) %>%
    ## This for changing the ordering of populations in plots.
    mutate(population=factor(population,levels = LTEE.pop.vec))

zitnik.single.KO.50K.resilience.df <- zitnik.single.KO.50K.data %>%
    filter(!(strain %in% zitnik.base.resilience.df$strain)) %>%
    mutate(Gene = str_split_fixed(strain, "_", n = 2)[,1]) %>%
    mutate(isEssentialKO = Gene %in% essential.genes$Gene) %>%
    mutate(Essentiality = ifelse(isEssentialKO, "Essential", "Non-essential")) %>%
    mutate(Dataset = "Zitnik") %>%
    ## This for changing the ordering of populations in plots.
    mutate(population=factor(population,levels = LTEE.pop.vec))


cong.single.KO.50K.resilience.df <- cong.single.KO.50K.data %>%
    filter(!(strain %in% cong.base.resilience.df$strain)) %>%
    mutate(Gene = str_split_fixed(strain, "_", n = 2)[,1]) %>%
    mutate(isEssentialKO = Gene %in% essential.genes$Gene) %>%
    mutate(Essentiality = ifelse(isEssentialKO, "Essential", "Non-essential")) %>%
    mutate(Dataset = "Cong") %>%
    ## This for changing the ordering of populations in plots.
    mutate(population=factor(population,levels = LTEE.pop.vec))



calc.delta.resilience <- function(pop.df) {
    ## This function calculates the differences between
    ## the mean of the resilience of essential KOs and
    ## the mean of the resilience of non-essential KOs.
    essential.mean <- mean(filter(pop.df, isEssentialKO == TRUE)$resilience)
    nonessential.mean <- mean(filter(pop.df, isEssentialKO == FALSE)$resilience)
    delta.mean <- essential.mean - nonessential.mean
    pop <- unique(pop.df$population)
    return.df <- data.frame(essential.mean.resilience = essential.mean,
                            nonessential.mean.resilience = nonessential.mean,
                            delta = delta.mean, population = pop)
    return(return.df)
}

## calculate the mean difference between essential and non-essential KOs on resilience.
zitnik.essential.KO.50K.delta.df <- zitnik.single.KO.50K.resilience.df %>%
    split(.$population) %>%
    map_dfr(.f = calc.delta.resilience)
## calculate the mean difference between essential and non-essential KOs on resilience.
cong.essential.KO.50K.delta.df <- cong.single.KO.50K.resilience.df %>%
    split(.$population) %>%
    map_dfr(.f = calc.delta.resilience)

## for both dataset, 10/12 populations have a negative difference, indicating
## that knockouts of essential genes, overall have a more negative effect
## on interactome resilience. However, the magnitude of the difference is small--
## < 0.001 in all cases.

## calculate Wilcoxon signed-rank test to do a nonparametric paired sample test.
wilcox.test(zitnik.essential.KO.50K.delta.df$essential.mean.resilience,
            zitnik.essential.KO.50K.delta.df$nonessential.mean.resilience,
            paired = TRUE, alternative = "less")

wilcox.test(cong.essential.KO.50K.delta.df$essential.mean.resilience,
            cong.essential.KO.50K.delta.df$nonessential.mean.resilience,
            paired = TRUE, alternative = "less")


S1Fig.title <- ggdraw() + 
    draw_label(
        "Single-gene disruptions of the 50,000 generation PPI networks",
        fontface = 'bold',
        x = 0,
        hjust = 0
    ) +
    theme(
        ## add margin on the left of the drawing canvas,
        ## so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
    )

S1Fig.panelA <- ggplot(zitnik.single.KO.50K.resilience.df,
                       aes(x = Essentiality, y = resilience)) +
    theme_classic() +
    ylab("PPI network resilience") +
    xlab("Type of disrupted gene") +
    facet_wrap(.~population,nrow=4) +
    geom_jitter(alpha = 0.1, height = 0, color = "moccasin") +
    geom_hline(data = zitnik.base.resilience.df,
               aes(yintercept = resilience),
               color = 'red', linetype = "dotted") +
    ggtitle("Zitnik PPI dataset")

S1Fig.panelB <- ggplot(cong.single.KO.50K.resilience.df,
                       aes(x = Essentiality, y = resilience)) +
    theme_classic() +
    ylab("PPI network resilience") +
    xlab("Type of disrupted gene") +
    facet_wrap(.~population,nrow=4) +
    geom_jitter(alpha = 0.1, height = 0, color = "lightskyblue") +
    geom_hline(data = cong.base.resilience.df,
               aes(yintercept = resilience),
               color = 'red', linetype = "dotted") +
    ggtitle("Cong PPI dataset")

S1Fig.panels <- plot_grid(S1Fig.panelA, S1Fig.panelB, ncol = 2)

S1Fig <- plot_grid(
    S1Fig.title, S1Fig.panels,
    ncol = 1,
    ## rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
)
    
ggsave("../results/resilience/figures/S1Fig.pdf", S1Fig, height = 6, width = 12)

#########################################################################################
## examine resilience in the MAE genomes.
## not enough evidence to make any definite conclusions in comparison to the LTEE.
## might need a 60,000 generation MAE experiment to do so!

cong.MAE.resilience.df <- read.csv("../results/resilience/resilience-analysis-runs/Cong_PPI_MAE_resilience_Rep999.csv") %>%
    filter(strain != "REL606") %>%
    mutate(log.resilience=log(resilience)) %>%
    ## annotate the MAE genomes.
    mutate(time = 13750) %>%
    mutate(Generation=time/10000) %>%
    mutate(Treatment = "actual data") %>%
    mutate(mutator_status = "non-mutator") %>%
    mutate(dataset = "Cong") %>%
    mutate(run_type = "MAE_genome_resilience") %>%
    mutate(population = "MAE")

cong.LTEE.10to20K.resilience.df <- big.cong.resilience.df %>%
    filter(run_type == "LTEE_genome_resilience") %>%
    filter(Generation >= 1) %>%
    filter(Generation <= 2)

cong.LTEE.to.MAE.comparison.df <- cong.MAE.resilience.df %>%
    full_join(cong.LTEE.10to20K.resilience.df)

## make the comparison plot.
cong.MAE.LTEE.comp.plot <- ggplot(cong.LTEE.to.MAE.comparison.df,
                                  aes(x = Generation,
                                      y = resilience,
                                      color = population,
                                      shape = mutator_status)) +
    theme_classic() +
    ylab("PPI network resilience") +
    xlab("Time (x 10,000 generations)") +
    geom_jitter(size = 0.5)


zitnik.MAE.resilience.df <- read.csv("../results/resilience/resilience-analysis-runs/Zitnik_PPI_MAE_resilience_Rep999.csv") %>%
    filter(strain != "REL606") %>%
    mutate(log.resilience=log(resilience)) %>%
    ## annotate the MAE genomes.
    mutate(time = 13750) %>%
    mutate(Generation=time/10000) %>%
    mutate(Treatment = "actual data") %>%
    mutate(mutator_status = "non-mutator") %>%
    mutate(dataset = "Zitnik") %>%
    mutate(run_type = "MAE_genome_resilience") %>%
    mutate(population = "MAE")


zitnik.LTEE.10to20K.resilience.df <- big.zitnik.resilience.df %>%
    filter(run_type == "LTEE_genome_resilience") %>%
    filter(Generation >= 1) %>%
    filter(Generation <= 2)

zitnik.LTEE.to.MAE.comparison.df <- zitnik.MAE.resilience.df %>%
    full_join(zitnik.LTEE.10to20K.resilience.df)


## make the comparison plot.
zitnik.MAE.LTEE.comp.plot <- ggplot(zitnik.LTEE.to.MAE.comparison.df,
                                  aes(x = Generation,
                                      y = resilience,
                                      color = population,
                                      shape = mutator_status)) +
    theme_classic() +
    ylab("PPI network resilience") +
    xlab("Time (x 10,000 generations)") +
    geom_jitter(size = 0.5)
