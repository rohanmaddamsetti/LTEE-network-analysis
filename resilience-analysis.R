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
