## resilience-analysis.R by Rohan Maddamsetti.

## CRITICAL TODO: use the REL606 resilience value for the 14 genomes without
## any KO mutations, then re-calculate statistics.

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


zitnik.network.resilience.df <- read.csv("../results/resilience/Zitnik_PPI_LTEE_genome_resilience.csv") %>% full_join(LTEE.genomes.KO.metadata) %>%
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

cong.network.resilience.df <- read.csv("../results/resilience/Cong_PPI_LTEE_genome_resilience.csv") %>% full_join(LTEE.genomes.KO.metadata) %>%
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
