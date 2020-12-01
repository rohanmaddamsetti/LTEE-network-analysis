## metagenomics-library.R by Rohan Maddamsetti
## This file is a library of common functions used across these three papers.

library(tidyverse)
library(cowplot)

##########################################################################
## FUNCTIONS FOR DATA ANALYSIS
##########################################################################

## function for plotting better y-axis labels.
## see solution here for nice scientific notation on axes.
## https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
fancy_scientific <- function(x) {
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

## look at accumulation of stars over time
## in other words, look at the rates at which the mutations occur over time.
## To normalize, we need to supply the number of sites at risk
## (such as sum of gene length).
## If plot.to.end is TRUE, then add one final row.
calc.cumulative.muts <- function(d, d.metadata, manual.pop.levels.vec=NA,
                                 plot.to.end=TRUE) {

    cumsum.per.pop.helper.func <- function(pop) {
        finalgen <- 6.3 ## this is outside of the data collection
        ## for nice plotting (final generation in mutation.data is 6.275).

        ## This constant is to make sure that all pops are in the levels
        ## of the Population factor after mergers, etc.
        pop.levels <- manual.pop.levels.vec
        if (any(is.na(manual.pop.levels.vec))) {
            pop.levels <- c("Ara-5", "Ara-6", "Ara+1",
                            "Ara+2", "Ara+4", "Ara+5",
                            "Ara-1", "Ara-2", "Ara-3",
                            "Ara-4", "Ara+3", "Ara+6")
        }
        
        df <- d %>% filter(Population==pop)
        if (nrow(df) == 0) { ## if no mutations in this pop.
            almost.done.df <- tibble(Population = factor(pop, levels = pop.levels),
                                     Generation=finalgen,
                                     count=0,
                                     cs=0)
        } else {
            summary.df <- df %>%
                arrange(t0) %>%
                group_by(Population, Generation) %>%
                summarize(count=n(), .groups = "drop_last") %>%
                mutate(cs=cumsum(count)) %>%
                ungroup()
            ## if the final generation is not in ret.df,
            ## then add one final row (for nicer plots).
            final.row.df <- tibble(Population=factor(pop, levels = pop.levels),
                                   Generation=finalgen,
                                   count=max(summary.df$count),
                                   cs=max(summary.df$cs))
            if (plot.to.end) {
                almost.done.df <- bind_rows(summary.df, final.row.df)
            } else {
                almost.done.df <- summary.df
            }
        }
        ## add an row for Generation == 0 (for nicer plots).
        init.row.df <- tibble(
            Population = factor(pop, levels = pop.levels),
            Generation = 0,
            count = 0,
            cs = 0)
        
        ret.df <- bind_rows(init.row.df,almost.done.df)
        return(ret.df)
    }

    ## normalize by the total length of genes
    ## in the given module (in d.metadata).
    my.genes <- d.metadata %>% dplyr::select(Gene, gene_length) %>% distinct()
    normalization.constant <- sum(my.genes$gene_length)

    
    c.dat <- map_dfr(.x=levels(d$Population),
                     .f=cumsum.per.pop.helper.func) %>%
        mutate(normalized.cs=cs/normalization.constant) %>%
        ## remove any NA values.
        na.omit()
    
    return(c.dat)
}

###############################

## calculate the tail probabilities of the true cumulative mutation trajectory
## of a given vector of genes (a 'module'), based on resampling
## random sets of genes. Returns the upper tail of null distribution,
## or P(random trajectory >= the actual trajectory).
## Output: a dataframe with three columns: Population, count, p.val
calc.traj.pvals <- function(data, REL606.genes, gene.vec, N=10000) {

    ## check type for gene.vec (e.g., if factor, change to vanilla vector of strings)
    gene.vec <- as.character(gene.vec)
    
    ## each sample has the same cardinality as the gene.vec.
    subset.size <- length(gene.vec)

    ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations.
    ## IMPORTANT: this function depends on variables defined in
    ## calculate.trajectory.tail.probs.
    generate.cumulative.mut.subset <- function(idx) {
        rando.genes <- base::sample(unique(REL606.genes$Gene),subset.size)
        mut.subset <- filter(data,Gene %in% rando.genes)
        mut.subset.metadata <- filter(REL606.genes, Gene %in% rando.genes)
        c.mut.subset <- calc.cumulative.muts(mut.subset,
                                             mut.subset.metadata) %>%
            mutate(bootstrap_replicate=idx)
        return(c.mut.subset)
    }
    
    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    bootstrapped.trajectories <- map_dfr(.x=seq_len(N),
                                         .f=generate.cumulative.mut.subset)
    
    gene.vec.data <- data %>% filter(Gene %in% gene.vec)
    gene.vec.metadata <- REL606.genes %>% filter(Gene %in% gene.vec)
    data.trajectory <- calc.cumulative.muts(gene.vec.data, gene.vec.metadata)
    data.trajectory.summary <- data.trajectory %>%
        group_by(Population,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup()
    
    trajectory.summary <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate, Population,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup()
    
    trajectory.filter.helper <- function(pop.trajectories) {
        pop <- unique(pop.trajectories$Population)
        data.traj <- filter(data.trajectory.summary,Population == pop)
        final.data.norm.cs <- unique(data.traj$final.norm.cs)
        tail.trajectories <- filter(pop.trajectories, final.norm.cs >= final.data.norm.cs)
        return(tail.trajectories)
    }
    
    ## split by Population, then filter for bootstraps > data trajectory.
    uppertail.probs <- trajectory.summary %>%
        split(.$Population) %>%
        map_dfr(.f=trajectory.filter.helper) %>%
        group_by(Population,.drop=FALSE) %>%
        summarize(count=n()) %>%
        mutate(p.val=count/N)
    
    return(uppertail.probs)
}

## This plot visualizes a two-tailed test (alphaval = 0.05)
## against a bootstrapped null distribution.
## Throughout, plots use the minimum subsample size to subsample the null distribution,
## to increase the variance in order to make a conservative comparison.
plot.base.layer <- function(data, REL606.genes, subset.size=50, N=1000, alphaval = 0.05, manual.pop.levels.vec=NA, my.color="gray", plot.rows=4) {
    
    ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations.
    ## IMPORTANT: this function depends on variables defined in plot.base.layer.
    generate.cumulative.mut.subset <- function(idx) {
        rando.genes <- base::sample(unique(REL606.genes$Gene),subset.size)
        mut.subset <- filter(data,Gene %in% rando.genes)
        mut.subset.metadata <- filter(REL606.genes, Gene %in% rando.genes)
        c.mut.subset <- calc.cumulative.muts(mut.subset,
                                             mut.subset.metadata,
                                             manual.pop.levels.vec) %>%
            mutate(bootstrap_replicate=idx)
        return(c.mut.subset)
    }

    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    ## I want to plot the distribution of cumulative mutations over time for
    ## say, 1000 or 10000 random subsets of genes.

    bootstrapped.trajectories <- map_dfr(.x=seq_len(N),.f=generate.cumulative.mut.subset)

    ## filter out the top alphaval/2 and bottom alphaval/2 trajectories from each population,
    ## for a two-sided test. default is alphaval == 0.05.
    
    trajectory.summary <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate, Population,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup() 
    
    top.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        slice_max(prop=alphaval/2, order_by=final.norm.cs) %>%
        dplyr::select(-final.norm.cs) %>%
        mutate(in.top=TRUE)
    
    bottom.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        slice_min(prop=alphaval/2, order_by=final.norm.cs) %>%
        dplyr::select(-final.norm.cs) %>%
        mutate(in.bottom=TRUE)
    
    filtered.trajectories <- bootstrapped.trajectories %>%
        left_join(top.trajectories) %>%
        left_join(bottom.trajectories) %>%
        filter(is.na(in.top)) %>%
        filter(is.na(in.bottom)) %>%
        dplyr::select(-in.top,-in.bottom)

    ## hack to remove unwanted Population levels from the plot.
 ##   if (!any(is.na(manual.pop.levels.vec))) {
 ##       filtered.trajectories <- filtered.trajectories %>%
 ##           filter(Population %in% manual.pop.levels.vec) %>%
 ##           mutate(Population = factor(Population, levels=manual.pop.levels.vec))
 ##   }
    
    p <- ggplot(filtered.trajectories,aes(x=Generation,y=normalized.cs)) +
        ylab('Cumulative mutations (normalized)') +
        theme_classic() +
        geom_point(size=0.2, color=my.color) +
        facet_wrap(.~Population,scales='free',nrow=plot.rows) +
        xlab('Generations (x 10,000)') +
        xlim(0,6.3) +
        theme(axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              axis.text.x  = element_text(size=14),
              axis.text.y  = element_text(size=14)) +
        scale_y_continuous(labels=fancy_scientific,
                           breaks = scales::extended_breaks(n = 6),
                           limits = c(0, NA))
    return(p)
}

## take a ggplot object output by plot.cumulative.muts, and add an extra layer.
add.cumulative.mut.layer <- function(p, layer.df, my.color) {
    p <- p +
        geom_point(data=layer.df,
                   aes(x=Generation,y=normalized.cs),
                   color=my.color, size=0.2) +
        geom_step(data=layer.df, aes(x=Generation,y=normalized.cs),
                  size=0.2, color=my.color)
    return(p)
}

################################################################################
## These functions are used in my LTEE thermostability analysis.

## Examine the distribution of various classes of mutations across genes in the
## genomics or metagenomics data. Which genes are enriched? Which genes are depleted?
## Then, can look at the annotation of these genes in STRING.
calc.gene.mutation.density <- function(gene.mutation.data, mut_type_vec) {
    density.df <- gene.mutation.data %>%
        filter(Annotation %in% mut_type_vec) %>%
        filter(Gene!= "intergenic") %>%
        mutate(Gene=as.factor(Gene)) %>%
        group_by(Gene,gene_length) %>%
        summarize(mut.count=n()) %>%
        ungroup() %>%
        mutate(density=mut.count/gene_length) %>%
        arrange(desc(density))

    ## CRITICAL STEP: replace NAs with zeros.
    ## We need to keep track of genes that haven't been hit by any mutations
    density.df[is.na(density.df)] <- 0
    density.df <- tbl_df(density.df)

    
    return(density.df)
}

## Examine the gene mutation density by population.
pop.calc.gene.mutation.density <- function(gene.mutation.data, mut_type_vec) {
    density.df <- gene.mutation.data %>%
        filter(Annotation %in% mut_type_vec) %>%
        filter(Gene!= "intergenic") %>%
        mutate(Gene=as.factor(Gene)) %>%
        group_by(Population, Gene, gene_length) %>%
        summarize(mut.count=n()) %>%
        ungroup() %>%
        mutate(density=mut.count/gene_length) %>%
        arrange(desc(density))

    ## CRITICAL STEP: replace NAs with zeros.
    ## We need to keep track of genes that haven't been hit by any mutations
    density.df[is.na(density.df)] <- 0
    density.df <- tbl_df(density.df)
    
    return(density.df)
}
################################################################################
