## thermostability-analysis.R by Rohan Maddamsetti.

## get functions for dealing with LTEE metagenomics data.
source("metagenomics-library.R")
library(ggridges)
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

########################################################################
## calculate densities of mutations per gene.
## This is used for filtering data.
calc.gene.mutation.densities <- function(gene.mutation.data,use.pseudocount=FALSE) {
    all.mutation.density <- calc.gene.mutation.density(
        gene.mutation.data,
        c("missense", "sv", "synonymous", "noncoding", "indel", "nonsense"),
        add.pseudocount=use.pseudocount) %>%
        rename(all.mut.count = mut.count) %>%
        rename(all.mut.density = density)
    
    ## look at dN density.
    dN.density <- calc.gene.mutation.density(
        gene.mutation.data,c("missense", "nonsense"),
        add.pseudocount=use.pseudocount) %>%
        rename(dN.mut.count = mut.count) %>%
        rename(dN.mut.density = density)
    
    ## look at dS density.
    dS.density <- calc.gene.mutation.density(
        gene.mutation.data,c("synonymous"),
        add.pseudocount=use.pseudocount) %>%
        rename(dS.mut.count = mut.count) %>%
        rename(dS.mut.density = density)
    
    ## all except dS density.
    all.except.dS.density <- calc.gene.mutation.density(
        gene.mutation.data,c("sv", "indel", "nonsense", "missense"),
        add.pseudocount=use.pseudocount) %>%
        rename(all.except.dS.mut.count = mut.count) %>%
        rename(all.except.dS.mut.density = density)

    ## knockout (KO) density.
    KO.density <- calc.gene.mutation.density(
        gene.mutation.data,c("sv", "indel", "nonsense"),
        add.pseudocount=use.pseudocount) %>%
        rename(KO.mut.count = mut.count) %>%
        rename(KO.mut.density = density)
    
    ## combine these into one dataframe.
    gene.mutation.densities <- REL606.genes %>%
        full_join(all.mutation.density) %>%
        full_join(dN.density) %>%
        full_join(dS.density) %>%
        full_join(all.except.dS.density) %>%
        full_join(KO.density) %>%
        as_tibble()

    if (use.pseudocount) {
        gene.mutation.densities <- gene.mutation.densities %>%
            ## CRITICAL STEP: replace NAs with pseudocounts of 0.1
            ## We need to keep track of genes that haven't been hit by any mutations
            ## in a given mutation class (sv, indels, dN, etc.)
            replace_na(list(all.mut.count = 0.1,
                            dN.mut.count = 0.1,
                            dS.mut.count = 0.1,
                            all.except.dS.mut.count = 0.1,
                            KO.mut.count = 0.1)) %>%
            mutate(all.mut.density = all.mut.count/gene_length,
                            dN.mut.density = dN.mut.count/gene_length,
                            dS.mut.density = dS.mut.count/gene_length,
                            all.except.dS.mut.density = all.except.dS.mut.count/gene_length,
                            KO.mut.density = KO.mut.count/gene_length)
    } else {
        gene.mutation.densities <- gene.mutation.densities %>%
            ## CRITICAL STEP: replace NAs with zeros.
            ## We need to keep track of genes that haven't been hit by any mutations
            ## in a given mutation class (sv, indels, dN, etc.)
            replace_na(list(all.mut.count = 0, all.mut.density = 0,
                            dN.mut.count = 0, dN.mut.density = 0,
                            dS.mut.count = 0, dS.mut.density = 0,
                            all.except.dS.mut.count = 0, all.except.dS.mut.density = 0,
                            KO.mut.count = 0, KO.mut.density = 0))
    }
    
    
    return(gene.mutation.densities)
}

gene.mutation.densities <- calc.gene.mutation.densities(gene.mutation.data)

## Calculate mutation densities per gene separately for
## nonmutators and hypermutators.

nonmut.mutation.densities <- gene.mutation.data %>%
    filter(Population %in% nonmutator.pops) %>%
    calc.gene.mutation.densities()

hypermut.mutation.densities <- gene.mutation.data %>%
    filter(Population %in% hypermutator.pops) %>%
    calc.gene.mutation.densities()

## For the supplementary information, examine the correlations when genes
## with zero mutations are omitted.
nozero.nonmut.mutation.densities <- nonmut.mutation.densities %>%
    filter(all.mut.density > 0)

nozero.hypermut.mutation.densities <- hypermut.mutation.densities %>%
    filter(all.mut.density > 0)

########################################################################
## Do mRNA and protein abundances predict evolutionary rates in the LTEE and in nature?
########################################################################

## analysis of reformatted Tim cDNA macroarray expression data
## from Cooper, Rozen, Lenski 2003.
cooper.data <- read.csv("../results/thermostability/reformatted-Cooper2003-expression-array.csv",as.is=TRUE) %>% filter(!is.na(locus_tag))

## reshape and merge the mRNA and protein abundance datasets using tidyr.
tidy.cooper.data <- cooper.data %>%
    gather(`REL607_1`,`REL607_2`,`REL607_3`,`REL607_4`,
           `AraPlus1_20K_1`,`AraPlus1_20K_2`,`AraPlus1_20K_3`,`AraPlus1_20K_4`,
           `REL606_1`,`REL606_2`,`REL606_3`,`REL606_4`,
           `AraMinus1_20K_1`,`AraMinus1_20K_2`,`AraMinus1_20K_3`,`AraMinus1_20K_4`,
           key="RawSample",value="mRNA") %>%
    mutate(Sample = str_replace(RawSample, "_[:digit:]$", "")) %>%
    mutate(Replicate = str_match(RawSample, "[:digit:]$")) %>%
    select(locus_tag, Sample, Replicate, mRNA)

anc.cooper.data <- tidy.cooper.data %>%
    filter(Sample %in% c("REL606", "REL607")) %>%
    group_by(locus_tag, Sample) %>%
    ## only analyze genes with reproducible expression values.
    summarize(mean_mRNA = mean(mRNA, na.rm=TRUE)) %>%
    ungroup()

evol.cooper.data <- tidy.cooper.data %>%
    filter(Sample %in% c("AraPlus1_20K", "AraMinus1_20K")) %>%
    group_by(locus_tag, Sample) %>%
    ## only analyze genes with reproducible expression values.
    summarize(mean_mRNA = mean(mRNA, na.rm=TRUE)) %>%
    ungroup()

anc.cooper.with.nonmut.mutation.density <- anc.cooper.data %>%
    left_join(nonmut.mutation.densities)

anc.cooper.with.hypermut.mutation.density <- anc.cooper.data %>%
    left_join(hypermut.mutation.densities)

anc.cooper.nonmut.cor <- cor.test(anc.cooper.with.nonmut.mutation.density$mean_mRNA,
         anc.cooper.with.nonmut.mutation.density$all.mut.density)

anc.cooper.hypermut.cor <- cor.test(anc.cooper.with.hypermut.mutation.density$mean_mRNA,
                                    anc.cooper.with.hypermut.mutation.density$all.mut.density)
anc.cooper.hypermut.cor$p.value


evol.cooper.with.nonmut.mutation.density <- evol.cooper.data %>%
    left_join(nonmut.mutation.densities)

evol.cooper.with.hypermut.mutation.density <- evol.cooper.data %>%
    left_join(hypermut.mutation.densities)

evol.cooper.nonmut.cor <- cor.test(evol.cooper.with.nonmut.mutation.density$mean_mRNA,
         evol.cooper.with.nonmut.mutation.density$all.mut.density)

evol.cooper.hypermut.cor <- cor.test(evol.cooper.with.hypermut.mutation.density$mean_mRNA,
                                     evol.cooper.with.hypermut.mutation.density$all.mut.density)
evol.cooper.hypermut.cor$p.value


########################################################################
## Favate et al. (2021) analysis of ancestral clones and 11 evolved 50K LTEE clones.

Favate.TableS1 <- read.csv("../data/Favate2021_table_s1_read_counts.csv",
                           as.is=TRUE, header=TRUE) %>%
    select(repl, seqtype, line, target_id, est_counts, tpm) %>%
    ## now, we need to match the target_id field to genes or IDs in the REL606 genome.
    ## I'm not sure what the ERC_00XXX IDs are. Ignore those and the tRNAs for now.
    mutate(locus_tag = target_id) %>%
    filter(str_detect(locus_tag,"^ECB")) %>%
    ## filter for genes that were used in the metagenomics data.
    filter(locus_tag %in% REL606.genes$locus_tag)
    
## Two datasets to examine separately: ribosome profiling and RNAseq.
## restrict analysis to evolved strains.

## RNAseq data
rnaseq.Favate.data <- Favate.TableS1 %>%
    filter(seqtype == "rna") %>%
    filter(!(line %in% c("REL606", "REL607"))) %>%
    group_by(line, locus_tag) %>%
    summarize(avg_est_counts = mean(est_counts,na.rm=TRUE),
              avg_tpm = mean(tpm,na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(log_est_counts = log(1+avg_est_counts))

## Riboseq data
riboseq.Favate.data <- Favate.TableS1 %>%
    filter(seqtype == "ribo") %>%
    filter(!(line %in% c("REL606", "REL607"))) %>%
    group_by(line, locus_tag) %>%
    summarize(avg_est_counts = mean(est_counts,na.rm=TRUE),
              avg_tpm = mean(tpm,na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(log_est_counts = log(1+avg_est_counts))

## examine the correlations in the RNAseq data.
rnaseq.with.nonmut.mutation.density <- rnaseq.Favate.data %>%
    left_join(nonmut.mutation.densities)

rnaseq.with.hypermut.mutation.density <- rnaseq.Favate.data %>%
    left_join(hypermut.mutation.densities)

rnaseq.nonmut.cor <- cor.test(rnaseq.with.nonmut.mutation.density$log_est_counts,
         rnaseq.with.nonmut.mutation.density$all.mut.density)

rnaseq.hypermut.cor <- cor.test(rnaseq.with.hypermut.mutation.density$log_est_counts,
                                rnaseq.with.hypermut.mutation.density$all.mut.density)
rnaseq.hypermut.cor$p.value

## examine the correlations in the riboseq data.
riboseq.with.nonmut.mutation.density <- riboseq.Favate.data %>%
    left_join(nonmut.mutation.densities)

riboseq.with.hypermut.mutation.density <- riboseq.Favate.data %>%
    left_join(hypermut.mutation.densities)

riboseq.nonmut.cor <- cor.test(riboseq.with.nonmut.mutation.density$avg_est_counts,
                                   riboseq.with.nonmut.mutation.density$all.mut.density)
                                   
riboseq.hypermut.cor <- cor.test(riboseq.with.hypermut.mutation.density$avg_est_counts,
                                 riboseq.with.hypermut.mutation.density$all.mut.density)
riboseq.hypermut.cor$p.value

## now make a figure for this analysis.
make.mut.density.favate.panel <- function(favate.data,
                                          lbl.xpos, rlbl.ypos, plbl.ypos,
                                          my.color="gray",
                                          muts.to.plot="all") {
    ## This is a helper function that plots a single panel.

    ## annotate r and p-values on the panel.
    if (muts.to.plot == "dN") {
        favate.result <- cor.test(favate.data$log_est_counts, favate.data$dN.mut.density)
    } else if (muts.to.plot == "dS") {
        favate.result <- cor.test(favate.data$log_est_counts, favate.data$dS.mut.density)
    } else { ## default is to plot the density of all mutations.
        favate.result <- cor.test(favate.data$log_est_counts, favate.data$all.mut.density)
    }
    
    pearson.r <- signif(favate.result$estimate,digits=3)
    pearson.p.value <- signif(favate.result$p.value,digits=3)

    lbl.rval <- paste("r", "=", pearson.r)
    lbl.p.value <- paste("p", "=", pearson.p.value)

    if (muts.to.plot == "dN") {
        favate.panel <- favate.data %>%
            ggplot(aes(x = log_est_counts, y = dN.mut.density))
    } else if (muts.to.plot == "dS") {
        favate.panel <- favate.data %>%
            ggplot(aes(x = log_est_counts, y = dS.mut.density))
    } else { ## default is to plot the density of all mutations.
        favate.panel <- favate.data %>%
            ggplot(aes(x = log_est_counts, y = all.mut.density))

    }

    favate.panel <- favate.panel +
        geom_point(color = my.color, alpha = 0.2) +
        geom_smooth(method = 'lm', formula = y~x) +
        theme_classic() +
        ylab("Mutation density") +
        xlab("log(1 + estimated mRNA counts)") +
        ## annotate correlation on the figure.
        annotate("text", x = lbl.xpos, y = rlbl.ypos, size=3,
                 label = lbl.rval, fontface = 'bold.italic') +
        ## annotate p-value on the figure.
        annotate("text", x = lbl.xpos, y = plbl.ypos, size=3,
                 label = lbl.p.value, fontface = 'bold.italic')

    return(favate.panel)
}

make.mut.density.favate.panel(rnaseq.with.hypermut.mutation.density,
                              lbl.xpos=2,
                              rlbl.ypos=0.06,
                              plbl.ypos=0.07)

make.mut.density.favate.panel(rnaseq.with.hypermut.mutation.density,
                              lbl.xpos=2,
                              rlbl.ypos=0.06,
                              plbl.ypos=0.07,
                              muts.to.plot="dN")

make.mut.density.favate.panel(rnaseq.with.hypermut.mutation.density,
                              lbl.xpos=2,
                              rlbl.ypos=0.06,
                              plbl.ypos=0.07,
                              muts.to.plot="dS")

make.mut.density.favate.panel(rnaseq.with.nonmut.mutation.density,
                              lbl.xpos=2,
                              rlbl.ypos=0.06,
                              plbl.ypos=0.07)

#######################################
## Caglar et al. analysis of REL606 over time.
## Import RNA and protein abundance data from Caglar et al. (2017).
## make sure to check abundances during BOTH exponential and stationary phase.
Caglar.samples <- read.csv("../results/thermostability/glucose-Caglar2017-REL606-data/glucose-Caglar2017-S1.csv", as.is=TRUE, header=TRUE) %>%
    ## turn growthTime_hr into a number from string.
    mutate(growthTime_hr = as.numeric(growthTime_hr)) %>%
    ## drop columns which are the same for all samples,
    ## and ones I don't care about.
    select(-sampleNum, -experiment, -harvestDate,
           -carbonSource, -RNA_Data_Freq, -Protein_Data_Freq,
           -Mg_mM, -Na_mM, -Mg_mM_Levels, -Na_mM_Levels,
           -rSquared, -uniqueCondition, -uniqueCondition02,
           -cellTotal, -cellsPerTube)

Caglar.mRNA <- read.csv("../results/thermostability/glucose-Caglar2017-REL606-data/glucose-Caglar2017-S2.csv", as.is=TRUE, header=TRUE)
Caglar.protein <- read.csv("../results/thermostability/glucose-Caglar2017-REL606-data/glucose-Caglar2017-S3.csv", as.is=TRUE, header=TRUE)
## reshape and merge the mRNA and protein abundance datasets using tidyr.
tidy.mRNA <- Caglar.mRNA %>%
    gather(`MURI_016`,`MURI_017`,`MURI_018`,`MURI_019`,`MURI_020`,`MURI_021`,
           `MURI_022`,`MURI_023`,`MURI_024`,`MURI_025`,`MURI_026`,`MURI_027`,
           `MURI_028`,`MURI_029`,`MURI_030`,`MURI_031`,`MURI_032`,`MURI_033`,
           `MURI_097`,`MURI_098`,`MURI_099`,`MURI_100`,`MURI_101`,`MURI_102`,
           `MURI_103`,`MURI_104`,`MURI_105`,key="dataSet",value="mRNA")
tidy.protein <- Caglar.protein %>%
    gather(`MURI_016`,`MURI_017`,`MURI_018`,`MURI_019`,`MURI_020`,`MURI_021`,
           `MURI_022`,`MURI_023`,`MURI_024`,`MURI_025`,`MURI_026`,`MURI_027`,
           `MURI_028`,`MURI_029`,`MURI_030`,`MURI_031`,`MURI_032`,`MURI_033`,
           `MURI_097`,`MURI_098`,`MURI_099`,`MURI_100`,`MURI_101`,`MURI_102`,
           `MURI_103`,`MURI_104`,`MURI_105`,key="dataSet",value="Protein") %>%
    select(-old_refseq)

full.Caglar.data <- Caglar.samples %>%
    inner_join(tidy.mRNA) %>%
    inner_join(tidy.protein)

Caglar.summary <- full.Caglar.data %>%
    group_by(locus_tag,growthPhase,growthTime_hr) %>%
    summarise(mRNA.mean=mean(mRNA), mRNA.sd=sd(mRNA),
              Protein.mean=mean(Protein), Protein.sd=sd(Protein))

## Now compare mRNA and protein abundance at each timepoint to
## mutation density in nonmutators and hypermutators.
## use an inner join to exclude any genes with no protein/RNA data.
nonmut.density.Caglar <- inner_join(nonmut.mutation.densities, Caglar.summary)
hypermut.density.Caglar <- inner_join(hypermut.mutation.densities, Caglar.summary)

## For supplementary information, re-run when genes with no mutations are excluded.
nozero.nonmut.density.Caglar <- inner_join(nozero.nonmut.mutation.densities,
                                           Caglar.summary)
nozero.hypermut.density.Caglar <- inner_join(nozero.hypermut.mutation.densities,
                                             Caglar.summary)

correlate.mut.density.with.timepoint <- function(density.Caglar, my.method="pearson") {

    print.correlations.given.timepoint <- function(my.data) {
        my.t <- unique(my.data$growthTime_hr)
        mRNA.result <- cor.test(my.data$mRNA.mean, my.data$all.mut.density,
                                method=my.method)
        Protein.result <- cor.test(my.data$Protein.mean, my.data$all.mut.density,
                                   method=my.method)
        
        print(paste("TIME:", my.t, 'hrs'))
        print("mRNA abundance correlation with mut density:")
        print(mRNA.result)
        print("Protein abundance correlation with mut density:")
        print(Protein.result)
    }

    density.Caglar %>%
        split(.$growthTime_hr) %>%
        map(.f = print.correlations.given.timepoint)        
}

## positive correlation in mutation density and mRNA,
## but no correlation with mutation density and protein abundance
## for nonmutators.
correlate.mut.density.with.timepoint(nonmut.density.Caglar)
## HIGHLY SIGNIFICANT anti-correlations with mutation density with
## both mRNA and protein abundances in REL606 in all timepoints!
correlate.mut.density.with.timepoint(hypermut.density.Caglar)


## For supplement: re-run, excluding zeros.
correlate.mut.density.with.timepoint(nozero.nonmut.density.Caglar)
correlate.mut.density.with.timepoint(nozero.hypermut.density.Caglar)


plot.mut.density.mRNA.anticorrelation <- function(my.data, my.method="pearson",
                                                  muts.to.plot="all") {
    ## This is a helper function to plot an mRNA panel for the full figure.
    my.t <- unique(my.data$growthTime_hr)
    time.title <- paste0(as.character(my.t), "h")

    if (muts.to.plot == "dN") {
        mRNA.result <- cor.test(my.data$mRNA.mean, my.data$dN.mut.density, method=my.method)
    } else if (muts.to.plot == "dS") {
        mRNA.result <- cor.test(my.data$mRNA.mean, my.data$dS.mut.density, method=my.method)
    } else { ## default is to plot the density of all mutations.
        mRNA.result <- cor.test(my.data$mRNA.mean, my.data$all.mut.density, method=my.method)
    }

    ## if the Pearson correlation is not significant, then the regression line is gray.
    ## factor of (1/9) is a Bonferroni-correction for multiple tests.
    my.color <- ifelse(mRNA.result$p.value < 0.05*(1/9), "blue", "light gray")

    ## annotate r and p-values for Pearson correlations.
    pearson.r <- signif(mRNA.result$estimate,digits=3)
    pearson.p.value <- signif(mRNA.result$p.value,digits=3)

    lbl.rval <- paste("r", "=", pearson.r)
    lbl.p.value <- paste("p", "=", pearson.p.value)
    

    if (muts.to.plot == "dN") {
          mRNA.plot <- my.data %>%
              ggplot(aes(x = mRNA.mean, y = dN.mut.density)) +
              ylim(0,0.05)
    } else if (muts.to.plot == "dS") {
         mRNA.plot <- my.data %>%
             ggplot(aes(x = mRNA.mean, y = dS.mut.density)) +
             ylim(0,0.05)
    } else { ## default is to plot the density of all mutations.
         mRNA.plot <- my.data %>%
        ggplot(aes(x = mRNA.mean, y = all.mut.density))

    }
    
    mRNA.plot <- mRNA.plot +
        geom_point(color = "violet", alpha = 0.5) +
        geom_smooth(method = 'lm', formula = y~x, color = my.color) +
        theme_classic() +
        coord_cartesian(xlim = c(0, 15)) +
        ggtitle(time.title) +
        ## annotate correlation on the figure.
        annotate("text", x = 7, y = 0.03, size = 3,
                 label = lbl.rval, fontface = 'bold.italic') +
        ## annotate p-value on the figure.
        annotate("text", x = 7, y = 0.04, size = 3,
                 label = lbl.p.value, fontface = 'bold.italic')

    ## use some conditional statements to label x-axis for the
    ## bottom center panel (168 hr) and to label y-axis for the
    ## left center panel (6 hr).
    my.xlabel <- ifelse(my.t == 168, "mRNA","")
    my.ylabel <- ifelse(my.t == 6, "Mutation density","")
    
    mRNA.plot <- mRNA.plot +
        xlab(my.xlabel) +
        ylab(my.ylabel)
    return(mRNA.plot)
}

plot.mut.density.protein.anticorrelation <- function(my.data, my.method="pearson",
                                                     muts.to.plot="all") {
    ## This is a helper function to plot a protein panel for the full figure.
    my.t <- unique(my.data$growthTime_hr)
    time.title <- paste0(as.character(my.t), "h")

    if (muts.to.plot == "dN") {
        Protein.result <- cor.test(my.data$Protein.mean,
                                   my.data$dN.mut.density,method=my.method)
        
    } else if (muts.to.plot == "dS") {
        Protein.result <- cor.test(my.data$Protein.mean,
                                   my.data$dS.mut.density,method=my.method)
        
    } else { ## default is to plot the density of all mutations.
        Protein.result <- cor.test(my.data$Protein.mean,
                                   my.data$all.mut.density,method=my.method)  
    }
        
    ## if the Pearson correlation is not significant, then the regression line is gray.
    ## factor of (1/9) is a Bonferroni-correction for multiple tests.
    my.color <- ifelse(Protein.result$p.value < 0.05*(1/9), "blue", "light gray")

    ## annotate r and p-values for Pearson correlations.
    pearson.r <- signif(Protein.result$estimate,digits=3)
    pearson.p.value <- signif(Protein.result$p.value,digits=3)

    lbl.rval <- paste("r", "=", pearson.r)
    lbl.p.value <- paste("p", "=", pearson.p.value)


    if (muts.to.plot == "dN") {
          Protein.plot <- my.data %>%
              ggplot(aes(x = Protein.mean, y = dN.mut.density)) +
              ylim(0,0.05)
    } else if (muts.to.plot == "dS") {
         Protein.plot <- my.data %>%
             ggplot(aes(x = Protein.mean, y = dS.mut.density)) +
             ylim(0,0.05)
    } else { ## default is to plot the density of all mutations.
         Protein.plot <- my.data %>%
        ggplot(aes(x = Protein.mean, y = all.mut.density))

    }
   
    Protein.plot <- Protein.plot +
        geom_point(color = "light green", alpha = 0.5) +
        geom_smooth(method = 'lm', formula = y~x, color = my.color) +
        theme_classic() +
        coord_cartesian(xlim = c(0, 10)) +
        scale_x_continuous(breaks=c(0,5,10)) +
        ggtitle(time.title) +
        ## no need for ylabel since the mRNA subfigure will be on the left.
        ylab("") +
        ## annotate correlation on the figure.
        annotate("text", x = 5, y = 0.03, size=3,
                 label = lbl.rval, fontface = 'bold.italic') +
        ## annotate p-value on the figure.
        annotate("text", x = 5, y = 0.04, size=3,
                 label = lbl.p.value, fontface = 'bold.italic')
  

    ## use a conditional statement to label x-axis for the
    ## bottom center panel (168 hr).
    my.xlabel <- ifelse(my.t == 168, "Protein","")
    Protein.plot <- Protein.plot +
        xlab(my.xlabel)
    
    return(Protein.plot)
}

make.mut.density.RNA.protein.expression.figure <- function(density.Caglar,
                                                           method="pearson",
                                                           muts.for.plot="all") {
    ## makes a figure that combines all 18 panels.
    ## generate a list of smaller panels, then pass
    ## to cowplot::plot_grid using do.call().

    ## pass plotting parameters through the use of partial function application.
    mRNA.plot.func <- partial(.f = plot.mut.density.mRNA.anticorrelation,
                              my.method=method, muts.to.plot=muts.for.plot)

    protein.plot.func <- partial(.f = plot.mut.density.protein.anticorrelation,
                              my.method=method, muts.to.plot=muts.for.plot)
    
    mRNA.panels <- density.Caglar %>%
        split(.$growthTime_hr) %>%
        map(.f = mRNA.plot.func)
    
    protein.panels <- density.Caglar %>%
        split(.$growthTime_hr) %>%
        map(.f = protein.plot.func)

    ## fill in some parameters for plot_grid using partial function application,
    ## and save the specialized versions.
    plot_subfigure <- partial(.f = plot_grid, nrow = 3)
    ## use do.call to unroll the panels as arguments for plot_grid.

    ## make two 3x3 subfigures for the mRNA and protein correlations.
    mRNA.plots <- do.call(plot_subfigure, mRNA.panels)
    protein.plots <- do.call(plot_subfigure, protein.panels)

    big.fig <- plot_grid(mRNA.plots, protein.plots)
    return(big.fig)
}

Fig1 <- make.mut.density.RNA.protein.expression.figure(hypermut.density.Caglar)
S1Fig <- make.mut.density.RNA.protein.expression.figure(nonmut.density.Caglar)

ggsave("../results/thermostability/figures/Fig1.pdf", Fig1, height = 5, width = 9)
ggsave("../results/thermostability/figures/S1Fig.pdf", S1Fig, height = 5, width = 9)


dN.Fig1 <- make.mut.density.RNA.protein.expression.figure(hypermut.density.Caglar, muts.for.plot="dN")
dS.Fig1 <- make.mut.density.RNA.protein.expression.figure(hypermut.density.Caglar, muts.for.plot="dS")
ggsave("../results/thermostability/figures/dN-Fig1.pdf", dN.Fig1, height = 5, width = 9)
ggsave("../results/thermostability/figures/dS-Fig1.pdf", dS.Fig1, height = 5, width = 9)


spearmanFig1 <- make.mut.density.RNA.protein.expression.figure(hypermut.density.Caglar, method="spearman")
spearmanS1Fig <- make.mut.density.RNA.protein.expression.figure(nonmut.density.Caglar,method="spearman")
ggsave("../results/thermostability/figures/spearmanFig1.pdf", spearmanFig1, height = 5, width = 9)
ggsave("../results/thermostability/figures/spearmanS1Fig.pdf", spearmanS1Fig, height = 5, width = 9)

nozeroFig1 <- make.mut.density.RNA.protein.expression.figure(nozero.hypermut.density.Caglar)
nozeroS1Fig <- make.mut.density.RNA.protein.expression.figure(nozero.nonmut.density.Caglar)

ggsave("../results/thermostability/figures/nozeroFig1.pdf", nozeroFig1, height = 5, width = 9)
ggsave("../results/thermostability/figures/nozeroS1Fig.pdf", nozeroS1Fig, height = 5, width = 9)


dSnozeroFig1 <- make.mut.density.RNA.protein.expression.figure(nozero.hypermut.density.Caglar, muts.for.plot="dS")
ggsave("../results/thermostability/figures/dSnozeroFig1.pdf", dSnozeroFig1, height = 5, width = 9)

########################################################

## Compare "divergence" in terms of number of observed mutations
## per population to the correlation with protein abundance.

## Compare to the theoretical expectation in Fig. 4 of
## Protein Biophysics Explains Why Highly Abundant Proteins Evolve Slowly
## by Serohijos et al. (2012).

## IMPORTANT BUG TO FIX: MAKE SURE ZEROS IN EACH POPULATION ARE NOT THROWN OUT!

calc.pop.gene.mutation.densities <- function(gene.mutation.data) {
    all.mutation.density <- pop.calc.gene.mutation.density(
        gene.mutation.data,
        c("missense", "sv", "synonymous", "noncoding", "indel", "nonsense")) %>%
        rename(all.mut.count = mut.count) %>%
        rename(all.mut.density = density)
    
    ## look at dN density.
    dN.density <- pop.calc.gene.mutation.density(
        gene.mutation.data,c("missense", "nonsense")) %>%
        rename(dN.mut.count = mut.count) %>%
        rename(dN.mut.density = density)
    
    ## look at dS density.
    dS.density <- pop.calc.gene.mutation.density(
        gene.mutation.data,c("synonymous")) %>%
        rename(dS.mut.count = mut.count) %>%
        rename(dS.mut.density = density)
    
    ## combine these into one dataframe.
    pop.gene.mutation.densities <- REL606.genes %>%
        full_join(all.mutation.density) %>%
        full_join(dN.density) %>%
        full_join(dS.density) %>%
        as_tibble()
    
    pop.gene.mutation.densities <- pop.gene.mutation.densities %>%
        ## CRITICAL STEP: replace NAs with zeros.
        ## We need to keep track of genes that haven't been hit by any mutations
        ## in a given mutation class (sv, indels, dN, etc.)
        replace_na(list(all.mut.count = 0, all.mut.density = 0,
                        dN.mut.count = 0, dN.mut.density = 0,
                        dS.mut.count = 0, dS.mut.density = 0))
    
    return(pop.gene.mutation.densities)
}

calc.correlation.helper <- function(pop.density.Caglar.subdf,my.method="pearson") {
    my.test <- cor.test(pop.density.Caglar.subdf$Protein.mean,
             pop.density.Caglar.subdf$all.mut.density,
             method=my.method)
    my.df <- data.frame(Population = unique(pop.density.Caglar.subdf$Population),
                        growthTime_hr = unique(pop.density.Caglar.subdf$growthTime_hr),
                        correlation = my.test$estimate)
    return(my.df)
}

## count the number of observed mutations per population.
total.muts.per.population <- gene.mutation.data %>%
    group_by(Population) %>%
    summarize(total.observed.muts = n()) %>%
    ungroup()

######################################################################

## calculate the density of mutations per gene per population.
pop.density.Caglar <- calc.pop.gene.mutation.densities(gene.mutation.data) %>%
    ## use an inner join to exclude any genes with no protein/RNA data.
    inner_join(Caglar.summary)

## CRITICAL BUG TO FIX!!!!!!!!!!!! working here.
buggy.df <- calc.pop.gene.mutation.densities(gene.mutation.data)
filter(buggy.df, all.mut.count == 0)


## scratch space for debugging.

ara.minus5.gene.mutation.data <- filter(gene.mutation.data,Population=="Ara-5")

ara.minus.5.gene.mutation.densities <- calc.pop.gene.mutation.densities(
    ara.minus5.gene.mutation.data)


ara.minus.5.density.Caglar <- ara.minus.5.gene.mutation.densities %>%
    inner_join(Caglar.summary)




quit()
###################################

## for recalculating the correlations when genes with zero mutation density are
## excluded.
nozero.pop.density.Caglar <- pop.density.Caglar %>% filter(all.mut.density > 0)


## now calculate the correlation between abundance and density of mutations
## per population and per timepoint.
pop.mut.density.protein.abundance.correlation.df <- pop.density.Caglar %>%
    split(list(.$growthTime_hr, .$Population)) %>%
    map_dfr(.f = calc.correlation.helper) %>%
    ## merge with the number of observed mutations per population
    left_join(total.muts.per.population) %>%
    mutate(is.Hypermutator = Population %in% hypermutator.pops)

pop.mut.density.protein.abundance.correlation.plot <- ggplot(
    data = pop.mut.density.protein.abundance.correlation.df,
    aes(x=total.observed.muts, y = correlation,
        color=Population, shape = is.Hypermutator)) +
    geom_point() + theme_classic() + xlab("Total observed mutations per population") +
    ylab("Correlation between protein abundance and mutation density")

## NOW-- recalculate the correlations when genes with zero mutation density are
## excluded.
nozero.pop.mut.density.protein.abundance.correlation.df <- pop.density.Caglar %>%
    split(list(.$growthTime_hr, .$Population)) %>%
    map_dfr(.f = calc.correlation.helper) %>%
    ## merge with the number of observed mutations per population
    left_join(total.muts.per.population) %>%
    mutate(is.Hypermutator = Population %in% hypermutator.pops)

nozero.pop.mut.density.prot.abundance.plot <- ggplot(
    data=nozero.pop.mut.density.protein.abundance.correlation.df,
    aes(x=total.observed.muts, y = correlation,
        color=Population, shape = is.Hypermutator)) +
    geom_point() + theme_classic() + xlab("Total observed mutations per population") +
    ylab("Correlation between protein abundance and mutation density")

## CRITICAL BUG: DOUBLE-CHECK THESE CORRELATIONS!
## NON-MUTATOR CORRELATIONS ARE NOT CONSISTENT WITH PREVIOUS RESULTS!!!


########################################################

## PPI network statistics analysis.
## to generate these files, run: python snap-ppi-analysis.py 
zitnik.network.df <- read.csv("../results/thermostability/Zitnik_network_statistics.csv",as.is=TRUE,header=TRUE) %>% inner_join(REL606.genes)
    
cong.network.df <- read.csv("../results/thermostability/Cong_network_statistics.csv",as.is=TRUE,header=TRUE) %>% inner_join(REL606.genes)

nonmut.PPI.zitnik <- zitnik.network.df %>%
    left_join(nonmut.mutation.densities)

hypermut.PPI.zitnik <- zitnik.network.df %>%
left_join(hypermut.mutation.densities)

nonmut.PPI.cong <- cong.network.df %>%
    left_join(nonmut.mutation.densities)

hypermut.PPI.cong <- cong.network.df %>%
    left_join(hypermut.mutation.densities)

## look at Zitnik network.

## all are significant.
cor.test(nonmut.PPI.zitnik$Degree, nonmut.PPI.zitnik$all.mut.density)
cor.test(nonmut.PPI.zitnik$DegreeCentrality, nonmut.PPI.zitnik$all.mut.density)

cor.test(nonmut.PPI.zitnik$Pagerank, nonmut.PPI.zitnik$all.mut.density)
cor.test(nonmut.PPI.zitnik$HubScore, nonmut.PPI.zitnik$all.mut.density)
cor.test(nonmut.PPI.zitnik$AuthorityScore, nonmut.PPI.zitnik$all.mut.density)
cor.test(nonmut.PPI.zitnik$ClosenessCentrality, nonmut.PPI.zitnik$all.mut.density)
cor.test(nonmut.PPI.zitnik$BetweenessCentrality, nonmut.PPI.zitnik$all.mut.density)
cor.test(nonmut.PPI.zitnik$EigenvectorCentrality, nonmut.PPI.zitnik$all.mut.density)
cor.test(nonmut.PPI.zitnik$IsArticulationPoint, nonmut.PPI.zitnik$all.mut.density)

## These ones are significant.
cor.test(hypermut.PPI.zitnik$Degree, hypermut.PPI.zitnik$all.mut.density)
cor.test(hypermut.PPI.zitnik$DegreeCentrality, hypermut.PPI.zitnik$all.mut.density)

cor.test(hypermut.PPI.zitnik$Pagerank, hypermut.PPI.zitnik$all.mut.density)
cor.test(hypermut.PPI.zitnik$HubScore, hypermut.PPI.zitnik$all.mut.density)
cor.test(hypermut.PPI.zitnik$AuthorityScore, hypermut.PPI.zitnik$all.mut.density)
cor.test(hypermut.PPI.zitnik$ClosenessCentrality, hypermut.PPI.zitnik$all.mut.density)
cor.test(hypermut.PPI.zitnik$BetweenessCentrality, hypermut.PPI.zitnik$all.mut.density)
cor.test(hypermut.PPI.zitnik$EigenvectorCentrality, hypermut.PPI.zitnik$all.mut.density)
## This last one is not significant.
cor.test(hypermut.PPI.zitnik$IsArticulationPoint, hypermut.PPI.zitnik$all.mut.density)

## now look at Cong network.
cor.test(nonmut.PPI.cong$Degree, nonmut.PPI.cong$all.mut.density) ## significant
cor.test(nonmut.PPI.cong$DegreeCentrality, nonmut.PPI.cong$all.mut.density) ## significant

cor.test(nonmut.PPI.cong$Pagerank, nonmut.PPI.cong$all.mut.density) ## significant
cor.test(nonmut.PPI.cong$HubScore, nonmut.PPI.cong$all.mut.density) ## NS
cor.test(nonmut.PPI.cong$AuthorityScore, nonmut.PPI.cong$all.mut.density) ## NS
cor.test(nonmut.PPI.cong$ClosenessCentrality, nonmut.PPI.cong$all.mut.density) ## NS
cor.test(nonmut.PPI.cong$BetweenessCentrality, nonmut.PPI.cong$all.mut.density) ## NS
cor.test(nonmut.PPI.cong$EigenvectorCentrality, nonmut.PPI.cong$all.mut.density) ## NS
cor.test(nonmut.PPI.cong$IsArticulationPoint, nonmut.PPI.cong$all.mut.density) ## NS


cor.test(hypermut.PPI.cong$Degree, hypermut.PPI.cong$all.mut.density) ## significant
cor.test(hypermut.PPI.cong$DegreeCentrality, hypermut.PPI.cong$all.mut.density) ## ditto

cor.test(hypermut.PPI.cong$Pagerank, hypermut.PPI.cong$all.mut.density) ## NS
cor.test(hypermut.PPI.cong$HubScore, hypermut.PPI.cong$all.mut.density) ## significant
cor.test(hypermut.PPI.cong$AuthorityScore, hypermut.PPI.cong$all.mut.density) ## sig
cor.test(hypermut.PPI.cong$ClosenessCentrality, hypermut.PPI.cong$all.mut.density) #sig
cor.test(hypermut.PPI.cong$BetweenessCentrality, hypermut.PPI.cong$all.mut.density) ## NS
cor.test(hypermut.PPI.cong$EigenvectorCentrality, hypermut.PPI.cong$all.mut.density) ## sig
cor.test(hypermut.PPI.cong$IsArticulationPoint, hypermut.PPI.cong$all.mut.density) ## NS

## Summarize these results by plotting mutation density against degree distribution
## for both networks, and both nonmutators and hyper-mutators. These results are
## consistent throughout, and even the non-significant correlations show the same,
## consistent trends seen for degree distribution.

make.mut.density.PPI.degree.panel <- function(PPI.data,
                                              lbl.xpos, rlbl.ypos, plbl.ypos,
                                              my.color="gray") {
    ## This is a helper function that plots a single panel.

    ## annotate r and p-values on the panel.
    PPI.result <- cor.test(PPI.data$Degree,
                           PPI.data$all.mut.density)

    pearson.r <- signif(PPI.result$estimate,digits=3)
    pearson.p.value <- signif(PPI.result$p.value,digits=3)

    lbl.rval <- paste("r", "=", pearson.r)
    lbl.p.value <- paste("p", "=", pearson.p.value)
    
    PPI.panel <- PPI.data %>%
        ggplot(aes(x = Degree, y = all.mut.density)) +
        geom_point(color = my.color, alpha = 0.2) +
        geom_smooth(method = 'lm', formula = y~x) +
        theme_classic() +
        ylab("Mutation density") +
        xlab("PPI degree") +
        ## annotate correlation on the figure.
        annotate("text", x = lbl.xpos, y = rlbl.ypos, size=3,
                 label = lbl.rval, fontface = 'bold.italic') +
        ## annotate p-value on the figure.
        annotate("text", x = lbl.xpos, y = plbl.ypos, size=3,
                 label = lbl.p.value, fontface = 'bold.italic')

    return(PPI.panel)
}

make.mut.density.PPI.degree.figure <- function(nonmut.PPI.zitnik, nonmut.PPI.cong,
                                               hypermut.PPI.zitnik, hypermut.PPI.cong) {

    ## label positions for each figure:
    zitnik.xpos <- 100
    cong.xpos <- 10

    nonmut.r.ypos <- 0.01
    nonmut.p.ypos <- 0.012
    hypermut.r.ypos <- 0.03
    hypermut.p.ypos <- 0.035
    
    panelA <- make.mut.density.PPI.degree.panel(nonmut.PPI.zitnik,
                                                zitnik.xpos,
                                                nonmut.r.ypos,
                                                nonmut.p.ypos,
                                                "lightsteelblue") +
        ggtitle("Nonmutators")
    panelB <- make.mut.density.PPI.degree.panel(nonmut.PPI.cong,
                                                cong.xpos,
                                                nonmut.r.ypos,
                                                nonmut.p.ypos,
                                                "moccasin") +
        ggtitle("Nonmutators")
    panelC <- make.mut.density.PPI.degree.panel(hypermut.PPI.zitnik,
                                                zitnik.xpos,
                                                hypermut.r.ypos,
                                                hypermut.p.ypos,
                                                "lightsteelblue") +
        ggtitle("Hypermutators")
    panelD <- make.mut.density.PPI.degree.panel(hypermut.PPI.cong,
                                                cong.xpos,
                                                hypermut.r.ypos,
                                                hypermut.p.ypos,
                                                "moccasin") +
        ggtitle("Hypermutators")

    
    fig <- plot_grid(panelA, panelB, panelC, panelD, nrow = 2,
                        labels = c('A', 'B', 'C', 'D'))
    return(fig)
}

PPI.figure <- make.mut.density.PPI.degree.figure(nonmut.PPI.zitnik,
                                                 nonmut.PPI.cong,
                                                 hypermut.PPI.zitnik,
                                                 hypermut.PPI.cong)
ggsave("../results/thermostability/figures/PPI-figure.pdf",
       PPI.figure, width=4, height=4)

################################################################################
## Analyze E. coli data from the ProteomeVis database.
################################################################################

## python ProteomeVis-to-REL606.py produces the table we need to join with REL606.genes.
REL606.to.ProteomeVis.df <- read.csv("../results/thermostability/REL606-to-ProteomeVis.csv", as.is = TRUE, header = TRUE)

## Import relevant E. coli ProteomeVis data.
## See Razban et al. (2018) in Bioinformatics:
## "ProteomeVis: a web app for exploration of protein properties
## from structure to sequence evolution across organisms’ proteomes"

## IMPORTANT: the metadata in proteomevis_inspect comes from github:
## https://github.com/rrazban/proteomevis/tree/master/proteomevis
## HOWEVER, the data was downloaded from:
## http://proteomevis.chem.harvard.edu.
## DO NOT USE THE proteomevis_chain.csv TABLE CURRENTLY AVAILBLE ON GITHUB!
## Those data are log10 transformed from the actual values that I want!

## available protein data. 1262 E. coli proteins here.
proteome.vis.inspect.df <- read.csv("../results/thermostability/Ecoli-proteomevis_inspect.csv", as.is = TRUE, header = TRUE)
proteome.vis.download.df <- read.csv("../data/ProteomeVis-data/NODES_ecoli_TM_0.00-0.00_SID_0.00-0.00_20201104_2235.csv", as.is = TRUE, header = TRUE) %>%
    select(-id) ## drop the id column as this is just an index.

proteome.vis.df <- full_join(proteome.vis.inspect.df, proteome.vis.download.df) %>%
    left_join(REL606.to.ProteomeVis.df) %>%
    filter(!(is.na(Gene))) %>% ## 1259 genes pass this filter.
    ## remove outliers with few residues in contact.
    ## 1255 proteins left after filtering.
    filter(contact_density > 4) 

nonmut.proteome.vis.comp.df <- proteome.vis.df %>%
    left_join(nonmut.mutation.densities)
    
hypermut.proteome.vis.comp.df <- proteome.vis.df %>%
    left_join(hypermut.mutation.densities)

## positive correction with PPI when looking at nonmutator data.
cor.test(nonmut.proteome.vis.comp.df$PPI_degree,
         nonmut.proteome.vis.comp.df$all.mut.density)

## negative correlation with PPI when
## looking at hypermutator data.
cor.test(hypermut.proteome.vis.comp.df$PPI_degree,
         hypermut.proteome.vis.comp.df$all.mut.density)

## negative correlation with contact density in nonmutators.
cor.test(nonmut.proteome.vis.comp.df$contact_density,
         nonmut.proteome.vis.comp.df$all.mut.density)

##  negative correlation with contact density in hypermutators.
cor.test(hypermut.proteome.vis.comp.df$contact_density,
         hypermut.proteome.vis.comp.df$all.mut.density)

## Let's make a figure to summarize these findings.
make.ProteomeVis.contact.density.figure <- function(nonmut.proteome.vis.comp.df,
                                                 hypermut.proteome.vis.comp.df) {

    add.formatting <- function(my.plot) {
        formatted.plot <- my.plot +
            geom_point(color = "gray", alpha = 0.2) +
            geom_smooth(method = 'lm', formula = y~x) +
            theme_classic() +
            ylab("Mutation density") +
            xlab("Contact density") 
        return(formatted.plot)
    }
    
    nonmut.plot <- ggplot(nonmut.proteome.vis.comp.df,
                        aes(x = contact_density, y = all.mut.density)) %>%
        add.formatting() + 
        ggtitle("Nonmutators")
    
    hypermut.plot <- ggplot(hypermut.proteome.vis.comp.df,
                          aes(x = contact_density, y = all.mut.density)) %>%
        add.formatting() + 
        ggtitle("Hypermutators")

    fig <- plot_grid(nonmut.plot, hypermut.plot)
    return(fig)
}

ProteomeVisFig <- make.ProteomeVis.contact.density.figure(nonmut.proteome.vis.comp.df,
                                                       hypermut.proteome.vis.comp.df)
ggsave("../results/thermostability/figures/ProteomeVisFig.pdf", ProteomeVisFig, width=4, height=2)

############################################################################
## E. COLI MELTOME ATLAS DATA ANALYSIS

categorize.by.meltPoint <- function(meltome) {
    ## Categories of thermal stability from the Meltome Atlas paper:
    ## "For the purpose of the subsequent analysis,
    ## we defined four categories of thermal stability.
    ## We distinguish proteins with low Tm (AUC<0.35),
    ## medium Tm (AUC 0.35–0.65), high Tm (AUC>0.65),
    ## and nonmelters (no Tm)."

    meltPoint.quantile.to.category <- function(quantile.num) {
        if (quantile.num == 1) {
            return("low")
        } else if (quantile.num == 2) {
            return("medium")
        } else if (quantile.num == 3) {
            return("high")
        } else {
            return("NA")
        }
    }
    
    ## I sort and then divide the data with measured meltPoints
    ## into thirds to get the low Tm, medium Tm, and high Tm categories.
    data.with.meltPoint <- meltome %>%
        filter(!is.na(meltPoint)) %>%
        arrange(meltPoint) %>%
        mutate(meltPoint.quantile = ntile(meltPoint,3)) %>%
        mutate(Tm.category = sapply(meltPoint.quantile,
                                    meltPoint.quantile.to.category)) %>%
        select(-meltPoint.quantile)

    ## merge back into the original data including nonmelters (meltPoint == NA),
    ## and categorize the nonmelters.
    meltome.with.Tm.category <- meltome %>%
        left_join(data.with.meltPoint) %>%
        mutate(Tm.category = ifelse(is.na(meltPoint), "nonmelter", Tm.category)) %>%
        mutate(Tm.category = fct_relevel(Tm.category, c("low","medium","high","nonmelter")))
    return(meltome.with.Tm.category)
}

## These data were written out by filter-meltome-atlas.R.
## Run that script in order to generate '../results/thermostability/Ecoli-meltome.csv'.
Ecoli.meltome <- read.csv("../results/thermostability/Ecoli-meltome.csv") %>%
    rename(Gene = gene_name) %>%
    select(run_name, Gene, meltPoint) %>%
    distinct() %>%
    left_join(REL606.genes) %>%
    filter(!is.na(locus_tag)) %>%
    ## bin genes into low Tm, medium Tm, high Tm, and nonmelters.
    categorize.by.meltPoint()

## to check correctness, plot the distribution of meltPoints in each Tm category.
## This should be a supplementary figure, if this strategies is used for other plots.
meltome.tertile.plot <- ggplot(Ecoli.meltome,
                               aes(x = meltPoint, fill = Tm.category)) +
    xlab("Melting point") + ylab("Count") +
    geom_histogram(bins=100) + theme_classic() + guides(fill = FALSE)
ggsave("../results/thermostability/figures/meltome-tertiles.pdf",
       meltome.tertile.plot, height=3, width = 3)

## examine the nonmelter proteins.
nonmelters <- Ecoli.meltome %>% filter(is.na(meltPoint)) %>%
    select(-run_name) %>% distinct() %>%
    ungroup()

## do a STRING analysis of the nonmelters.
## print out to file.
write.csv(nonmelters,file="../results/thermostability/nonmelters.csv")
## STRING shows that oxidative phosphorylation, and ribosome are
## enriched KEGG pathways for non-melters.

## examine nonmutators and hypermutators separately, as usual.
nonmut.density.meltome <- Ecoli.meltome %>%
    left_join(nonmut.mutation.densities)

hypermut.density.meltome <- Ecoli.meltome %>%
    left_join(hypermut.mutation.densities)

## positive correlation between melting point and mutation density.
## these data exclude nonmelters (NA values).
cor.test(hypermut.density.meltome$meltPoint,
         hypermut.density.meltome$all.mut.density)

## no such correlation in nonmutators.
cor.test(nonmut.density.meltome$meltPoint,
         nonmut.density.meltome$all.mut.density)

## plot for the positive correlation between Tm and mutation density in hypermutators.
hypermut.meltome.plot <- ggplot(hypermut.density.meltome,
                                aes(x = all.mut.density,
                                    y = Tm.category)) +
    stat_density_ridges(quantile_lines = TRUE, color = "white", fill = "dark gray") +
    theme_classic() +
    xlab("Mutation density") +
    ylab(bquote(T[m]~" category"))

ggsave("../results/thermostability/figures/Tm-hypermut-mut-density.pdf",
       hypermut.meltome.plot, height=2.5, width=6)

## on average, nonmelters have a higher mutation density in the hypermutators.
## there is no such pattern in the nonmutators.
h1 <- filter(hypermut.density.meltome, Tm.category == "nonmelter")$all.mut.density
h2 <- filter(hypermut.density.meltome, Tm.category != "nonmelter")$all.mut.density
mean(h1)
mean(h2)
## This is statistically significant.
wilcox.test(h1,h2)

## exclude non-melters for now.
correlate.meltPoint.with.timepoint <- function(meltPoint.Caglar) {
    for (t in sort(unique(meltPoint.Caglar$growthTime_hr))) {
        my.data <- meltPoint.Caglar %>%
            filter(growthTime_hr == t)
        
        mRNA.result <- cor.test(my.data$mRNA.mean, my.data$meltPoint)
        Protein.result <- cor.test(my.data$Protein.mean, my.data$meltPoint)
        print(paste("TIME:",t,'hrs'))
        print("mRNA abundance correlation with meltPoint:")
        print(mRNA.result)
        print("Protein abundance correlation with meltPoint:")
        print(Protein.result)
    }
}

## Now compare mRNA and protein abundance at each timepoint to
## Tm.
meltome.with.abundance <- Ecoli.meltome %>%
    left_join(Caglar.summary)
## IMPORTANT RESULT: melting temperature, excluding non-melters,
## is negatively correlated with protein and RNA abundance,
## especially during growth phase, but also during starvation.
correlate.meltPoint.with.timepoint(meltome.with.abundance)

plot.meltPoint.mRNA.anticorrelation <- function(my.data) {
    ## This is a helper function to plot an mRNA panel for the full figure.
    my.t <- unique(my.data$growthTime_hr)
    time.title = paste0(as.character(my.t), "h")

    ## if the correlation is not significant, then the fill is gray.
    mRNA.result <- cor.test(my.data$mRNA.mean, my.data$meltPoint)
        ## factor of (1/9) is a Bonferroni-correction for multiple tests.
    my.fill <- ifelse(mRNA.result$p.value < 0.05*(1/9), "lightskyblue", "light gray")
    
    mRNA.plot <- my.data %>%
        ggplot(aes(x = mRNA.mean, y = Tm.category)) +
        stat_density_ridges(quantile_lines = TRUE, fill = my.fill) +
        theme_classic() +
        coord_cartesian(xlim = c(0, 15)) +
        ggtitle(time.title)

    ## use some conditional statements to label x-axis for the
    ## bottom center panel (168 hr) and to label y-axis for the
    ## left center panel (6 hr).
    my.xlabel <- ifelse(my.t == 168, "mRNA","")
    ## have to use the full if else statement when using bquote
    ## because ifelse() is a vectorized construct than can't take
    ## language objects.
    if(my.t == 6) {
        my.ylabel <- bquote(T[m]~" category")
    } else {
        my.ylabel <- ""
    }
    
    mRNA.plot <- mRNA.plot +
        xlab(my.xlabel) +
        ylab(my.ylabel)
    
    ## use more conditional statements to only label y-axis categories
    ## for plots on the left-hand side of the big figure.
    do.not.label.yaxis <- TRUE #ifelse(!(my.t %in% c(3,6,48)),TRUE,FALSE)
    if (do.not.label.yaxis) {
        mRNA.plot <- mRNA.plot + theme(axis.text.y = element_blank())
    }
    
    return(mRNA.plot)
}

plot.meltPoint.protein.anticorrelation <- function(my.data) {
    ## This is a helper function to plot a protein panel for the full figure.
    my.t <- unique(my.data$growthTime_hr)
    time.title = paste0(as.character(my.t), "h")

    ## if the correlation is not significant, then the fill is gray.
    Protein.result <- cor.test(my.data$Protein.mean, my.data$meltPoint)
    ## factor of (1/9) is a Bonferroni-correction for multiple tests.
    my.fill <- ifelse(Protein.result$p.value < 0.05*(1/9), "moccasin", "light gray")
    
    Protein.plot <- my.data %>%
        ggplot(aes(x = Protein.mean, y = Tm.category)) +
        stat_density_ridges(quantile_lines = TRUE, fill = my.fill) +
        theme_classic() +
        coord_cartesian(xlim = c(0, 10)) +
        scale_x_continuous(breaks=c(0,5,10)) +
        ggtitle(time.title) +
        ylab("") ## no need for ylabel since the mRNA subfigure will be on the left.

    ## use a conditional statement to label x-axis for the
    ## bottom center panel (168 hr).
    my.xlabel <- ifelse(my.t == 168, "Protein","")
    Protein.plot <- Protein.plot + xlab(my.xlabel)
    ## y-axis categories are labeled on the left-hand side of the big figure.
    Protein.plot <- Protein.plot + theme(axis.text.y = element_blank())
    
    return(Protein.plot)
}

make.meltPoint.RNA.protein.expression.figure <- function(meltPoint.Caglar) {
    ## makes a figure that combines all panels.
    ## generate a list of smaller panels, then pass
    ## to cowplot::plot_grid using do.call().
    
    mRNA.panels <- meltPoint.Caglar %>%
        split(.$growthTime_hr) %>%
        map(.f = plot.meltPoint.mRNA.anticorrelation)
    
    protein.panels <- meltPoint.Caglar %>%
        split(.$growthTime_hr) %>%
        map(.f = plot.meltPoint.protein.anticorrelation)

    ## fill in some parameters for plot_grid using partial function application,
    ## and save the specialized versions.
    plot_subfigure <- partial(.f = plot_grid, nrow = 3)
    ## use do.call to unroll the panels as arguments for plot_grid.

    ## make two 3x3 subfigures for the mRNA and protein correlations.
    mRNA.plots <- do.call(plot_subfigure, mRNA.panels)
    protein.plots <- do.call(plot_subfigure, protein.panels)

    big.fig <- plot_grid(mRNA.plots, protein.plots)
    return(big.fig)
}

meltPointOverTime.Fig <- make.meltPoint.RNA.protein.expression.figure(meltome.with.abundance)
ggsave("../results/thermostability/figures/meltPointOverTime.pdf",
       meltPointOverTime.Fig,
       height = 4, width = 9)

## Put both figures for the meltome analysis together for Figure 3.

Fig3 <- plot_grid(plot_grid(hypermut.meltome.plot,NULL),
                  meltPointOverTime.Fig,
                  labels = c('A','B'), ncol = 1,
                  rel_heights = c(1, 2))
                  
ggsave("../results/thermostability/figures/Fig3.pdf",
       Fig3, height = 6, width = 7)
