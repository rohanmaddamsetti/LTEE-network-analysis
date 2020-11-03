## thermostability-analysis.R by Rohan Maddamsetti.

## IMPORTANT TODO: See if I can recapitulate the correlations
## in Figure 1 of Drummond and Wilke (2008) using LTEE hypermutator data.

## get functions for dealing with LTEE metagenomics data.
source("metagenomics-library.R")

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
    mutate(Population=factor(Population,levels=c(nonmutator.pops,hypermutator.pops))) %>%
    inner_join(REL606.genes) %>%
    filter(Gene!='intergenic')

########################################################################
## calculate densities of mutations per gene.
## This is used for filtering data.
calc.gene.mutation.densities <- function(gene.mutation.data) {
    all.mutation.density <- calc.gene.mutation.density(
        gene.mutation.data,
        c("missense", "sv", "synonymous", "noncoding", "indel", "nonsense")) %>%
        rename(all.mut.count = mut.count) %>%
        rename(all.mut.density = density)
    
    ## look at dN density.
    dN.density <- calc.gene.mutation.density(
        gene.mutation.data,c("nonsynonymous")) %>%
        rename(dN.mut.count = mut.count) %>%
        rename(dN.mut.density = density)
    
    ## look at dS density.
    dS.density <- calc.gene.mutation.density(
        gene.mutation.data,c("synonymous")) %>%
        rename(dS.mut.count = mut.count) %>%
        rename(dS.mut.density = density)
    
    ## all except dS density.
    all.except.dS.density <- calc.gene.mutation.density(
        gene.mutation.data,c("sv", "indel", "nonsense", "missense")) %>%
        rename(all.except.dS.mut.count = mut.count) %>%
        rename(all.except.dS.mut.density = density)

    ## knockout (KO) density.
    KO.density <- calc.gene.mutation.density(
        gene.mutation.data,c("sv", "indel", "nonsense")) %>%
        rename(KO.mut.count = mut.count) %>%
        rename(KO.mut.density = density)
    
    ## combine these into one dataframe.
    gene.mutation.densities <- REL606.genes %>%
        full_join(all.mutation.density) %>%
        full_join(dN.density) %>%
        full_join(dS.density) %>%
        full_join(all.except.dS.density) %>%
        full_join(KO.density) %>%
        tbl_df() %>%
        ## CRITICAL STEP: replace NAs with zeros.
        ## We need to keep track of genes that haven't been hit by any mutations
        ## in a given mutation class (sv, indels, dN, etc.)
        replace_na(list(all.mut.count = 0, all.mut.density = 0,
                        dN.mut.count = 0, dN.mut.density = 0,
                        dS.mut.count = 0, dS.mut.density = 0,
                        all.except.dS.mut.count = 0, all.except.dS.mut.density = 0,
                        KO.mut.count = 0, KO.density = 0))
    
    return(gene.mutation.densities)
}

gene.mutation.densities <- calc.gene.mutation.densities(gene.mutation.data)

## calculate mutation densities for genes, but separately for
## non-mutators and hyper-mutators.

nonmut.mutation.densities <- gene.mutation.data %>%
    filter(Population %in% nonmutator.pops) %>%
    calc.gene.mutation.densities()

hypermut.mutation.densities <- gene.mutation.data %>%
    filter(Population %in% hypermutator.pops) %>%
    calc.gene.mutation.densities()

########################################################################
## Get essential and near-essential genes reported in
## Supplementary Table 1 of Couce et al. 2017.
## I manually fixed the names of a couple genes in this dataset.
## The original names are in the "Name" column, and updated names
## are in the "Gene" column.
essential.genes <- read.csv("../data/Couce2017-LTEE-essential.csv") %>%
    inner_join(REL606.genes) %>% filter(!(is.na(locus_tag)))

###########################################################################
## Use LTEE Genomics data downloaded from barricklab.org/shiny/LTEE-Ecoli
## for extra filtering for genes with no mutations at all.
LTEE.genomics.muts <- read.csv("../data/LTEE-Ecoli-data 2020-02-24 14_35_41.csv") %>%
    ## ignore MOB, DEL, etc. that are annotated as intergenic.
    ## this is to prevent false negatives, i.e. removing genes in the gene_list
    ## that weren't directly affected by the intergenic mutation.
    filter(!(str_detect(html_mutation_annotation,"intergenic")))
    
## remove genes from 'no.mutation.genes' that have mutations in LTEE.genomics.muts.
## make a big string of all mutated genes in the LTEE-genomics data, and search this for
## matches.
LTEE.genomics.mutated.genestr <- str_c(unique(LTEE.genomics.muts$gene_list),collapse = ",")

no.mutation.genes <- REL606.genes %>%
    ## no hits allowed in the metagenomics.
    filter(!(Gene %in% gene.mutation.data$Gene)) %>%
    ## and no hits (other than dS, amps, and intergenic) allowed in the genomics.
    filter(!(str_detect(LTEE.genomics.mutated.genestr,Gene))) %>%
    arrange(desc(gene_length))

## genes that are only affected by synonymous mutations.
only.dS.allowed.genes <- gene.mutation.densities %>%
    filter(all.except.dS.mut.count == 0) %>%
    ## no hits (other than dS, amps, and intergenic) allowed in the LTEE genomics
    ## to remove false positives.
    filter(!(str_detect(LTEE.genomics.mutated.genestr,Gene))) %>%
    arrange(desc(gene_length))

########################################################################
## Do mRNA and protein abundances predict evolutionary rates in the LTEE and in nature?
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
## mutation density in non-mutators and hypermutators.
## use an inner join to exclude any genes with no protein/RNA data.
nonmut.density.Caglar <- inner_join(nonmut.mutation.densities, Caglar.summary)
hypermut.density.Caglar <- inner_join(hypermut.mutation.densities, Caglar.summary)

correlate.mut.density.with.timepoint <- function(density.Caglar) {

    print.correlations.given.timepoint <- function(my.data) {
        my.t <- unique(my.data$growthTime_hr)
        mRNA.result <- cor.test(my.data$mRNA.mean, my.data$all.mut.density)
        Protein.result <- cor.test(my.data$Protein.mean, my.data$all.mut.density)
        
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
## but no correlation with mutation density and  protein abundance
## for non-mutators.
correlate.mut.density.with.timepoint(nonmut.density.Caglar)
## HIGHLY SIGNIFICANT anti-correlations with mutation density with
## both mRNA and protein abundances in REL606 in all timepoints!
correlate.mut.density.with.timepoint(hypermut.density.Caglar)

plot.mut.density.mRNA.anticorrelation <- function(my.data) {
    ## This is a helper function to plot a single small panel for the full figure.
    my.t <- unique(my.data$growthTime_hr)
    time.title = paste0(as.character(my.t), "h")

    ## if the correlation is not significant, then the regression line is gray.
    mRNA.result <- cor.test(my.data$mRNA.mean, my.data$all.mut.density)
    my.color <- ifelse(mRNA.result$p.value < 0.05, "blue", "light gray")
    
    mRNA.plot <- my.data %>%
        ggplot(aes(x = mRNA.mean, y = all.mut.density)) +
        geom_point(color = "violet", alpha = 0.5) +
        geom_smooth(method = 'lm', formula = y~x, color = my.color) +
        theme_classic() +
        coord_cartesian(xlim = c(0, 15)) +
        ggtitle(time.title)

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

plot.mut.density.protein.anticorrelation <- function(my.data) {
    ## This is a helper function to plot a single panel for larger figures.
    my.t <- unique(my.data$growthTime_hr)
    time.title = paste0(as.character(my.t), "h")

    ## if the correlation is not significant, then the regression line is gray.
    Protein.result <- cor.test(my.data$Protein.mean, my.data$all.mut.density)
    my.color <- ifelse(Protein.result$p.value < 0.05, "blue", "light gray")
    
    Protein.plot <- my.data %>%
        ggplot(aes(x = Protein.mean, y = all.mut.density)) +
        geom_point(color = "light green", alpha = 0.5) +
        geom_smooth(method = 'lm', formula = y~x, color = my.color) +
        theme_classic() +
        coord_cartesian(xlim = c(0, 15)) +
        ggtitle(time.title) +
        ylab("") ## no need for ylabel since the mRNA subfigure will be on the left.

    ## use a conditional statement to label x-axis for the
    ## bottom center panel (168 hr).
    my.xlabel <- ifelse(my.t == 168, "Protein","")
    Protein.plot <- Protein.plot +
        xlab(my.xlabel)
    
    return(Protein.plot)
}

make.mut.density.RNA.protein.expression.figure <- function(density.Caglar) {
    ## makes a figure that combines all 18 panels.
    ## generate a list of smaller panels, then pass
    ## to cowplot::plot_grid using do.call().
    
    mRNA.panels <- density.Caglar %>%
        split(.$growthTime_hr) %>%
        map(.f = plot.mut.density.mRNA.anticorrelation)
    
    protein.panels <- density.Caglar %>%
        split(.$growthTime_hr) %>%
        map(.f = plot.mut.density.protein.anticorrelation)

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

ggsave("../results/thermostability/figures/Fig1.pdf", Fig1, height = 4, width = 9)
ggsave("../results/thermostability/figures/S1Fig.pdf", S1Fig, height = 4, width = 9)

########################################################

## PPI network statistics analysis.
## to generate these files, run: python snap-ppi-analysis.py 
zitnik.network.df <- read.csv("../results/thermostability/Zitnik_network_statistics.csv",as.is=TRUE,header=TRUE) %>% inner_join(REL606.genes)
    
cong.network.df <- read.csv("../results/thermostability/Cong_network_statistics.csv",as.is=TRUE,header=TRUE) %>% inner_join(REL606.genes)

nonmut.PPI.zitnik <- zitnik.network.df %>% left_join(nonmut.density) %>%
    ## turn NAs to 0s.
    mutate(mut.count=ifelse(is.na(mut.count),0,mut.count)) %>%
    mutate(density=ifelse(is.na(density),0,density))

hypermut.PPI.zitnik <- zitnik.network.df %>% left_join(hypermut.density) %>%
    ## turn NAs to 0s.
    mutate(mut.count=ifelse(is.na(mut.count),0,mut.count)) %>%
    mutate(density=ifelse(is.na(density),0,density))

nonmut.PPI.cong <- cong.network.df %>% left_join(nonmut.density) %>%
    ## turn NAs to 0s.
    mutate(mut.count=ifelse(is.na(mut.count),0,mut.count)) %>%
    mutate(density=ifelse(is.na(density),0,density))

hypermut.PPI.cong <- cong.network.df %>% left_join(hypermut.density) %>%
    ## turn NAs to 0s.
    mutate(mut.count=ifelse(is.na(mut.count),0,mut.count)) %>%
    mutate(density=ifelse(is.na(density),0,density))

## look at Zitnik network.

## all are significant.
cor.test(nonmut.PPI.zitnik$Pagerank, nonmut.PPI.zitnik$density)
cor.test(nonmut.PPI.zitnik$HubScore, nonmut.PPI.zitnik$density)
cor.test(nonmut.PPI.zitnik$AuthorityScore, nonmut.PPI.zitnik$density)
cor.test(nonmut.PPI.zitnik$ClosenessCentrality, nonmut.PPI.zitnik$density)
cor.test(nonmut.PPI.zitnik$BetweenessCentrality, nonmut.PPI.zitnik$density)
cor.test(nonmut.PPI.zitnik$EigenvectorCentrality, nonmut.PPI.zitnik$density)
cor.test(nonmut.PPI.zitnik$Degree, nonmut.PPI.zitnik$density)
cor.test(nonmut.PPI.zitnik$DegreeCentrality, nonmut.PPI.zitnik$density)
cor.test(nonmut.PPI.zitnik$IsArticulationPoint, nonmut.PPI.zitnik$density)

## These ones are significant.
cor.test(hypermut.PPI.zitnik$Pagerank, hypermut.PPI.zitnik$density)
cor.test(hypermut.PPI.zitnik$HubScore, hypermut.PPI.zitnik$density)
cor.test(hypermut.PPI.zitnik$AuthorityScore, hypermut.PPI.zitnik$density)
cor.test(hypermut.PPI.zitnik$ClosenessCentrality, hypermut.PPI.zitnik$density)
cor.test(hypermut.PPI.zitnik$BetweenessCentrality, hypermut.PPI.zitnik$density)
cor.test(hypermut.PPI.zitnik$EigenvectorCentrality, hypermut.PPI.zitnik$density)
cor.test(hypermut.PPI.zitnik$Degree, hypermut.PPI.zitnik$density)
cor.test(hypermut.PPI.zitnik$DegreeCentrality, hypermut.PPI.zitnik$density)
## This last one is not significant.
cor.test(hypermut.PPI.zitnik$IsArticulationPoint, hypermut.PPI.zitnik$density)

## now look at Cong network.
cor.test(nonmut.PPI.cong$Pagerank, nonmut.PPI.cong$density) ## significant
cor.test(nonmut.PPI.cong$HubScore, nonmut.PPI.cong$density) ## NS
cor.test(nonmut.PPI.cong$AuthorityScore, nonmut.PPI.cong$density) ## NS
cor.test(nonmut.PPI.cong$ClosenessCentrality, nonmut.PPI.cong$density) ## NS
cor.test(nonmut.PPI.cong$BetweenessCentrality, nonmut.PPI.cong$density) ## NS
cor.test(nonmut.PPI.cong$EigenvectorCentrality, nonmut.PPI.cong$density) ## NS
cor.test(nonmut.PPI.cong$Degree, nonmut.PPI.cong$density) ## significant
cor.test(nonmut.PPI.cong$DegreeCentrality, nonmut.PPI.cong$density) ## significant
cor.test(nonmut.PPI.cong$IsArticulationPoint, nonmut.PPI.cong$density) ## NS


cor.test(hypermut.PPI.cong$Pagerank, hypermut.PPI.cong$density) ## NS
cor.test(hypermut.PPI.cong$HubScore, hypermut.PPI.cong$density) ## significant
cor.test(hypermut.PPI.cong$AuthorityScore, hypermut.PPI.cong$density) ## sig
cor.test(hypermut.PPI.cong$ClosenessCentrality, hypermut.PPI.cong$density) #sig
cor.test(hypermut.PPI.cong$BetweenessCentrality, hypermut.PPI.cong$density) ## NS
cor.test(hypermut.PPI.cong$EigenvectorCentrality, hypermut.PPI.cong$density) ## sig
cor.test(hypermut.PPI.cong$Degree, hypermut.PPI.cong$density) ## significant
cor.test(hypermut.PPI.cong$DegreeCentrality, hypermut.PPI.cong$density) ## significant

cor.test(hypermut.PPI.cong$IsArticulationPoint, hypermut.PPI.cong$density) ## NS

####### Let's specifically ask about knockout mutations.
nonmut.KO.PPI.zitnik <- zitnik.network.df %>% left_join(nonmut.KO.density) %>%
    ## turn NAs to 0s.
    mutate(mut.count=ifelse(is.na(mut.count),0,mut.count)) %>%
    mutate(density=ifelse(is.na(density),0,density))

hypermut.KO.PPI.zitnik <- zitnik.network.df %>% left_join(hypermut.KO.density) %>%
    ## turn NAs to 0s.
    mutate(mut.count=ifelse(is.na(mut.count),0,mut.count)) %>%
    mutate(density=ifelse(is.na(density),0,density))

nonmut.KO.PPI.cong <- cong.network.df %>% left_join(nonmut.KO.density) %>%
    ## turn NAs to 0s.
    mutate(mut.count=ifelse(is.na(mut.count),0,mut.count)) %>%
    mutate(density=ifelse(is.na(density),0,density))

hypermut.KO.PPI.cong <- cong.network.df %>% left_join(hypermut.KO.density) %>%
    ## turn NAs to 0s.
    mutate(mut.count=ifelse(is.na(mut.count),0,mut.count)) %>%
    mutate(density=ifelse(is.na(density),0,density))

######## test for KO mutation density with Zitnik PPI data.
## None are significant for non-mutators.
cor.test(nonmut.KO.PPI.zitnik$Pagerank, nonmut.KO.PPI.zitnik$density)
cor.test(nonmut.KO.PPI.zitnik$HubScore, nonmut.KO.PPI.zitnik$density)
cor.test(nonmut.KO.PPI.zitnik$AuthorityScore, nonmut.KO.PPI.zitnik$density)
cor.test(nonmut.KO.PPI.zitnik$ClosenessCentrality, nonmut.KO.PPI.zitnik$density)
cor.test(nonmut.KO.PPI.zitnik$BetweenessCentrality, nonmut.KO.PPI.zitnik$density)
cor.test(nonmut.KO.PPI.zitnik$EigenvectorCentrality, nonmut.KO.PPI.zitnik$density)
cor.test(nonmut.KO.PPI.zitnik$Degree, nonmut.KO.PPI.zitnik$density)
cor.test(nonmut.KO.PPI.zitnik$DegreeCentrality, nonmut.KO.PPI.zitnik$density)
cor.test(nonmut.KO.PPI.zitnik$IsArticulationPoint, nonmut.KO.PPI.zitnik$density)

## These are all significant for hypermutators.
cor.test(hypermut.KO.PPI.zitnik$Pagerank, hypermut.KO.PPI.zitnik$density)
cor.test(hypermut.KO.PPI.zitnik$HubScore, hypermut.KO.PPI.zitnik$density)
cor.test(hypermut.KO.PPI.zitnik$AuthorityScore, hypermut.KO.PPI.zitnik$density)
cor.test(hypermut.KO.PPI.zitnik$ClosenessCentrality, hypermut.KO.PPI.zitnik$density)
cor.test(hypermut.KO.PPI.zitnik$BetweenessCentrality, hypermut.KO.PPI.zitnik$density)
cor.test(hypermut.KO.PPI.zitnik$EigenvectorCentrality, hypermut.KO.PPI.zitnik$density)
cor.test(hypermut.KO.PPI.zitnik$Degree, hypermut.KO.PPI.zitnik$density)
cor.test(hypermut.KO.PPI.zitnik$DegreeCentrality, hypermut.KO.PPI.zitnik$density)
## not significant
cor.test(hypermut.KO.PPI.zitnik$IsArticulationPoint, hypermut.KO.PPI.zitnik$density)

## now look at KOs in Cong network.
## significant negative correlation for closeness centrality,
## NS for the rest.
cor.test(nonmut.KO.PPI.cong$Pagerank, nonmut.KO.PPI.cong$density)
cor.test(nonmut.KO.PPI.cong$HubScore, nonmut.KO.PPI.cong$density)
cor.test(nonmut.KO.PPI.cong$AuthorityScore, nonmut.KO.PPI.cong$density) 
cor.test(nonmut.KO.PPI.cong$ClosenessCentrality, nonmut.KO.PPI.cong$density) 
cor.test(nonmut.KO.PPI.cong$BetweenessCentrality, nonmut.KO.PPI.cong$density) 
cor.test(nonmut.KO.PPI.cong$EigenvectorCentrality, nonmut.KO.PPI.cong$density) 
cor.test(nonmut.KO.PPI.cong$Degree, nonmut.KO.PPI.cong$density) 
cor.test(nonmut.KO.PPI.cong$DegreeCentrality, nonmut.KO.PPI.cong$density)
cor.test(nonmut.KO.PPI.cong$IsArticulationPoint, nonmut.KO.PPI.cong$density)

## significant for some but not all.
cor.test(hypermut.KO.PPI.cong$Pagerank, hypermut.KO.PPI.cong$density)
cor.test(hypermut.KO.PPI.cong$HubScore, hypermut.KO.PPI.cong$density)
cor.test(hypermut.KO.PPI.cong$AuthorityScore, hypermut.KO.PPI.cong$density)
cor.test(hypermut.KO.PPI.cong$ClosenessCentrality, hypermut.KO.PPI.cong$density)
cor.test(hypermut.KO.PPI.cong$BetweenessCentrality, hypermut.KO.PPI.cong$density)
cor.test(hypermut.KO.PPI.cong$EigenvectorCentrality, hypermut.KO.PPI.cong$density)
cor.test(hypermut.KO.PPI.cong$Degree, hypermut.KO.PPI.cong$density)
cor.test(hypermut.KO.PPI.cong$DegreeCentrality, hypermut.KO.PPI.cong$density)

cor.test(hypermut.KO.PPI.cong$IsArticulationPoint, hypermut.KO.PPI.cong$density)

##################################
## what are the essential genes that are not even in the
## PPI networks?

## 541 essential genes in REL606, total.
essential.not.in.cong <- essential.genes %>%
    filter(!(Gene %in% cong.network.df$Gene))
## 250 essential genes not in Cong network.
essential.not.in.zitnik <- essential.genes %>%
    filter(!(Gene %in% zitnik.network.df$Gene))
## 56 essential genes not in Zitnik network.

## these still show evidence of purifying selection in the LTEE.
essential.not.cong.data <- gene.mutation.data %>%
    filter(Gene %in% essential.not.in.cong$Gene)
c.essential.not.cong <- calc.cumulative.muts(essential.not.cong.data)
essential.not.cong.base.layer <- plot.base.layer(gene.mutation.data,
                                        subset.size=length(unique(essential.not.in.cong$Gene)))
essential.not.cong.STIMS.fig <- essential.not.cong.base.layer %>%
    add.cumulative.mut.layer(c.essential.not.cong, my.color="black")
ggsave("../results/thermostability/essential-not-cong-STIMS.pdf")


## these still show evidence of purifying selection in the LTEE.
essential.not.zitnik.data <- gene.mutation.data %>%
    filter(Gene %in% essential.not.in.zitnik$Gene)
c.essential.not.zitnik <- calc.cumulative.muts(essential.not.zitnik.data)
essential.not.zitnik.base.layer <- plot.base.layer(gene.mutation.data,
                                        subset.size=length(unique(essential.not.in.zitnik$Gene)))
essential.not.zitnik.STIMS.fig <- essential.not.zitnik.base.layer %>%
    add.cumulative.mut.layer(c.essential.not.zitnik, my.color="black")
ggsave("../results/thermostability/essential-not-zitnik-STIMS.pdf")

## IMPORTANT TODO: Use STIMS approach, to resample essential genes within
## the PPI network. HYPOTHESIS: the essential genes OUTSIDE of the PPI network
## are actually under significantly stronger purifying selection that those
## essential genes within the PPI network!

################################################################################
## Analyze data from the ProteomeVis database (the data specific to E. coli).
##########################################################################
## python ProteomeVis-to-REL606.py produces the table we need to join with REL606.genes.
REL606.to.ProteomeVis.df <- read.csv("../results/thermostability/REL606-to-ProteomeVis.csv",
                                     as.is=TRUE, header=TRUE)

## Import relevant E. coli ProteomeVis data. See Razban et al. (2018) in Bioinformatics:
## "ProteomeVis: a web app for exploration of protein properties from structure to sequence evolution across organisms’ proteomes"

## available protein data. 1262 E. coli proteins here.
proteome.vis.inspect.df <- read.csv("../results/thermostability/Ecoli-ProteomeVis-data/Ecoli-proteomevis_inspect.csv", as.is=TRUE, header=TRUE)
proteome.vis.chain.df <- read.csv("../results/thermostability/Ecoli-ProteomeVis-data/Ecoli-proteomevis_chain.csv", as.is=TRUE, header=TRUE)
proteome.vis.df <- full_join(proteome.vis.chain.df,proteome.vis.inspect.df) %>%
    left_join(REL606.to.ProteomeVis.df) %>%
    filter(!(is.na(Gene))) ## 1237 genes pass filters. 

## just looking at nonmut data, skipping dS.
LTEE.nonmut.proteome.vis.comp.df <- inner_join(proteome.vis.df,nonmut.density)

## significant positive correlation with abundant proteins.
cor.test(LTEE.nonmut.proteome.vis.comp.df$abundance,
         LTEE.nonmut.proteome.vis.comp.df$density)

### significant negative correlation between evolutionary rates in nature
### and rates in the LTEE, as reported in Maddamsetti et al. (2017).
cor.test(LTEE.nonmut.proteome.vis.comp.df$evolutionary_rate,
         LTEE.nonmut.proteome.vis.comp.df$density)

## negative correlation here with contact density, but why?
cor.test(LTEE.nonmut.proteome.vis.comp.df$contact_density,
         LTEE.nonmut.proteome.vis.comp.df$density)

## no correlation with PPI when only looking at non-mutator data. only when
## looking at hypermutator data.
cor.test(LTEE.nonmut.proteome.vis.comp.df$PPI_degree,
         LTEE.nonmut.proteome.vis.comp.df$density)
####################
## just looking at hypermut data, skipping dS.
LTEE.hypermut.proteome.vis.comp.df <- inner_join(proteome.vis.df,hypermut.density)

## no correlation with abundant proteins.
cor.test(LTEE.hypermut.proteome.vis.comp.df$abundance,
         LTEE.hypermut.proteome.vis.comp.df$density)

### no correlation between evolutionary rates in nature
### and rates in hypermutators.
cor.test(LTEE.hypermut.proteome.vis.comp.df$evolutionary_rate,
         LTEE.hypermut.proteome.vis.comp.df$density)

## again, a negative correlation here with contact density. but why?
cor.test(LTEE.hypermut.proteome.vis.comp.df$contact_density,
         LTEE.hypermut.proteome.vis.comp.df$density)
## see PNAS paper: Protein misinteraction avoidance causes
## highly expressed proteins to evolve slowly


## super strong negative correlation with PPI when
## looking at hypermutator data.
cor.test(LTEE.hypermut.proteome.vis.comp.df$PPI_degree,
         LTEE.hypermut.proteome.vis.comp.df$density)

####################
## using all LTEE data.
LTEE.proteome.vis.comp.df <- inner_join(proteome.vis.df, gene.mutation.densities)

## significant correlation between rate and abundance as reported in Razban et al. 2018
## and other papers.
cor.test(LTEE.proteome.vis.comp.df$abundance,
         LTEE.proteome.vis.comp.df$evolutionary_rate)

#######################################################################
## METABOLIC ENZYME ANALYSIS.
########################################################

## Only Ara-1 and Ara+6 show evidence of purifying selection on
## superessential metabolic enzymes. Why only these two? No idea why.

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
## look at generalist and specialist enzymes in Nam et al. (2012)
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

#############################
## ask about melting temperature of essential genes from Couce data set,
## versus those with no mutations.


## Import thermostability data for the E. coli proteome, from Leuenberger et al. (2017)
## in Science. 757 proteins here.
LeuenbergerS3.df <- read.csv("../results/thermostability/Ecoli-Leuenberger-data/Ecoli-Leuenberger_Table-S3.csv",as.is=TRUE, header=TRUE) %>%
    select(Tm.Protein,Length,Protein_ID,Protinfo,is.disordered,Protein.Abundance,
           Essential,Synthesis.Rate,Synthesis.Rate..Min.,T.90..Unfolded) %>%
    distinct() %>% ## examine proteins rather than peptides/domains.
    ## set uniprot column to be compatible with ProteomeViz
    rename(uniprot = Protein_ID) %>%
    left_join(proteome.vis.df)

## IMPORTANT TODO: NOT ALL OF THESE GENES MAP TO ProteomeViz!
## I will have to find a separate solution to match these to REL606 genes.
Leuentest <- LeuenbergerS3.df %>% filter(!(is.na(Gene)))
## only 442 out of 757 are in ProteomeViz.

## as another solution, use the supplementary data from Razban (2019).
## get 577 genes with abundance, Tm, and evolutionary rates.
Razban2019.df <- read.csv("../results/thermostability/Ecoli-Razban2019.csv") %>%
    left_join(REL606.genes)

nonmut.thermo.df <- Razban2019.df %>% inner_join(nonmut.density)
hypermut.thermo.df <- Razban2019.df %>% inner_join(hypermut.density)

## significant negative correlation with non-mutators.
cor.test(nonmut.thermo.df$evolutionary_rate_seq_identity,
         nonmut.thermo.df$density)
## really significant POSITIVE correlation with hypermutators!
cor.test(hypermut.thermo.df$evolutionary_rate_seq_identity,
         hypermut.thermo.df$density)

## positive correlation with abundance in non-mutators.
cor.test(nonmut.thermo.df$abundance_absolute_counts,
         nonmut.thermo.df$density)
## no correlation with abundance in hypermutators.
cor.test(hypermut.thermo.df$abundance_absolute_counts,
         hypermut.thermo.df$density)

## positive correlation with melting temperature in non-mutators.
cor.test(nonmut.thermo.df$melting_temperature_Celsius,
         nonmut.thermo.df$density)
## no correlation with melting temperature in hypermutators.
cor.test(hypermut.thermo.df$melting_temperature_Celsius,
         hypermut.thermo.df$density)

############################################################
## compare melting temperature and abundance for
## essential genes in REL606 and genes with no dS in LTEE.

## IMPORTANT TODO: ask whether results depend on including or excluding
## essential genes that are under selection.

nonmut.essential.thermo.df <- right_join(essential.genes,nonmut.thermo.df) %>%
    ## turn NAs to 0s.
    mutate(mut.count=ifelse(is.na(mut.count),0,mut.count)) %>%
    mutate(density=ifelse(is.na(density),0,density))

hypermut.essential.thermo.df <- right_join(essential.genes,hypermut.thermo.df) %>%
        ## turn NAs to 0s.
    mutate(mut.count=ifelse(is.na(mut.count),0,mut.count)) %>%
    mutate(density=ifelse(is.na(density),0,density))

## negative correlation between conservation and mutation density.
cor.test(nonmut.essential.thermo.df$evolutionary_rate_seq_identity,
         nonmut.essential.thermo.df$density)

## negative correlation between conservation and mutation density.
cor.test(hypermut.essential.thermo.df$evolutionary_rate_seq_identity,
         hypermut.essential.thermo.df$density)

## positive correlation between abundance and mutation density
cor.test(nonmut.essential.thermo.df$abundance_absolute_counts,
         nonmut.essential.thermo.df$density)

## no correlation between abundance and mutation density
cor.test(hypermut.thermo.df$abundance_absolute_counts,
         hypermut.thermo.df$density)

## positive correlation between melting temperature and mutation density
cor.test(nonmut.thermo.df$melting_temperature_Celsius,
         nonmut.thermo.df$density)

## no correlation between melting temperature and mutation density
cor.test(hypermut.thermo.df$melting_temperature_Celsius,
         hypermut.thermo.df$density)

############################################################################
## E. COLI MELTOME ATLAS DATA ANALYSIS

## TODO: put this stuff in a more sensible place in this script.

## These data were written out by filter-meltome-atlas.R.
## Run that script in order to generate '../results/thermostability/Ecoli-meltome.csv'
Ecoli.meltome <- read.csv("../results/thermostability/Ecoli-meltome.csv") %>%
    rename(Gene = gene_name) %>%
    select(run_name, Gene, meltPoint) %>%
    distinct() %>%
    left_join(REL606.genes) %>%
    filter(!is.na(locus_tag))

nonmut.density.meltome <- Ecoli.meltome %>%
    left_join(nonmut.density) %>%
    ## turn NAs into zeros.
    mutate(mut.count = replace(mut.count, is.na(mut.count), 0)) %>%
    mutate(density = replace(density, is.na(density), 0))

hypermut.density.meltome <- Ecoli.meltome %>%
    left_join(hypermut.density) %>%
    ## turn NAs into zeros.
    mutate(mut.count = replace(mut.count, is.na(mut.count), 0)) %>%
    mutate(density = replace(density, is.na(density), 0))

nonmut.meltome.plot <- ggplot(nonmut.density.meltome,
                                 aes(x=density,y=meltPoint)) +
    geom_point() +
    theme_classic() +
    facet_wrap(run_name ~ .) +
    geom_smooth()

hypermut.meltome.plot <- ggplot(hypermut.density.meltome,
                                 aes(x=density,y=meltPoint)) +
    geom_point() +
    theme_classic() +
    facet_wrap(run_name ~ .) +
    geom_smooth()

nonmut.meltome.plot
hypermut.meltome.plot

nonmelters <- Ecoli.meltome %>% filter(is.na(meltPoint)) %>%
    select(-run_name) %>% distinct() %>%
    ungroup()
## are non-melters under purifying selection?
## let's use STIMS.
nonmelter.mut.data <- gene.mutation.data %>%
    filter(Gene %in% nonmelters$Gene)
c.nonmelters <- calc.cumulative.muts(nonmelter.mut.data)
nonmelter.base.layer <- plot.base.layer(gene.mutation.data,
                                        subset.size=length(unique(nonmelters$Gene)))

nonmelter.STIMS.fig <- nonmelter.base.layer %>%
    add.cumulative.mut.layer(c.nonmelters, my.color="black")
ggsave("../results/thermostability/nonmelters-STIMS.pdf")

## unlike what I expected: more mutations than average in Ara+3,
## and in Ara-1?

hypermut.density.meltome2 <- hypermut.density.meltome %>%
    mutate(is.nonmelter = ifelse(Gene %in% nonmelters$Gene,TRUE, FALSE)) %>%
    mutate(is.essential = ifelse(Gene %in% essential.genes$Gene,TRUE, FALSE))

nonmelter.vs.melter.hist <- ggplot(hypermut.density.meltome2,
                                   aes(x=density,fill=is.nonmelter)) +
    geom_histogram()

g1 <- filter(hypermut.density.meltome2,is.nonmelter==TRUE)$density
g2 <- filter(hypermut.density.meltome2,is.nonmelter==FALSE)$density
mean(g1) ## on average, nonmelters have a higher mean.
mean(g2)
## statistically significant.
wilcox.test(g1,g2)

## Is there an association between essentiality and being a nonmelter?
## let's use a contingency test.
essentiality.nonmelter.data <- REL606.genes %>%
    mutate(is.nonmelter = ifelse(Gene %in% nonmelters$Gene,TRUE, FALSE)) %>%
    mutate(is.essential = ifelse(Gene %in% essential.genes$Gene,TRUE, FALSE))
x1 <- nrow(filter(essentiality.nonmelter.data,is.nonmelter & is.essential))
x2 <- nrow(filter(essentiality.nonmelter.data,is.nonmelter & !is.essential))
x3 <- nrow(filter(essentiality.nonmelter.data,!is.nonmelter & is.essential))
x4 <- nrow(filter(essentiality.nonmelter.data,!is.nonmelter & !is.essential))

## p = 0.0028. positive association between being a nonmelter and essentiality.
fisher.test(matrix(c(x1,x2,x3,x4),2))

## what if we remove ribosomal proteins?
essen.nonmelt2 <- essentiality.nonmelter.data %>%
    filter(!str_detect(product,"ribosomal"))

y1 <- nrow(filter(essen.nonmelt2,is.nonmelter & is.essential))
y2 <- nrow(filter(essen.nonmelt2,is.nonmelter & !is.essential))
y3 <- nrow(filter(essen.nonmelt2,!is.nonmelter & is.essential))
y4 <- nrow(filter(essen.nonmelt2,!is.nonmelter & !is.essential))
## p = 0.1. association goes away when we remove ribosomal proteins.
fisher.test(matrix(c(y1,y2,y3,y4),2))

## what if we remove hypothetical proteins?
essen.nonmelt3 <- essen.nonmelt2 %>%
    filter(!str_detect(product,"hypothetical"))

z1 <- nrow(filter(essen.nonmelt3,is.nonmelter & is.essential))
z2 <- nrow(filter(essen.nonmelt3,is.nonmelter & !is.essential))
z3 <- nrow(filter(essen.nonmelt3,!is.nonmelter & is.essential))
z4 <- nrow(filter(essen.nonmelt3,!is.nonmelter & !is.essential))
## p = 0.0376. significant. seems like hypothetical/unstructured proteins
## in here.
fisher.test(matrix(c(z1,z2,z3,z4),2))

## do a STRING analysis of the nonmelters.
## print out to file.
write.csv(nonmelters,file="../results/thermostability/nonmelters.csv")
## STRING shows that oxidative phosphorylation, and ribosome are
## enriched KEGG pathways for non-melters.

## plot essentiality against melting temperature.
meltome.with.essentiality <- Ecoli.meltome %>%
    ## for plotting, give non-melters a Tm = 100.
    mutate(Tm = ifelse(is.na(meltPoint), 100, meltPoint)) %>%
    mutate(is.essential = ifelse(Gene %in% essential.genes$Gene,TRUE, FALSE))

## This plot is interesting: excluding non-melters,
## essential genes have a slightly lower Tm, on average.
essen.Tm.plot <- ggplot(meltome.with.essentiality,
                        aes(x=Tm,fill=is.essential)) +
    geom_histogram()

## Now compare mRNA and protein abundance at each timepoint to
## Tm.
meltome.with.abundance <- meltome.with.essentiality %>%
    left_join(Caglar.summary)

## exclude non-melters for now.
correlate.meltPoint.with.timepoint <- function(Tm.Caglar) {
    for (t in sort(unique(Tm.Caglar$growthTime_hr))) {
        my.data <- Tm.Caglar %>%
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

## IMPORTANT RESULT: melting temperature, excluding non-melters,
## is negatively correlated with protein and RNA abundance,
## especially during growth phase, but also during starvation.
correlate.meltPoint.with.timepoint(meltome.with.abundance)
## IMPORTANT TODO: make a figure for this fascinating result!

