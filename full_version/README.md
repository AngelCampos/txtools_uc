txtools use cases
================
Miguel A. García-Campos
2024-02-04

- [Background](#background)
- [Setup](#setup)
- [Use case \#1. Pseudouridine in yeast rRNA (Carlile et al.,
  2014)](#use-case-1-pseudouridine-in-yeast-rrna-carlile-et-al-2014)
- [Use case \#2. m6A miCLIP2 in mESC (Körtel et al.,
  2021)](#use-case-2-m6a-miclip2-in-mesc-körtel-et-al-2021)
- [Use case \#3. ac4C in the archea T. kodakarensis (Sas-Chen et al.,
  2020)](#use-case-3-ac4c-in-the-archea-t-kodakarensis-sas-chen-et-al-2020)
- [Supplementary analysis - Bridging paired-end
  reads](#supplementary-analysis---bridging-paired-end-reads)
- [References](#references)
- [Session Info](#session-info)

# Background

[***txtools***](https://github.com/AngelCampos/txtools) is an R package
that enables the processing, analysis, and visualization of RNA-seq data
at the nucleotide-level resolution, seamlessly integrating alignments to
the genome with transcriptomic representation. txtools’ main inputs are
BAM files and a transcriptome annotation, and the main output is a
table, capturing mismatches, deletions, and the number of reads
beginning and ending at each nucleotide in the transcriptomic space.
txtools further facilitates downstream visualization and analyses. We
showcase, using examples from the epitranscriptomic field, how a few
calls to txtools functions can yield insightful and ready-to-publish
results. txtools is of broad utility also in the context of structural
mapping and RNA:protein interaction mapping. By providing a simple and
intuitive framework, we believe that txtools will be a useful and
convenient tool and pave the path for future discovery. txtools is
available for installation from its GitHub repository at:
<https://github.com/AngelCampos/txtools>.

To show the functionality of ***txtools*** in real life scenarios we
will use external data as use cases. We use the following studies’s data
and show how to reach some of the main results the original authors did.

| RNAmod        | Method  | Paper                  | GSE       | Target | Organism                         |     |
|---------------|---------|------------------------|-----------|--------|----------------------------------|-----|
| Pseudouridine | CMC     | (Carlile et al. 2014)  | GSE58200  | rRNA   | *Saccharomyces cerevisiae* (Sk1) |     |
| m6A           | miCLIP  | (Körtel et al. 2021)   | GSE163500 | mRNA   | *Mus musculus*                   |     |
| ac4C          | NaCNBH3 | (Sas-Chen et al. 2020) | GSE135826 | rRNA   | *Thermococcus kodakarensis*      |     |

# Setup

## Packages

``` r
library("txtools")
# Supplementary packages
library("GEOfastq") # Obtain FASTQ files and metadata from SRA/ENA databases
library("Rsubread") # Align FASTQ files to reference genome
library("magrittr") # Pipe operator
library("data.table") # data.table class
library("ggplot2") # Additional plot tweaking
library("gridExtra") # Arrange plots in one panel
library("dplyr") # Data.frame manipulation
library("stringr") # Regex capture
```

## Global variables

``` r
NCORES <- 10 # Number of cores to be used in all multi-core processes
```

## Local functions

``` r
# Shortcut function: Get GSE tables using GEOfastq
get_GSEtable <- function(GSE){
    require(GEOfastq)
    require(magrittr)
    crawl_gse(GSE) %>% extract_gsms() %>% crawl_gsms()
}
```

## Directories

``` r
nDirs <- c("./figs", "./results", "./omicRefs", "./data")  
sapply(nDirs[!nDirs %in% list.dirs()], "dir.create") %>% invisible()
```

## Download omic references

``` r
download.file(url = "https://zenodo.org/records/8278410/files/uc1_geneAnnot.bed", destfile = "omicRefs/uc1_geneAnnot.bed")
download.file(url = "https://zenodo.org/records/8278410/files/uc1_genome.fa", destfile = "omicRefs/uc1_genome.fa")
download.file(url = "https://zenodo.org/record/8278045/files/uc1_RNAmods.rds", destfile = "omicRefs/uc1_RNAmods.rds")
download.file(url = "https://zenodo.org/records/8278410/files/uc2_geneAnnot.bed", destfile = "omicRefs/uc2.geneAnnot.bed")
download.file(url = "https://zenodo.org/records/8278410/files/uc3_geneAnnot.bed", destfile = "omicRefs/uc3_geneAnnot.bed")
download.file(url = "https://zenodo.org/records/8278410/files/uc3_genome.fa", destfile = "omicRefs/uc3_genome.fa")
pathTomm9Genome = NULL # Please add here a path to a local copy of the mm9 reference genome. Not included in downloads due to size and potential redundancy.
if(is.null(pathTomm9Genome)){stop("Path to a FASTA file of the mm9 genome is required")}
```

# Use case \#1. Pseudouridine in yeast rRNA (Carlile et al., 2014)

## Download data

For the first use case, we analyzed pseudouridine mapping data, acquired
via Pseudo-seq, a method for genome-wide, single-nucleotide resolution
identification of pseudouridine (Carlile et al. 2014). Carlile and
collaborators employed N-cyclohexyl-N′-(2-morpholinoethyl)-carbodiimide
metho-p-toluenesulphonate (CMC) which reacts with pseudouridines
creating an adduct that blocks cDNA’s reverse transcription one
nucleotide downstream to the pseudouridylated site.

The FASTQ data of a CMC-treated sample and a control sample were
downloaded from GEO: GSE58200, and aligned to the yeast ribosome
reference sequence (Taoka et al. 2016).

``` r
# Omics references
sc_faGenome <- "omicRefs/uc1_genome.fa"
sc_geneAnno <- "omicRefs/uc1_geneAnnot.bed"
sc_TaokaRNAmods <- "omicRefs/uc1_RNAmods.rds"  # Taoka et al., 2016
# Get Metadata
GSEtab_GSE58200 <- get_GSEtable("GSE58200") %>% subset(run %in% c("SRR1327248", "SRR1327249"))
META_GSE58200 <- data.table(run = GSEtab_GSE58200$run,
                            title = GSEtab_GSE58200$title,
                            FASTQ = paste0(GSEtab_GSE58200$run, ".fastq.gz"),
                            BAM = paste0(GSEtab_GSE58200$run, ".bam"))
# Download FASTQ files
FASTQDir <- "data/fastq/GSE58200"; dir.create(FASTQDir, recursive = TRUE)
BAMDir <- "data/bam/GSE58200"; dir.create(BAMDir, recursive = TRUE)
```

``` r
get_fastqs(GSEtab_GSE58200, "data/fastq/GSE58200")
```

## Read mapping

Reads are mapped to the reference sequence using Rsubread (Liao, Smyth,
and Shi 2019).

``` r
# Align to reference
buildindex(basename = "sc_rib", reference = sc_faGenome)
lapply(META_GSE58200$FASTQ, function(FASTQ){
    align(index = "sc_rib",
          readfile1 = file.path(FASTQDir, FASTQ),
          output_file = file.path(BAMDir, gsub(x = FASTQ, pattern = "fastq.gz", replacement = "bam")),
          phredOffset = 64,
          nthreads = NCORES)
})
file.remove(list.files(pattern = "sc_rib.*")) # Remove index files
```

## txtools processing

Below are the steps and txtools functions used to analyze the data.

1.  Load the reference genome with tx_load_genome(), and the gene
    annotation with tx_load_bed().
2.  Loop through the resulting BAM file for loading with tx_load_bam();
    processing into transcriptomic reads with tx_reads(); create a
    summarized txDT with coverage, read-start, and read-ends data;
    calculate the start-ratio 1bp down-stream (SR_1bpDS) with
    tx_add_startRatio1bpDS(); and finally add a position identifier that
    merges gene and transcriptomic position with tx_add_pos().
3.  Join the SR_1bpDSs results into one txDT and subtract the one of the
    control from that of the CMC-treated sample, yielding the start-rate
    difference 1bp down-stream (SRD_1bpDS).
4.  Join with the Taoka RNA modifications table.
5.  Visualize the processed data with tx_plot_staEndCov() and observe
    the effect of CMC treatment on RT premature stoppage at a known
    pseudouridylated site, which manifests as an abrupt increase of
    read-starts compared to the control sample (Figure 1, top, and
    middle panel). a metric that we can use to detect pseudouridines in
    CMC-treatment data.
6.  Plot all the results using the tx_plot_numeric() function to plot
    numeric variables in a txDT along the same windown of the known
    pseudouridylated sites (Figure 1, bottom panel).
7.  Using the resulting txDT and a common call to ggplot2’s scatter
    plot, the results for all four ribosomal transcripts are shown in
    Figure 2.

``` r
# Load omic references
GENOME <- tx_load_genome(sc_faGenome)
TXOME <- tx_load_bed(sc_geneAnno)

sc_txDTL <- lapply(file.path(BAMDir, META_GSE58200$BAM), function(bam){
    tx_load_bam(bam, pairedEnd = FALSE, verbose = FALSE) %>% 
        tx_reads(geneAnnot = TXOME, minReads = 1, nCores = NCORES, verbose = FALSE) %>% 
        tx_makeDT_coverage(geneAnnot = TXOME, genome = GENOME, nCores = NCORES) %>% 
        tx_add_startRatio1bpDS(minCov = 50) %>% 
        tx_add_pos()
})
names(sc_txDTL) <- c("control", "CMC")
```

``` r
# Join control and treatment results in one txDT for easier plotting
sc_txDTL$RES <- sc_txDTL[[1]][,1:7]
sc_txDTL$RES$SR_1bpDS_CMC <- sc_txDTL$CMC$startRatio1bpDS
sc_txDTL$RES$SR_1bpDS_ctrl <- sc_txDTL$control$startRatio1bpDS
sc_txDTL$RES$SRD_1bpDS <- sc_txDTL$RES$SR_1bpDS_CMC - sc_txDTL$RES$SR_1bpDS_ctrl
# Joining Taoka's data set (RNA modifications identity)
sc_RNAmods <- readRDS(sc_TaokaRNAmods)
sc_txDTL$RES <- dplyr::left_join(sc_txDTL$RES, dplyr::select(sc_RNAmods, pos, nuc), by = "pos")
sc_txDTL$RES$nucleotide <- "N"
sc_txDTL$RES$nucleotide[sc_txDTL$RES$nuc == "Y"] <- "Y"
sc_txDTL$RES$nucleotide[sc_txDTL$RES$nuc == "Ym"] <- "Ym"
```

## Results and plots

``` r
gg_Y_1 <- tx_plot_staEndCov(sc_txDTL$control, gene = "28s", show_yLabels = F,
                            txRange = window_around(2973, 15), bar_border = F) +
    ggtitle("rRNA:28s - Control") + theme(legend.position="top") + xlab(NULL)
gg_Y_2 <- tx_plot_staEndCov(sc_txDTL$CMC, gene = "28s", show_yLabels = F,
                            txRange = window_around(2973, 15), bar_border = F) +
    ggtitle("rRNA:28s - CMC") + theme(legend.position="none") + xlab(NULL)
gg_Y_3 <- tx_plot_numeric(DT = sc_txDTL$RES, gene = "28s", txRange = window_around(2973, 15),
                          colVars = c("SRD_1bpDS", "SR_1bpDS_CMC", "SR_1bpDS_ctrl"),
                          plot_type = "lineplot") + ggtitle("rRNA:28s - StartRatios 1bp-downstream")
ggPsiU_1 <- gridExtra::grid.arrange(gg_Y_1, gg_Y_2, gg_Y_3, layout_matrix = matrix(c(1,1,1,2,2,3,3,3), ncol = 1))

ggsave("figs/pseudoU_1.png", ggPsiU_1, width = 6, height = 9, dpi = 200, bg = "gray98")
```

![](figs/pseudoU_1.png)

**Figure 1.** Use case \#1 - rRNA pseudouridylation. Top. txtools
Starts-Ends-Coverage plot for the mock-treatment sample centered at
position 2973 of the yeast 28s rRNA, a known pseudouridine site. Middle.
Same plot as above but for the CMC-treated sample. A dramatic difference
in read-starts is evident at 1 bp downstream from the pseudouridylated
site. Bottom. Lineplot showing SR_1bpDS for both CMC and control
treatment and the resulting difference, SRD_1bpDS, at 28s:2973 and
surrounding nucleotides.

Using txtools´ txDT and `ggplot2` we can easily visualize the data per
nucleotide and RNA modification identity according to MS-based
detection.

``` r
# Plot results
ggPsiU_2 <- ggplot(sc_txDTL$RES, aes(x = txcoor, y = SRD_1bpDS, colour = nucleotide)) + 
    geom_point(size = 2, alpha = 0.5) + facet_grid(.~gene, scales = "free_x") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom") +
    xlab("Transcriptome coordinate")
ggsave("figs/pseudoU_2.png", ggPsiU_2, width = 6, height = 3.5, dpi = 200, bg = "gray99")
```

![](figs/pseudoU_2.png)

**Figure 2.** Use case \#1 - Scatterplots of full rRNA transcripts
showing the SRD_1bpDS per nucleotide. Marked in green are the known
pseudouridylated sites and marked in blue is the sole 2’-O-methylated
pseudouridine.

The resulting SRD_1bpDS metric shows that a simple threshold can
discriminate between pseudouridine harboring and non-harboring sites on
rRNA.

# Use case \#2. m6A miCLIP2 in mESC (Körtel et al., 2021)

## Download data

For this case study we analyzed the data of the study that presented
miCLIP2 (Körtel et al. 2021), a miCLIP enhancement allowing
antibody-based single-nucleotide resolution mapping of m6A sites,
relying on crosslinking of an antibody to methylated sites (Linder et
al. 2015). Similarly to its predecessor, miCLIP2 relies on premature
termination of reverse transcription during cDNA synthesis at the
cross-linked residue. Data was downloaded from GEO: GSE163500, selecting
for samples of wild-type (WT) and methyltransferase- like 3 (Mettl3)
knockout in mouse embryonic stem cells (mESC).

``` r
# Omic references
mm_faGenome <- pathTomm9Genome
mm_geneAnno <- "omicRefs/uc2_geneAnnot.bed"
# GEO table - selecting samples
GSEtab_GSE163500 <- get_GSEtable("GSE163500") %>% 
    subset(grepl(pattern = "miCLIP mESC", title)) %>% 
    subset(!grepl(pattern = "technical rep2", title))
# Experimental design
META_mm9_eCLIP <- data.table(run = GSEtab_GSE163500$run,
                            title = GSEtab_GSE163500$title,
                            FASTQ = paste0(GSEtab_GSE163500$run, ".fastq.gz"),
                            BAM = paste0(GSEtab_GSE163500$run, ".bam"))
META_mm9_eCLIP$genotype <- stringr::str_extract(GSEtab_GSE163500$title, pattern = "(KO|WT)")
META_mm9_eCLIP$replicate <- stringr::str_extract(GSEtab_GSE163500$title, pattern = "rep([:digit:])")
META_mm9_eCLIP$id <- paste(META_mm9_eCLIP$genotype, META_mm9_eCLIP$replicate, sep = "_")

# Download FASTQ files
FASTQDir <- "data/fastq/GSE163500"; dir.create(FASTQDir, recursive = TRUE)
BAMDir <- "data/bam/GSE163500"; dir.create(BAMDir, recursive = TRUE)
```

``` r
get_fastqs(GSEtab_GSE163500, "data/fastq/GSE163500") # download FASTQS
```

## Read mapping

Reads are mapped to the reference sequence using Rsubread (Liao, Smyth,
and Shi 2019).

``` r
# Align to reference 
buildindex(basename = "mm9_refGenome", reference = mm_faGenome)
lapply(list.files(FASTQDir), function(FASTQ){
    align(index = "mm9_refGenome",
          readfile1 = file.path(FASTQDir, FASTQ),
          output_file = file.path(BAMDir, gsub(x = FASTQ, pattern = "fastq.gz", replacement = "bam")),
          phredOffset = 64,
          nthreads = NCORES)
    cat(FASTQ, "aligned.")
})
file.remove(list.files(pattern = "mm9_refGenome.*")) # Remove index files
```

## txtools processing

Below are the steps and txtools functions used to analyze the data.

1.  Load the reference genome with tx_load_genome(), and the gene
    annotation with tx_load_bed().
2.  Process all BAM files by looping through each with:
    1.  tx_load_bam(): Loading the BAM file
    2.  tx_reads(): Processing into transcriptomic reads
    3.  tx_makeDT_coverage(): Processing into a table with summarized
        data on coverage, read-start, and read-ends.
    4.  tx_add_startRatio1bpDS(): Adding the start ratio 1bp down-stream
3.  Unify all the resulting tables with tx_unifyTxDTL(), this is to have
    them share the same transcriptomic coordinates, by selecting the
    intersection of genes in all data tables.
4.  Perform t-tests using the txtools’ inbuilt ‘genefilter’ wrapper
    (Gentleman et al. 2019) tx_test_ttests(), testing for differences in
    startRatio_1bpDS between WT and KO samples both with the
    immunoprecipitation and crosslinking treatment.
5.  Add DRACH motif presence with tx_add_motifPresence().
6.  Plot a volcano plot using a call to ggplot2’s scatterplot,
    color-coding for DRACH presence, and select putative m6A sites at
    mean difference \> 0.05 and p-value \< 0.01 (Figure 3).
7.  Plot a metagene profile using tx_plot_metageneRegions(). (Figure 4)

``` r
# Load omic references
mm9_geneAnnot <- tx_load_bed(mm_geneAnno)
mm9_genome <- tx_load_genome(mm_faGenome)
# Processing loop
txDTL <- lapply(META_mm9_eCLIP$BAM, function(BAM){
    tx_load_bam(BAM, pairedEnd = FALSE, loadSeq = FALSE) %>% 
        tx_reads(geneAnnot = mm9_geneAnnot, minReads = 20, withSeq = FALSE, nCores = NCORES) %>% 
        tx_makeDT_coverage(geneAnnot = mm9_geneAnnot, genome = mm9_genome, nCores = NCORES) %>% 
        tx_add_startRatio1bpDS(minCov = 20)
})
txDTL <- tx_unifyTxDTL(txDTL, type = "intersection", nCores = NCORES) # Unify txDTs to have the same genes
names(txDTL) <- META_mm9_eCLIP$id
saveRDS(txDTL, file = "results/useCase_2.rds")
rm(mm9_genome); gc() # free memory

txRES <- tx_test_ttest(DTL = txDTL,       # Calculate pValues
                       cont_var = "startRatio1bpDS", 
                       test_groups = factor(META_mm9_eCLIP$genotype, levels = c("WT", "KO")),
                       test_na.rm = FALSE)
# Calling putative sites based on arbitrary thresholds
txRES$putativem6A <- txRES$dm > 0.05 & -log10(txRES$p.value) > -log10(0.01) & txRES$refSeq == "A"
txRES$putativem6A[is.na(txRES$putativem6A)] <- FALSE
# Annotating DRACH motif 
txRES <- tx_add_motifPresence(txRES, motif = "DRACH", nucPositions = 3, nCores = NCORES)
# Adding gene region annotation for faster metagene plotting
txRES <- tx_add_geneRegion(txRES, geneAnnot = mm9_geneAnnot, nCores = NCORES)
saveRDS(txRES, "results/useCase_2_RES.rds")
```

## Results and plots

``` r
# Volcano plot
gg_m6A1 <- ggplot(na.omit(txRES),
                  aes(x = dm, y = -log10(p.value), colour = DRACH_motif_3)) + 
    geom_point() + 
    geom_vline(xintercept = 0.05, colour = "black") + 
    geom_vline(xintercept = -0.05, colour = "black") + 
    geom_hline(yintercept = -log10(0.01), colour = "black") + 
    xlab("Mean StartRate Difference") + 
    theme_minimal() + theme(legend.position = "bottom") + xlim(c(-0.6, 0.6))
ggplot2::ggsave(gg_m6A1, filename = "figs/m6A_plot1.png", width = 5, height = 5, bg = "white")
```

![](figs/m6A_plot1.png)

**Figure 3.** Case study \#2 - m6A epitranscriptome in mESC using
miCLIP2. Volcano plot showing the 1bp-down-stream mean start-rate
difference at each queried position of the transcriptome (x-axis) and
the -log10 p-value, calculated using a t-test (y-axis). Colored in blue
are all sites that are centered in a DRACH motif.

``` r
ggtx_MGR1 <- tx_plot_metageneRegions(
    txDT = txRES, geneAnnot = mm9_geneAnnot, 
    colVars = c("putativem6A", "DRACH_motif_3"), nBins_5UTR = 20,
    nBins_CDS = 50, nBins_3UTR = 50, spar = 0.5, nCores = NCORES, normalize = TRUE) 
ggtx_MGR1 <- ggtx_MGR1 + ggtitle("m6A putative sites vs. DRACH motif distribution")
ggsave("figs/ggtx_MGR1.png", plot = ggtx_MGR1, width = 6, height = 4, bg = "gray98")
```

![](figs/ggtx_MGR1.png)

**Figure 4.** Case study \#2 - m6A epitranscriptome in mESC using
miCLIP2. Metagene plot aligned at the end of CDS. Showing the relative
abundance of putative m6A sites in blue, compared to the baseline
presence of the DRACH motif in red. Putative sites thresholds used were
startRatio_1bpDS difference \> 0.05, p-value \< 0.01, only considering
adenines.

``` r
GG_m6A_1 <- tx_plot_metageneAtCDS(txDT = txRES, geneAnnot = mm9_geneAnnot, 
                                  colVars = c("putativem6A", "DRACH_motif_3"), 
                                  CDS_align = "start", upFlank = 500, doFlank = 500, spar = 0.5, normalize = TRUE) + 
    theme(legend.position = "top")

GG_m6A_2 <- tx_plot_metageneAtCDS(txDT = txRES, geneAnnot = mm9_geneAnnot, 
                                  colVars = c("putativem6A", "DRACH_motif_3"), 
                                  CDS_align = "end", upFlank = 500, doFlank = 500, spar = 0.5, normalize = TRUE) +
    theme(legend.position = "none")

GG_m6A_3 <- tx_plot_metageneAtCDS(txDT = txRES, geneAnnot = mm9_geneAnnot, 
                                  colVars = c("putativem6A", "DRACH_motif_3"), 
                                  CDS_align = "spliceSite", upFlank = 500, doFlank = 500, spar = 0.5, normalize = TRUE) +
    theme(legend.position = "none")

ggsave(filename = "figs/m6A_plot_3.png", gridExtra::grid.arrange(GG_m6A_1, GG_m6A_2, GG_m6A_3, ncol = 3), width = 15, height = 5)
```

![](figs/m6A_plot_3.png)

**Figure 5.** Case study \#2 - m6A epitranscriptome in mESC using
miCLIP2. Additional plots of metagene analysis centered at the
CDS-start, CDS-end, and splice sites positions. These plots show data
not summarizing bins along the transcriptome but by overlapping the
genes at the reference position and softening the resulting line with a
spline.

# Use case \#3. ac4C in the archea T. kodakarensis (Sas-Chen et al., 2020)

For this case study we analyzed RNA acetylation data, acquired via
ac4C-seq, a chemical method for the transcriptome-wide quantitative
mapping of N4-acetylcytidine (ac4C) at single-nucleotide resolution
(Sas-Chen et al. 2020). Sas-Chen and collaborators employed the reaction
of ac4C with sodium cyanoborohydride (NaCNBH3) under acidic conditions,
which leads to C-\>T mutations at acetylated positions (Thomas et al.
2018). A key result in this study was the discovery of an abundance of
acetylation sites on rRNA derived from the hyperthermophilic archaea T.
kodakarensis. Data from GEO: GSE135826 was downloaded and aligned to the
T. kodakarensis genome.

## Download data

``` r
# Omic references
tk_faGenome <- "omicRefs/uc3_genome.fa"
tk_geneAnno <- "omicRefs/uc3_geneAnnot.bed"

# T. Kodakarensis
GSEtab_GSE135826 <- get_GSEtable("GSE135826") %>% 
    subset(organism_ch1 == "Thermococcus kodakarensis") %>% 
    subset(characteristics_ch1 == "genetic background: WT") %>%
    subset(characteristics_ch1.1 != "treatment: deacetylation") %>% 
    subset(grepl(pattern = "deg_", title)) %>% 
    subset(instrument_model == "Illumina NextSeq 500")
GSEtab_GSE135826 <- GSEtab_GSE135826[c(1:8, 13:14), ]

FASTQDir <- "data/fastq/GSE135826"; dir.create(FASTQDir, recursive = T)
BAMDir <- "data/bam/GSE135826"; dir.create(BAMDir, recursive = T)
```

``` r
get_fastqs(GSEtab_GSE135826, "data/fastq/GSE135826")
```

``` r
# Format meta data from experimental design
META_GSE135826 <- data.table(run = GSEtab_GSE135826$run,
                            title = GSEtab_GSE135826$title,
                            FASTQ = paste0(GSEtab_GSE135826$run, ".fastq.gz"),
                            BAM = paste0(GSEtab_GSE135826$run, ".bam"))
META_GSE135826 <- META_GSE135826[!META_GSE135826$run %in% c("SRR11178346", "SRR11178347"),]
META_GSE135826$temp <- str_extract(META_GSE135826$title, pattern = "..deg")
META_GSE135826$treat <- str_extract(META_GSE135826$title, pattern = "(NaCNBH3|mock)")
META_GSE135826$id <- with(META_GSE135826, paste(treat, temp, sep = "_"))
META_GSE135826$BAM <- paste0(file.path(BAMDir, META_GSE135826$run), ".bam")
META_GSE135826$RDS <- gsub(".bam$", ".txDT.rds", META_GSE135826$BAM)
META_GSE135826 <- META_GSE135826[order(META_GSE135826$temp, META_GSE135826$treat),]
```

## Read mapping

Reads are mapped to the reference sequence using Rsubread (Liao, Smyth,
and Shi 2019).

``` r
# Align to reference with Rsubread
buildindex(basename = "tk_RsubRef", reference = tk_faGenome)
r1Files <- grep(list.files(FASTQDir), pattern = "_1.fastq.gz", value = TRUE)
tk_alignReport <- lapply(r1Files, function(FASTQ){
    align(index = "tk_RsubRef",
          readfile1 = file.path(FASTQDir, FASTQ),
          readfile2 = file.path(FASTQDir, gsub(FASTQ, pattern = "_1.fastq.gz", replacement = "_2.fastq.gz")),
          output_file = file.path(BAMDir, gsub(x = FASTQ, pattern = "_1.fastq.gz", replacement = ".bam")),
          nthreads = NCORES, minFragLength = 23, sortReadsByCoordinates = TRUE)
})
file.remove(list.files(pattern = "tk_RsubRef.*")) # Remove index files
```

## txtools processing

Below are the steps and txtools functions used to analyze the data.

1.  Load the reference genome with tx_load_genome(), and the gene
    annotation with tx_load_bed().
2.  Process the BAM files with the bam2txDT.R script provided along the
    txtools installation, with the parameters “-p TRUE -d covNuc -r 300
    -m 0” to generate a summarized count data table with coverage,
    read-starts, read-ends, nucleotide frequency, and deletion frequency
    information.
3.  Plot a nucleotide frequency plot using the tx_plot_nucFreq()
    function. To observe the high levels of misincorporations of
    cytidines in place of thymines exclusively in the NaCNBH3-treated
    samples (Figure 6).
4.  Calculate the C to T misincorporation rate across every position of
    the reference transcriptome with tx_add_CtoTMR(). and subtract the C
    to T misincorporation between the control and the NaCNBH3 treatment
    samples for each growth temperature to obtain the rate C to T
    misincorporation rate difference (MRD_CtoT).
5.  Call putative ac4C sites if a position in any of the samples
    surpassed a threshold of MRD_CtoT \> 0.005.
6.  Plot a boxplot showing the calculated MRD_CtoT across temperatures
    (Figure 7). The box plot shows how the stoichiometry of ac4C sites
    increases as the temperature of growth increases, reaching the
    highest levels at 95 °C.
7.  Plot a sequence logo with the ggseqlogo wrapper tx_plot_ggseqlogo().
    (Figure 8)

``` r
path_bam2txDT <- system.file("bam2txDT.R", package = "txtools") # bam2txDT.R script path in local txtools installation
lapply(GSEtab_GSE135826$BAM, function(bamFile){
    system(paste("Rscript", path_bam2txDT, "-i ", bamFile, "-p TRUE -g", 
                 tk_geneAnno, "-f", tk_faGenome, "-d covNuc -r 300 -m 0 -n 10 -v TRUE"))
})

txDTL_Tk <- lapply(META_GSE135826$RDS, function(x){
                     tmpDT <- tx_load_rdsDT(x) %>% 
                         tx_add_misincRateNucSpec(refNuc = "C", misNuc = "T", minNucReads = 20) %>%
                         tx_add_misincCount()
                     return(tmpDT)}) %>% set_names(META_GSE135826$id)
saveRDS(txDTL_Tk, file = "results/useCase3_dat.rds")
txRES_tk <- txDTL_Tk[[1]][,1:6]
txRES_tk$'55deg' <- txDTL_Tk$NaCNBH3_55deg$MR_CtoT - txDTL_Tk$mock_55deg$MR_CtoT
txRES_tk$'65deg' <- txDTL_Tk$NaCNBH3_65deg$MR_CtoT - txDTL_Tk$mock_65deg$MR_CtoT
txRES_tk$'75deg' <- txDTL_Tk$NaCNBH3_75deg$MR_CtoT - txDTL_Tk$mock_75deg$MR_CtoT
txRES_tk$'85deg' <- txDTL_Tk$NaCNBH3_85deg$MR_CtoT - txDTL_Tk$mock_85deg$MR_CtoT
txRES_tk$'95deg' <- txDTL_Tk$NaCNBH3_95deg$MR_CtoT - txDTL_Tk$mock_95deg$MR_CtoT
### De novo calling ac4C sites
THR_ac4C <- 0.05 # Threshold for calling ac4C sites
txRES_tk$putative_ac4C <- rowSums(txRES_tk[,7:11] >= THR_ac4C, na.rm = TRUE) > 0
txRES_tk$putative_ac4C_fct <- factor(ifelse(txRES_tk$putative_ac4C,
                                                "putative_ac4C", "background"))
saveRDS(txRES_tk, file = "results/useCase3_res.rds")
```

## Plotting and results

``` r
gene_i <- "Tk_rRNA_23S"; pos_i <- 1462
gg_Tk1 <- tx_plot_nucFreq(txDTL_Tk$NaCNBH3_85deg, 
                          gene = gene_i, txRange = window_around(pos_i, 10),
                          bar_border = F, show_yLabels = F) +
    theme(legend.position = "none") + ggtitle("NaCNBH3-treatment")
gg_Tk2 <- tx_plot_nucFreq(txDTL_Tk$mock_85deg, 
                          gene = gene_i, txRange = window_around(pos_i, 10),
                          bar_border = F) + ggtitle("Control")

ggsave(filename = "figs/tk_plot_1.png",
       gridExtra::grid.arrange(gg_Tk1, gg_Tk2, ncol = 1), height = 6, width = 5)
```

![](figs/tk_plot_1.png)

**Figure 6.** Case study \# 3 - Dynamic RNA acetylation in T.
kodakarensis’ across a temperature gradient. txtools’ nucleotide
frequency plot.

``` r
longDT <- tidyr::pivot_longer(txRES_tk, cols = 7:11, names_to = "group", 
                              values_to = "MRD_CtoT")
gg_tk3 <- ggplot(longDT) + 
    geom_boxplot(aes(x = group, y = MRD_CtoT, colour = putative_ac4C_fct)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    theme(legend.position = "bottom") 
ggsave(filename = "figs/boxplot_tk.png", gg_tk3, height = 5, width = 4, bg = "white")
```

![](figs/boxplot_tk.png)

**Figure 7.** Case study \# 3 - Dynamic RNA acetylation in T.
kodakarensis’ across a temperature gradient. Boxplot of MRD_CtoT at
putative ac4C sites at increasing growth temperatures.

``` r
gg_tk4 <- tx_plot_ggseqlogo(txRES_tk, "putative_ac4C", upFlank = 5, doFlank = 5)
ggsave(filename = "figs/seqlogo_tk.png", gg_tk4, height = 3, width = 6, bg = "white")
```

![](figs/seqlogo_tk.png)

**Figure 8.** Case study \# 3 - Dynamic RNA acetylation in T.
kodakarensis’ across a temperature gradient. Sequence logo at detected
putative acetylated sites

# Supplementary analysis - Bridging paired-end reads

> The code to execute this analysis is provided in the Rscript
> `supp_PairedvsSingle.R`

To showcase the importance of bridging paired-end reads, which is
implemented as an option by txtools but not, to our knowledge, by other
tools, we analyzed m6A-seq data in yeast. In m6A-seq (Dominissini et al.
2012), RNA is first fragmented and then subjected to immunoprecipitation
using an anti-m6A antibody, resulting in selective capturing of
methylated fragments. These fragments are then sequenced from both ends.
The typical analytic pipeline consists of peak calling, based on
coverage, in the immunoprecipitated data. While m6A-seq is not
considered, inherently, a single-nucleotide resolution methodology, in
an idealized scenario (infinite coverage, no sources of noise), the
signal should peak precisely over the methylated adenosine in the DRAC
motif.

Given that the size of the immunoprecipitated fragments is around
100-150 nt and oftentimes only ~30 nt are sequenced from each end,
relying only on the sequenced ends results only in partial coverage of
the insert. txtools offers the possibility of computationally bridging
the reads, allowing to restore the full-length fragment. To assess to
what extent the lack of complete coverage resulted in inaccurate peak
calling, we analyzed existing m6A-seq data in yeast (Schwartz et al.
2013) using txtools, in two modes: either in paired-end mode (in which
paired ends are bridged) or in single-end mode (in which each read is
considered a separate entity). We observed a roughly 8% increase in the
number of peaks detected precisely centered on a DRAC motif (330 out of
1475 peaks were centered on the motif in single-end processing vs 429
out of 1405 peaks in paired-end processing) (Fig. S1).

<figure>
<img src="figs/readsProc_BP2.png" alt="Figure S1" />
<figcaption aria-hidden="true">Figure S1</figcaption>
</figure>

This was accompanied by a mild increase in the number of peaks called in
single-end mode, in comparison to paired-end mode (Fig. S2).

<figure>
<img src="figs/readsProc_BP1.png" alt="Figure S2" />
<figcaption aria-hidden="true">Figure S2</figcaption>
</figure>

Thus, without bridging paired ends reads, more peaks are called but of
poorer quality. An investigation into the source of these miscalled
peaks revealed that oftentimes a single peak - centered within a DRAC
motif - in paired-end mode, got split into multiple smaller peaks in
single-end mode, resulting in calling of multiple, erroneous sites
instead of a single correct one (Fig. S3). This analysis thus highlights
the added value of paired-end information for the accurate calling of
sites.

<figure>
<img src="figs/readsProc_covPeaks.png" alt="Figure S3" />
<figcaption aria-hidden="true">Figure S3</figcaption>
</figure>

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Carlile2014-ve" class="csl-entry">

Carlile, Thomas M, Maria F Rojas-Duran, Boris Zinshteyn, Hakyung Shin,
Kristen M Bartoli, and Wendy V Gilbert. 2014. “Pseudouridine Profiling
Reveals Regulated <span class="nocase">mRNA</span> Pseudouridylation in
Yeast and Human Cells.” *Nature* 515 (7525): 143–46.

</div>

<div id="ref-dominissini2012topology" class="csl-entry">

Dominissini, Dan, Sharon Moshitch-Moshkovitz, Schraga Schwartz, Mali
Salmon-Divon, Lior Ungar, Sivan Osenberg, Karen Cesarkas, et al. 2012.
“Topology of the Human and Mouse m6A RNA Methylomes Revealed by
m6A-Seq.” *Nature* 485 (7397): 201–6.

</div>

<div id="ref-Gentleman2019-ch" class="csl-entry">

Gentleman, Robert, V Carey, Wolfgang Huber, Florian Hahne, and
Maintainer Bioconductor Package Maintainer. 2019. “Package
‘Genefilter’.”

</div>

<div id="ref-Kortel2021-nn" class="csl-entry">

Körtel, Nadine, Cornelia Rücklé, You Zhou, Anke Busch, Peter Hoch-Kraft,
F X Reymond Sutandy, Jacob Haase, et al. 2021. “Deep and Accurate
Detection of m6A RNA Modifications Using
<span class="nocase">miCLIP2</span> and m6Aboost Machine Learning.”
*Nucleic Acids Res.* 49 (16): e92.

</div>

<div id="ref-Liao2019-oo" class="csl-entry">

Liao, Yang, Gordon K Smyth, and Wei Shi. 2019. “The R Package Rsubread
Is Easier, Faster, Cheaper and Better for Alignment and Quantification
of RNA Sequencing Reads.” *Nucleic Acids Res.* 47 (8): e47.

</div>

<div id="ref-Linder2015-yi" class="csl-entry">

Linder, Bastian, Anya V Grozhik, Anthony O Olarerin-George, Cem Meydan,
Christopher E Mason, and Samie R Jaffrey. 2015.
“Single-Nucleotide-Resolution Mapping of m6A and m6Am Throughout the
Transcriptome.” *Nat. Methods* 12 (8): 767–72.

</div>

<div id="ref-Sas-Chen2020-px" class="csl-entry">

Sas-Chen, Aldema, Justin M Thomas, Donna Matzov, Masato Taoka, Kellie D
Nance, Ronit Nir, Keri M Bryson, et al. 2020. “Dynamic RNA Acetylation
Revealed by Quantitative Cross-Evolutionary Mapping.” *Nature* 583
(7817): 638–43.

</div>

<div id="ref-schwartz2013high" class="csl-entry">

Schwartz, Schraga, Sudeep D Agarwala, Maxwell R Mumbach, Marko
Jovanovic, Philipp Mertins, Alexander Shishkin, Yuval Tabach, et al.
2013. “High-Resolution Mapping Reveals a Conserved, Widespread, Dynamic
mRNA Methylation Program in Yeast Meiosis.” *Cell* 155 (6): 1409–21.

</div>

<div id="ref-Taoka2016-hj" class="csl-entry">

Taoka, Masato, Yuko Nobe, Yuka Yamaki, Yoshio Yamauchi, Hideaki
Ishikawa, Nobuhiro Takahashi, Hiroshi Nakayama, and Toshiaki Isobe.
2016. “The Complete Chemical Structure of Saccharomyces Cerevisiae
<span class="nocase">rRNA</span>: Partial Pseudouridylation of U2345 in
25S <span class="nocase">rRNA</span> by
<span class="nocase">snoRNA</span> snR9.” *Nucleic Acids Res.* 44 (18):
8951–61.

</div>

<div id="ref-Thomas2018-dh" class="csl-entry">

Thomas, Justin M, Chloe A Briney, Kellie D Nance, Jeffrey E Lopez,
Abigail L Thorpe, Stephen D Fox, Marie-Line Bortolin-Cavaille, et al.
2018. “A Chemical Signature for Cytidine Acetylation in RNA.” *J. Am.
Chem. Soc.* 140 (40): 12667–70.

</div>

</div>

# Session Info

``` r
sessionInfo()
```
