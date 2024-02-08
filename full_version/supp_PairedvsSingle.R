# Setup ####
library("txtools")
## Supplementary packages ####
library("GEOfastq") # Obtain FASTQ files and metadata from SRA/ENA databases
library("Rsubread") # Align FASTQ files to reference genome
library("magrittr") # Pipe operator
library("data.table") # data.table class full functionality
library("gridExtra") # Arrange plots in one panel
library("genefilter") # Fast t-tests
library("parallel") # Multi-thread functionalization
library("tidyverse") # Collection of data-science packages

## Values ####
NCORES <- 10 # Number of cores to be used in all multi-core processes

## Functions ####
get_GSEtable <- function(GSE){
    require(GEOfastq)
    require(magrittr)
    crawl_gse(GSE) %>% extract_gsms() %>% crawl_gsms()
}

labelRegions <- function(x, prefix = "m6AseqRegion", minRegionSize = 20){
    tmpLabel <- character(length(x))
    factor_num <- 1
    for(i in seq(length(x))){
        if(!is.na(x[i])){
            if(x[i] == FALSE){
                tmpLabel[i] <- NA
            }else if(x[i] == TRUE){
                tmpLabel[i] <- paste(prefix, factor_num, sep = "_")
                if(x[i + 1] == FALSE | is.na(x[i + 1])){
                    factor_num <- factor_num + 1
                }
            }
        }else{
            tmpLabel[i] <- NA
        }
    }
    numLvls <- length(unique(na.omit(tmpLabel)))
    tmpFct <- factor(tmpLabel, levels = paste(prefix, seq(numLvls), sep = "_"))
    regCounts <- table(tmpFct)
    shortRegions <- names(regCounts)[regCounts < minRegionSize]
    tmpFct[which(tmpFct %in% shortRegions)] <- NA
    tmpFct <- fct_drop(tmpFct)
    levels(tmpFct) <- paste(prefix, seq(length(unique(na.omit(as.character(tmpFct))))), sep = "_")
    tmpFct
}

USdistToMotif <- function(sites, motifs, genes){
    tmpS <- split(sites, genes)
    tmpT <- split(motifs, genes)
    tmpV <- lapply(seq(length(tmpS)), function(i){
        tmpU <- rep(NA, length(tmpS[[i]])) %>% as.integer()
        if(sum(!is.na(tmpS[[i]])) > 0){
            dist <- sapply(which(tmpS[[i]]), function(a){
                tmpDiff <- a - which(tmpT[[i]])
                selDist <- tmpDiff[tmpDiff >= 0]
                if(length(selDist) > 0){
                    min(abs(selDist))
                }else{
                    NA
                }
            })
            tmpU[which(tmpS[[i]])] <- dist
            return(tmpU)
        }else{
            return(tmpU)
        }
    }) 
    unlist(tmpV)
}

DSdistToMotif <- function(sites, motifs, genes){
    tmpS <- split(sites, genes)
    tmpT <- split(motifs, genes)
    tmpV <- lapply(seq(length(tmpS)), function(i){
        tmpU <- rep(NA, length(tmpS[[i]])) %>% as.integer()
        if(sum(!is.na(tmpS[[i]])) > 0){
            dist <- sapply(which(tmpS[[i]]), function(a){
                tmpDiff <- a - which(tmpT[[i]])
                selDist <- tmpDiff[tmpDiff <= 0]
                if(length(selDist) > 0){
                    min(abs(selDist))
                }else{
                    NA
                }
            })
            tmpU[which(tmpS[[i]])] <- dist
            return(tmpU)
        }else{
            return(tmpU)
        }
    }) 
    unlist(tmpV)
}

# Rolling/moving mean
roll_mean <- function(x, n = 10){as.numeric(stats::filter(x, rep(1 / n, n), sides = 2))}

# Winscore between txDTs
calculate_winScore <- function(txDT_IP, txDT_IN, winSize = 10, min_IP_cov = 20){
    tmp_IP <- split(txDT_IP$cov, txDT_IP$gene) %>% lapply(function(x){
        tmpA <- roll_mean(x, n = winSize) 
        tmpA[tmpA < min_IP_cov] <- NA
        tmpA / median(x, na.rm = T)
    }) %>% do.call(what = "c")
    tmp_Input <- split(txDT_IN$cov, txDT_IN$gene) %>% lapply(function(x){
        roll_mean(x, n = winSize) / median(x, na.rm = T)
    }) %>% do.call(what = "c")
    out <- log2(tmp_IP / tmp_Input)
    out[is.infinite(out)] <- NA
    out
}

# Minimum of pairs of two vectors
minByPairs <- Vectorize(function(a, b){min(a, b)}, vectorize.args = c("a", "b"))

# Aggregate coverage, calculate winScore, de novo detect m6A sites WT.IP coverage
process_m6A <- function(txDTL, META, winScore_thr = 2, nCores){
    # Aggregate data from replicates
    sumDTL <- lapply(levels(META$group), function(iGroup){
        tmp <- txDTL[META$group == iGroup]
        baseDT <- tmp[[1]][, names(tmp[[1]]) %in% txCoreCols_refSeq, with = FALSE]
        colsToAdd <- c("start_5p", "end_3p", "cov")
        tmpM <- lapply(colsToAdd, function(iCol){
            Reduce("+", lapply(tmp, function(x) x[[iCol]]))
        }) %>% do.call(what = cbind) %>% data.table::data.table() %>% magrittr::set_colnames(colsToAdd)
        cbind(baseDT, tmpM)
    }) %>% set_names(levels(META$group))
    
    # WinScores
    resDT <- sumDTL[[1]][, ..txCoreCols_refSeq]
    resDT$winScore_WT <- calculate_winScore(sumDTL$`WT:IP`, sumDTL$`WT:input`)
    resDT$winScore_DL <- calculate_winScore(sumDTL$`IME4del:IP`, sumDTL$`IME4del:input`)
    resDT$WT_IP_cov <- sumDTL$`WT:IP`$cov
    resDT$winScore_diff <- resDT$winScore_WT - resDT$winScore_DL
    resDT$putRegion <- resDT$winScore_diff > winScore_thr &
        resDT$winScore_DL >= 0
    resDT$putRegion <- labelRegions(resDT$putRegion)
    
    # putative site - coverage
    resDT$putSite_cov <- FALSE
    for(iReg in levels(resDT$putRegion)){
        iTxRange <- which(resDT$putRegion == iReg)
        resDT$putSite_cov[iTxRange[which.max(resDT$WT_IP_cov[iTxRange])]] <- TRUE
    }
    
    # DRAC site distance
    resDT <- tx_add_motifPresence(resDT, motif = "DRAC", nucPositions = 3, nCores = nCores)
    resDT$USdist_DRAC <- USdistToMotif(sites = resDT$putSite_cov, motifs = resDT$DRAC_motif_3, genes = resDT$gene)
    resDT$DSdist_DRAC <- DSdistToMotif(sites = resDT$putSite_cov, motifs = resDT$DRAC_motif_3, genes = resDT$gene)
    resDT$ABdist_DRAC <- minByPairs(resDT$USdist_DRAC, resDT$DSdist_DRAC)
    
    resDT
}

# Download data from Zenodo ####
download.file(url = "https://zenodo.org/records/10612727/files/supp_geneAnnot.bed", destfile = "omicRefs/supp_geneAnnot.bed")
download.file(url = "https://zenodo.org/records/10612727/files/supp_genome.fa", destfile = "omicRefs/supp_genome.bed")
download.file(url = "https://zenodo.org/records/10612727/files/supp_resPE.rds", destfile = "results/supp_resPE.rds")
download.file(url = "https://zenodo.org/records/10612727/files/supp_resSE.rds", destfile = "results/supp_resSE.rds")

# Omic references 
sk1_faGenome <- "omicRefs/supp_genome.fa"
sk1_geneAnno <- "omicRefs/supp_geneAnnot.bed"

# Schwartz et al. 2013 m6A-seq yeast ####
# Experiment design
GSEtab_GSE51583 <- get_GSEtable("GSE51583")
GSEtab_GSE51583 <- GSEtab_GSE51583[str_detect(GSEtab_GSE51583$source_name_ch1, pattern = "Replicate"),]
GSEtab_GSE51583$rep <- str_extract(GSEtab_GSE51583$source_name_ch1, pattern = "Replicate [[:digit:]]")
GSEtab_GSE51583$bioTreat <- str_extract(GSEtab_GSE51583$source_name_ch1, pattern = "(WT|IME4 deletion)")
GSEtab_GSE51583 <- GSEtab_GSE51583[GSEtab_GSE51583$source_name_ch1 %>% order(),]

META_GSE51583 <- data.table(run = GSEtab_GSE51583$run,
                            title = GSEtab_GSE51583$title,
                            rep = GSEtab_GSE51583$rep,
                            bioTreat = GSEtab_GSE51583$bioTreat,
                            FASTQ = paste0(file.path(FASTQDir, GSEtab_GSE51583$run), "_1.fastq.gz"),
                            BAM = paste0(file.path(BAMDir, GSEtab_GSE51583$run), ".bam"),
                            RDS = paste0(file.path(BAMDir, GSEtab_GSE51583$run), ".txDT.rds"))
META_GSE51583$libTreat <- str_extract(META_GSE51583$title, pattern = "(IP|input)")
META_GSE51583$rep <- META_GSE51583$rep |> str_replace("Replicate", "rep") |> str_replace(" ", "")
META_GSE51583$bioTreat <- META_GSE51583$bioTreat |> str_replace(" deletion", "del")
META_GSE51583$id <- paste(META_GSE51583$bioTreat, META_GSE51583$libTreat, META_GSE51583$rep, sep = "_")
META_GSE51583$group <- factor(paste(META_GSE51583$bioTreat, META_GSE51583$libTreat, sep = ":"))
META_GSE51583 <- META_GSE51583[-9, ] # low quality of WT_IP_rep2

# FASTQ downloading 
FASTQDir <- "data/fastq/GSE51583"; dir.create(FASTQDir, recursive = TRUE)
BAMDir <- "data/bam/GSE51583"; dir.create(BAMDir, recursive = TRUE)
get_fastqs(GSEtab_GSE51583, FASTQDir) # download FASTQ data

# Mapping reads to genome
buildindex(basename = "sk1_RsubRef", reference = sk1_faGenome)
r1Files <- "SRR1019460_1.fastq.gz"
r1Files <- grep(list.files(FASTQDir), pattern = "_1.fastq.gz", value = TRUE)
sk1_alignReport <- lapply(r1Files, function(FASTQ){
    align(index = "sk1_RsubRef",
          readfile1 = file.path(FASTQDir, FASTQ),
          readfile2 = file.path(FASTQDir, gsub(FASTQ, pattern = "_1.fastq.gz", replacement = "_2.fastq.gz")),
          output_file = file.path(BAMDir, gsub(x = FASTQ, pattern = "_1.fastq.gz", replacement = ".bam")),
          nthreads = NCORES, minFragLength = 23, sortReadsByCoordinates = TRUE)
})
file.remove(list.files(pattern = "sk1_RsubRef.*")) # Remove index files

# Process BAM files to txDT using the bam2txDT.R script included in txtools
path_bam2txDT <- system.file("bam2txDT.R", package = "txtools")

## Paired-end txDTs
walk(META_GSE51583$BAM, function(bamFile){
    system(paste("nohup Rscript", path_bam2txDT, "-i", bamFile, "-p TRUE -g", 
                 sk1_geneAnno, "-f", sk1_faGenome,
                 "-d covNuc -r 1000 -S FALSE -m 200 -n 12 -v TRUE >", 
                 str_replace(bamFile, pattern = ".bam$", ".log"), "&"))
})

## Single-end txDTs
walk(META_GSE51583$BAM, function(bamFile){
    system(paste("nohup Rscript", path_bam2txDT, "-i", bamFile, "-o", 
                 str_replace(bamFile, pattern = ".bam$", replacement = ".SE.txDT.rds"),
                 "-p FALSE -g", sk1_geneAnno, "-f", sk1_faGenome, 
                 "-d covNuc -r 1000 -S TRUE -m 200 -n 12 -v TRUE >", 
                 str_replace(bamFile, pattern = ".bam$", ".SE.log"), "&"))
})

# Analysis ####
# Load data processes as paired-end reads (PE)
txDTL_PE <- lapply(META_GSE51583$RDS, function(rds){
    readRDS(rds)
}) %>% set_names(META_GSE51583$id)
txDTL_PE <- tx_unifyTxDTL(txDTL_PE, nCores = NCORES)

# Load data processes as single-end reads (SE)
SE_RDS <- str_replace(META_GSE51583$RDS, pattern = ".txDT.rds", 
                      replacement = ".SE.txDT.rds")
txDTL_SE <- lapply(SE_RDS , function(rds){
    readRDS(rds)
}) %>% set_names(META_GSE51583$id)
txDTL_SE[[1]] <- txDTL_SE[[1]][gene %in% txDTL_PE[[1]]$gene,]
txDTL_SE <- tx_unifyTxDTL(txDTL_SE, nCores = NCORES)

# m6A-sites detection
# Briefly process_m6A() aggregates coverage data of all replicates for each 
# group (WT_Input, WT_IP, IME4KO_Input, IME4KO_IP). Then WinScores are calculated
# according to Dominissini et al., 2012, using the calculate_winScore() function. 
# Interest regions are selected as those where the difference of Winscore(WT) - Winscore(IME4_deletion), 
# is greater than 2, and the Winscore(IME4_deletion) is equal or greater than 0. 
# Across each interest region the position with the greatest coverage in the 
# WT_IP sample is selected as the putative m6A site.

RES_PE <- process_m6A(txDTL = txDTL_PE, META_GSE51583, nCores = NCORES) # Paired-end processed data
RES_SE <- process_m6A(txDTL = txDTL_SE, META_GSE51583, nCores = NCORES) # Single-end processed data

# # Uncomment for loading results (downloaded at lines 161-162)
# RES_PE <- readRDS("results/supp_resPE.rds")
# RES_SE <- readRDS("results/supp_resSE.rds")

# Results ####
tmpPE <- select(RES_PE, c("gene", "ABdist_DRAC", "WT_IP_cov")) %>% na.omit()
tmpPE$readType <- "PE"
tmpSE <- select(RES_SE, c("gene", "ABdist_DRAC", "WT_IP_cov")) %>% na.omit()
tmpSE$readType <- "SE"
tmpBO <- rbind(tmpPE, tmpSE)
tmpDF <- as.data.frame.table(table(tmpBO$readType)) %>% set_colnames(c("readType", "freq"))

# Plots ####
# Number of putative m6A sites detected at exactly a DRAC motif.
gg_BP1 <- ggplot(tmpDF, aes(x = readType, y = freq)) + 
    geom_bar(stat = "identity", width = 0.8, position = "dodge") + theme_bw() +
    theme(legend.position = "bottom") + 
    labs(x = "Read processing type", y = "Putative m6A sites detected", fill = "read-type\n")

# Proportion of detected m6A sites that fall in a DRAC motif with binomial confidence interval
summAtMotif <- lapply(c("SE", "PE"), function(rT){
    data.table(readType = rT, 
               succ = sum(tmpBO[readType == rT, ]$ABdist_DRAC == 0),
               n = nrow(tmpBO[readType == rT, ])) %>% 
        mutate(prop = succ/n, 
               ci_min = prop.test(succ, n, p = 0.95, correct = FALSE)$conf.int[1],
               ci_max = prop.test(succ, n, p = 0.95, correct = FALSE)$conf.int[2])
}) %>% do.call(what = rbind)
gg_BP2 <- ggplot(summAtMotif, aes(x = readType, fill = readType)) +
    geom_bar(aes(y = prop), stat = "identity") +
    geom_errorbar(aes(ymin = ci_min, ymax = ci_max), width=0.2, alpha=0.9) + 
    theme_bw() + labs(y = "Proportion of detected m6A sites \nat DRAC sites",
                      x = "Reads processing", fill = "read-processing\n") +
    theme(legend.position = "none")

# Individual peaks PE vs SE
# Selected from inspection of sites that fall in a DRAC motif in PE and not present in SE (336 out of 429)
# which(RES_PE$ABdist_DRAC == 0 & RES_SE$putSite_cov == FALSE) <- this would show similar results
SITES <- c(306109, 5129446, 5931841, 1291758) 

boF <- 200 # Number of bp to both sides of center
ggl_covPeaks <- lapply(SITES, function(site_i){
    # site_i <- which(RES_PE$putSite_cov == TRUE & RES_SE$putSite_cov == FALSE & RES_PE$ABdist_DRAC < 2) %>%
    #     sample(1)
    gene_i <- RES_PE[site_i,]$gene
    txco_i <- RES_PE[site_i,]$txcoor
    xcoor_1 <- paste(txco_i, RES_PE[gene == gene_i & txcoor == txco_i]$refSeq, sep = "-")
    # xcoor_1 <- paste(RES_PE[gene == gene_i & putSite_cov == TRUE]$txcoor, 
    #                  RES_PE[gene == gene_i & putSite_cov == TRUE]$refSeq, sep = "-")
    tmpGG1 <- tx_plot_numeric(RES_PE, colVars = c("WT_IP_cov"), gene = gene_i, 
                              txRange = window_around(txco_i, boF), plot_type = "barplot") + 
        ggplot2::scale_fill_manual(values = c("gray", "black")) + 
        ggplot2::geom_vline(xintercept = xcoor_1, color = "red") +
        labs(title = paste(gene_i, txco_i, sep = ":"), y = "Coverage", x = NULL, subtitle = "Paired-End processing") +
        theme(legend.position = "none")
    xcoor_2 <- paste(RES_SE[gene == gene_i & putSite_cov == TRUE]$txcoor, 
                     RES_SE[gene == gene_i & putSite_cov == TRUE]$refSeq, sep = "-")
    tmpGG2 <- tx_plot_numeric(RES_SE, colVars = c("WT_IP_cov"), gene = gene_i, 
                              txRange = window_around(txco_i, boF), plot_type = "barplot") + 
        ggplot2::scale_fill_manual(values = c("gray", "black")) +
        ggplot2::geom_vline(xintercept = xcoor_2, color = "blue") + 
        labs(title = NULL, y = "Coverage", subtitle = "Single-End processing") + theme(legend.position = "none")
    # grid.arrange(grobs = list(tmpGG1, tmpGG2), ncol = 1)
    list(tmpGG1, tmpGG2)
})
gg_covPeaks <- grid.arrange(grobs = unlist(ggl_covPeaks, recursive = FALSE), 
                            layout_matrix = matrix(seq(length(ggl_covPeaks)*2), ncol = 2))

# saving plots
ggsave("figs/readsProc_BP1.png", gg_BP1, width = 3, height = 4)
ggsave("figs/readsProc_BP2.png", gg_BP2, width = 3, height = 4)
ggsave("figs/readsProc_covPeaks.png", plot = gg_covPeaks, width = 8, height = 10)
