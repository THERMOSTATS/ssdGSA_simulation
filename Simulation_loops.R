library(tidyverse)
library(limma)
library(edgeR)
library(gridExtra)
library(org.Hs.eg.db)
library(GSVA)
library(ggpubr)
library(GSVAdata)
data(c2BroadSets)
library(ssdGSA)
library(nlme)
library(lme4)
library(emmeans)
library(hgu133plus2.db)
library(genefilter)
library(Hmisc)
library(annotate)
library(biomaRt)
library("GEOquery")




#below public data are gene expression matrix for `Data` object
#Brodalumab data download and analysis
d <- GEOquery::getGEO("GSE117468")[[1]]
d@annotation <- "hgu133plus2.db" # revise annotation name GPL570 to acutal
df<- nsFilter(d, require.entrez=TRUE, remove.dupEntrez=TRUE,
              var.func=IQR, var.filter=TRUE, var.cutoff=0.1, filterByQuantile=TRUE,
              feature.exclude="^AFFX") #filter the probes with small variance and remove probes with duplicated entrez id or no entrez id

data <- df$eset@assayData$exprs %>% as.data.frame %>% tibble::rownames_to_column(var = "ID")
meta <- df$eset@phenoData@data %>%
  as_tibble(rownames = "SAMPLE_ID") %>%
  dplyr::select(SAMPLE_ID,
                PATIENT_ID = `patientid:ch1`,
                Treatment = `treatment:ch1`,
                Visit = `visit:ch1`,
                Tissue = `tissue:ch1`) %>%
  mutate(Tissue = recode(Tissue, `lesional skin` = "LS",
                         `non-lesional skin` = "NL"),
         Treatment = recode(Treatment, `Brodalumab 140 mg Q2W` = "Brodalumab_140mg",
                            `Brodalumab 210 mg Q2W` = "Brodalumab_210mg"))


### Choose those are LS + BASELINE
#EXTRACT probe annotation information from GPL570
anno <- read_delim("GPL570-55999.txt",delim = "\t", skip = 16) %>%
  dplyr::select(ID, SYMBOL = `Gene Symbol`, ENTREZID = ENTREZ_GENE_ID) %>%
  mutate(SYMBOL = word(SYMBOL),
         ENTREZID = word(ENTREZID))

dat <- data %>%
  left_join(anno %>% dplyr::select(ID, ENTREZID)) %>%
  dplyr::select(-ID) %>%
  distinct(ENTREZID, .keep_all = TRUE) %>%
  filter(!is.na(ENTREZID)) %>%
  dplyr::select(ENTREZID,meta$SAMPLE_ID) %>%
  column_to_rownames("ENTREZID") %>%
  as.matrix()

Sample_ID_used <- meta %>%
  filter(Visit == "BL", Tissue == "LS") %>%
  dplyr::select(SAMPLE_ID) %>%
  pull() %>%
  as.vector()

dat <- dat[,colnames(dat) %in% Sample_ID_used]


# Withdraw KEGG pathway information
gene.set <- c2BroadSets[grep("KEGG_", names(c2BroadSets))] %>% geneIds
gene.set <- gene.set[lengths(gene.set)>20]




##################################################################
### Next, we need to add treatment effects
####################################################################

# Randomly select ten pathways
nGenesets = 10
B = 1200
Pval.hypo.ssdGSA <- matrix(nrow = B, ncol = nGenesets)
Pval.hypo.ssGSA <- matrix(nrow = B, ncol = nGenesets)
for (b in 1:B){
  skip_to_next <- FALSE
  tryCatch({
    GeneSetID <- sample(1:length(gene.set), size = nGenesets, replace = FALSE)
    gene.set2 <- gene.set[GeneSetID]

    # Remove the common elements in all the lists
    gene.set2 <- enframe(gene.set2) %>%
      unnest(value) %>%
      group_by(value) %>%
      filter(n_distinct(name) == 1) %>%
      with(., split(value, name))

    # According to genes in the gene sets, we form the gene expression matrix
    genes_in_gene_sets <- unlist(gene.set2) %>% as.vector

    dat.mat <- dat %>%
      as.data.frame() %>%
      rownames_to_column(var = "gene") %>%
      filter(gene %in% genes_in_gene_sets) %>%
      column_to_rownames(var = "gene")

    genes_in_datamat <- rownames(dat.mat) %>% as.vector

    ## Deal with genes in gene set that are not in the gene expression matrix
    gene.set3 <- enframe(gene.set2) %>%
      unnest(value) %>%
      filter(value %in% genes_in_datamat) %>%
      with(., split(value, name))
    # sum(lengths(gene.set3))
    GenB <- function(mean, sd = 0.05){
      x <- rnorm(n = length(mean) * ncol(dat.mat),
                 mean = rep(mean, each = ncol(dat.mat)), sd = sd)
      beta <- t(matrix(x, nrow = ncol(dat.mat)))
      return(beta)
    }

    bMat <- c()
    ESmeanVec <- c()
    for (k in 1:nGenesets){
      #############################################################################
      ## Simulate disease effect sizes on each gene, e.g. 80% of genes in pathways are 
      ## positively correlated with the disease   
    
      es.pov <- runif(round(lengths(gene.set3)[k]*0.8), min = 0, max = 2)
      es.neg <- runif((lengths(gene.set3)[k]-length(es.pov)), min =-2, max = 0)
      esMean <- c(es.pov, es.neg)
      ESmeanVec <- c(ESmeanVec, esMean)
      
      #############################################################################
      ## According to disease effect sizes, simulate corresponding treatment effect means
      mean.p <- runif(round(lengths(gene.set3)[k]*0.8), min = -1, max = 0)
      mean.n <- runif((lengths(gene.set3)[k]-length(mean.p)), min = 0, max = 1)
      mean <- c(mean.p, mean.n)
      
      #############################################################################
      ## Simulate treatment effect for each sample
      bMat.s <- GenB(mean = mean, sd = 0.05)
      bMat <- rbind(bMat, bMat.s)
    }

    ############################################
    # Add treatment effect #####################
    ## Align the order in gene sets
    gene.list <- unlist(gene.set3) %>% as.vector()
    dat.mat2 <- dat.mat[gene.list,]
    TreatGroup <- bMat + dat.mat2


    ############################################
    # Format direction matrix
    DirectionMat <- cbind(gene.list,ESmeanVec)
    colnames(DirectionMat) <- c("Gene", "ES")
    data.dis <- DirectionMat %>%
      as_tibble() %>%
      transform(ES = as.numeric(ES))

    #################################################
    ## Run Gene Set Analysis

    ssgsea_avgES_ssdGSA <- ssdGSA(Data = as.matrix(TreatGroup),
                                  Gene_sets = gene.set3,
                                  Direction_matrix = data.dis,
                                  GSA_weight = "group_weighted",
                                  GSA_weighted_by = "avg.ES",
                                  GSA_method = "ssgsea",
                                  min.sz = 1,
                                  max.sz = 2000,
                                  mx.diff = TRUE)

    ssgsea_avgES_ssdGSA_beforeTreat <- ssdGSA(Data = as.matrix(dat.mat2),
                                              Gene_sets = gene.set3,
                                              Direction_matrix = data.dis,
                                              GSA_weight = "group_weighted",
                                              GSA_weighted_by = "avg.ES",
                                              GSA_method = "ssgsea",
                                              min.sz = 1,
                                              max.sz = 2000,
                                              mx.diff = TRUE)

    ssgsea_avgES_ssGSA <- ssdGSA(Data = as.matrix(TreatGroup),
                                 Gene_sets = gene.set3,
                                 Direction_matrix = NULL,
                                 GSA_weight = "group_weighted",
                                 GSA_weighted_by = "avg.ES",
                                 GSA_method = "ssgsea",
                                 min.sz = 1,
                                 max.sz = 2000,
                                 mx.diff = TRUE)

    ssgsea_avgES_ssGSA_beforeTreat <- ssdGSA(Data = as.matrix(dat.mat2),
                                             Gene_sets = gene.set3,
                                             Direction_matrix = NULL,
                                             GSA_weight = "group_weighted",
                                             GSA_weighted_by = "avg.ES",
                                             GSA_method = "ssgsea",
                                             min.sz = 1,
                                             max.sz = 2000,
                                             mx.diff = TRUE)

    ## Compare scores from ssdGSA before and after treatment

    for (k in 1:nGenesets){
      Pval.hypo.ssdGSA[b,k] <- t.test(ssgsea_avgES_ssdGSA[k,], ssgsea_avgES_ssdGSA_beforeTreat[k,])$p.value
      Pval.hypo.ssGSA[b,k] <- t.test(ssgsea_avgES_ssGSA[k,], ssgsea_avgES_ssGSA_beforeTreat[k,])$p.value
    }

  }, error = function(e){skip_to_next <- TRUE})
  if(skip_to_next){next}

}

## Compare scores from ssGSA before and after treatment

Pevaluation.ssdGSA.all <- Pevaluation.ssGSA.all <- c()
for (k in 1:nGenesets){
  Pevaluation.ssdGSA.all <- c(Pevaluation.ssdGSA.all, sum(na.omit(Pval.hypo.ssdGSA[,k])[1:1000]>0.05))
  Pevaluation.ssGSA.all <- c(Pevaluation.ssGSA.all, sum(na.omit(Pval.hypo.ssGSA[,k])[1:1000]>0.05))
}
sum(Pevaluation.ssdGSA.all)/10000
sum(Pevaluation.ssGSA.all)/10000
