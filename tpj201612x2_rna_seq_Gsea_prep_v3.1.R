#biomart tutorial
#http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.R
#tpj201612x2_rna_seq_Gsea_prep_v1.R

###See version 3 for why the first version had different genes due to possibly Biomart updating names
#I'm using v3 (based on v1) 
#v3.1 is the same as v3 but I'm going to add further annotation, download the annotation file, and remove extraneous code

rm(list=ls())
setwd("~/Desktop/projects/POST/data/CHORI/CAP/expression/theusch_2016/")




## ----biomaRt----------------------------------------------------------------------------------------------------------
library("biomaRt")
#listMarts()
## ----ensembl1---------------------------------------------------------------------------------------------------------
#ensembl=useMart("ensembl")


## ----ensembl3---------------------------------------------------------------------------------------------------------
#if you run this too quickly before the library has loaded you will get an extra error.
#you may need to try doing it a few times

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

## ----attributes-------------------------------------------------------------------------------------------------------
attributes = listAttributes(ensembl)
attributes

#Example of submitting a query
head(getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol','hgnc_id' ), mart = ensembl))


annot <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol' ), mart = ensembl)


####Written out on October 12, 2018 (do no write over)
###write.csv(annot, file = "ensembl_gene_id_hgnc_symbol_biomart_2018.10.12.output_v3.1.csv")

dim(annot)#64932



#Downloaded data from https://www.nature.com/articles/tpj201612#supplementary-information
#Need to convert characters to numeric

#Import and skip the first and last lines which were merged with some info about the data
#First line: "Table S1. Statin-responsiveness of 13,931 human and 73 EBV LCL genes."
#Last line : "The average delta was calculated by subtracting the library-size adjusted
#and variance stabilized estimates of the control from the corresponding statin libraries 
#and averaging across all LCLs. The q value is the FDR-adjusted p-value, and entries in 
#bold were considered significantly statin responsive (q<0.0001)."
getwd()
rna_xls <- as.data.frame(readxl::read_xls(path = "/Users/jemiller2/Desktop/projects/POST/data/CHORI/CAP/expression/theusch_2016/tpj201612x2.xls", 
                            skip = 1,n_max = 14004 , col_names = TRUE, col_types = "text"))

#Convert to numeric
colnames(rna_xls) <- c("Gene", "AvgDelta1", "p1", "q1")
rna_xls$AvgDelta <- as.numeric(rna_xls$AvgDelta1)
rna_xls$p <- as.numeric(rna_xls$p1)
rna_xls$q <- as.numeric(rna_xls$q1)
rna_seq <- rna_xls[,c(1,5:7)]

#Just taking a look at the data to check 
head(rna_seq)
tail(rna_seq)
dim(rna_seq)
str(rna_seq)
subset(rna_xls, AvgDelta > 1)
subset(rna_seq, AvgDelta > 1)

#Merge RNA-seq data with HGNC symbols so that it can be used with custom gmt file in GSEA
rna_an <- merge(annot, rna_seq, by.x = "ensembl_gene_id", "Gene")
head(rna_an)
dim(rna_an)#13552
head(rna_an[order(rna_an$hgnc_symbol, decreasing = FALSE),])

#Beth et al used a cutoff of q < 0.0001
rna_sig <- subset(rna_an, q < 0.0001)
dim(rna_sig)#5170
head(rna_sig)
rna_sig <- rna_sig[,c(2,3)]

#There were many genes with no hgnc symbol so I will replace with NA then remove
head(rna_sig[order(rna_sig$hgnc_symbol, decreasing = FALSE),])
rna_sig$GeneName <- as.character(rna_sig$hgnc_symbol)
head(rna_sig)
head(rna_sig[order(rna_sig$hgnc_symbol, decreasing = FALSE),])


#Add NA
rna_sig$GeneName[rna_sig$GeneName==""] <- "NA"
head(rna_sig[order(rna_sig$hgnc_symbol, decreasing = FALSE),])
#Remove NA rows
rna_sig_noNA <- subset(rna_sig, GeneName != "NA")
head(rna_sig_noNA[order(rna_sig_noNA$hgnc_symbol, decreasing = FALSE),])
dim(rna_sig_noNA)#5109
#Order
rna_sig_noNA_order <- rna_sig_noNA[order(rna_sig_noNA$AvgDelta, decreasing = TRUE),]



max(rna_sig_noNA_order$AvgDelta)
subset(rna_sig_noNA_order, AvgDelta > .8)
hist(rna_sig_noNA_order$AvgDelta)
dim(rna_sig_noNA_order)
dim(na.omit(rna_sig_noNA_order))

#Select gene and delta 
rna_sig_noNA_order_fi <- rna_sig_noNA_order[,c(3,2)]
head(rna_sig_noNA_order_fi)
dim(rna_sig_noNA_order_fi)

dim(rna_sig_noNA_order_fi)
#FInal version No need to rewrite, written on 2018.10.12
##write.table(rna_sig_noNA_order_fi, file = "tpj201612x2_rna_seq_order_sig_v3.1_2018.10.12_gsea_hgnc_input.txt", quote = FALSE,sep = "\t", row.names = FALSE)




