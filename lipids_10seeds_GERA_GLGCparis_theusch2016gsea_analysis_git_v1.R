#TG_GLGCparis_theusc2016gsea_analysis_v3.R
#v1based on File =  TG_GLGCparis_theusc2016gsea_analysis_v1.R
#v2 lipids_10seeds... I'm add the gera part
#from tc_analysis_v1
#2018 analysis


rm(list = ls())


conflicts()
library(plyr)
library(tidyverse)

#######################
#Count SNPs
#GERA p-value threshold of 5x10^-8
#sqrtHDL <- 10,916/11,196,891 loci
#logTG  <- 12,202/11,196,891 loci
#LDL <-9,078/11,196,885 loci
#TC <- 13,552/11,196,884 loci

mean(11196891, 11196891, 11196885, 11196884)
#Mean = 11,196,891


#GLGC p-value threshold of 5x10^-8
#HDL =         3,518/2,447,441 loci
#TG =           3,249/2,439,432 loci
#LDL =         3,077 /2,437,751 loci
#TC =           4,169/2,446,981 loci


mean(2447441, 2439432, 2437751, 2446981)
#Mean = 2,447,441




setwd("~/Desktop/projects/POST/lipids_GLGC_GERA/output/paris_results/biofilter_2.4.2/glgc/version4/")

tg.glgc.1 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed1.paris-summary_formatted", header = TRUE)
tg.glgc.2 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed2.paris-summary_formatted", header = TRUE)
tg.glgc.3 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed3.paris-summary_formatted", header = TRUE)
tg.glgc.4 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed4.paris-summary_formatted", header = TRUE)
tg.glgc.5 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed5.paris-summary_formatted", header = TRUE)
tg.glgc.6 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed6.paris-summary_formatted", header = TRUE)
tg.glgc.7 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed7.paris-summary_formatted", header = TRUE)
tg.glgc.8 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed8.paris-summary_formatted", header = TRUE)
tg.glgc.9 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed9.paris-summary_formatted", header = TRUE)
tg.glgc.10 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed10.paris-summary_formatted", header = TRUE)



process_paris <- function(x, y) {
  
  colnames(x) <- c("pathway" ,"description" ,"genes" ,"features", "simple", "simple.sig.", "complex" ,"complex.sig.pval")
  z <- cbind(y,x)
  z$p_vals2 <- gsub( "<", "", as.character(z$complex.sig.pval) )
  z$p_vals <- as.numeric(gsub( ">=", "", as.character(z$p_vals2) ))
  z
  return(z)
  #paste(colnames(x), y, sep = ".")
}
#sed 's/<//g' | sed 's/>=//g'
tg.glgc.1_pp <- process_paris(tg.glgc.1, "tg.glgc.1")
head(tg.glgc.1_pp)
tg.glgc.2_pp <- process_paris(tg.glgc.2, "tg.glgc.2")
tg.glgc.3_pp <- process_paris(tg.glgc.3, "tg.glgc.3")
tg.glgc.4_pp <- process_paris(tg.glgc.4, "tg.glgc.4")
tg.glgc.5_pp <- process_paris(tg.glgc.5, "tg.glgc.5")
tg.glgc.6_pp <- process_paris(tg.glgc.6, "tg.glgc.6")
tg.glgc.7_pp <- process_paris(tg.glgc.7, "tg.glgc.7")
tg.glgc.8_pp <- process_paris(tg.glgc.8, "tg.glgc.8")
tg.glgc.9_pp <- process_paris(tg.glgc.9, "tg.glgc.9")
tg.glgc.10_pp <- process_paris(tg.glgc.10, "tg.glgc.10")


tg.glgc.1to10_pp <- rbind.fill(tg.glgc.1_pp,tg.glgc.2_pp,tg.glgc.3_pp,tg.glgc.4_pp,tg.glgc.5_pp,tg.glgc.6_pp,tg.glgc.7_pp,tg.glgc.8_pp,tg.glgc.9_pp,tg.glgc.10_pp)
dim(tg.glgc.1to10_pp)
head(tg.glgc.1to10_pp)

tg.glgc.1to10_pp %>% filter(p_vals < 0.01)
tg.glgc.1.2.3v2_ol <- as.data.frame(tg.glgc.1to10_pp %>% filter(p_vals < 0.01) %>% group_by(pathway) %>% summarize(n = n() ))
tg.glgc.1.2.3v2_ol
dim(tg.glgc.1.2.3v2_ol)#Total significant pathways (31)
##############################################################################################################

#inesrt GSEA results

setwd("/Users/jemiller2/Desktop/projects/POST/data/POST/theusch_2016/gsea_tpj201612x2_analysis/output/tpj201612x2_rna_seq_order_sig_gsea_hgnc_2018lokiKEGG_plot50.GseaPreranked.1520864978540/")



###############Functions
gsea_pos_ol <- function(x) {
  #positive expression overlap
  gsea_pos_fx <- read.csv( "gsea_report_for_na_pos_1520864978540.csv", header =TRUE, sep = ",")
  gsea_pos_fx <- gsea_pos_fx[,c(1,8)]
  gsea_pos_sig25_fx <- subset(gsea_pos_fx, FDR.q.val < 0.25)
  x$NAME <- toupper(x$pathway)
  gsea_pos_sig25_fx$NAME <- gsub(" ", "_", gsea_pos_sig25_fx$NAME, fixed=TRUE)
  gsea_tg.glgc_paris_kegg_fx<- merge(gsea_pos_sig25_fx, x,  by = "NAME")
  #negative expression overlap
  return(subset(gsea_tg.glgc_paris_kegg_fx, n > 7))
}

gsea_pos_ol(tg.glgc.1.2.3v2_ol)#Overlap with up-regulated genes (4)
dim(gsea_pos_ol(tg.glgc.1.2.3v2_ol))#pathways with up-regulated genes

gsea_neg_ol <- function(x) {
  #positive expression overlap
  gsea_neg_fx <- read.csv( "gsea_report_for_na_neg_1520864978540.csv", header =TRUE, sep = ",")
  gsea_neg_fx <- gsea_neg_fx[,c(1,8)]
  gsea_neg_sig25_fx <- subset(gsea_neg_fx, FDR.q.val < 0.25)
  x$NAME <- toupper(x$pathway)
  gsea_neg_sig25_fx$NAME <- gsub(" ", "_", gsea_neg_sig25_fx$NAME, fixed=TRUE)
  gsea_tg.glgc_paris_kegg_neg_fx<- merge(gsea_neg_sig25_fx, x,  by = "NAME")
  #negative expression overlap
  return(subset(gsea_tg.glgc_paris_kegg_neg_fx, n > 7))
}

gsea_neg_ol(tg.glgc.1.2.3v2_ol)#Overlap with down-regulated genes
dim(gsea_neg_ol(tg.glgc.1.2.3v2_ol))#Pathways with down-regulated genes



##############################################################################################################
##############################################################################################################
##############################################################################################################

setwd("~/Desktop/projects/POST/lipids_GLGC_GERA/output/paris_results/biofilter_2.4.2/glgc/version4/")


hdl.glgc.1 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed1.paris-summary_formatted", header = TRUE)
hdl.glgc.2 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed2.paris-summary_formatted", header = TRUE)
hdl.glgc.3 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed3.paris-summary_formatted", header = TRUE)
hdl.glgc.4 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed4.paris-summary_formatted", header = TRUE)
hdl.glgc.5 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed5.paris-summary_formatted", header = TRUE)
hdl.glgc.6 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed6.paris-summary_formatted", header = TRUE)
hdl.glgc.7 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed7.paris-summary_formatted", header = TRUE)
hdl.glgc.8 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed8.paris-summary_formatted", header = TRUE)
hdl.glgc.9 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed9.paris-summary_formatted", header = TRUE)
hdl.glgc.10 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed10.paris-summary_formatted", header = TRUE)


hdl.glgc.1_pp <- process_paris(hdl.glgc.1, "hdl.glgc.1")
hdl.glgc.2_pp <- process_paris(hdl.glgc.2, "hdl.glgc.2")
hdl.glgc.3_pp <- process_paris(hdl.glgc.3, "hdl.glgc.3")
hdl.glgc.4_pp <- process_paris(hdl.glgc.4, "hdl.glgc.4")
hdl.glgc.5_pp <- process_paris(hdl.glgc.5, "hdl.glgc.5")
hdl.glgc.6_pp <- process_paris(hdl.glgc.6, "hdl.glgc.6")
hdl.glgc.7_pp <- process_paris(hdl.glgc.7, "hdl.glgc.7")
hdl.glgc.8_pp <- process_paris(hdl.glgc.8, "hdl.glgc.8")
hdl.glgc.9_pp <- process_paris(hdl.glgc.9, "hdl.glgc.9")
hdl.glgc.10_pp <- process_paris(hdl.glgc.10, "hdl.glgc.10")


hdl.glgc.1to10_pp <- rbind.fill(hdl.glgc.1_pp,hdl.glgc.2_pp,hdl.glgc.3_pp,hdl.glgc.4_pp,hdl.glgc.5_pp,hdl.glgc.6_pp,hdl.glgc.7_pp,hdl.glgc.8_pp,hdl.glgc.9_pp,hdl.glgc.10_pp)
dim(hdl.glgc.1to10_pp)
head(hdl.glgc.1to10_pp)
dim(hdl.glgc.1to10_pp %>% filter(p_vals < 0.01))
dim(hdl.glgc.1to10_pp %>% filter(p_vals <= 0.01))

hdl.glgc.1.2.3v2_ol <- as.data.frame(hdl.glgc.1to10_pp %>% filter(p_vals < 0.01) %>% group_by(pathway) %>% summarize(n = n() ))

hdl.glgc.1.2.3v2_ol
dim(hdl.glgc.1.2.3v2_ol)#Total significant pathways (18)
##############################################################################################################

#inesrt GSEA results

setwd("/Users/jemiller2/Desktop/projects/POST/data/POST/theusch_2016/gsea_tpj201612x2_analysis/output/tpj201612x2_rna_seq_order_sig_gsea_hgnc_2018lokiKEGG_plot50.GseaPreranked.1520864978540/")



###############Functions

gsea_pos_ol(hdl.glgc.1.2.3v2_ol)#Overlap with up-regulated genes
dim(gsea_pos_ol(hdl.glgc.1.2.3v2_ol))#pathways with up-regulated genes (6)


gsea_neg_ol(hdl.glgc.1.2.3v2_ol)#Overlap with down-regulated genes
dim(gsea_neg_ol(hdl.glgc.1.2.3v2_ol))#Pathways with down-regulated genes

##############################################################################################################
##############################################################################################################
##############################################################################################################

setwd("~/Desktop/projects/POST/lipids_GLGC_GERA/output/paris_results/biofilter_2.4.2/glgc/version4/")


ldl.glgc.1 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed1.paris-summary_formatted", header = TRUE)
ldl.glgc.2 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed2.paris-summary_formatted", header = TRUE)
ldl.glgc.3 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed3.paris-summary_formatted", header = TRUE)
ldl.glgc.4 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed4.paris-summary_formatted", header = TRUE)
ldl.glgc.5 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed5.paris-summary_formatted", header = TRUE)
ldl.glgc.6 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed6.paris-summary_formatted", header = TRUE)
ldl.glgc.7 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed7.paris-summary_formatted", header = TRUE)
ldl.glgc.8 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed8.paris-summary_formatted", header = TRUE)
ldl.glgc.9 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed9.paris-summary_formatted", header = TRUE)
ldl.glgc.10 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed10.paris-summary_formatted", header = TRUE)


ldl.glgc.1_pp <- process_paris(ldl.glgc.1, "ldl.glgc.1")
ldl.glgc.2_pp <- process_paris(ldl.glgc.2, "ldl.glgc.2")
ldl.glgc.3_pp <- process_paris(ldl.glgc.3, "ldl.glgc.3")
ldl.glgc.4_pp <- process_paris(ldl.glgc.4, "ldl.glgc.4")
ldl.glgc.5_pp <- process_paris(ldl.glgc.5, "ldl.glgc.5")
ldl.glgc.6_pp <- process_paris(ldl.glgc.6, "ldl.glgc.6")
ldl.glgc.7_pp <- process_paris(ldl.glgc.7, "ldl.glgc.7")
ldl.glgc.8_pp <- process_paris(ldl.glgc.8, "ldl.glgc.8")
ldl.glgc.9_pp <- process_paris(ldl.glgc.9, "ldl.glgc.9")
ldl.glgc.10_pp <- process_paris(ldl.glgc.10, "ldl.glgc.10")

head(ldl.glgc.1_pp)

ldl.glgc.1to10_pp <- rbind.fill(ldl.glgc.1_pp,ldl.glgc.2_pp,ldl.glgc.3_pp,ldl.glgc.4_pp,ldl.glgc.5_pp,ldl.glgc.6_pp,ldl.glgc.7_pp,ldl.glgc.8_pp,ldl.glgc.9_pp,ldl.glgc.10_pp)
dim(ldl.glgc.1to10_pp)

ldl.glgc.1.2.3v2_ol <- as.data.frame(ldl.glgc.1to10_pp %>% filter(p_vals < 0.01) %>% group_by(pathway) %>% summarize(n = n() ))
ldl.glgc.1.2.3v2_ol
dim(ldl.glgc.1.2.3v2_ol)#Total significant pathways  (26)
##############################################################################################################

#inesrt GSEA results

setwd("/Users/jemiller2/Desktop/projects/POST/data/POST/theusch_2016/gsea_tpj201612x2_analysis/output/tpj201612x2_rna_seq_order_sig_gsea_hgnc_2018lokiKEGG_plot50.GseaPreranked.1520864978540/")



###############Functions

gsea_pos_ol(ldl.glgc.1.2.3v2_ol)#Overlap with up-regulated genes
dim(gsea_pos_ol(ldl.glgc.1.2.3v2_ol))#pathways with up-regulated genes (9)


gsea_neg_ol(ldl.glgc.1.2.3v2_ol)#Overlap with down-regulated genes
dim(gsea_neg_ol(ldl.glgc.1.2.3v2_ol))#Pathways with down-regulated genes

##############################################################################################################
##############################################################################################################
##############################################################################################################


setwd("~/Desktop/projects/POST/lipids_GLGC_GERA/output/paris_results/biofilter_2.4.2/glgc/version4/")


tc.glgc.1 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_tc.txt_kegg_prep_v4.1_loop_v4.1_seed1.paris-summary_formatted", header = TRUE)
tc.glgc.2 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_tc.txt_kegg_prep_v4.1_loop_v4.1_seed2.paris-summary_formatted", header = TRUE)
tc.glgc.3 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_tc.txt_kegg_prep_v4.1_loop_v4.1_seed3.paris-summary_formatted", header = TRUE)
tc.glgc.4 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_tc.txt_kegg_prep_v4.1_loop_v4.1_seed4.paris-summary_formatted", header = TRUE)
tc.glgc.5 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_tc.txt_kegg_prep_v4.1_loop_v4.1_seed5.paris-summary_formatted", header = TRUE)
tc.glgc.6 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_tc.txt_kegg_prep_v4.1_loop_v4.1_seed6.paris-summary_formatted", header = TRUE)
tc.glgc.7 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_tc.txt_kegg_prep_v4.1_loop_v4.1_seed7.paris-summary_formatted", header = TRUE)
tc.glgc.8 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_tc.txt_kegg_prep_v4.1_loop_v4.1_seed8.paris-summary_formatted", header = TRUE)
tc.glgc.9 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_tc.txt_kegg_prep_v4.1_loop_v4.1_seed9.paris-summary_formatted", header = TRUE)
tc.glgc.10 <- read.table(file = "paris_hg19_pval_input_glgc_jointGwasMc_tc.txt_kegg_prep_v4.1_loop_v4.1_seed10.paris-summary_formatted", header = TRUE)


tc.glgc.1_pp <- process_paris(tc.glgc.1, "tc.glgc.1")
tc.glgc.2_pp <- process_paris(tc.glgc.2, "tc.glgc.2")
tc.glgc.3_pp <- process_paris(tc.glgc.3, "tc.glgc.3")
tc.glgc.4_pp <- process_paris(tc.glgc.4, "tc.glgc.4")
tc.glgc.5_pp <- process_paris(tc.glgc.5, "tc.glgc.5")
tc.glgc.6_pp <- process_paris(tc.glgc.6, "tc.glgc.6")
tc.glgc.7_pp <- process_paris(tc.glgc.7, "tc.glgc.7")
tc.glgc.8_pp <- process_paris(tc.glgc.8, "tc.glgc.8")
tc.glgc.9_pp <- process_paris(tc.glgc.9, "tc.glgc.9")
tc.glgc.10_pp <- process_paris(tc.glgc.10, "tc.glgc.10")


tc.glgc.1to10_pp <- rbind.fill(tc.glgc.1_pp,tc.glgc.2_pp,tc.glgc.3_pp,tc.glgc.4_pp,tc.glgc.5_pp,tc.glgc.6_pp,tc.glgc.7_pp,tc.glgc.8_pp,tc.glgc.9_pp,tc.glgc.10_pp)
dim(tc.glgc.1to10_pp)

tc.glgc.1.2.3v2_ol <- as.data.frame(tc.glgc.1to10_pp %>% filter(p_vals < 0.01) %>% group_by(pathway) %>% summarize(n = n() ))
tc.glgc.1.2.3v2_ol
dim(tc.glgc.1.2.3v2_ol)#Total significant pathways 
##############################################################################################################

#inesrt GSEA results

setwd("/Users/jemiller2/Desktop/projects/POST/data/POST/theusch_2016/gsea_tpj201612x2_analysis/output/tpj201612x2_rna_seq_order_sig_gsea_hgnc_2018lokiKEGG_plot50.GseaPreranked.1520864978540/")



###############Functions

gsea_pos_ol(tc.glgc.1.2.3v2_ol)#Overlap with up-regulated genes
dim(gsea_pos_ol(tc.glgc.1.2.3v2_ol))#pathways with up-regulated genes (11)


gsea_neg_ol(tc.glgc.1.2.3v2_ol)#Overlap with down-regulated genes
dim(gsea_neg_ol(tc.glgc.1.2.3v2_ol))#Pathways with down-regulated genes

##############################################################################################################
##############################################################################################################
##############################################################################################################

#########   ###########   ##########      ####
##         ##             ##      ##    ##   ##
##         ##             ########     ##     ##
##    ###  ##########     ##     ##    #########
##     ##  ##             ##     ##   ##      ##
#########  ############   ##     ##   ##      ##

##############################################################################################################
##############################################################################################################
##############################################################################################################


setwd("~/Desktop/projects/POST/lipids_GLGC_GERA/output/paris_results/biofilter_2.4.2/gera/version4/")

tg.gera.1 <- read.table(file = "paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed1.paris-summary_formatted", header = TRUE)
tg.gera.2 <- read.table(file = "paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed2.paris-summary_formatted", header = TRUE)
tg.gera.3 <- read.table(file = "paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed3.paris-summary_formatted", header = TRUE)
tg.gera.4 <- read.table(file = "paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed4.paris-summary_formatted", header = TRUE)
tg.gera.5 <- read.table(file = "paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed5.paris-summary_formatted", header = TRUE)
tg.gera.6 <- read.table(file = "paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed6.paris-summary_formatted", header = TRUE)
tg.gera.7 <- read.table(file = "paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed7.paris-summary_formatted", header = TRUE)
tg.gera.8 <- read.table(file = "paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed8.paris-summary_formatted", header = TRUE)
tg.gera.9 <- read.table(file = "paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed9.paris-summary_formatted", header = TRUE)
tg.gera.10 <- read.table(file = "paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed10.paris-summary_formatted", header = TRUE)



process_paris <- function(x, y) {
  
  colnames(x) <- c("pathway" ,"description" ,"genes" ,"features", "simple", "simple.sig.", "complex" ,"complex.sig.pval")
  z <- cbind(y,x)
  z$p_vals2 <- gsub( "<", "", as.character(z$complex.sig.pval) )
  z$p_vals <- as.numeric(gsub( ">=", "", as.character(z$p_vals2) ))
  z
  return(z)
  #paste(colnames(x), y, sep = ".")
}
#sed 's/<//g' | sed 's/>=//g'
tg.gera.1_pp <- process_paris(tg.gera.1, "tg.gera.1")
head(tg.gera.1_pp)
tg.gera.2_pp <- process_paris(tg.gera.2, "tg.gera.2")
tg.gera.3_pp <- process_paris(tg.gera.3, "tg.gera.3")
tg.gera.4_pp <- process_paris(tg.gera.4, "tg.gera.4")
tg.gera.5_pp <- process_paris(tg.gera.5, "tg.gera.5")
tg.gera.6_pp <- process_paris(tg.gera.6, "tg.gera.6")
tg.gera.7_pp <- process_paris(tg.gera.7, "tg.gera.7")
tg.gera.8_pp <- process_paris(tg.gera.8, "tg.gera.8")
tg.gera.9_pp <- process_paris(tg.gera.9, "tg.gera.9")
tg.gera.10_pp <- process_paris(tg.gera.10, "tg.gera.10")


tg.gera.1to10_pp <- rbind.fill(tg.gera.1_pp,tg.gera.2_pp,tg.gera.3_pp,tg.gera.4_pp,tg.gera.5_pp,tg.gera.6_pp,tg.gera.7_pp,tg.gera.8_pp,tg.gera.9_pp,tg.gera.10_pp)
dim(tg.gera.1to10_pp)
head(tg.gera.1to10_pp)

tg.gera.1.2.3v2_ol <- as.data.frame(tg.gera.1to10_pp %>% filter(p_vals < 0.01) %>% group_by(pathway) %>% summarize(n = n() ))
tg.gera.1.2.3v2_ol
dim(tg.gera.1.2.3v2_ol)#Total significant pathways (47)
##############################################################################################################

#inesrt GSEA results

setwd("/Users/jemiller2/Desktop/projects/POST/data/POST/theusch_2016/gsea_tpj201612x2_analysis/output/tpj201612x2_rna_seq_order_sig_gsea_hgnc_2018lokiKEGG_plot50.GseaPreranked.1520864978540/")



###############Functions
gsea_pos_ol <- function(x) {
  #positive expression overlap
  gsea_pos_fx <- read.csv( "gsea_report_for_na_pos_1520864978540.csv", header =TRUE, sep = ",")
  gsea_pos_fx <- gsea_pos_fx[,c(1,8)]
  gsea_pos_sig25_fx <- subset(gsea_pos_fx, FDR.q.val < 0.25)
  x$NAME <- toupper(x$pathway)
  gsea_pos_sig25_fx$NAME <- gsub(" ", "_", gsea_pos_sig25_fx$NAME, fixed=TRUE)
  gsea_tg.gera_paris_kegg_fx<- merge(gsea_pos_sig25_fx, x,  by = "NAME")
  #positive expression overlap
  return(subset(gsea_tg.gera_paris_kegg_fx, n > 7))
}

gsea_pos_ol(tg.gera.1.2.3v2_ol)#Overlap with up-regulated genes 
dim(gsea_pos_ol(tg.gera.1.2.3v2_ol))#pathways with up-regulated genes (7)

gsea_neg_ol <- function(x) {
  #positive expression overlap
  gsea_neg_fx <- read.csv( "gsea_report_for_na_neg_1520864978540.csv", header =TRUE, sep = ",")
  gsea_neg_fx <- gsea_neg_fx[,c(1,8)]
  gsea_neg_sig25_fx <- subset(gsea_neg_fx, FDR.q.val < 0.25)
  x$NAME <- toupper(x$pathway)
  gsea_neg_sig25_fx$NAME <- gsub(" ", "_", gsea_neg_sig25_fx$NAME, fixed=TRUE)
  gsea_tg.gera_paris_kegg_neg_fx<- merge(gsea_neg_sig25_fx, x,  by = "NAME")
  #negative expression overlap
  return(subset(gsea_tg.gera_paris_kegg_neg_fx, n > 7))
}

gsea_neg_ol(tg.gera.1.2.3v2_ol)#Overlap with down-regulated genes
dim(gsea_neg_ol(tg.gera.1.2.3v2_ol))#Pathways with down-regulated genes (3)

##############################################################################################################
##############################################################################################################
##############################################################################################################

setwd("~/Desktop/projects/POST/lipids_GLGC_GERA/output/paris_results/biofilter_2.4.2/gera/version4/")


hdl.gera.1 <- read.table(file = "paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed1.paris-summary_formatted", header = TRUE)
hdl.gera.2 <- read.table(file = "paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed2.paris-summary_formatted", header = TRUE)
hdl.gera.3 <- read.table(file = "paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed3.paris-summary_formatted", header = TRUE)
hdl.gera.4 <- read.table(file = "paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed4.paris-summary_formatted", header = TRUE)
hdl.gera.5 <- read.table(file = "paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed5.paris-summary_formatted", header = TRUE)
hdl.gera.6 <- read.table(file = "paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed6.paris-summary_formatted", header = TRUE)
hdl.gera.7 <- read.table(file = "paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed7.paris-summary_formatted", header = TRUE)
hdl.gera.8 <- read.table(file = "paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed8.paris-summary_formatted", header = TRUE)
hdl.gera.9 <- read.table(file = "paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed9.paris-summary_formatted", header = TRUE)
hdl.gera.10 <- read.table(file = "paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed10.paris-summary_formatted", header = TRUE)


hdl.gera.1_pp <- process_paris(hdl.gera.1, "hdl.gera.1")
hdl.gera.2_pp <- process_paris(hdl.gera.2, "hdl.gera.2")
hdl.gera.3_pp <- process_paris(hdl.gera.3, "hdl.gera.3")
hdl.gera.4_pp <- process_paris(hdl.gera.4, "hdl.gera.4")
hdl.gera.5_pp <- process_paris(hdl.gera.5, "hdl.gera.5")
hdl.gera.6_pp <- process_paris(hdl.gera.6, "hdl.gera.6")
hdl.gera.7_pp <- process_paris(hdl.gera.7, "hdl.gera.7")
hdl.gera.8_pp <- process_paris(hdl.gera.8, "hdl.gera.8")
hdl.gera.9_pp <- process_paris(hdl.gera.9, "hdl.gera.9")
hdl.gera.10_pp <- process_paris(hdl.gera.10, "hdl.gera.10")


hdl.gera.1to10_pp <- rbind.fill(hdl.gera.1_pp,hdl.gera.2_pp,hdl.gera.3_pp,hdl.gera.4_pp,hdl.gera.5_pp,hdl.gera.6_pp,hdl.gera.7_pp,hdl.gera.8_pp,hdl.gera.9_pp,hdl.gera.10_pp)
dim(hdl.gera.1to10_pp)
head(hdl.gera.1to10_pp)
dim(hdl.gera.1to10_pp %>% filter(p_vals < 0.01))
dim(hdl.gera.1to10_pp %>% filter(p_vals <= 0.01))

hdl.gera.1.2.3v2_ol <- as.data.frame(hdl.gera.1to10_pp %>% filter(p_vals < 0.01) %>% group_by(pathway) %>% summarize(n = n() ))

hdl.gera.1.2.3v2_ol

dim(hdl.gera.1.2.3v2_ol)#Total significant pathways (38)
##############################################################################################################

#inesrt GSEA results

setwd("/Users/jemiller2/Desktop/projects/POST/data/POST/theusch_2016/gsea_tpj201612x2_analysis/output/tpj201612x2_rna_seq_order_sig_gsea_hgnc_2018lokiKEGG_plot50.GseaPreranked.1520864978540/")



###############Functions

gsea_pos_ol(hdl.gera.1.2.3v2_ol)#Overlap with up-regulated genes(9)
dim(gsea_pos_ol(hdl.gera.1.2.3v2_ol))#pathways with up-regulated genes


gsea_neg_ol(hdl.gera.1.2.3v2_ol)#Overlap with down-regulated genes
dim(gsea_neg_ol(hdl.gera.1.2.3v2_ol))#Pathways with down-regulated genes

##############################################################################################################
##############################################################################################################
##############################################################################################################

setwd("~/Desktop/projects/POST/lipids_GLGC_GERA/output/paris_results/biofilter_2.4.2/gera/version4/")


ldl.gera.1 <- read.table(file = "paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed1.paris-summary_formatted", header = TRUE)
ldl.gera.2 <- read.table(file = "paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed2.paris-summary_formatted", header = TRUE)
ldl.gera.3 <- read.table(file = "paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed3.paris-summary_formatted", header = TRUE)
ldl.gera.4 <- read.table(file = "paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed4.paris-summary_formatted", header = TRUE)
ldl.gera.5 <- read.table(file = "paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed5.paris-summary_formatted", header = TRUE)
ldl.gera.6 <- read.table(file = "paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed6.paris-summary_formatted", header = TRUE)
ldl.gera.7 <- read.table(file = "paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed7.paris-summary_formatted", header = TRUE)
ldl.gera.8 <- read.table(file = "paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed8.paris-summary_formatted", header = TRUE)
ldl.gera.9 <- read.table(file = "paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed9.paris-summary_formatted", header = TRUE)
ldl.gera.10 <- read.table(file = "paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed10.paris-summary_formatted", header = TRUE)

ldl.gera.1_pp <- process_paris(ldl.gera.1, "ldl.gera.1")
ldl.gera.2_pp <- process_paris(ldl.gera.2, "ldl.gera.2")
ldl.gera.3_pp <- process_paris(ldl.gera.3, "ldl.gera.3")
ldl.gera.4_pp <- process_paris(ldl.gera.4, "ldl.gera.4")
ldl.gera.5_pp <- process_paris(ldl.gera.5, "ldl.gera.5")
ldl.gera.6_pp <- process_paris(ldl.gera.6, "ldl.gera.6")
ldl.gera.7_pp <- process_paris(ldl.gera.7, "ldl.gera.7")
ldl.gera.8_pp <- process_paris(ldl.gera.8, "ldl.gera.8")
ldl.gera.9_pp <- process_paris(ldl.gera.9, "ldl.gera.9")
ldl.gera.10_pp <- process_paris(ldl.gera.10, "ldl.gera.10")

head(ldl.gera.1_pp)

ldl.gera.1to10_pp <- rbind.fill(ldl.gera.1_pp,ldl.gera.2_pp,ldl.gera.3_pp,ldl.gera.4_pp,ldl.gera.5_pp,ldl.gera.6_pp,ldl.gera.7_pp,ldl.gera.8_pp,ldl.gera.9_pp,ldl.gera.10_pp)
dim(ldl.gera.1to10_pp)

ldl.gera.1.2.3v2_ol <- as.data.frame(ldl.gera.1to10_pp %>% filter(p_vals < 0.01) %>% group_by(pathway) %>% summarize(n = n() ))
ldl.gera.1.2.3v2_ol
dim(ldl.gera.1.2.3v2_ol)#Total significant pathways  (60)
##############################################################################################################

#inesrt GSEA results

setwd("/Users/jemiller2/Desktop/projects/POST/data/POST/theusch_2016/gsea_tpj201612x2_analysis/output/tpj201612x2_rna_seq_order_sig_gsea_hgnc_2018lokiKEGG_plot50.GseaPreranked.1520864978540/")



###############Functions

gsea_pos_ol(ldl.gera.1.2.3v2_ol)#Overlap with up-regulated genes
dim(gsea_pos_ol(ldl.gera.1.2.3v2_ol))#pathways with up-regulated genes (12)


gsea_neg_ol(ldl.gera.1.2.3v2_ol)#Overlap with down-regulated genes
dim(gsea_neg_ol(ldl.gera.1.2.3v2_ol))#Pathways with down-regulated genes

##############################################################################################################
##############################################################################################################
##############################################################################################################


setwd("~/Desktop/projects/POST/lipids_GLGC_GERA/output/paris_results/biofilter_2.4.2/gera/version4/")


tc.gera.1 <- read.table(file = "paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed1.paris-summary_formatted", header = TRUE)
tc.gera.2 <- read.table(file = "paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed2.paris-summary_formatted", header = TRUE)
tc.gera.3 <- read.table(file = "paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed3.paris-summary_formatted", header = TRUE)
tc.gera.4 <- read.table(file = "paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed4.paris-summary_formatted", header = TRUE)
tc.gera.5 <- read.table(file = "paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed5.paris-summary_formatted", header = TRUE)
tc.gera.6 <- read.table(file = "paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed6.paris-summary_formatted", header = TRUE)
tc.gera.7 <- read.table(file = "paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed7.paris-summary_formatted", header = TRUE)
tc.gera.8 <- read.table(file = "paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed8.paris-summary_formatted", header = TRUE)
tc.gera.9 <- read.table(file = "paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed9.paris-summary_formatted", header = TRUE)
tc.gera.10 <- read.table(file = "paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed10.paris-summary_formatted", header = TRUE)




tc.gera.1_pp <- process_paris(tc.gera.1, "tc.gera.1")
tc.gera.2_pp <- process_paris(tc.gera.2, "tc.gera.2")
tc.gera.3_pp <- process_paris(tc.gera.3, "tc.gera.3")
tc.gera.4_pp <- process_paris(tc.gera.4, "tc.gera.4")
tc.gera.5_pp <- process_paris(tc.gera.5, "tc.gera.5")
tc.gera.6_pp <- process_paris(tc.gera.6, "tc.gera.6")
tc.gera.7_pp <- process_paris(tc.gera.7, "tc.gera.7")
tc.gera.8_pp <- process_paris(tc.gera.8, "tc.gera.8")
tc.gera.9_pp <- process_paris(tc.gera.9, "tc.gera.9")
tc.gera.10_pp <- process_paris(tc.gera.10, "tc.gera.10")


tc.gera.1to10_pp <- rbind.fill(tc.gera.1_pp,tc.gera.2_pp,tc.gera.3_pp,tc.gera.4_pp,tc.gera.5_pp,tc.gera.6_pp,tc.gera.7_pp,tc.gera.8_pp,tc.gera.9_pp,tc.gera.10_pp)
dim(tc.gera.1to10_pp)

tc.gera.1.2.3v2_ol <- as.data.frame(tc.gera.1to10_pp %>% filter(p_vals < 0.01) %>% group_by(pathway) %>% summarize(n = n() ))
tc.gera.1.2.3v2_ol
dim(tc.gera.1.2.3v2_ol)#Total significant pathways (55)
##############################################################################################################

#inesrt GSEA results

setwd("/Users/jemiller2/Desktop/projects/POST/data/POST/theusch_2016/gsea_tpj201612x2_analysis/output/tpj201612x2_rna_seq_order_sig_gsea_hgnc_2018lokiKEGG_plot50.GseaPreranked.1520864978540/")



###############Functions

gsea_pos_ol(tc.gera.1.2.3v2_ol)#Overlap with up-regulated genes
dim(gsea_pos_ol(tc.gera.1.2.3v2_ol))#pathways with up-regulated genes (13)


gsea_neg_ol(tc.gera.1.2.3v2_ol)#Overlap with down-regulated genes
dim(gsea_neg_ol(tc.gera.1.2.3v2_ol))#Pathways with down-regulated genes


##############################################################################################################
##############################################################################################################
##############################################################################################################
##
##   Comparing GERA to GLGC 
##
##############################################################################################################
##############################################################################################################
##############################################################################################################


#Total Cholesterol

#GERA
dim(tc.gera.1.2.3v2_ol)
dim(subset(tc.gera.1.2.3v2_ol, n > 7))
tc.gera.1.2.3v2_ol7 <- subset(tc.gera.1.2.3v2_ol, n > 7)


#GLGC
dim(tc.glgc.1.2.3v2_ol)
dim(subset(tc.glgc.1.2.3v2_ol, n > 7))
tc.glgc.1.2.3v2_ol7 <- subset(tc.glgc.1.2.3v2_ol, n > 7)

dim(merge(tc.gera.1.2.3v2_ol7 , tc.glgc.1.2.3v2_ol7, by = "pathway"))#OL = 37
tc.gera.glgc.ol7 <- merge(tc.gera.1.2.3v2_ol7 , tc.glgc.1.2.3v2_ol7, by = "pathway", all = TRUE)
sum(is.na(tc.gera.glgc.ol7$n.x))#Not in gera, in GLGC (9)
sum(is.na(tc.gera.glgc.ol7$n.y))#Not in glgc, in GERA (15)

con.table <- matrix(c(37,9,15,264),2,2)
con.table
chisq.test((con.table))
fisher.test((con.table))#2.2e-16


#GERA
tc.gera.1.2.3v2_ol7 <- subset(tc.gera.1.2.3v2_ol, n > 7)
#GLGC
tc.glgc.1.2.3v2_ol7 <- subset(tc.glgc.1.2.3v2_ol, n > 7)

continTT <- dim(merge(tc.gera.1.2.3v2_ol7 , tc.glgc.1.2.3v2_ol7, by = "pathway"))[1]#OL = 37
merged_data <- merge(tc.gera.1.2.3v2_ol7 , tc.glgc.1.2.3v2_ol7, by = "pathway", all = TRUE)
continTF <- sum(is.na(tc.gera.glgc.ol7$n.x))#Not in gera, in GLGC (9)
continFT <- sum(is.na(tc.gera.glgc.ol7$n.y))#Not in glgc, in GERA (15)
continFF <- 325-continTT-continTF-continFT

con.table <- matrix(c(continTT,continTF,continFT,continFF),2,2)
con.table
chisq.test((con.table))
fisher.test((con.table))#2.2e-16


gera_glgc_stats <- function(gera_all, glgc_all){
  #GERA
  gera7 <- subset(gera_all, n > 7)
  #GLGC
  glgc7 <- subset(glgc_all, n > 7)

  continTT <- dim(merge(gera7 , glgc7, by = "pathway"))[1]#OL = 37
  merged_data <- merge(gera7 , glgc7, by = "pathway", all = TRUE)
  continTF <- sum(is.na(merged_data$n.x))#Not in gera, in GLGC (9)
  continFT <- sum(is.na(merged_data$n.y))#Not in glgc, in GERA (15)
  continFF <- 325-continTT-continTF-continFT

  con.table <- matrix(c(continTT,continTF,continFT,continFF),2,2)
  con.table
  list(continTT,fisher.test((con.table)),chisq.test((con.table)) )
  

}


gera_glgc_stats(tc.gera.1.2.3v2_ol,tc.glgc.1.2.3v2_ol)

#TC
gera_glgc_stats(tc.gera.1.2.3v2_ol,tc.glgc.1.2.3v2_ol)

#LDL
gera_glgc_stats(ldl.gera.1.2.3v2_ol,ldl.glgc.1.2.3v2_ol)

#HDL
gera_glgc_stats(hdl.gera.1.2.3v2_ol,hdl.glgc.1.2.3v2_ol)

#Triglycerides
gera_glgc_stats(tg.gera.1.2.3v2_ol,tg.glgc.1.2.3v2_ol)



####################################
#sanity check
#GERA
dim(ldl.gera.1.2.3v2_ol)
dim(subset(ldl.gera.1.2.3v2_ol, n > 7))
ldl.gera.1.2.3v2_ol7 <- subset(ldl.gera.1.2.3v2_ol, n > 7)
#GLGC
dim(ldl.glgc.1.2.3v2_ol)
dim(subset(ldl.glgc.1.2.3v2_ol, n > 7))
ldl.glgc.1.2.3v2_ol7 <- subset(ldl.glgc.1.2.3v2_ol, n > 7)

dim(merge(ldl.gera.1.2.3v2_ol7 , ldl.glgc.1.2.3v2_ol7, by = "pathway"))#OL = 23
####################################

dim(subset(tc.gera.1.2.3v2_ol, n > 7))
dim(subset(tc.glgc.1.2.3v2_ol, n > 7))

dim(subset(ldl.gera.1.2.3v2_ol, n > 7))
dim(subset(ldl.glgc.1.2.3v2_ol, n > 7))

dim(subset(hdl.gera.1.2.3v2_ol, n > 7))
dim(subset(hdl.glgc.1.2.3v2_ol, n > 7))

dim(subset(tg.gera.1.2.3v2_ol, n > 7))
dim(subset(tg.glgc.1.2.3v2_ol, n > 7))


##############################################################################################################
##############################################################################################################
##############################################################################################################
##
##   GERA and GLGC overlap with up-reg RNA-seq data
##
##############################################################################################################
##############################################################################################################
##############################################################################################################

gsea_pos_ol(tc.glgc.1.2.3v2_ol)[,3:4]
gsea_pos_ol(tc.gera.1.2.3v2_ol)

head(tc.gera.1.2.3v2_ol)



gera_glgc_gsea_stats <- function(gera_all, glgc_all){

  gera7.gs <- gsea_pos_ol(gera_all)[,3:4]
  #GLGC
  glgc7.gs <- gsea_pos_ol(glgc_all)[,3:4]

  continTT.gs <- dim(merge(gera7.gs , glgc7.gs, by = "pathway"))[1]#OL = 10
  merged_data.gs <- merge(gera7.gs , glgc7.gs, by = "pathway", all = TRUE)
  continTF.gs <- sum(is.na(merged_data.gs$n.x))#Not in gera, in GLGC (1)
  continFT.gs <- sum(is.na(merged_data.gs$n.y))#Not in glgc, in GERA (15)
  continFF.gs <- 325-continTT.gs-continTF.gs-continFT.gs

  con.table.gs <- matrix(c(continTT.gs,continTF.gs,continFT.gs,continFF.gs),2,2)
  con.table.gs
  list(continTT.gs,fisher.test(con.table.gs),chisq.test(con.table.gs) )

}

gera_glgc_gsea_stats(tc.gera.1.2.3v2_ol,tc.glgc.1.2.3v2_ol)
merge(gsea_pos_ol(tc.gera.1.2.3v2_ol)[,3:4], gsea_pos_ol(tc.glgc.1.2.3v2_ol)[,3:4], by ="pathway")

gera_glgc_gsea_stats(ldl.gera.1.2.3v2_ol,ldl.glgc.1.2.3v2_ol)
merge(gsea_pos_ol(ldl.gera.1.2.3v2_ol)[,3:4], gsea_pos_ol(ldl.glgc.1.2.3v2_ol)[,3:4], by ="pathway")

gera_glgc_gsea_stats(hdl.gera.1.2.3v2_ol,hdl.glgc.1.2.3v2_ol)
merge(gsea_pos_ol(hdl.gera.1.2.3v2_ol)[,3:4], gsea_pos_ol(hdl.glgc.1.2.3v2_ol)[,3:4], by ="pathway")

gera_glgc_gsea_stats(tg.gera.1.2.3v2_ol,tg.glgc.1.2.3v2_ol)
merge(gsea_pos_ol(tg.gera.1.2.3v2_ol)[,3:4], gsea_pos_ol(tg.glgc.1.2.3v2_ol)[,3:4], by ="pathway")


##############################################################################################################
##############################################################################################################
##############################################################################################################
##
##   GERA gene-based results
##
##############################################################################################################
##############################################################################################################
##############################################################################################################

setwd("~/Desktop/projects/POST/lipids_GLGC_GERA/output/paris_results/biofilter_2.4.2/gera/version4/detailed_output/")
dir()

#paris_hg19 rw/ pathway_genes_paris_hg19
#paris-summary_formatted rw/ paris-detail
#gera. rw gera.genes


tc.gera.genes.1  <- read.table(file = "pathway_genes_paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed1.paris-detail", header = FALSE)
tc.gera.genes.2  <- read.table(file = "pathway_genes_paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed2.paris-detail", header = FALSE)
tc.gera.genes.3  <- read.table(file = "pathway_genes_paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed3.paris-detail", header = FALSE)
tc.gera.genes.4  <- read.table(file = "pathway_genes_paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed4.paris-detail", header = FALSE)
tc.gera.genes.5  <- read.table(file = "pathway_genes_paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed5.paris-detail", header = FALSE)
tc.gera.genes.6  <- read.table(file = "pathway_genes_paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed6.paris-detail", header = FALSE)
tc.gera.genes.7  <- read.table(file = "pathway_genes_paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed7.paris-detail", header = FALSE)
tc.gera.genes.8  <- read.table(file = "pathway_genes_paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed8.paris-detail", header = FALSE)
tc.gera.genes.9  <- read.table(file = "pathway_genes_paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed9.paris-detail", header = FALSE)
tc.gera.genes.10 <- read.table(file = "pathway_genes_paris_hg19_pval_input_TC-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed10.paris-detail", header = FALSE)

head(tc.gera.genes.1)

tc.gera.genes.1$seed <- 1
tc.gera.genes.2$seed <- 2
tc.gera.genes.3$seed <- 3
tc.gera.genes.4$seed <- 4
tc.gera.genes.5$seed <- 5
tc.gera.genes.6$seed <- 6
tc.gera.genes.7$seed <- 7
tc.gera.genes.8$seed <- 8
tc.gera.genes.9$seed <- 9
tc.gera.genes.10$seed <- 10


tc.gera.genes.1to10_pp <- rbind.fill(tc.gera.genes.1, 
                                     tc.gera.genes.2, 
                                     tc.gera.genes.3, 
                                     tc.gera.genes.4, 
                                     tc.gera.genes.5, 
                                     tc.gera.genes.6, 
                                     tc.gera.genes.7, 
                                     tc.gera.genes.8, 
                                     tc.gera.genes.9, 
                                     tc.gera.genes.10)
dim(tc.gera.genes.1to10_pp)
colnames(tc.gera.genes.1to10_pp) <- c("pathway","gene","p","seed")
head(tc.gera.genes.1to10_pp)

#tc.gera.genes.1to10_ol <- 

#Count 
tc.gera.genes.1to10_pp %>% 
  filter( p < 0.01 , pathway %in%  c("Hepatitis_C" ,"Cholesterol_metabolism" )) %>% 
  group_by(pathway, gene) %>% summarize(n = n())

#LDLR in HC and CM

########################################################################################################################################################
########################################################################################################################################################

TG.gera.genes.1  <- read.table(file = "pathway_genes_paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed1.paris-detail", header = FALSE)
TG.gera.genes.2  <- read.table(file = "pathway_genes_paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed2.paris-detail", header = FALSE)
TG.gera.genes.3  <- read.table(file = "pathway_genes_paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed3.paris-detail", header = FALSE)
TG.gera.genes.4  <- read.table(file = "pathway_genes_paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed4.paris-detail", header = FALSE)
TG.gera.genes.5  <- read.table(file = "pathway_genes_paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed5.paris-detail", header = FALSE)
TG.gera.genes.6  <- read.table(file = "pathway_genes_paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed6.paris-detail", header = FALSE)
TG.gera.genes.7  <- read.table(file = "pathway_genes_paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed7.paris-detail", header = FALSE)
TG.gera.genes.8  <- read.table(file = "pathway_genes_paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed8.paris-detail", header = FALSE)
TG.gera.genes.9  <- read.table(file = "pathway_genes_paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed9.paris-detail", header = FALSE)
TG.gera.genes.10 <- read.table(file = "pathway_genes_paris_hg19_pval_input_logTG-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed10.paris-detail", header = FALSE)

head(TG.gera.genes.1)

TG.gera.genes.1$seed <- 1
TG.gera.genes.2$seed <- 2
TG.gera.genes.3$seed <- 3
TG.gera.genes.4$seed <- 4
TG.gera.genes.5$seed <- 5
TG.gera.genes.6$seed <- 6
TG.gera.genes.7$seed <- 7
TG.gera.genes.8$seed <- 8
TG.gera.genes.9$seed <- 9
TG.gera.genes.10$seed <- 10


TG.gera.genes.1to10_pp <- rbind.fill(TG.gera.genes.1, 
                                     TG.gera.genes.2, 
                                     TG.gera.genes.3, 
                                     TG.gera.genes.4, 
                                     TG.gera.genes.5, 
                                     TG.gera.genes.6, 
                                     TG.gera.genes.7, 
                                     TG.gera.genes.8, 
                                     TG.gera.genes.9, 
                                     TG.gera.genes.10)
dim(TG.gera.genes.1to10_pp)
colnames(TG.gera.genes.1to10_pp) <- c("pathway","gene","p","seed")

#TG.gera.genes.1to10_ol <- 

#Count 
TG.gera.genes.1to10_pp %>% 
  filter( p < 0.01 , pathway %in%  c("Hepatitis_C" ,"Cholesterol_metabolism" )) %>% 
  group_by(pathway, gene) %>% summarize(n = n())

########################################################################################################################################################
########################################################################################################################################################

LDL.gera.genes.1  <- read.table(file = "pathway_genes_paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed1.paris-detail", header = FALSE)
LDL.gera.genes.2  <- read.table(file = "pathway_genes_paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed2.paris-detail", header = FALSE)
LDL.gera.genes.3  <- read.table(file = "pathway_genes_paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed3.paris-detail", header = FALSE)
LDL.gera.genes.4  <- read.table(file = "pathway_genes_paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed4.paris-detail", header = FALSE)
LDL.gera.genes.5  <- read.table(file = "pathway_genes_paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed5.paris-detail", header = FALSE)
LDL.gera.genes.6  <- read.table(file = "pathway_genes_paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed6.paris-detail", header = FALSE)
LDL.gera.genes.7  <- read.table(file = "pathway_genes_paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed7.paris-detail", header = FALSE)
LDL.gera.genes.8  <- read.table(file = "pathway_genes_paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed8.paris-detail", header = FALSE)
LDL.gera.genes.9  <- read.table(file = "pathway_genes_paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed9.paris-detail", header = FALSE)
LDL.gera.genes.10 <- read.table(file = "pathway_genes_paris_hg19_pval_input_LDLf-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed10.paris-detail", header = FALSE)

head(LDL.gera.genes.1)

LDL.gera.genes.1$seed <- 1
LDL.gera.genes.2$seed <- 2
LDL.gera.genes.3$seed <- 3
LDL.gera.genes.4$seed <- 4
LDL.gera.genes.5$seed <- 5
LDL.gera.genes.6$seed <- 6
LDL.gera.genes.7$seed <- 7
LDL.gera.genes.8$seed <- 8
LDL.gera.genes.9$seed <- 9
LDL.gera.genes.10$seed <- 10


LDL.gera.genes.1to10_pp <- rbind.fill(LDL.gera.genes.1, 
                                     LDL.gera.genes.2, 
                                     LDL.gera.genes.3, 
                                     LDL.gera.genes.4, 
                                     LDL.gera.genes.5, 
                                     LDL.gera.genes.6, 
                                     LDL.gera.genes.7, 
                                     LDL.gera.genes.8, 
                                     LDL.gera.genes.9, 
                                     LDL.gera.genes.10)
dim(LDL.gera.genes.1to10_pp)
colnames(LDL.gera.genes.1to10_pp) <- c("pathway","gene","p","seed")

#LDL.gera.genes.1to10_ol <- 

#Count 
LDL.gera.genes.1to10_pp %>% 
  filter( p < 0.01 , pathway %in%  c("Hepatitis_C" ,"Cholesterol_metabolism" )) %>% 
  group_by(pathway, gene) %>% summarize(n = n())

#LDLR in HC and CM

########################################################################################################################################################
########################################################################################################################################################

HDL.gera.genes.1  <- read.table(file = "pathway_genes_paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed1.paris-detail", header = FALSE)
HDL.gera.genes.2  <- read.table(file = "pathway_genes_paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed2.paris-detail", header = FALSE)
HDL.gera.genes.3  <- read.table(file = "pathway_genes_paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed3.paris-detail", header = FALSE)
HDL.gera.genes.4  <- read.table(file = "pathway_genes_paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed4.paris-detail", header = FALSE)
HDL.gera.genes.5  <- read.table(file = "pathway_genes_paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed5.paris-detail", header = FALSE)
HDL.gera.genes.6  <- read.table(file = "pathway_genes_paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed6.paris-detail", header = FALSE)
HDL.gera.genes.7  <- read.table(file = "pathway_genes_paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed7.paris-detail", header = FALSE)
HDL.gera.genes.8  <- read.table(file = "pathway_genes_paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed8.paris-detail", header = FALSE)
HDL.gera.genes.9  <- read.table(file = "pathway_genes_paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed9.paris-detail", header = FALSE)
HDL.gera.genes.10 <- read.table(file = "pathway_genes_paris_hg19_pval_input_sqrtHDL-EUR.tsv_kegg_prep_v4.1_loop_v4.1_seed10.paris-detail", header = FALSE)

head(HDL.gera.genes.1)

HDL.gera.genes.1$seed <- 1
HDL.gera.genes.2$seed <- 2
HDL.gera.genes.3$seed <- 3
HDL.gera.genes.4$seed <- 4
HDL.gera.genes.5$seed <- 5
HDL.gera.genes.6$seed <- 6
HDL.gera.genes.7$seed <- 7
HDL.gera.genes.8$seed <- 8
HDL.gera.genes.9$seed <- 9
HDL.gera.genes.10$seed <- 10


HDL.gera.genes.1to10_pp <- rbind.fill(HDL.gera.genes.1, 
                                      HDL.gera.genes.2, 
                                      HDL.gera.genes.3, 
                                      HDL.gera.genes.4, 
                                      HDL.gera.genes.5, 
                                      HDL.gera.genes.6, 
                                      HDL.gera.genes.7, 
                                      HDL.gera.genes.8, 
                                      HDL.gera.genes.9, 
                                      HDL.gera.genes.10)
dim(HDL.gera.genes.1to10_pp)
colnames(HDL.gera.genes.1to10_pp) <- c("pathway","gene","p","seed")

#HDL.gera.genes.1to10_ol <- 

#Count 
HDL.gera.genes.1to10_pp %>% 
  filter( p < 0.01 , pathway %in%  c("Hepatitis_C" ,"Cholesterol_metabolism" )) %>% 
  group_by(pathway, gene) %>% summarize(n = n())

#HDL SCARB1 in HC and CM


##############################################################################################################
##############################################################################################################
##############################################################################################################
##
##   GLGC gene-based analysis
##
##############################################################################################################
##############################################################################################################
##############################################################################################################

setwd("~/Desktop/projects/POST/lipids_GLGC_GERA/output/paris_results/biofilter_2.4.2/glgc/version4/detailed_output/")
dir()

#paris_hg19 rw/ pathway_genes_paris_hg19
#paris-summary_formatted rw/ paris-detail
#gera. rw gera.genes

#replcae TC-EUR.tsv w glgc_jointGwasMc_TC.txt
tc.glgc.genes.1  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TC.txt_kegg_prep_v4.1_loop_v4.1_seed1.paris-detail", header = FALSE)
tc.glgc.genes.2  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TC.txt_kegg_prep_v4.1_loop_v4.1_seed2.paris-detail", header = FALSE)
tc.glgc.genes.3  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TC.txt_kegg_prep_v4.1_loop_v4.1_seed3.paris-detail", header = FALSE)
tc.glgc.genes.4  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TC.txt_kegg_prep_v4.1_loop_v4.1_seed4.paris-detail", header = FALSE)
tc.glgc.genes.5  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TC.txt_kegg_prep_v4.1_loop_v4.1_seed5.paris-detail", header = FALSE)
tc.glgc.genes.6  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TC.txt_kegg_prep_v4.1_loop_v4.1_seed6.paris-detail", header = FALSE)
tc.glgc.genes.7  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TC.txt_kegg_prep_v4.1_loop_v4.1_seed7.paris-detail", header = FALSE)
tc.glgc.genes.8  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TC.txt_kegg_prep_v4.1_loop_v4.1_seed8.paris-detail", header = FALSE)
tc.glgc.genes.9  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TC.txt_kegg_prep_v4.1_loop_v4.1_seed9.paris-detail", header = FALSE)
tc.glgc.genes.10 <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TC.txt_kegg_prep_v4.1_loop_v4.1_seed10.paris-detail", header = FALSE)

head(tc.glgc.genes.1)

tc.glgc.genes.1$seed <- 1
tc.glgc.genes.2$seed <- 2
tc.glgc.genes.3$seed <- 3
tc.glgc.genes.4$seed <- 4
tc.glgc.genes.5$seed <- 5
tc.glgc.genes.6$seed <- 6
tc.glgc.genes.7$seed <- 7
tc.glgc.genes.8$seed <- 8
tc.glgc.genes.9$seed <- 9
tc.glgc.genes.10$seed <- 10


tc.glgc.genes.1to10_pp <- rbind.fill(tc.glgc.genes.1, 
                                     tc.glgc.genes.2, 
                                     tc.glgc.genes.3, 
                                     tc.glgc.genes.4, 
                                     tc.glgc.genes.5, 
                                     tc.glgc.genes.6, 
                                     tc.glgc.genes.7, 
                                     tc.glgc.genes.8, 
                                     tc.glgc.genes.9, 
                                     tc.glgc.genes.10)
dim(tc.glgc.genes.1to10_pp)
colnames(tc.glgc.genes.1to10_pp) <- c("pathway","gene","p","seed")
head(tc.glgc.genes.1to10_pp)

#tc.glgc.genes.1to10_ol <- 

#Count 
tc.glgc.genes.1to10_pp %>% 
  filter( p < 0.01 , pathway %in%  c("Hepatitis_C" ,"Cholesterol_metabolism" )) %>% 
  group_by(pathway, gene) %>% summarize(n = n())

#Total cholesterol
#LDLR in HC and CM

########################################################################################################################################################
########################################################################################################################################################

TG.glgc.genes.1  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed1.paris-detail", header = FALSE)
TG.glgc.genes.2  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed2.paris-detail", header = FALSE)
TG.glgc.genes.3  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed3.paris-detail", header = FALSE)
TG.glgc.genes.4  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed4.paris-detail", header = FALSE)
TG.glgc.genes.5  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed5.paris-detail", header = FALSE)
TG.glgc.genes.6  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed6.paris-detail", header = FALSE)
TG.glgc.genes.7  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed7.paris-detail", header = FALSE)
TG.glgc.genes.8  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed8.paris-detail", header = FALSE)
TG.glgc.genes.9  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed9.paris-detail", header = FALSE)
TG.glgc.genes.10 <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_TG.txt_kegg_prep_v4.1_loop_v4.1_seed10.paris-detail", header = FALSE)

head(TG.glgc.genes.1)

TG.glgc.genes.1$seed <- 1
TG.glgc.genes.2$seed <- 2
TG.glgc.genes.3$seed <- 3
TG.glgc.genes.4$seed <- 4
TG.glgc.genes.5$seed <- 5
TG.glgc.genes.6$seed <- 6
TG.glgc.genes.7$seed <- 7
TG.glgc.genes.8$seed <- 8
TG.glgc.genes.9$seed <- 9
TG.glgc.genes.10$seed <- 10


TG.glgc.genes.1to10_pp <- rbind.fill(TG.glgc.genes.1, 
                                     TG.glgc.genes.2, 
                                     TG.glgc.genes.3, 
                                     TG.glgc.genes.4, 
                                     TG.glgc.genes.5, 
                                     TG.glgc.genes.6, 
                                     TG.glgc.genes.7, 
                                     TG.glgc.genes.8, 
                                     TG.glgc.genes.9, 
                                     TG.glgc.genes.10)
dim(TG.glgc.genes.1to10_pp)
colnames(TG.glgc.genes.1to10_pp) <- c("pathway","gene","p","seed")

#TG.glgc.genes.1to10_ol <- 

#Count 
TG.glgc.genes.1to10_pp %>% 
  filter( p < 0.01 , pathway %in%  c("Hepatitis_C" ,"Cholesterol_metabolism" )) %>% 
  group_by(pathway, gene) %>% summarize(n = n())

########################################################################################################################################################
########################################################################################################################################################

LDL.glgc.genes.1  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed1.paris-detail", header = FALSE)
LDL.glgc.genes.2  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed2.paris-detail", header = FALSE)
LDL.glgc.genes.3  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed3.paris-detail", header = FALSE)
LDL.glgc.genes.4  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed4.paris-detail", header = FALSE)
LDL.glgc.genes.5  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed5.paris-detail", header = FALSE)
LDL.glgc.genes.6  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed6.paris-detail", header = FALSE)
LDL.glgc.genes.7  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed7.paris-detail", header = FALSE)
LDL.glgc.genes.8  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed8.paris-detail", header = FALSE)
LDL.glgc.genes.9  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed9.paris-detail", header = FALSE)
LDL.glgc.genes.10 <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_LDL.txt_kegg_prep_v4.1_loop_v4.1_seed10.paris-detail", header = FALSE)

head(LDL.glgc.genes.1)

LDL.glgc.genes.1$seed <- 1
LDL.glgc.genes.2$seed <- 2
LDL.glgc.genes.3$seed <- 3
LDL.glgc.genes.4$seed <- 4
LDL.glgc.genes.5$seed <- 5
LDL.glgc.genes.6$seed <- 6
LDL.glgc.genes.7$seed <- 7
LDL.glgc.genes.8$seed <- 8
LDL.glgc.genes.9$seed <- 9
LDL.glgc.genes.10$seed <- 10


LDL.glgc.genes.1to10_pp <- rbind.fill(LDL.glgc.genes.1, 
                                      LDL.glgc.genes.2, 
                                      LDL.glgc.genes.3, 
                                      LDL.glgc.genes.4, 
                                      LDL.glgc.genes.5, 
                                      LDL.glgc.genes.6, 
                                      LDL.glgc.genes.7, 
                                      LDL.glgc.genes.8, 
                                      LDL.glgc.genes.9, 
                                      LDL.glgc.genes.10)
dim(LDL.glgc.genes.1to10_pp)
colnames(LDL.glgc.genes.1to10_pp) <- c("pathway","gene","p","seed")

#LDL.glgc.genes.1to10_ol <- 

#Count 
LDL.glgc.genes.1to10_pp %>% 
  filter( p < 0.01 , pathway %in%  c("Hepatitis_C" ,"Cholesterol_metabolism" )) %>% 
  group_by(pathway, gene) %>% summarize(n = n())

#LDL gwas
#LDLR in HC and CM

########################################################################################################################################################
########################################################################################################################################################

HDL.glgc.genes.1  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed1.paris-detail", header = FALSE)
HDL.glgc.genes.2  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed2.paris-detail", header = FALSE)
HDL.glgc.genes.3  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed3.paris-detail", header = FALSE)
HDL.glgc.genes.4  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed4.paris-detail", header = FALSE)
HDL.glgc.genes.5  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed5.paris-detail", header = FALSE)
HDL.glgc.genes.6  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed6.paris-detail", header = FALSE)
HDL.glgc.genes.7  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed7.paris-detail", header = FALSE)
HDL.glgc.genes.8  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed8.paris-detail", header = FALSE)
HDL.glgc.genes.9  <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed9.paris-detail", header = FALSE)
HDL.glgc.genes.10 <- read.table(file = "pathway_genes_paris_hg19_pval_input_glgc_jointGwasMc_HDL.txt_kegg_prep_v4.1_loop_v4.1_seed10.paris-detail", header = FALSE)

head(HDL.glgc.genes.1)

HDL.glgc.genes.1$seed <- 1
HDL.glgc.genes.2$seed <- 2
HDL.glgc.genes.3$seed <- 3
HDL.glgc.genes.4$seed <- 4
HDL.glgc.genes.5$seed <- 5
HDL.glgc.genes.6$seed <- 6
HDL.glgc.genes.7$seed <- 7
HDL.glgc.genes.8$seed <- 8
HDL.glgc.genes.9$seed <- 9
HDL.glgc.genes.10$seed <- 10


HDL.glgc.genes.1to10_pp <- rbind.fill(HDL.glgc.genes.1, 
                                      HDL.glgc.genes.2, 
                                      HDL.glgc.genes.3, 
                                      HDL.glgc.genes.4, 
                                      HDL.glgc.genes.5, 
                                      HDL.glgc.genes.6, 
                                      HDL.glgc.genes.7, 
                                      HDL.glgc.genes.8, 
                                      HDL.glgc.genes.9, 
                                      HDL.glgc.genes.10)
dim(HDL.glgc.genes.1to10_pp)
colnames(HDL.glgc.genes.1to10_pp) <- c("pathway","gene","p","seed")

#HDL.glgc.genes.1to10_ol <- 

#Count 
HDL.glgc.genes.1to10_pp %>% 
  filter( p < 0.01 , pathway %in%  c("Hepatitis_C" ,"Cholesterol_metabolism" )) %>% 
  group_by(pathway, gene) %>% summarize(n = n())

#Scarb1 in CM and HC



##############################################################################################################
##############################################################################################################
##############################################################################################################
##
##   Creating input for pathview
##  
##
##############################################################################################################
##############################################################################################################
##############################################################################################################
#creat file with gene, gera, glgc, and GSEA results using CM and HC genes

#,"Cholesterol_metabolism"

#select all genes with p < 0.01 in hep C pathway
kegg_hcv <- function (x) {
  xdf <- as.data.frame(x %>% 
    filter( p < 0.01 , pathway %in%  c("Hepatitis_C"  )) %>% 
    group_by(pathway, gene) %>% summarize(n = n()))
    return(xdf[,c(-1)])
}
##################################################################
###Genes that are seen across all the data > 64 seeds
head(rbind(HDL.glgc.genes.1to10_pp, HDL.gera.genes.1to10_pp,LDL.glgc.genes.1to10_pp, LDL.gera.genes.1to10_pp,
           TG.glgc.genes.1to10_pp, TG.gera.genes.1to10_pp,tc.glgc.genes.1to10_pp, tc.gera.genes.1to10_pp))

all_glgc_gera_genes <- rbind(HDL.glgc.genes.1to10_pp, HDL.gera.genes.1to10_pp,LDL.glgc.genes.1to10_pp, LDL.gera.genes.1to10_pp,
                             TG.glgc.genes.1to10_pp, TG.gera.genes.1to10_pp,tc.glgc.genes.1to10_pp, tc.gera.genes.1to10_pp) %>% 
                              filter(p < 0.01)%>%  group_by(pathway, gene)  %>% summarise(count = n()) %>% filter(count > 64)
dim(all_glgc_gera_genes)
##################################################################

hdl_glgc_hcv_genes <- kegg_hcv(HDL.glgc.genes.1to10_pp)
ldl_glgc_hcv_genes <- kegg_hcv(LDL.glgc.genes.1to10_pp)
tg_glgc_hcv_genes <- kegg_hcv(TG.glgc.genes.1to10_pp)
tc_glgc_hcv_genes <- kegg_hcv(tc.glgc.genes.1to10_pp)

hdl_gera_hcv_genes <- kegg_hcv(HDL.gera.genes.1to10_pp)
ldl_gera_hcv_genes <- kegg_hcv(LDL.gera.genes.1to10_pp)
tg_gera_hcv_genes <- kegg_hcv(TG.gera.genes.1to10_pp)
tc_gera_hcv_genes <- kegg_hcv(tc.gera.genes.1to10_pp)


hcv_gera_glgc_hdl <- merge(hdl_gera_hcv_genes, hdl_glgc_hcv_genes, by = "gene", all = TRUE)
colnames(hcv_gera_glgc_hdl) <- c("gene","HDL.gera.count","HDL.glgc.count" )

hcv_gera_glgc_ldl <- merge(ldl_gera_hcv_genes, ldl_glgc_hcv_genes, by = "gene", all = TRUE)
colnames(hcv_gera_glgc_ldl) <- c("gene", "LDL.gera.count","LDL.glgc.count" )

hcv_gera_glgc_tc <-  merge(tc_gera_hcv_genes,  tc_glgc_hcv_genes, by = "gene", all = TRUE)
colnames(hcv_gera_glgc_tc) <- c("gene", "TC.gera.count", "TC.glgc.count" )

hcv_gera_glgc_tg <-  merge(tg_gera_hcv_genes,  tg_glgc_hcv_genes, by = "gene", all = TRUE)
colnames(hcv_gera_glgc_tg) <- c("gene", "TG.gera.count", "TG.glgc.count" )

hcv_gera_glgc_hdl.ldl <- merge(hcv_gera_glgc_hdl, hcv_gera_glgc_ldl, by = "gene", all = TRUE)
hcv_gera_glgc_hdl.ldl
hcv_gera_glgc_hdl.ldl.tc <- merge(hcv_gera_glgc_hdl.ldl, hcv_gera_glgc_tc, by = "gene", all = TRUE)
hcv_gera_glgc_hdl.ldl.tc
hcv_gera_glgc_all <- merge(hcv_gera_glgc_hdl.ldl.tc, hcv_gera_glgc_tg, by = "gene", all = TRUE)
hcv_gera_glgc_all

hcv_gera_glgc_all[is.na(hcv_gera_glgc_all)] <- "0"
hcv_gera_glgc_all

setwd("~/Desktop/projects/POST/data/POST/theusch_2016/gsea_tpj201612x2_analysis/output/tpj201612x2_rna_seq_order_sig_gsea_hgnc_2018lokiKEGG_plot50.GseaPreranked.1520864978540/")


gsea_hepc <- read.table(file = "HEPATITIS_C_gsea_output.txt", sep = "\t", header =TRUE)
head(gsea_hepc)


hcv_gera_glgc_gsea <- merge(hcv_gera_glgc_all, gsea_hepc, by.x = "gene", by.y = "PROBE", all = TRUE)
hcv_gera_glgc_gsea





colnames(hcv_gera_glgc_gsea)
hcv_gera_glgc_gsea %>%  select(1:9,16)
hcv_gera_glgc_gsea2 <- hcv_gera_glgc_gsea %>%  select(1:9,16)
hcv_gera_glgc_gsea2

setwd("~/Desktop/projects/POST/lipids_GLGC_GERA/output/gsea_paris_results/")

write.table(hcv_gera_glgc_gsea2, file = "hcv_gera_glgc_up_reg_rna_gsea_all.txt",sep = "\t",quote = FALSE,row.names = FALSE)


kegg_cm <- function (x) {
  xdf <- as.data.frame(x %>% 
                         filter( p < 0.01 , pathway %in%  c("Cholesterol_metabolism"  )) %>% 
                         group_by(pathway, gene) %>% summarize(n = n()))
  return(xdf[,c(-1)])
}


make_pathview_table <- function(HDL.glgc.genes.1to10_pp, LDL.glgc.genes.1to10_pp, TG.glgc.genes.1to10_pp,tc.glgc.genes.1to10_pp,
                                HDL.gera.genes.1to10_pp, LDL.gera.genes.1to10_pp, TG.gera.genes.1to10_pp,tc.gera.genes.1to10_pp){
  
  
  hdl_glgc_path_genes <- kegg_cm(HDL.glgc.genes.1to10_pp)
  ldl_glgc_path_genes <- kegg_cm(LDL.glgc.genes.1to10_pp)
  tg_glgc_path_genes <- kegg_cm(TG.glgc.genes.1to10_pp)
  tc_glgc_path_genes <- kegg_cm(tc.glgc.genes.1to10_pp)
  
  hdl_gera_path_genes <- kegg_cm(HDL.gera.genes.1to10_pp)
  ldl_gera_path_genes <- kegg_cm(LDL.gera.genes.1to10_pp)
  tg_gera_path_genes <- kegg_cm(TG.gera.genes.1to10_pp)
  tc_gera_path_genes <- kegg_cm(tc.gera.genes.1to10_pp)
 
  
  path_gera_glgc_hdl <- merge(hdl_gera_path_genes, hdl_glgc_path_genes, by = "gene", all = TRUE)
  colnames(path_gera_glgc_hdl) <- c("gene","HDL.gera.count","HDL.glgc.count" )
  
  path_gera_glgc_ldl <- merge(ldl_gera_path_genes, ldl_glgc_path_genes, by = "gene", all = TRUE)
  colnames(path_gera_glgc_ldl) <- c("gene", "LDL.gera.count","LDL.glgc.count" )
  
  path_gera_glgc_tc <-  merge(tc_gera_path_genes,  tc_glgc_path_genes, by = "gene", all = TRUE)
  colnames(path_gera_glgc_tc) <- c("gene", "TC.gera.count", "TC.glgc.count" )
  
  path_gera_glgc_tg <-  merge(tg_gera_path_genes,  tg_glgc_path_genes, by = "gene", all = TRUE)
  colnames(path_gera_glgc_tg) <- c("gene", "TG.gera.count", "TG.glgc.count" )
  
  path_gera_glgc_hdl.ldl <- merge(path_gera_glgc_hdl, path_gera_glgc_ldl, by = "gene", all = TRUE)
  
  path_gera_glgc_hdl.ldl.tc <- merge(path_gera_glgc_hdl.ldl, path_gera_glgc_tc, by = "gene", all = TRUE)
  
  path_gera_glgc_all <- merge(path_gera_glgc_hdl.ldl.tc, path_gera_glgc_tg, by = "gene", all = TRUE)
  
  
  path_gera_glgc_all[is.na(path_gera_glgc_all)] <- "0"
  
  setwd("~/Desktop/projects/POST/data/POST/theusch_2016/gsea_tpj201612x2_analysis/output/tpj201612x2_rna_seq_order_sig_gsea_hgnc_2018lokiKEGG_plot50.GseaPreranked.1520864978540/")
  
  gsea_path <- read.table(file = "CHOLESTEROL_METABOLISM_gsea_output.txt", sep = "\t", header =TRUE)
  
  
  path_gera_glgc_gsea <- merge(path_gera_glgc_all, gsea_path, by.x = "gene", by.y = "PROBE", all = TRUE)
  
  #path_gera_glgc_gsea$Core_Enrichment <- ifelse(path_gera_glgc_gsea$CORE.ENRICHMENT == "Yes", 10 , as.character("NA"))
  path_gera_glgc_gsea2 <- path_gera_glgc_gsea %>%  select(1:9,16)
  return(path_gera_glgc_gsea2)
}
  


cm_gera_glgc_gsea2 <- make_pathview_table(HDL.glgc.genes.1to10_pp, LDL.glgc.genes.1to10_pp, TG.glgc.genes.1to10_pp,tc.glgc.genes.1to10_pp,
                    HDL.gera.genes.1to10_pp, LDL.gera.genes.1to10_pp, TG.gera.genes.1to10_pp,tc.gera.genes.1to10_pp)




setwd("~/Desktop/projects/POST/lipids_GLGC_GERA/output/gsea_paris_results/")



write.table(cm_gera_glgc_gsea2, file = "cm_gera_glgc_up_reg_rna_gsea_all.txt",sep = "\t",quote = FALSE,row.names = FALSE)




#######################################################################################################################


##############################################################################################################
##############################################################################################################
##############################################################################################################
##
##   Heatmap Viz
##  
##
##############################################################################################################
##############################################################################################################
##############################################################################################################

################################################################
################################################################
#Heat map cluster

tg.glgc.gera_merge <- merge(tg.glgc.1.2.3v2_ol, tg.gera.1.2.3v2_ol, by = "pathway" , all = TRUE )
colnames(tg.glgc.gera_merge) <- c("pathway","TG_GLGC","TG_GERA")

tc.glgc.gera_merge <- merge(tc.glgc.1.2.3v2_ol, tc.gera.1.2.3v2_ol, by = "pathway" , all = TRUE )
colnames(tc.glgc.gera_merge) <- c("pathway","TC_GLGC","TC_GERA")

ldl.glgc.gera_merge <- merge(ldl.glgc.1.2.3v2_ol, ldl.gera.1.2.3v2_ol, by = "pathway" , all = TRUE )
colnames(ldl.glgc.gera_merge) <- c("pathway","LDL_GLGC","LDL_GERA")

hdl.glgc.gera_merge <- merge(hdl.glgc.1.2.3v2_ol, hdl.gera.1.2.3v2_ol, by = "pathway" , all = TRUE )
colnames(hdl.glgc.gera_merge) <- c("pathway","HDL_GLGC","HDL_GERA")


tg_tc_glgc.gera_merge <- merge(tg.glgc.gera_merge, tc.glgc.gera_merge, by = "pathway", all = TRUE)
tg_tc_ldl_glgc.gera_merge <- merge(tg_tc_glgc.gera_merge, ldl.glgc.gera_merge, by = "pathway", all = TRUE)
all_counts_glgc.gera_merge <- merge(tg_tc_ldl_glgc.gera_merge, hdl.glgc.gera_merge, by = "pathway", all = TRUE)


head(all_counts_glgc.gera_merge)

################################################################
################################################################
up_all_counts_glgc.gera_merge <- all_counts_glgc.gera_merge
up_all_counts_glgc.gera_merge$PATHWAY <- toupper(up_all_counts_glgc.gera_merge$pathway)
head(up_all_counts_glgc.gera_merge)

###PLOT with GSEA
library("gplots")
library("devtools")

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

setwd("/Users/jemiller2/Desktop/projects/POST/data/POST/theusch_2016/gsea_tpj201612x2_analysis/output/tpj201612x2_rna_seq_order_sig_gsea_hgnc_2018lokiKEGG_plot50.GseaPreranked.1520864978540/")

gsea_pos <- as.data.frame(read.csv( "gsea_report_for_na_pos_1520864978540.csv", header =TRUE, sep = ","))
gsea_neg <- as.data.frame(read.csv( "gsea_report_for_na_neg_1520864978540.csv", header =TRUE, sep = ","))
colnames(gsea_pos)
colnames(gsea_neg)
gsea_all <- rbind(gsea_pos,gsea_neg)
dim(gsea_all)
head(gsea_all)
#FDR.q.val < 0.25

library(R.utils)
gsea_all <- as.data.frame(gsea_all %>% select(NAME, FDR.q.val )) 
head(gsea_all)

gsea_all$Name2 <- gsub( " ", "_", as.character(gsea_all$NAME) )

#all_counts_glgc.gera_merge.df <- as.data.frame(all_counts_glgc.gera_merge)
#class(gsea_all)
#class(all_counts_glgc.gera_merge.df)
up_all_counts_glgc.gera_merge[is.na(up_all_counts_glgc.gera_merge)] <- 0
up_all_counts_glgc.gera.gsea_merge <- merge(up_all_counts_glgc.gera_merge, gsea_all, by.x = "PATHWAY", by.y = "Name2", all.x = TRUE )
dim(up_all_counts_glgc.gera.gsea_merge)

head(up_all_counts_glgc.gera.gsea_merge)
up_all_counts_glgc.gera.gsea_merge$FDR.q.val

colnames(up_all_counts_glgc.gera.gsea_merge)
up_all_counts_glgc.gera.gsea_merge_input <- up_all_counts_glgc.gera.gsea_merge[,-c(2,11)]
head(up_all_counts_glgc.gera.gsea_merge_input)

rownames(up_all_counts_glgc.gera.gsea_merge_input) <- up_all_counts_glgc.gera.gsea_merge_input[,1]
up_all_counts_glgc.gera.gsea_merge_input <- up_all_counts_glgc.gera.gsea_merge_input[,-1]
head(up_all_counts_glgc.gera.gsea_merge_input)
blue_brew2 <- c('#f0f9e8','#bae4bc','#7bccc4','#43a2ca','#0868ac')
#pdf(file = "all_counts_glgc.gera_merge_v3.heatmap.pdf", width = 20, height = 34)
heatmap.2(as.matrix(up_all_counts_glgc.gera.gsea_merge_input),na.rm = TRUE, trace = "none",margins=c(16,50),
          dendrogram = "both",col = blue_brew2, sepcolor = "black",colsep = 1:9, rowsep = 1:98,
          sepwidth =c(0.01,0.01), cexRow = 1.75,cexCol = 3,na.color = "grey",
          keysize=0.75, key.par = list(cex=0.5))

#dev.off()


mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

heatmap.3(as.matrix(up_all_counts_glgc.gera.gsea_merge_input), hclustfun=myclust, distfun=mydist, na.rm = FALSE, scale="none", 
          dendrogram="both", margins=c(13,12),
          Rowv=TRUE, Colv=TRUE, RowSideColors = rlab_cmM, symbreaks=FALSE, key=TRUE, symkey=FALSE,
          density.info="none", trace="none", main=main_title, labRow=cm_genes2_names,labCol = cm_genes2_colnames2 , 
          cexRow=1.4,cexCol = 1.7, col=light_blue_col,
          RowSideColorsSize=2, KeyValueName="", sepcolor = "black",
          sepwidth =c(0.01,0.01),colsep = 1:8, rowsep = 1:27)
legend(.03,.79,legend=c("Yes","No","NA"),
       fill=c("#f4a582","black","grey"), border=FALSE, bty="n")



################################################################
################################################################







row.names(all_counts_glgc.gera_merge) <- all_counts_glgc.gera_merge[,1]
all_counts_glgc.gera_merge <- as.matrix(all_counts_glgc.gera_merge[,-1])
head(all_counts_glgc.gera_merge)

library(gplots)
all_counts_glgc.gera_merge[is.na(all_counts_glgc.gera_merge)] <- 0

#my_palette <- colorRampPalette(c("white", "light_blue", "blue"))(n = 1000)
dim(all_counts_glgc.gera_merge)


#pdf(file = "all_counts_glgc.gera_merge.heatmap.pdf", width = 20, height = 11)
heatmap.2(all_counts_glgc.gera_merge,na.rm = TRUE, trace = "none",margins=c(8,20),
          Colv="NA",dendrogram = "row",col = "bluered", sepcolor = "white",colsep = 1:8, rowsep = 1:98,
          sepwidth =c(0.001,0.0001), cexRow = .6,
           keysize=0.75, key.par = list(cex=0.5))
#dev.off()

blue_brew <- c('#ffffcc','#a1dab4','#41b6c4','#2c7fb8','#253494')
green_brew <- c('#ffffcc','#c2e699','#78c679','#31a354','#006837')
violet_brew <- c('#feebe2','#fbb4b9','#f768a1','#c51b8a','#7a0177')
dark_violet_brew <- c('#f1eef6','#d7b5d8','#df65b0','#dd1c77','#980043')
blue_brew2 <- c('#f0f9e8','#bae4bc','#7bccc4','#43a2ca','#0868ac')

##pdf(file = "all_counts_glgc.gera_merge_v2.heatmap.pdf", width = 20, height = 34)
heatmap.2(all_counts_glgc.gera_merge,na.rm = TRUE, trace = "none",margins=c(8,50),
          Colv="NA",dendrogram = "row",col = blue_brew2, sepcolor = "black",colsep = 1:8, rowsep = 1:98,
          sepwidth =c(0.01,0.01), cexRow = 1.75,
          keysize=0.75, key.par = list(cex=0.5))

##dev.off()


#Final no RNA
#pdf(file = "all_counts_glgc.gera_merge_v3.heatmap.pdf", width = 20, height = 34)
heatmap.2(all_counts_glgc.gera_merge,na.rm = TRUE, trace = "none",margins=c(16,50),
          dendrogram = "both",col = blue_brew2, sepcolor = "black",colsep = 1:8, rowsep = 1:98,
          sepwidth =c(0.01,0.01), cexRow = 1.75,cexCol = 3,
          keysize=0.75, key.par = list(cex=0.5))

#dev.off()


heatmap.2(all_counts_glgc.gera_merge,na.rm = FALSE, trace = "none",
          Colv="NA",dendrogram = "row",col = blue_brew2, sepcolor = "black",colsep = 1:8, rowsep = 1:98,
          sepwidth =c(0.001,0.001), cexRow = 1.75,
          keysize=0.75, key.par = list(cex=0.5))



heatmap.2(t(all_counts_glgc.gera_merge),na.rm = TRUE, trace = "none",margins=c(20,8),
          Rowv="NA",dendrogram = "col",col = "bluered", sepcolor = "white",colsep = 1:98, rowsep = 1:8,
          sepwidth =c(0.0001,0.001), cexCol = .8,
          keysize=0.75, key.par = list(cex=0.5))



heatmap.2(t(all_counts_glgc.gera_merge),na.rm = TRUE, trace = "none",margins=c(4,8) )



#install.packages("plotly")
#install.packages("heatmaply")

library(heatmaply)
heatmaply(scale(mtcars), k_row = 3, k_col = 2)


heatmaply(all_counts_glgc.gera_merge, k_row = 3, k_col = 2)
#Margin (bottom,left,top, right)
heatmaply(t(all_counts_glgc.gera_merge), k_row = 2, k_col = 3,margins = c(400,150,0,10),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", 
        midpoint = 5, limits = c(0, 10)), column_text_angle = 65) 

dir.create("folder")
heatmaply(t(all_counts_glgc.gera_merge), k_row = 2, k_col = 3,margins = c(400,150,0,10),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", 
           midpoint = 5, limits = c(0, 10)), column_text_angle = 65,
          file = "folder/all_counts_glgc.gera_merge_heatmaply_plot.html")
#browseURL("folder/all_counts_glgc.gera_merge_heatmaply_plot.html")
pdf(file = "all_counts_glgc.gera_merge.heatmaply.pdf", width = 15, height = 10)
heatmaply(t(all_counts_glgc.gera_merge), k_row = 2, k_col = 3,margins = c(400,150,0,10),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", 
                                                                  midpoint = 5, limits = c(0, 10)), column_text_angle = 65)

dev.off()
getwd()
#install.packages("KEGGREST")
#library(KEGGREST)
#listDatabases()


##########GGPLOT

tg.glgc.1.2.3v2_ol$id <- "TG.GLGC"
tg.gera.1.2.3v2_ol$id <- "TG.GERA"
tc.glgc.1.2.3v2_ol$id <- "TC.GLGC"
tc.gera.1.2.3v2_ol$id <- "TC.GERA"
ldl.glgc.1.2.3v2_ol$id <- "LDL.GLGC"
ldl.gera.1.2.3v2_ol$id <- "LDL.GERA"
hdl.glgc.1.2.3v2_ol$id <- "HDL.GLGC"
hdl.gera.1.2.3v2_ol$id <- "HDL.GERA"


lipids_rbind <- rbind(tg.glgc.1.2.3v2_ol,tg.gera.1.2.3v2_ol,tc.glgc.1.2.3v2_ol,tc.gera.1.2.3v2_ol,
                    ldl.glgc.1.2.3v2_ol,ldl.gera.1.2.3v2_ol,hdl.glgc.1.2.3v2_ol,hdl.gera.1.2.3v2_ol)
library(reshape2)               
dim(lipids_rbind)
head(lipids_rbind)
head(melt(lipids_rbind))
lipids_rbind_melt <- melt(lipids_rbind)[,c(1,2,4)]
head(lipids_rbind_melt)
dim(lipids_rbind_melt)

ggplot(lipids_rbind, aes(pathway, id )) +
  geom_tile(aes(fill = n), colour = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  xlab("List of Pathways") +
  ylab("Dataset") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "times p-value < 0.001")


ggplot(lipids_rbind, aes(id,pathway )) +
  geom_tile(aes(fill = n), color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ylab("List of Pathways") +
  xlab("Dataset") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "times p-value < 0.001")
setwd("/Users/jemiller2/Desktop/projects/POST/lipids_GLGC_GERA/output/gsea_paris_results/")
ggsave(filename = "lipids_heatmap_ggplot.pdf",device = "pdf")




