#TG_GLGCparis_theusc2016gsea_analysis_v3.R
#v1based on File =  TG_GLGCparis_theusc2016gsea_analysis_v1.R
#v2 lipids_10seeds... I'm add the gera part
#from tc_analysis_v1
#2018 analysis


rm(list = ls())


conflicts()
library(plyr)
library(tidyverse)



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
  return(gsea_tg.glgc_paris_kegg_fx)
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
  return(gsea_tg.glgc_paris_kegg_neg_fx)
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
dim(gsea_pos_ol(hdl.glgc.1.2.3v2_ol))#pathways with up-regulated genes


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
dim(gsea_pos_ol(ldl.glgc.1.2.3v2_ol))#pathways with up-regulated genes


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
dim(gsea_pos_ol(tc.glgc.1.2.3v2_ol))#pathways with up-regulated genes


gsea_neg_ol(tc.glgc.1.2.3v2_ol)#Overlap with down-regulated genes
dim(gsea_neg_ol(tc.glgc.1.2.3v2_ol))#Pathways with down-regulated genes

##############################################################################################################
##############################################################################################################
##############################################################################################################

