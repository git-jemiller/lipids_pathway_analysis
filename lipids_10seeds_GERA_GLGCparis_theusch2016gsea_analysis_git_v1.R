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
  #negative expression overlap
  return(gsea_tg.gera_paris_kegg_fx)
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
  return(gsea_tg.gera_paris_kegg_neg_fx)
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

gsea_pos_ol(hdl.gera.1.2.3v2_ol)#Overlap with up-regulated genes
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
dim(gsea_pos_ol(ldl.gera.1.2.3v2_ol))#pathways with up-regulated genes


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
dim(gsea_pos_ol(tc.gera.1.2.3v2_ol))#pathways with up-regulated genes


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



