library(ggplot2)
library(gridExtra)
library(grid)
library(Hmisc)
library(stringr)
library(ggpubr)
library("tidyverse")
# df <- readr::read_table("/Users/ejong19/repos/DIMSdb/DIMSdb/outlist_identified_combined_Zscore_slice_IS.txt")

df <- read.table("/Users/ejong19/repos/DIMSdb/DIMSdb/outlist_identified_combined_Zscore_slice_IS.txt", sep="\t", header=TRUE, quote='"')

mzmed_input = 126.136408975421
#patient <- "P725.4_Zscore"
patient <- "Single_P439.1_Zscore"

pivot_longer(df, cols = ends_with("Zscore"), values_to= "zscore", names_to="sample")

dfs_sample <-  subset(df, mzmed.pgrp==mzmed_input)
dfs_controls <-  dfs_sample[,grepl("C", names(dfs_sample))]
dfs_patients <-  dfs_sample[,!grepl("C|modus|mzmed|assi", names(dfs_sample))]
dfs_patients_exclude_patient <-  dfs_patients[,!grepl(patient, names(dfs_patients))]
dfs_patient <-  dfs_sample[,grepl(patient, names(dfs_sample))]
dfs_controls_t <- t(dfs_controls)
dfs_controls_t_mzmed <- transform(dfs_controls_t, mzmed = mzmed_input)
colnames(dfs_controls_t_mzmed)[1] <-"Zscore"
dfs_patient_mzmed <- transform(dfs_patient, mzmed = mzmed_input)
colnames(dfs_patient_mzmed)[1] <-"Zscore"
dfs_patients_exclude_patient_t <- t(dfs_patients_exclude_patient)
colnames(dfs_patients_exclude_patient_t)[1] <-"Zscore"
dfs_patients_exclude_patient_t_mzmed <- transform(dfs_patients_exclude_patient_t, mzmed = mzmed_input)


plotRepeat <- function(df, df_patient, cohort){
  p <-  ggplot(df, aes(x = mzmed, y = Zscore)) +
    geom_violin() + 
    geom_point() +
    scale_y_continuous(limits=c(-5,20), breaks=seq(-5,20,1)) + 
    ggtitle( paste(mzmed_input, "\n", cohort)) +
    geom_hline(yintercept = 2.0, linetype = "dashed", color = "red") + 
    geom_hline(yintercept = -2.5, linetype = "dashed", color = "red") + 
    theme(
      legend.position="none",
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      plot.title = element_text(color="black", size=12, hjust = 0.5, face="bold")
    ) +
    theme_bw() +  
    geom_point(data=df_patient, aes(x = mzmed, y = Zscore), color='red')
  return(p)
}

grid.arrange(
  plotRepeat(dfs_controls_t_mzmed, dfs_patient_mzmed, "controls"),
  plotRepeat(dfs_patients_exclude_patient_t_mzmed, dfs_patient_mzmed, "patients"),
  ncol=2
)











#  
# p_controls <-  ggplot(dfs_controls_t_mzmed, aes(x = mzmed, y = Zscore)) +
#   geom_violin() +
#   geom_point() +
#   scale_y_continuous(limits=c(-5,20), breaks=seq(-5,20,1)) +
#   ggtitle( paste(mzmed_input, "\nControles")) +
#   geom_hline(yintercept = 2.0, linetype = "dashed", color = "red") +
#   geom_hline(yintercept = -2.5, linetype = "dashed", color = "red") +
#   theme(
#     legend.position="none",
#     axis.title.x=element_blank(),
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank(),
#     plot.title = element_text(color="black", size=12, hjust = 0.5, face="bold")
#   ) +
#   theme_bw() +
#   #label y "Zscore" +
#   geom_point(data=dfs_patient_mzmed, aes(x = mzmed, y = Zscore), color='red')
# p_controls
# 
# dfs_patients_exclude_patient_t <- t(dfs_patients_exclude_patient)
# colnames(dfs_patients_exclude_patient_t)[1] <-"Zscore"
# dfs_patients_exclude_patient_t_mzmed <- transform(dfs_patients_exclude_patient_t, mzmed = mzmed_input)
# 
# p_patients <-  ggplot(dfs_patients_exclude_patient_t_mzmed, aes(x = mzmed, y = Zscore)) +
#   geom_violin() +
#   geom_point() +
#   scale_y_continuous(limits=c(-5,20), breaks=seq(-5,20,1)) +
#   ggtitle( paste(mzmed_input, "\nPatients")) +
#   geom_hline(yintercept = 2.0, linetype = "dashed", color = "red") +
#   geom_hline(yintercept = -2.5, linetype = "dashed", color = "red") +
#   theme(
#     legend.position="none",
#     axis.title.x=element_blank(),
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank(),
#     plot.title = element_text(color="black", size=12, hjust = 0.5, face="bold")
#   ) +
#   theme_bw() +
#   geom_point(data=dfs_patient_mzmed, aes(x = mzmed, y = Zscore), color='red')
# p_patients
# 
# grid.arrange(
#   p_controls,
#   p_patients,
#   ncol=2
# )
