#############################################################
#
# Ref to the ARTICLE
# 
# Ridhdhi Rathore
# Ridhdhi.Rathore@itcarlow.ie
#
# Davide Bulgarelli
# d.bulgarelli@dundee.ac.uk
# 
# Revison July 2017
# 
# script to reproduce ANCOM calculations presented in the manuscript
# 
# Disclaimer: the manuscript is currently submitted, the script might be subjected to changes derived from the revision process
#
#############################################################

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#Requirements for ANCOM package (installed from zip file)
#source("http://bioconductor.org/biocLite.R")
#biocLite("doParallel")
#biocLite("foreach")
#biocLite("ggplot2")
#biocLite("Rcpp")
#biocLite("shiny")
#biocLite("exactRankTests")
#biocLite("openxlsx")
#biocLite("DT")
#biocLite("coin")
#biocLite("grid")
#biocLite("futile.logger")


#Requirements for ANCOM package (installed from zip file)
library("doParallel")
library("foreach")
library("ggplot2")
library("Rcpp")
library("shiny")
library("exactRankTests")
library("openxlsx")
library("DT")
library("coin")
library("ancom.R")


#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()


#set working directory
#setwd("C:/Revised R script_manuscript_050517")
#setwd("C:/R script_Chap.2/ANCOM")

#davide cpu
#setwd("C:/Users/DB42008/Box Sync/Davide_lab_manuscripts/Ridhdhi_OSR_2017/OSR_revision_190717/R script")
#############################################################
#                     ANCOM analysis                        #
#############################################################


#############################################################
# Import OTU table in 'wide format', generated manually in Excel.
OTU_table_transposed <- read.delim("Ridhdhi_OTU_table_transposed.txt")

# ANCOM calculation
ancom_OSR_phyla <- ANCOM(OTU_table_transposed, multcorr = 2, sig = 0.05)

# Generate data frame with output of the calculation
ancom_result_dataframe <- as.data.frame(ancom_OSR_phyla$detected)

plot_ancom(ancom_OSR_phyla)

# Save result of ANCOM analysis as data frame
#write.table(ancom_result_dataframe, file="ANCOM.txt", sep="\t")

# End