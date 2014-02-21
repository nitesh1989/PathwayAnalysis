######################
# Set my local path

my.path = "~/TestRun/Git_Projects/PathwayAnalysis"
setwd(my.path)

#################
#Source the pathway functions
source("pathway.analysis01.R")
source("pathway.plot02.R")
source("pathway_barplot03.R")
#################


# From source pathway.analysis01.R

intersectPathways(dataSet="~/TestRun/TCGA_Data_Analysis/CESC_DNA_Methylation/bumps_CESC.rda",dataName="CESC",length.DMR=3,savePath = "TestRun/TCGA_Data_Analysis/CESC_DNA_Methylation/")

dataSet="~/TestRun/TCGA_Data_Analysis/CESC_DNA_Methylation/bumps_CESC.rda";dataName="CESC";length.DMR=3;savePath = "TestRun/TCGA_Data_Analysis/CESC_DNA_Methylation/"
