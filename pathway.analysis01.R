# Nitesh Turaga
# Pathway analysis 01
# Jan 28th 2014

setwd("~/TestRun/Gene_Lists/Gene_Lists_CSV/")
library(plyr)

##########################################
# 01: Process the pathways files
##########################################

processPathways = function(pathwaysFile) { 
    
    a.pathway = read.csv(pathwaysFile)
    
    if(grep("PCR.Array.Catalog...",colnames(a.pathway))) {
        a.pathway$PCR.Array.Catalog... = NULL
    }
    
    gene.symb.index = which(a.pathway[1,] == "Symbol" | colnames(a.pathway) == "Symbol")
    pathway.genes = a.pathway[,gene.symb.index]
    pathway.genes = data.frame(as.character(pathway.genes))
    
    if (pathway.genes[1,1] == "Symbol") {
        pathway.genes = data.frame(Symbol = pathway.genes[-1,])
    }
    pathway.genes = data.frame(name = pathway.genes[which(pathway.genes$Symbol != ""),])
    return(pathway.genes)
    
}


##########################################
# 02 3: Compare with Data set with DMRS
##########################################
# Input : Dataset with a column descriptor which has the gene symbols
# It also requires the path to location of the pathways

# Arguments:
# dataSet = name of bumps file
# dataName = All the intersection files will have this name


intersectPathways = function(dataSet,dataName,length.DMR=3,savePath) {
    
    # Default path to pathways
    pathToPathways = "~/TestRun/Gene_Lists/Gene_Lists_CSV"
    
    pathways = list.files(pathToPathways,pattern = ".csv",full.names=T)
    if(grep(".csv",dataSet)) {
        tab = read.csv(dataSet)
        data = tab
    }
        
    data = data[data$L>length.DMR,]
    
    for (i in 1:length(pathways)) {
        
        genes = processPathways(pathways[i])        
        genes_present = intersect(genes$name,as.character(data$name))
        if (length(genes_present)<1){next}
        
        toWrite = join(data.frame(name = genes_present),data,by = "name",type = "inner",match = "all")                
        #make the file name
        pathwayName = gsub(".+/","",pathways[i])
        fileName = paste0(savePath,strsplit(pathwayName,split=".csv")[1],"_",dataName,".csv")      
        write.csv(toWrite,file = fileName)
        
    }
}


dataSet="~/TestRun/TCGA_Data_Analysis/Lung_Cancer_AD_SC/LUNG_Cancer_SCC/bumps_analysis.rda";dataName="LUSC";length.DMR=3;savePath = "~/TestRun/TCGA_Data_Analysis/Lung_Cancer_AD_SC/"

intersectPathways(dataSet="~/TestRun/TCGA_Data_Analysis/Lung_Cancer_AD_SC/LUNG_Cancer_Ad/bumps_64v32.rda",dataName="LUAD",length.DMR=3,savePath = "~/TestRun/TCGA_Data_Analysis/Lung_Cancer_AD_SC/LUNG_Cancer_Ad/")
intersectPathways(dataSet="~/TestRun/TCGA_Data_Analysis/stomach/STAD_bumps.rda",dataName="STAD",length.DMR=3,savePath = "~/TestRun/TCGA_Data_Analysis/stomach/")



dataSet="~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/BLCA/11_bumps.csv";dataName="BLCA";length.DMR=3;savePath = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/BLCA/"
intersectPathways(dataSet,dataName,length.DMR,savePath)



# Function calls on Feb 21st.

##############################################
# TCGA - LUAD function call
##############################################


dir.create("~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUAD/pathways")
dataSet="~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUAD/11_bumps.csv";dataName="LUAD";length.DMR=3;savePath = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUAD/pathways/"
intersectPathways(dataSet,dataName,length.DMR,savePath)


##############################################
# TCGA - LUSC function call
##############################################


dir.create("~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC/pathways")
dataSet="~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC/11_bumps.csv";dataName="LUSC";length.DMR=3;savePath = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC/pathways/"
intersectPathways(dataSet,dataName,length.DMR,savePath)


##############################################
# Combined LUNG CANCER = union(LUAD + LUSC) function call
##############################################

dir.create("~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/pathways")
dataSet="~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/bumps_combined_LUAD_LUSC.csv";dataName="combined_LU";length.DMR=3;savePath = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/pathways/"
intersectPathways(dataSet,dataName,length.DMR,savePath)


dir.create("~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/pathways")
dataSet="~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/bumps_ind_luad.csv";dataName="ind_LUAD";length.DMR=3;savePath = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/pathways/"
intersectPathways(dataSet,dataName,length.DMR,savePath)

dataSet="~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/bumps_ind_lusc.csv";dataName="ind_LUSC";length.DMR=3;savePath = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/pathways/"
intersectPathways(dataSet,dataName,length.DMR,savePath)

##################################################
# TCGA function calls Feb 24th 2014
##################################################

dataSet="~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/CESC/bumps_CESC.csv";dataName="CESC";length.DMR=3;savePath = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/CESC/pathways/"
intersectPathways(dataSet,dataName,length.DMR,savePath)

dataSet="~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LIHC/11_bumps.csv";dataName="LIHC";length.DMR=3;savePath = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LIHC/pathways/"
intersectPathways(dataSet,dataName,length.DMR,savePath)

#SKCM -- NO BUMPS
dataSet="~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/SKCM/11_bumps.csv";dataName="SKCM";length.DMR=3;savePath = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/SKCM/pathways/"
intersectPathways(dataSet,dataName,length.DMR,savePath)

# TCHA
name = "TCHA"
name = "BRCA"
dataSet=paste0("~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/",name,"/11_bumps.csv");dataName=name;length.DMR=3;savePath = paste0("~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/",name,"/pathways/")
intersectPathways(dataSet,dataName,length.DMR,savePath)


