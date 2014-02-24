# Nitesh Turaga
# nturaga1@jhmi.edu, Jan 25th 2014
# Make plots for DMR for pathway analysis
##################################################################################

#setwd("~/Google Drive/Deliverables/Gastric_Cancer_Pathway_Analysis/1_25_2104/pathwayPlots")

#######################
# Load the required files
####################### 
# data.sources = list.files("AllGastritis_vs_Cancer_pathways",pattern="*.rda",full.names=T)
# sapply(data.sources,load,.GlobalEnv)



source("../../HNSCC_Local/FUNCTIONS_2012.R")

#Model Matrix


# library(charm)
library(minfi)
library(limma)




#Path here is path to file.
pathwayPlot = function(path,T1,T2) {
    files.to.plot = list.files(path,pattern = ".csv",full.names=T)    
#     pd=pData(object)
#     keep=pd$Phenotype%in%c(T1,T2)
#     tt=factor(pd$Phenotype[keep],c(T1,T2))
#     X=model.matrix(~tt)
#     design=model.matrix(~tt)
#     pos = start(object)
#     
    
    
    for (i in 1:length(files.to.plot)) { 
#                 i=1
        tab = read.csv(files.to.plot[i])        
        if (nrow(tab) == 1){
            next;
        }
        tab = tab[tab$L>4,]
        if(nrow(tab)<1){
            next
        }
        y=M
        plotName = paste0(strsplit(files.to.plot[i],".csv"),"_dmrPlot.pdf")
        pdf(file=plotName)
        for(i in 1:nrow(tab)){
            #                     i=1    
            par(1,1)
            Index =tab$indexStart[i]:tab$indexEnd[i]
            #                     Index
            x=pos[Index]
            yy=y[Index,]
            
            #                     x
            #                     yy
            matplot(jitter(x),ilogit2(yy),col=as.fumeric(tt),
                    pch=16,cex=0.75, main=paste(tab$region[i], sep=" : ", tab$name[i]),xlab=paste("location on",tab$chr[i]),ylab="Beta")
            
            tmpIndexes=split(1:ncol(yy),tt)
            for(j in seq(along=tmpIndexes)){
                yyy=rowMeans(yy[,tmpIndexes[[j]]])
                lfit=loess(yyy~x,span=0.75)
                lines(x,ilogit2(lfit$fitted),col=j)
            }
            legend("bottomleft",levels(tt),col=1:2,lty=1)
        }
        dev.off()
        
    }
    
}

getwd()


# Function call with Cancer
data.sources = list.files("AllGastritis_vs_Cancer_pathways",pattern="*.rda",full.names=T)
sapply(data.sources,load,.GlobalEnv)
mypath = "AllGastritis_vs_Cancer_pathways"
pathwayPlot(mypath,"gastritis","cancer")

#Function call 2
data.sources = list.files("Gastritis_vs_Cancer_pathways",pattern="*.rda",full.names=T)
sapply(data.sources,load,.GlobalEnv)
mypath = "Gastritis_vs_Cancer_pathways"
pathwayPlot(mypath,"gastritis","cancer")

#Function call 3
data.sources = list.files("MisGastritis_vs_Cancer_pathways",pattern="*.rda",full.names=T)
sapply(data.sources,load,.GlobalEnv)
pd=pData(object)
T1 = "mis_gastritis";T2 = "cancer"
keep=pd$Phenotype%in%c(T1,T2)
tt=factor(pd$Phenotype[keep],c(T1,T2))
X=model.matrix(~tt)
design=model.matrix(~tt)
pos = start(object)

mypath = "MisGastritis_vs_Cancer_pathways"
pathwayPlot(mypath,"mis_gastritis","cancer")

#Function call 4
data.sources = list.files("Gastritis_vs_Misgastritis_pathways",pattern="*.rda",full.names=T)
sapply(data.sources,load,.GlobalEnv)
pd=pData(object)
T1 = "gastritis";T2 = "mis_gastritis"
keep=pd$Phenotype%in%c(T1,T2)
tt=factor(pd$Phenotype[keep],c(T1,T2))
X=model.matrix(~tt)
design=model.matrix(~tt)
pos = start(object)

mypath = "Gastritis_vs_Misgastritis_pathways"
pathwayPlot(mypath,"gastritis","mis_gastritis")


###############################################################################

#Function call 5
data.sources = list.files("~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC/",pattern="*.rda",full.names=T)
sapply(data.sources,load,.GlobalEnv)
#Model Matrix
pd=pData(object)
T1="normal";T2="cancer"
keep=pd$phenotype%in%c(T1,T2)
tt=factor(pd$phenotype[keep],c(T1,T2))
X=model.matrix(~tt)
design=model.matrix(~tt)
pos = start(object)

mypath = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC/pathways"
path = mypath
pathwayPlot(mypath,T1,T2)




#Function call 6
data.sources = list.files("~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUAD/",pattern="*.rda",full.names=T)
sapply(data.sources,load,.GlobalEnv)
#Model Matrix
pd=pData(object)
T1="normal";T2="cancer"
keep=pd$phenotype%in%c(T1,T2)
tt=factor(pd$phenotype[keep],c(T1,T2))
X=model.matrix(~tt)
design=model.matrix(~tt)
pos = start(object)

mypath = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUAD/pathways"
path = mypath
pathwayPlot(mypath,T1,T2)



###############################################################################
# TCGA data set populate
# Feb 24th 2014
###############################################################################


# BLCA

data.sources = list.files("~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/BLCA/",pattern="*.rda",full.names=T)
sapply(data.sources,load,.GlobalEnv)
#Model Matrix
pd=pData(object)
T1="normal";T2="cancer"
keep=pd$phenotype%in%c(T1,T2)
tt=factor(pd$phenotype[keep],c(T1,T2))
X=model.matrix(~tt)
design=model.matrix(~tt)
pos = start(object)

mypath = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/BLCA/pathways"
path = mypath
pathwayPlot(mypath,T1,T2)

rm(object,design,tab,x,M)

# CESC

data.sources = list.files("~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/CESC/",pattern="*.rda",full.names=T)
sapply(data.sources,load,.GlobalEnv)
#Model Matrix
pd=pData(object)
T1="normal";T2="cancer"
keep=pd$phenotype%in%c(T1,T2)
tt=factor(pd$phenotype[keep],c(T1,T2))
X=model.matrix(~tt)
design=model.matrix(~tt)
pos = start(object)

mypath = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/CESC/pathways"
path = mypath
pathwayPlot(mypath,T1,T2)


# Function calls in a loop

name = c("BRCA")
for (i in 1:length(name)){
    
    mypath = file.path("~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014",name)
    
    data.sources = list.files(mypath,pattern="*.rda",full.names=T)
    sapply(data.sources,load,.GlobalEnv)
    #Model Matrix
    pd=pData(object)
    T1="normal";T2="cancer"
    keep=pd$phenotype%in%c(T1,T2)
    tt=factor(pd$phenotype[keep],c(T1,T2))
    X=model.matrix(~tt)
    design=model.matrix(~tt)
    pos = start(object)
    
    mypath = file.path(mypath,"pathways")
    path = mypath
    pathwayPlot(mypath,T1,T2)
    gc()
}
