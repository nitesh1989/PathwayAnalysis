# Nitesh Turaga
# Intersect all TCGA bumps in R

# Feb 28th 2014


# Import Libraries
library(plyr)

# set path
setwd("~/Google Drive/Deliverables/TCGA-450K Data Sets/")


# List files
files = list.files(".",pattern = "11_bumps.csv",recursive =T,full.names =T)

# Make a data frame with all the files
files = data.frame(path = files, length = (sapply( files, function(f) nrow(read.csv(f)))) )

files = files[order(files$length,decreasing =TRUE),]

rownames(files) = gsub("/11_bumps.csv","",files$path)
rownames(files) = gsub("./","",rownames(files))


# Intersect all files
x = read.csv(as.character(files$path[1]), stringsAsFactors = FALSE, header = T)

x = x$name
files$Intersection_size = 0
files$Intersection_size[1] = length(x)
for (i in 2:length(x)){   
    f = read.csv(as.character(files$path[i]),stringsAsFactors =FALSE,header =T)
    x = intersect(x,f$name)
    
    print(length(x))
    files$Intersection_size[i] = length(x)
    if(length(x)<20 && length(x)>10){
        a = x    
        print(x)
    }
    if(length(x)<10 && length(x)>1){
        b = x    
        print(x)
    }
}

# Remove magic variables a, and b
common_genes = a
rm(a)


# Match with all the bumps, from 1 -> 16 in files df
# match for STAD only 6, in row 17

for(i in 1:17) {
    
    bumps = read.csv(as.character(files$path[i]),stringsAsFactors = FALSE,header =T)
    
    common_bumps = bumps[which( bumps$name %in% common_genes),]
    #common_bumps = common_bumps[common_bumps$L>4,]
    filename = sub("11_bumps","common_genes",files$path[i])
    if(length(unique(common_bumps$name)) == 17){
        print(i);
        print("PASS");
    }
    # STAD fails as the number of intersected genes in STAD is only 6
    else{
        print("FAIL")
    }
    write.csv(common_bumps, file = filename)
}

##############################################################################
# Plot the DMRs for the intersected genes which are found to be common 
# in all bumps in TCGA cancer data sets
##############################################################################

# Set path
setwd("~/Google Drive/Deliverables/TCGA-450K Data Sets/")

# List files which have the common_genes
files = list.files(".",pattern = "common_genes.csv",recursive =T,full.names =T)

source("~/Documents/TestRun/HNSCC_Local/FUNCTIONS_2012.R")

# library(charm)
library(minfi)
library(limma)


# Path here is path to file.
path = "~/Google Drive/Deliverables/TCGA-450K Data Sets"

# Plot all DRMS for the common genes

#dmrPlot = function(path) {
    
    files.to.plot = list.files(path,pattern = "common_genes.csv",recursive =T, full.names=T)    
    
    
    for (i in 1:length(files.to.plot)) { 
       
        
        name = files.to.plot[i]
        name = sub("C:/Users/nturaga1/Google Drive/Deliverables/TCGA-450K Data Sets/","",name)
        name = sub("/common_genes.csv","",name)
        name = sub(".+_","",name)
        print(name)
            
        objs_path = "~/Documents/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014"
        objs_data = file.path(objs_path,name)
        objs_data
        
        data.sources = list.files(objs_data,pattern="*.rda",full.names=T)
        data.sources
        sapply(data.sources,load,.GlobalEnv)
        
        T1 = "normal";T2 = "cancer";
        pd=pData(object)
        keep=pd$phenotype%in%c(T1,T2)
        tt=factor(pd$phenotype[keep],c(T1,T2))
        X=model.matrix(~tt)
        design=model.matrix(~tt)
        pos = start(object)
        
        tab = read.csv(files.to.plot[i])        
        if (nrow(tab) == 1){
            print("File unable to print");
            print(name);
            next;
        }
        
        tab = tab[tab$L>4,]
        if(nrow(tab)<1){
            print("File unable to print")
            print(name)
            next;
        }
        
        y=M
        plotName = paste0(strsplit(files.to.plot[i],".csv"),"_dmrPlot.pdf")
        pdf(file=plotName)
        
        for(i in 1:nrow(tab)){
            
            
            par(1,1)
            Index =tab$indexStart[i]:tab$indexEnd[i]
        
            
            x=pos[Index]
            yy=y[Index,]
            
        
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
    
#}







