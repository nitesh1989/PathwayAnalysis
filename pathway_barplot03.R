# Nitesh Turaga
# Make bar plots with frequency of pathways.

# Workflow: 
# 1. Collect pathways
# 2. unlist frequencies
# 3. make barplots
# Commit changes and push

# Barplot for pathways in Gastric cancer.

library(reshape)
library(ggplot2)


################################################################################
# Collect pathways in a format for barplots
################################################################################
collect.pathways = function(status,path) {
    files = list.files(path,full.names=T,pattern = ".csv")
    new_df = data.frame()
    for (i in 1:length(files)) {     
        
        x = read.csv(files[i])
        pathwayName = gsub(pattern = paste0("_",status,".csv"),replacement="",x=(gsub(".+/","",files[i])))
        # CHANGE HERE TO 7:22 for original TSS peak files
        x = x[,c(4:7,1:3,7:20)]
        knock_df = data.frame(x,pathway = pathwayName,pathwayGene = paste0(x$name,"-",pathwayName))
        new_df = rbind(new_df,knock_df)
    }
    if("L" %in% colnames(new_df)){
    new_df = new_df[new_df$L>3,]
    }
    return(new_df)
}


################################################################################
# Unlist the frequency of tables
################################################################################
unlist.frequency = function(df) {
    x = which(df$Freq >= 1)
    n.df = data.frame()
    for (i in x){
        i
        freq = as.numeric(df$Freq[i])
        c.df = data.frame(Pathway = rep(df$Pathway[i],times=freq),Freq = rep(1,times = freq),Gene = rep(df$Gene[i],times = freq))
        n.df = rbind(n.df,c.df)
    }
    df = df[-(x),]
    df = rbind(df,n.df)
    return(df)
}


################################################################################
# Make bar plots for pathways
################################################################################
make.plot = function(collected.pathway,plotTitle){
    data = as.data.frame(table(collected.pathway$name,collected.pathway$pathway))
    data = cast(data=data,Var2~Var1)
    dat.lng = melt(data=data,value.name=pathways)
    dat.lng = dat.lng[dat.lng$value!=0,]
    
    row.names(dat.lng) = NULL
    colnames(dat.lng) = c("Pathway","Freq","Gene")
    
    # Custom function to unlist the frequency tables
    dat.lng = unlist.frequency(dat.lng)
       
    p = ggplot(dat.lng,aes(Pathway,fill =Gene)) + geom_bar(position = "stack",color = "black",width=.5) + coord_flip()
    p + ggtitle(plotTitle) +  guides(fill=guide_legend(ncol=4)) + scale_fill_discrete("Genes")
}

################################################################################



################################################################################
# GASTRIC COMPARISIONS
################################################################################

#Test Gastric comparisions
setwd("~/TestRun/Gastric_Comparision/Gastrtic_local_comparisions_Jan242014/pathwayPlots/")
# Function CALL
allG.c = collect.pathways(status = "AllGastric_vs_Cancer_pathways","AllGastritis_vs_Cancer_pathways")[,c("name","description","pathway")]
gas.c = collect.pathways(status="Gastritis_vs_Cancer",path = "Gastritis_vs_Cancer_pathways")[,c("name","description","pathway")]
mis.c = collect.pathways(status="MisGastric__vs_Cancer",path = "MisGastritis_vs_Cancer_pathways")[,c("name","description","pathway")]
gas.mis = collect.pathways(status="Gastritis_vs_Misgastritis",path = "Gastritis_vs_Misgastritis_pathways")[,c("name","description","pathway")]


# CHANGE PATHWAY NAME HERE

local_gastric = collect.pathways(status = "Gastritis_vs_Cancer_pathways",path = "~/TestRun/TCGA_Data_Analysis/STAD/Gastritis_vs_cancer_local_pathways")[,c("name","description","pathway")]
make.plot(local_gastric,"") + ggtitle(expression(atop("Peru Gastric Cancer Samples", atop(italic(" Pathway analysis"), ""))))
ggsave(file = "~/TestRun/TCGA_Data_Analysis/STAD/Gastritis_vs_cancer_local_pathways/peru_gastric_cancer_pathways.png",width = 20,height = 8)


make.plot(allG.c,"All_Gastritis_vs_Cancer")
ggsave(file = "Allgastritis_vs_cancer_pathways.png",width = 20,height = 8)
write.csv(allG.c,file = "allgastritis_vs_cancer_pathway_table.csv")

collected.pathway = gas.c;plotTitle = "test"
make.plot(gas.c,"Gastritis_vs_Cancer")
ggsave(file = "gastritis_vs_cancer_pathways.png",width = 20,height = 8)
write.csv(gas.c,file = "gastritis_vs_cancer_pathway_table.csv")

make.plot(gas.mis,"Gastritis_vs_misgastritis")
ggsave(file = "gastritis_vs_misgastritis_pathways.png",width = 20,height = 8)
write.csv(gas.mis,file = "gastritis_vs_misgastritis_pathway_table.csv")

make.plot(mis.c,"misgastritis_vs_cancer")
ggsave(file = "misgastritis_vs_cancer_pathways.png",width = 20,height = 8)
write.csv(mis.c,file = "misgastritis_vs_cancer_pathway_table.csv")


################################################################################
# TCGA Data sets
################################################################################

#Function Call
cesc = collect.pathways(status = "CESC",path = "~/TestRun/TCGA_Data_Analysis/CESC_DNA_Methylation/pathways")[,c("name","description","pathway")]
luad = collect.pathways(status = "LUAD",path = "~/TestRun/TCGA_Data_Analysis/Lung_Cancer_AD_SC/LUNG_Cancer_Ad/pathways")[,c("name","description","pathway")]
lusc = collect.pathways(status = "LUSC",path = "~/TestRun/TCGA_Data_Analysis/Lung_Cancer_AD_SC/LUNG_Cancer_SCC/pathways")[,c("name","description","pathway")]


make.plot(luad,"LUAD")
make.plot(lusc,"LUSC")
make.plot(cesc,"CESC")

# FULL PIPELINE FUNCTION CALL
stad = collect.pathways(status = "STAD",path = "~/TestRun/TCGA_Data_Analysis/stomach/pathways")
make.plot(stad,"Stomach Adenocarcinoma")
ggsave(file = "stad_pathways.png",width = 20,height = 8)




################################################################################
# Biswal pathways before TSS change
################################################################################

knockin = collect.pathways(status = "H460_knock",path = "~/TestRun/ShyamBiswal/circos_plot_test/Pathways_knock/")[,c("name","description","pathway")]
make.plot(knockin,"H460 knock in")


parent = collect.pathways(status = "H460_parent",path = "~/TestRun/ShyamBiswal/circos_plot_test/Pathways_parent")[,c("name","description","pathway")]
make.plot(parent,"H460 Parent")

################################################################################
# Biswal pathways after TSS change.
################################################################################

knockin_change = collect.pathways(status = "H460_knock",path = "~/TestRun/ShyamBiswal/circos_plot_test/Pathways_knock_TSS_change/")[,c("name","description","pathway")]
make.plot(collected.pathway=knockin_change,plotTitle="H460 knock-in pathway-gene frequency barplot")
ggsave(file = "~/TestRun/ShyamBiswal/circos_plot_test/collected_pathways/h460knockin_pathways.png",width = 25,height = 12)


parent_change = collect.pathways(status = "H460_parent",path = "~/TestRun/ShyamBiswal/circos_plot_test/Pathways_parent_TSS_change/")[,c("name","description","pathway")]
make.plot(parent_change,"H460 Parent pathway-gene frequency barplot")
ggsave(file = "~/TestRun/ShyamBiswal/circos_plot_test/collected_pathways/h460parent_pathways.png",width = 25,height = 12)



# Function call on feb 21st
################################################################################
# LUSC bar plot AFTER Jan30 2014.
################################################################################

lusc = collect.pathways(status = "LUSC",path = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSc/pathways")
lusc = lusc[lusc$L>10,]
lusc = lusc[,c("name","description","pathway")]
make.plot(lusc,"Lung Cancer Squamous Cell Carcinoma pathway-gene frequency barplot")
ggsave(file = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC/LUSC_pathwayBarplot_L_over_10.png",width = 25,height = 15)


################################################################################
# LUAD bar plot AFTER Jan30 2014.
################################################################################

luad = collect.pathways(status = "LUAD",path = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUAD/pathways")
luad = luad[luad$L>10,]
luad = luad[,c("name","description","pathway")]
make.plot(luad,"Lung Cancer Adenocarcinoma pathway-gene frequency barplot")
ggsave(file = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUAD/LUAD_pathwayBarplot_L_over_10.png",width = 25,height = 15)

################################################################################
# LUAD + LUSC bar plot AFTER Jan30 2014.
################################################################################

comb = collect.pathways(status = "combined_LU",path = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/pathways")
comb = comb[comb$L>10,]
comb = comb[,c("name","description","pathway")]
make.plot(comb,"Lung Cancer Adenocarcinoma union Lung Cancer Squamous pathway-gene frequency barplot")
ggsave(file = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/LUSC_LUAD_pathwayBarplot.png",width = 25,height = 15)


ind_luad = collect.pathways(status = "ind_LUAD",path = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/pathways/ind_LUAD")
ind_luad = (ind_luad[ind_luad$L>10,])[,c("name","description","pathway")]
make.plot(ind_luad,"Independent LUAD pathway-gene frequency barplot")
ggsave(file = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/ind_LUAD.png",width = 25,height = 15)


ind_lusc = collect.pathways(status = "ind_LUSC",path = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/pathways/ind_LUSC")
ind_lusc = (ind_lusc[ind_lusc$L>10,])[,c("name","description","pathway")]
make.plot(ind_lusc,"Independent LUSC pathway-gene frequency barplot")
ggsave(file = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/LUSC_LUAD_Combined/ind_LUSC.png",width = 25,height = 15)


# pathway plots for TCGA datasets on Feb 24th 2014

x = collect.pathways(status = "BLCA",path = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/BLCA/pathways/")
x = (x[x$L>10,])[,c("name","description","pathway")]
make.plot(x,"BLCA pathway-gene frequency barplot")
ggsave(file = "~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/BLCA/BLCA_barplot.png",width = 25,height = 15)


name = c("BRCA")
for(i in 1:length(name)){
    x = collect.pathways(status = name[i],path = file.path("~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014",name[i],"pathways"))
    #x = (x[x$L>10,])[,c("name","description","pathway")]
    x = x[,c("name","description","pathway")]
    make.plot(x,paste(name[i]," pathway-gene frequency barplot"))
    ggsave(file = paste0("~/TestRun/TCGA_Data_Analysis/TCGA-Analysed-After-Jan302014/"
                         ,name[i],"/",name[i],"_barplot.png"),width = 25,height = 15)
    
}