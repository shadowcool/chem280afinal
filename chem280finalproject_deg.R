# Script for Chem 280 Finals
# TYPE 2 Diabetes

source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")


require(GEOquery)

#File for differential expression:
#GDS4909_full.soft.gz

library(Biobase)
library(GEOquery)

gset=getGEO(filename='GDS3782_full.soft.gz',GSEMatrix = TRUE)

Meta(gset)$sample_count
Meta(gset)$sample_organism
colnames(Table(gset))            # Get column Name

numbergenes=Meta(gset)$feature_count

#Filter by genes and expression levels 
table=Table(gset)[1:numbergenes,2:30]

#Controls: 
# "GSM524151","GSM524152","GSM524153","GSM524154","GSM524155","GSM524156","GSM524157","GSM524158","GSM524159","GSM524160"
#Cases:
# "GSM524161","GSM524162","GSM524163","GSM524164","GSM524165","GSM524166","GSM524167","GSM524168","GSM524169","GSM524170"

table2=subset(table,select=c("GSM524161","GSM524162","GSM524163","GSM524164","GSM524165","GSM524166","GSM524167","GSM524168","GSM524169","GSM524170","GSM524151","GSM524152","GSM524153","GSM524154","GSM524155","GSM524156","GSM524157","GSM524158","GSM524159","GSM524160","Gene ID"))
  
# Make a list containing old names
list_old=c("GSM524161","GSM524162","GSM524163","GSM524164","GSM524165","GSM524166","GSM524167","GSM524168","GSM524169","GSM524170","GSM524151","GSM524152","GSM524153","GSM524154","GSM524155","GSM524156","GSM524157","GSM524158","GSM524159","GSM524160")

# Make a list containing new names
list_new=c("Case1","Case2","Case3","Case4","Case5","Case6","Case7","Case8","Case9","Case10","Control1","Control2","Control3","Control4","Control5","Control6","Control7","Control8","Control9","Control10")
  
# Update Names, loop through using for loop

for (i in 1:length(list_old)){
  names(table2)[names(table2) == list_old[i]] <- list_new[i]
}

#Make sure to remove genes that have no names
table3=table2[!(is.na(table2$`Gene ID`) | table2$`Gene ID`==""), ]

# Prepare for DEG analysis:

library("gplots")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")



library("DESeq2")

head(table3)     #Look at head of file
dim(table3)      #Look at dimensions of file

#Remove duplicated columns:
table4 = table3[!duplicated(table3$`Gene ID`),]

#Designate row name:
rownames(table4) <- table4$`Gene ID`

#Subset so that only expression data is within the table

table5<-subset(table4,select=c("Case1","Case2","Case3","Case4","Case5","Case6","Case7","Case8","Case9","Case10","Control1","Control2","Control3","Control4","Control5","Control6","Control7","Control8","Control9","Control10"))

#Remove all NA's

table6=table5[complete.cases(table5),]

table7=table6

table7$Case1=as.integer(table6$Case1)
table7$Case2=as.integer(table6$Case2)
table7$Case3=as.integer(table6$Case3)
table7$Case4=as.integer(table6$Case4)
table7$Case5=as.integer(table6$Case5)
table7$Case6=as.integer(table6$Case6)
table7$Case7=as.integer(table6$Case7)
table7$Case8=as.integer(table6$Case8)
table7$Case9=as.integer(table6$Case9)
table7$Case10=as.integer(table6$Case10)

table7$Control1=as.integer(table6$Control1)
table7$Control2=as.integer(table6$Control2)
table7$Control3=as.integer(table6$Control3)
table7$Control4=as.integer(table6$Control4)
table7$Control5=as.integer(table6$Control5)
table7$Control6=as.integer(table6$Control6)
table7$Control7=as.integer(table6$Control7)
table7$Control8=as.integer(table6$Control8)
table7$Control9=as.integer(table6$Control9)
table7$Control10=as.integer(table6$Control10)



gene_counts<-as.matrix(table7)

write.table(gene_counts,"t2d_cluster.txt")


condition=read.table("final_condition.txt")


# Command to start the differential analysis:
dds<-DESeq(DESeqDataSetFromMatrix(countData=gene_counts
                                  ,colData=condition,design=~condition))

deg=results(dds)    #results of analysis
deg

#Significant DEGS
subset(deg,padj<0.1) #Significant degs with pvalues <0.1

subset(deg,padj<0.05) #Significant degs with pvalues <0.05

subset(deg,abs(log2FoldChange)>1) #degs with log 2 fold change >1 

subset(deg,(abs(log2FoldChange)>1) 
       & (padj<0.05) ) # degs with log 2 fold change > 1 and p value < 0.05

# Output results:
write.table(as.data.frame(subset(deg,(abs(log2FoldChange)>1) & (padj<0.05))),"DEG.txt",quote=F, sep="\t")

#Output table:
table8=table7[c("5349","1490","8529","285025","389541","162963","135112","79623","12","1908","140462","9899","5973","1358","5967","563","100507064","307","23643","26034","8908","11082","221395","29969","140701","51365","290","10343","1266","29842","771","10125","91582","9750","6425","131096","643911","5646","7113","10263","27429","9023","1840","101929550","1164","10930","28513","10924","64881","55796","6514","54810","8459","3486","29951","163782","26577","8842","6357","5105","84677","54947","5104","1057","10066","1672","1809","4837","729","116496","23566","57007","100506699","57402","84937","1208","602","486","144501","57393","4852"),]
write.table(table8,"top100genes.txt")
