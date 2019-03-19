# Name: Steven Cao
# ID: A12083782

#Gene expression and clustering (First 400)

install.packages("gplots")   # Install Package g plots

library("gplots")      # Load library g plots

gene_counts=read.table("t2d_cluster.txt")   # Load gene expression files

#First 400 genes
gene_counts2=gene_counts[1:400,]

# Hierarchal clustering

head(gene_counts2)       # Look at top elements of file
gene_counts3<-as.matrix(gene_counts2)       # Load expression data as matrix


myCol<-colorRampPalette(c("blue","yellow","red"))     #select color Scheme
heatmap.2(gene_counts3,trace="none",scale="row",col=myCol,       
          hclustfun=function(x)hclust(x,method = "complete"),
          cexCol=0.5,cexRow=0.1, distfun=function(x) dist(x,"euclidean"))


# k means clustering

t_gc<-t(gene_counts)
cl<-kmeans(t_gc,2)
cl$cluster

write.table(cl$cluster,"first400kmeans.txt")




#Gene expression and clustering (DEG genes)


gene_counts=read.table("top100genes.txt")   # Load gene expression files


#Top 100 genes
gene_counts2=gene_counts

# Hierarchal clustering

head(gene_counts2)       # Look at top elements of file
gene_counts3<-as.matrix(gene_counts2)       # Load expression data as matrix


myCol<-colorRampPalette(c("blue","yellow","red"))     #select color Scheme
heatmap.2(gene_counts3,trace="none",scale="row",col=myCol,       
          hclustfun=function(x)hclust(x,method = "complete"),
          cexCol=0.5,cexRow=0.1, distfun=function(x) dist(x,"euclidean"))


# k means clustering

t_gc<-t(gene_counts2)
cl<-kmeans(t_gc,2)
cl$cluster

write.table(cl$cluster,"deg_k_means.txt")

