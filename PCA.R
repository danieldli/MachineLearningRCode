#Generate a fake dataset.
#Generate a matrix of data with 10 samples and 100 genes in each sample.
data.matrix <- matrix(nrow=100, ncol=10)
#Name the samples.
colnames(data.matrix) <- c(paste("wt", 1:5, sep=""),
                           paste("ko", 1:5, sep=""))
rownames(data.matrix) <- paste("gene", 1:100, sep="")
#Give fake genes fake read counts.
for (i in 1:100) {
  wt.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  ko.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  
  data.matrix[i,] <- c(wt.values, ko.values)
}
#Show first few rows of the data.
head(data.matrix)
#Use PCA analysis. Transpose the matrix using t() function.
#scale=TRUE center and scale the data.
pca <- prcomp(t(data.matrix), scale=TRUE)
#prcomp() returns three things: 1) x 2) sdev 3)rotation.
#x contains the PCs for drawing a graph.
plot(pca$x[,1], pca$x[,2])
#sdev "standard deviation" calculates how much variation each PC accounts for.
pca.var <- pca$sdev^2
#Calculate the percentages of each pc acconts for.
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
#Plotting the percentages.
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component",
        ylab="Percent Variation")
#Use gglpot2 to make a fancy PCA plot.
library(ggplot2)
#Format the data the way ggplot2 likes it: SampleID, X, Y.
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) + 
  geom_text() +                                            #Plot labels.
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +     #paste() func to combine
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +      #the percentage and text.
  theme_bw() +                                             #Makes the background white.
  ggtitle("My PCA Graph")                                  #Add a title.
#Loading scores.
loading_scores <- pca$rotation[,1]
genes_scores <- abs(loading_scores)
genes_score_ranked <- sort(genes_scores, decreasing=TRUE)
top_10_genes <- names(genes_score_ranked[1:10])
top_10_genes
##Show the scores (and +/- sign)
pca$rotation[top_10_genes, 1]
pca$rotation[top_10_genes, ]

