library("FactoMineR")
library("factoextra")

data <- read.table("./pca.eigenvec",header = F)

#data <- read.table("pca.eigenvec",header = F)
data2 <- as.data.frame(data)
pdf("SNP_pca.pdf")
ggplot(data2, aes(x=V4,y=V5,color=V2)) + geom_point(shape=18,size=3) + scale_shape_manual(values=c(1:8)) + labs(title="PCA",x="PC1",y="PC2") + theme(legend.title = element_blank()) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
dev.off()

