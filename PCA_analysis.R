# This script was used to plot the PCA of the environmental variables, as well as draw the correlation bubble plot. 
library(Hmisc)
install.packages("corrplot")
library(corrplot)
library(devtools)
library(vegan)
library(ggbiplot)
library(FactoMineR)
library(factoextra)
library(grid)
install_github("vqv/ggbiplot")
library(ggbiplot)

#import the table with environmental variables
table=read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")

#select variable that are continuous
env.bio <- table[,c(1:21)]
env <- table[,23]

#standardized the variables so that they are comparable with each other
chemistry_norm<-decostand(env.bio, "stand") 

#Calculate the pearson correlations between the variables
res_pearson <- rcorr(as.matrix(chemistry_norm),type = c("pearson"))

#Plot the correlation bubble plot
flattenCorrMatrix(res2$r, res2$P)
corrplot(res_pearson$r, type="upper", order="hclust", 
         p.mat = res_pearson$P, sig.level = 0.01, insig = "n")

##########################################################################################################################
#PCA analysis

#plot the PCA variable graph without the sample point
prc = PCA(chemistry_norm)


#Plot PCA of samples without the environmental variable arrows

manual = (values=c("Benin" = "blue", "Botswana" = "red2", "Cote d_Ivoire" ="green", "Kenya" = "skyblue4", "Mozambique" = "cyan", "Namibia" = "darkmagenta", "South Africa" = "deeppink", "Zambia" = "orange", "Zimbabwe" = "palevioletred"))

g <- ggbiplot(prc, groups = table$Country,ellipse=T, var.scale = 1, var.axes = F,circle = F) + theme_bw()
g2 <- g 
g3 <- g2 + scale_color_manual(values = manual)
g5 = g3 + theme_bw()
p6 = g5 + theme(axis.title.x = element_text(size=12, face="bold", margin = margin(t = 20,r = 20, b = 20,l = 0)),axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 20)))
p7 = p6 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + annotation_custom(my_grob)

#add this text after calculatin the sifnificance and R-squared of the sample distribution in the PCA (see scripts bellow)
my_text <- "R2 = 0.6, p-value > 0.01"
my_grob = grid.text(my_text, x=0.3,  y=0.8, gp=gpar(col="firebrick", fontsize=12, fontface="bold"))

#Calculate the significane and R-suared of sample distribution. 
#First calculate PCA distance matric
prC <- prcomp(chemistry_norm, center=TRUE, scale=TRUE) 

#Extrace the values for the first and second PCA axis from the distance matrix (PC1 and PC2)
dtp <- data.frame('Site' = as.character(table$Country),
                  prC$x[,1:2])

#Create a table with the PC1 and PC2 values
PCAcoords <- dtp[,c("PC1","PC2")]

#Calculate the euclidean distance between samples based on the PC1 and PC2 values
PCAdist <- dist(PCAcoords)

#Use PERMANOVA to calculate significance and R2 of the sample PC1 and PC2 distribution according to country. 
adonis(PCAdist  ~ Site, data = dtp, permutations = 1000)

