library(phyloseq)
library(ggplot2)
library(igraph)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(svglite)
library(qiime2R)
library(ggpubr)
library(microbiomeSeq)
devtools::install_github("jbisanz/qiime2R")
install.packages("remotes")
install.packages("adespatial")
 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")
BiocManager::install("preprocessCore")
BiocManager::install("GO.db")

#The first stage would be to use phyloseq to calculate both alpha and beta diversity. This can be replaced by vegan or other packages that do the job.
# First we need to import the tables

ASV_QZA = read_qza(choose.files())
otu_table_p = as.matrix(ASV_QZA$data)
write.csv(otu_table_p, "ITS_table_810_samples.csv")
otu_ITS = read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")

taxonomy_p = read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")
metadata_p = read.table(file.choose(), header = TRUE, row.names = 1, sep = ",")
metadata_p2 = read.table(file.choose(), header = TRUE, row.names = 1, sep = ",")
#Then we convert these tables into a phyloseq object (the tree object is not essential if we are not using unifrac)
otu_table_matrix = as.matrix(otu_ITS)
taxonomy_matrix = as.matrix (taxonomy_p)
AfSM_otu = otu_table (otu_table_matrix, taxa_are_rows = TRUE)
AfSM_tax = tax_table(taxonomy_matrix)
AfSM_sampledata = sample_data(metadata_p)

physeq = phyloseq(AfSM_otu, AfSM_tax, AfSM_sampledata)

#if we want to filter out any OTUS with less than 10 reads(I tend to do this manually)
physeq_cutoff = filter_taxa(physeq, function(x) sum(x) > 10, TRUE)

# I normally delete all the samples that do not show adequate depth and then run the following code - however, we can also replace the sample.size term with a hard number
physeq_ra = rarefy_even_depth(physeq_cutoff, rngseed = 1, sample.size = 11400, replace = F)



###################################################
#Alpha-diversity calculations  
#I don't actually agree with rarefying the data before calculating alpha-diversity, as it will mask the "true" differences between samples. However, I understand that this is the most accepted methodology. 
alpha_diversity = estimate_richness(physeq_ra)
#convert to a table if not already
alpha_diversity = as.matrix(alpha_diversity)
#At this point, I export the table and merge it manually with the groupings of the samples (countries, aridity, etc...). I more clever person could probably figure out a way of doing it in R without having to export/import again
write.csv(alpha_diversity, "alpha_diversity_ITS.csv")
# And then import the merged table 
table_alpha = read.table (file.choose(), header = TRUE, row.names = 1 , sep = ",")

#I normally choose a couple of metrics to analyse rather than all alpha-diversity. Observed,Shannon, and InvSimpson always make sense to me. First I test for normal distribution of that metric.
shapiro.test(table_alpha$InvSimpson)

#If it is normally distributed I wiil calculate significance using anova and Tukey. The term site can be substituted for whatever grouping we want to test (aridity, country, etc)
aov.observed.b = aov(Observed ~ Site, table_alpha)
summary(aov.observed.b)
TukeyHSD(aov.observed.b)
#If not normally distributed
kruskal.test(InvSimpson ~ table_alpha$pH_general ,data=table_alpha)
pairwise.wilcox.test(table_alpha$InvSimpson, table_alpha$pH_general, p.adjust.method="fdr")

#Afterward, I plot the significant results as a boxplot and draw the significance brackets manually after exporting the graph as a metafile
observed_compare = ggplot (table_alpha, aes( x= Country, y = InvSimpson, fill = Country)) + geom_boxplot() +theme_classic() + scale_fill_manual(values = c("blue", "red2", "green", "skyblue4", "cyan", "gold", "deeppink", "orange",  "darkgreen"))
observed_compare1 = ggplot (table_alpha, aes( x= Aridity_class_2000, y = Shannon, fill = Aridity_class_2000)) + geom_boxplot() +theme_classic() 
observed_compare2 = ggplot (table_alpha, aes( x= pH_general, y = InvSimpson, fill = pH_general)) + geom_boxplot() +theme_classic() 

###################################################
#Beta-diversity calculations 
#The first step in this process is to transform the rarefied data in order to have more normal distribution. I use the log function, but the best method will depend on the data. 
physeq_transformed = transform_sample_counts(physeq, function (x) log(x+1))

write.csv(otu_table(physeq_transformed), "otu_ITS_transformed.csv")


erie_PCoA= phyloseq::distance(physeq_transformed, method = "bray")



sampledf <- data.frame(sample_data(physeq_transformed))

#we then calculate the significance of the dissimilarity with PERMANOVA
adonis(erie_PCoA  ~ sampledf$Country, data = sampledf, permutations = 1000)
#In addition, I also calculate the beta-dispersivity, which tells me if the groups are highly variable within each other. If this is significant, we cannot trust the PERMANOVA results. 
beta <- betadisper(erie_PCoA, sampledf$Aridity_class_2000)
permutest(beta)

manual = (values=c("Benin" = "blue", "Botswana" = "red2", "Cote d_Ivoire" ="green", "Kenya" = "skyblue4", "Mozambique" = "cyan", "Namibia" = "darkmagenta", "South Africa" = "deeppink", "Zambia" = "orange", "Zimbabwe" = "palevioletred"))


erie_PCoA_ordination <- ordinate(physeq = physeq_transformed, method = "PCoA", distance = erie_PCoA_bact)

PCoA_country = plot_ordination(physeq = physeq_transformed, ordination = erie_PCoA_ordination , color = "Country", axes = 1:2)  +scale_color_manual(values = manual) + theme_bw()

PCoA2_country = PCoA_country  + geom_point(size = 3) + theme(axis.text.x = element_text(face="bold", size=10),axis.text.y = element_text(face="bold", size=10), legend.title=element_text(size=16), legend.text=element_text(size=14))
PCoA3_country  = PCoA2_country + theme(axis.title.x = element_text(size=12, face="bold", margin = margin(t = 20,r = 20, b = 20,l = 0)),axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 20)))
PCoA4_country  = PCoA3_country +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +stat_ellipse()
PCoA4_country

#export the resulting figure as a metagile or svg. 


###################################################
#RDA analysis
#This script was taken from Surendra's analysis.
## Which variable(s) do best explain the variation in prokaryotic distruibution 

col.gr <- c("blue",  "red2", "green", "skyblue4", "cyan", "darkmagenta", "deeppink", "orange", "palevioletred")


#import transformed otu table from phyloseq
otu_table_b_rda=otu_table(physeq_transformed)

#import environmental table with country groupings on the last column
env_table=read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")

#table with just continuous variables
env.bio <- env_table[,c(1:21)]
#standardize table
env.z.L5 <- decostand(env.bio, method="standardize")

#check if CCA or RDA should be used (for Bacteria the result indicated  CCA)
DCA = decorana(t(otu_table_b_rda), iweigh = 1, ira = 0, iresc = 100)
DCA 

#First model
spe.rda <- cca(t(otu_table_b_rda)~., data=env.z.L5)
summary(spe.rda)
#anova of model
anova.cca(spe.rda, step=1000)

#to check for outliers. Any dots outside of the standard deviations for CCA1 or CCA2 should be outliers. 
Boxplot(spe.rda$CCA$wa)
#You can get the values for these as a table. then you can manually check the values for any outliers to remove from the otu table before rerunning the rda script. 
CCA_table =spe.rda$CCA$wa

#ordistep with upper and lower levels
ordiR2step(cca(t(otu_table_b_rda)~1, data=env.z.L5), scope= formula(spe.rda), direction= "forward", R2scope=TRUE, pstep=1000)

#subset the env table with just the significant factors from the ordistep results 
data.env.subs.L5 <- env.z.L5[,c("Longitude", "Latitude",  "Isothermality.", "Mean_Ann_Temp", "Precip_Driest_Month", "Mean_Ann_Precip",  "Precip_Seasonality")]

#final model
spe.rda.signif <- cca(t(otu_table_b_rda)~., data=data.env.subs.L5)

summary(spe.rda.signif,  display=NULL)

#check if variable are co-linear (in this case they weren't). If they are, we would run the model again after removing the co-linear variables with the highest VIF values, in a step-wise procedure. 
vif.cca(spe.rda.signif)

#calculated adjusted R-square
R2adj <- RsquareAdj(spe.rda.signif)$adj.r.squared

#Anova of final model
anova.cca(spe.rda.signif, step=1000)

#ploting the RDA
site <- factor(env_table$Country)

gr.use <- factor(site)

plot(spe.rda.signif, scaling=1, cex.lab =1.3, cex.axis=1.3, 
     font.lab = 1, display = "sites", 
     xlim=c(-3,3), ylim=c(-3,3),
     type="n", xlab=c("CCA1 (18.7%)"), ylab=c("CCA2 (18.4%)")) # these are the explanatory percentages for CCA1 and CCA2 extracted from summary(spe.rda.signif,  display=NULL))
points(scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=1),
       pch = 21, col = "black", bg = col.gr[gr.use], cex=1, lwd = 1.3)
legend(0.6, 0.7, legend = levels(gr.use), col= "black", pt.bg = col.gr, pch=c(21), cex=0.75, pt.cex = 1.5)

text(spe.rda.signif, display = "bp", col="black", cex=1)

#partition variation analysis using the same data, after calculating the explanatory variables with the ordistep modelling.Variables should be clustered into 2 to 4 groups
mod=varpart(t(otu_table_b_rda),   ~ Isothermality. + Mean_Ann_Temp + Precip_Driest_Month + Mean_Ann_Precip + Precip_Seasonality , ~ Longitude + Latitude , data=data.env.subs.L5,scale = FALSE)
#plot as venn diagram
plot(mod, cutoff = -Inf, cex = 0.7, bg=2:5)
# We then exported the plot as a metafile or svg for further formatting. 



