library(phyloseq)
library(ggplot2)
library(igraph)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(svglite)
library(ggordiplots)

`#The first stage would be to use phyloseq to calculate both alpha and beta diversity. This can be replaced by vegan or other packages that do the job.
# First we need to import the tables
otu_table_p=read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")
taxonomy_p = read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")
metadata_p = read.table(file.choose(), header = TRUE, row.names = 1, sep = ",")
metadata_p2 = read.table(file.choose(), header = TRUE, row.names = 1, sep = ",")
metadata_f2 = read.table(file.choose(), header = TRUE, row.names = 1, sep = ",")
#Then we convert these tables into a phyloseq object (the tree object is not essential if we are not using unifrac)
otu_table_matrix = as.matrix(otu_table_p)
taxonomy_matrix = as.matrix (taxonomy_p)
AfSM_otu = otu_table (otu_table_matrix, taxa_are_rows = TRUE)
AfSM_tax = tax_table(taxonomy_matrix)
AfSM_sampledata = sample_data(metadata_p)

physeq = phyloseq(AfSM_otu, AfSM_sampledata)

#for prokaryotic data, I do an initial filtering step to remove any mislabelled taxa. Additionally, we can also implement a filtering step to remove any OTUs that are less than a specific number of counts. 
physeq_Prokaryotes = subset_taxa( physeq,tax_table(physeq) != "Eukaryota")

physeq_species = tax_glom(physeq_Prokaryotes, taxrank = "Phylum", NArm = F)
physeq_Bacteria = subset_taxa (physeq_species, Domain == "Bacteria") 
physeq_ACE_bac_RA = transform_sample_counts(physeq_Bacteria, function (x) x/sum(x))

write.csv(otu_table(physeq_species), "otu_table_S16_class.csv")
write.csv(otu_table(physeq_ACE_bac_RA), "otu_table_bacteria_phylotype_RA.csv")
write.csv(tax_table(physeq_species), "tax_table_S16_class.csv")
write.csv(sample_data(physeq_Archea_phylotype), "metadata_table_archaea_phylotype.csv")


physeq_bact_filter = filter_taxa(physeq_Bacteria, function(x) count(x) > (0.5*length(x)), TRUE)

#if we want to filter out any OTUS with less than 10 reads(I tend to do this manually)
physeq_cutoff = filter_taxa(physeq_Prokaryotes, function(x) sum(x) > 10, TRUE)

#The next step would be to rarefy the data to an even depth
rarecurve(t(otu_table(physeq_cutoff)), step=50, cex=0.5)

# I normally delete all the samples that do not show adequate depth and then run the following code - however, we can also replace the sample.size term with a hard number. 
physeq_ra = rarefy_even_depth(physeq_cutoff, rngseed = 1, sample.size =27000, replace = F)

#subdividing the phyloseq into different Domains
physeq_Archea = subset_taxa(physeq_ra, Domain == "Archaea")
physeq_Bacteria = subset_taxa (physeq_ra, Domain == "Bacteria")

#we needed to prune the Archaeal phyloseq object again to remove ASVs that have 0 counts due to rarefraction step
physeq_Archea_pruned = prune_samples(sample_sums(physeq_Archea)>=1, physeq_Archea)


###################################################
#Alpha-diversity calculations (example script that was used for both bacterial and archaeal phyloseq objects) 
#I don't actually agree with rarefying the data before calculating alpha-diversity, as it will mask the "true" differences between samples. However, I understand that this is the most accepted methodology. 
alpha_diversity = estimate_richness(physeq_Bacteria)
#convert to a table if not already
alpha_diversity_matrix = as.matrix(alpha_diversity)
#At this point, I export the table and merge it manually with the groupings of the samples (countries, aridity, etc...). I more clever person could probably figure out a way of doing it in R without having to export/import again
write.csv(alpha_diversity_matrix, "alpha_diversity.csv")
# And then import the merged table 
table_alpha = read.table (file.choose(), header = TRUE, row.names = 1 , sep = ",")

#I normally choose a couple of metrics to analyse rather than all alpha-diversity. Observed,Shannon, and InvSimpson always make sense to me. First I test for normal distribution of that metric.
shapiro.test(table_alpha$Observed)

#If it is normally distributed I wiil calculate significance using anova and Tukey. The term site can be substituted for whatever grouping we want to test (aridity, country, etc)
aov.observed.b = aov(Observed ~ Site, table_alpha)
summary(aov.observed.b)
TukeyHSD(aov.observed.b)
#If not normally distributed
kruskal.test(Observed ~ Site, data=table_alpha)
pairwise.wilcox.test(table_alpha$Observed, table_alpha$Site, p.adjust.method="fdr")

#Afterward, I plot the significant results as a boxplot and draw the significance brackets manually after exporting the graph as a metafile
manual = (values=c("Benin" = "blue", "Botswana" = "red2", "Cote d_Ivoire" ="green", "Kenya" = "skyblue4", "Mozambique" = "cyan", "Namibia" = "darkmagenta", "South Africa" = "deeppink", "Zambia" = "orange", "Zimbabwe" = "palevioletred"))

observed_compare = ggplot (table_f_alpha, aes( x= Site, y = Observed, fill = Site)) + geom_boxplot() + geom_point (aes(fill = Site), size = 3, position = position_jitterdodge())+theme_classic()
observed_compare2 = observed_compare + scale_fill_manual(values = manual)

###################################################
#Beta-diversity calculations  (example script, was used for both bacterial and archaeal phyloseqs)
#The first step in this process is to transform the rarefied data in order to have more normal distribution. I use the log function, but the best method will depend on the data. 
physeq_transformed_Bact = transform_sample_counts(physeq_Bacteria, function (x) log(x+1))

#Afterward, the beta-diversity is calculated using the bray-curtis method. This can be changed depending on the metric we use, and we should agree on a specific method to use. 
erie_PCoA_bact = phyloseq::distance(physeq_transformed_Bact, method = "bray")

sampledf_bact = data.frame(sample_data(physeq_transformed_Bact))

#we then calculate the significance of the dissimilarity with PERMANOVA
adonis(erie_PCoA_bact  ~ sampledf_bact$Country, data = sampledf_bact, permutations = 1000)

#In addition, I also calculate the beta-dispersivity, which tells me if the groups are highly variable within each other. If this is significant, we cannot trust the PERMANOVA results. 
beta <- betadisper(erie_PCoA_bact , sampledf_bact$Country)
permutest(beta)


erie_PCoA_ordination <- ordinate(physeq = physeq_transformed_Bact, method = "PCoA", distance = erie_PCoA_bact)


PCoA_country = plot_ordination(physeq = physeq_transformed_Bact, ordination = erie_PCoA_ordination , color = "Country", axes = 1:2)  +scale_color_manual(values = manual) + theme_bw()

PCoA2_country = PCoA_country  + geom_point(size = 3) + theme(axis.text.x = element_text(face="bold", size=10),axis.text.y = element_text(face="bold", size=10), legend.title=element_text(size=16), legend.text=element_text(size=14))
PCoA3_country  = PCoA2_country + theme(axis.title.x = element_text(size=12, face="bold", margin = margin(t = 20,r = 20, b = 20,l = 0)),axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 20)))
PCoA4_country  = PCoA3_country +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +stat_ellipse()
PCoA4_country

#export the resulting figure as a metagile or svg. 


###################################################
#RDA analysis (ran for both Bacteria and Archaea separately)
#This script was taken from Surendra's analysis.
## Which variable(s) do best explain the variation in prokaryotic distruibution 

col.gr <- c("blue",  "red2", "green", "skyblue4", "cyan", "darkmagenta", "deeppink", "orange", "palevioletred")


#import transformed otu table from phyloseq
otu_table_b_rda=otu_table(physeq_transformed_Bact)

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
data.env.subs.L5 <- env.z.L5[,c("Longitude", "Latitude",  "Isothermality.", "Mean_Ann_Temp", "Precip_Driest_Month", "Mean_Ann_Precip",  "Precip_Seasonality", "EVI2", "pH", "Ca", "N_Percent", "C_Percent", "Na", "K", "Mg", "Sand_Percent")]

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
     type="n", xlab=c("CCA1 (14.5%)"), ylab=c("CCA2 (10.1%)")) # these are the explanatory percentages for CCA1 and CCA2 extracted from summary(spe.rda.signif,  display=NULL))
points(scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=1),
       pch = 21, col = "black", bg = col.gr[gr.use], cex=1, lwd = 1.3)
legend(0.6, 0.7, legend = levels(gr.use), col= "black", pt.bg = col.gr, pch=c(21), cex=0.75, pt.cex = 1.5)

text(spe.rda.signif, display = "bp", col="black", cex=1)

#partition variation analysis using the same data, after calculating the explanatory variables with the ordistep modelling.Variables should be clustered into 2 to 4 groups, for instance, for the bacterial and archaeal data, environmental data was clustered into soil chemistry, soil texture, climate, and distance
mod=varpart(t(otu_table_b_rda),  ~ pH + Ca + N_Percent + C_Percent  + Na + K + Mg  , ~ Sand_Percent + EVI2 , ~ Isothermality. + Mean_Ann_Temp + Precip_Driest_Month + Mean_Ann_Precip + Precip_Seasonality , ~ Longitude + Latitude , data=data.env.subs.L5,scale = FALSE)
#plot as venn diagram
plot(mod, cutoff = -Inf, cex = 0.7, bg=2:5)
# We then exported the plot as a metafile or svg for further formatting. 


