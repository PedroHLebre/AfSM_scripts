#This script was used to calculate the environmental factor with the highest correlation to each dominant phylotypes 
# Firstly, using the rarefied phyloseq objects created for Bacteria, Archaea, and Fungi, we clustered the otu_tables accorting to species level, which for our study we cibsudeed as phylotype. The example script below was performed on the three phyloseq objects seperately
Physeq_bact_phylotype = tax_glom(physeq_Bacteria, taxrank = "Species", NArm = F)

#we then export the ASV and taxonomy tables, which we will then use to manually select the dominant phylotypes, i.e. those that represent the top 10% of reads in 95% or more of samples. 
write.csv(otu_table(Physeq_bact_phylotype), "asv_table_phylotype_bacteria.csv")
write.csv(tax_table(Physeq_bact_phylotype), "tax_table_phylotype_bacteria.csv")

#After manual curation, the dominant phylotypes for Bacteria, Archaea, and ITS were combined in a single table (with rarefied counts)

################################################################################################################################################
#Randon Forest analysis
#After this, we performed random-forest analysis to identify the phylotypes that were affected by environmental variables

library(randomForest)
library(ppcor)
library(rfPermute)

#import a table with the environmental variables and  a column with the abundance (in rarefied counts) of a dominant phylotype (named "response". Due to the output format on the random forest analysis, a new table was loaded for each phylotype. 

metadata_f2 = read.table(file.choose(), header = TRUE, row.names = 1, sep = ",")

#run random forest on the variation in abundance of the phylotype.
set.seed(2)
erie.classify <- randomForest(response~., data = metadata_f2, ntree = 1000,do.trace=50, importance=TRUE)

#These outputs will tell you wheather a phylotype has a specific envirnomental preference (only considered if the percentage of explained variation of the environemtnal variables is above 30%)
print(erie.classify)
varImpPlot(erie.classify)
round(importance(erie.classify), 2)

#The phylotypes that shou environmental preference are subsequently concatenated into a count table and imported again for the semi-partial correlation analysis. This analysis is performed to identify the environmental factor that is most strongly and significantly correlated to the dominant phylotype.

############################################################################################
#Best-subset correlation analysis 

#import data 
#table of environemtal varibles
table1 = read.table(file.choose(), header = TRUE, row.names = 1, sep = ",")
#table(rarefied counts) of phylotypes that were identified as having environmental preferences
table2 = read.table(file.choose(), header = TRUE, row.names = 1, sep = ",")

#construct tables for results
table_e <- as.data.frame(matrix(0,23,ncol(table2))) #the numer of rows for these two tables is equal to the number of environmental variables in table1
colnames(table_e) <- colnames(table2)

table_p <- as.data.frame(matrix(0,23,ncol(table2)))
colnames(table_p) <- colnames(table2)

#perform best-subset correlation with Spearman correlations
for(i in 1:length(colnames(table2))) {
  #i <- 1
  #print((colnames(table1)[i]))
  table2_sub <- table2[,i]
  table3 <- cbind(table1, table2_sub)
  colnames(table3)[23] <- colnames(table2)[i]
  table = spcor(table3, method =c("spearman"))
  
  table_e[,i] <- table$estimate[,23]
  
  table_p[,i] <- table$p.value[,23]
}

row.names(table_e) <- row.names(table$estimate)
row.names(table_p) <- row.names(table$p.value)

#export the tables with the correlation estimates and corresponding p-values for the phylotypes. These tables were then used to manually select the environmental predictor with the highest correlation and lowest p-vale for each phylotype. 
write.csv(table_e, "ITS_pcor_table_estimate.csv")
write.csv(table_p, "ITS_pcor_table_pvalue.csv")
