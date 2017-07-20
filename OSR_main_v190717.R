#############################################################
#
# Ref to the ARTICLE
# 
# Ridhdhi Rathore
# Ridhdhi.Rathore@itcarlow.ie
#
# Davide Bulgarelli
# d.bulgarelli@dundee.ac.uk
# 
# script to reproduce calculations and figures presented in the manuscript
# 
#
# Revision July 2017
# Disclaimer: the manuscript is currently submitted, the script might be subjected to changes derived from the revision process
#
#############################################################

#############################################################
# Libraries and functions required
#############################################################

#These initial commands are required to clean-up the memory and start a new session
rm(list=ls())
dev.off()


#set working directory 
#setwd("C:/Revised R script_manuscript_050517/")

#davide cpu
#setwd("C:/Users/DB42008/Box Sync/Davide_lab_manuscripts/Ridhdhi_OSR_2017/OSR_revision_190717/R script")


#1st time Phyloseq installation
#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
#biocLite("DESeq2")
#biocLite("VennDiagram")
#biocLite("PMCMR")

#load the required packages (all included in Phyloseq, but is necessary to invoke them)
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("vegan")
library ("ape")
library("VennDiagram")
library("PMCMR")

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()
#############################################################

#############################################################
#import the count matrix and the desing file

#OTU table this file has been generated using USEARCH v8 and QIIME 1.9.0. In the OTU ids, OTU abundance information has been removed
dat_info <- read.delim("Ridhdhi_OTU_table_R2.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#inspect the file using the command dim
dim(dat_info)

#extract the total number of reads clustered at OTU 97% identiy (the number underneath the Bn identifier represents the total number of reads clustered for that sample) 
OTU_97_reads <- sort(colSums(dat_info[, 1:24]))
OTU_97_reads
#total reads
OTU_97_reads_sum <- sum(colSums(dat_info[, 1:24]))
OTU_97_reads_sum

#design file
design <- read.delim("Ridhdhi_Mapping_file_r090616.txt", sep = "\t", header=TRUE, row.names=1)
design

#remove chloroplast e mitochondria OTUs from the original dataset
Chloroplast <- dat_info[grepl("Chloroplast", dat_info$taxonomy), ]
dim(Chloroplast)
Chloroplast[1:10, ]

mitochondria <- dat_info[grepl("mitochondria", dat_info$taxonomy), ]
dim(mitochondria)
mitochondria[1:3, ]

#Filter plant-derived OTUs from the OTU table
noPlants <- setdiff(rownames(dat_info), c(rownames(Chloroplast), rownames(mitochondria)))
#inspect the results
length(rownames(dat_info))
length(noPlants)

#save this piece of information for qiime purposes. We need to create a a new OTU table in QIIME to generate the taxa tables
#write(noPlants, "Ridhdhi_noPlant_OTUs_id.txt")

#Generate a new OTU table depleted with chloroplast and mitochondria OTUs
dat_info_noPlants <- dat_info[noPlants, ]

#create a new count matrix without OTUs assigned to Choloplast and Mitochondria
dat_count <- dat_info[, rownames(design)]
dat_count_noplants <- dat_info[noPlants, rownames(design)]
dim(dat_count_noplants)

#and a new taxa table
dat_tax_noPlants <- as.data.frame(dat_info[rownames(dat_count_noplants), 25])
rownames(dat_tax_noPlants) <- rownames(dat_count_noplants)
#Save the above file and in excel generate a tax table where column represents a taxonomic rank
#write.table(dat_tax_noPlants, file="Ridhdhi_dat_tax_noPlants.txt", sep="\t")

#check the effect of mitochondria/chloroplast depletion on the new OTU table
#with plant-derived sequences
dim(dat_count)

#w/o plant-derived sequences
dim(dat_count_noplants)

#total number of reads w/o plant sequences
OTU_97_reads_noPlants <- colSums(dat_count_noplants)
OTU_97_reads_noPlants

#now sort per sample
sort(OTU_97_reads_noPlants)

#total number of reads
OTU_97_reads_noPlants_sum <- sum(OTU_97_reads_noPlants)
OTU_97_reads_noPlants_sum 

#now determine the proportion of non-plant reads in the original dataset 
useful_reads <- (OTU_97_reads_noPlants_sum/OTU_97_reads_sum)*100
useful_reads

#create a dataset to visualise the proportion of reads per sample before and after removing OTUs assigned to chloroplast and mitochondria
OTU_97_reads_noPlants <- as.data.frame(OTU_97_reads_noPlants)
OTU_97_microbial_reads_proportion <- as.data.frame(colSums(dat_count_noplants)/colSums(dat_count))*100

#rename the columns in the generated datasets
colnames(OTU_97_reads_noPlants) <- c("reads")
colnames(OTU_97_microbial_reads_proportion) <- c("Microbial_OTUs_reads")

#combine these datasets with the design file
design_info <- cbind(design, OTU_97_reads_noPlants)
design_info_2 <- cbind(design_info, OTU_97_microbial_reads_proportion)

##########
#Figure S1
##########

#re-order the factors
design_info_2$Compartments <- ordered(design_info_2$Compartments, levels=c("Bulk", "Rhizosphere", "Root", "Shoot"))
design_info_2$Tillage <- ordered(design_info_2$Tillage, levels=c("Conventional", "Conservation"))                                                                                                 

#Figure S1A
with(design_info_2, boxplot(Microbial_OTUs_reads ~ Compartments * Tillage, xlab = "Samples", ylab = "%  sequencing reads",   main = "Reads assigned to Microbial OTUs"))

#Figure S1B
with(design_info_2, boxplot(reads ~ Compartments * Tillage, xlab = "Samples", ylab = " sequencing reads",   main = "Reads assigned to Microbial OTUs"))

#Calculate the max, min and mean number of reads for the dataset
mean(design_info_2$reads)
max(design_info_2$reads)
min(design_info_2$reads)

##########
#Figure 1
##########

#This table was generated in QIIME w/o plant-derived OTUs and biological replicates are already averaged according to the levels of the factor "Description" in the mapping file
dat_info_taxa_Phylum <- read.delim("Description_otu_table_L2.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
colnames(dat_info_taxa_Phylum)

# transform the data in % for visualisation
dat_count_taxa_Phylum <- dat_info_taxa_Phylum
dat_norm_Phylum <- ((dat_count_taxa_Phylum)/colSums(dat_count_taxa_Phylum,na=T)) * 100 
dat_norm_Phylum[1:5, ]

# determine the average % of reads for each phylum
Phylum_mean_sorted <- dat_norm_Phylum[(order(-rowSums(dat_norm_Phylum))), ] 

#Calculate the contribution of the top 10 Phyla to the total dataset
Phylum_mean_topRank <- Phylum_mean_sorted[1:10, ]
Phylum_mean_topRank 
colSums(Phylum_mean_topRank)

#overall
mean(colSums(Phylum_mean_topRank))

Phylum_mean_topRank_10 <- as.matrix(Phylum_mean_topRank[1:10, ]) 
#first we need to arrange the samples in a coherent way (i.e. according to the experiment)
colnames(Phylum_mean_topRank_10)


Phylum_mean_topRank_10_samples <- c("BkConventional", "RzConventional", "RtConventional", "ShConventional",
                                    "BkConservation", "RzConservation", "RtConservation", "ShConservation")

#now we will use the order of samples we have generated above to sort the count matrix
Phylum_mean_topRank_10_ordered <- Phylum_mean_topRank_10[ ,Phylum_mean_topRank_10_samples] 
colnames(Phylum_mean_topRank_10_ordered)

#Inspect the generated files
Phylum_mean_topRank_10_ordered[1:5, ]
Phylum_mean_topRank_10[1:5, ]

#Now we can plot them
# Stacked Bar Plot with Colors and Legend
barplot(Phylum_mean_topRank_10_ordered, main="Phylum Distribution",
        xlab="Samples", ylab = "Reads %", ylim = c(0,100), col=c("darkblue","red", "green", "cyan", "orange", "brown", "black", "magenta", "grey", "yellow"), beside=FALSE,   legend = rownames(Phylum_mean_topRank_10_ordered))

 
#############################################################
#generate a phyloseq object
################################
#a) The OTU able counts
Ridhdhi_OTU <- otu_table(dat_count_noplants, taxa_are_rows=TRUE)

#b) The taxonomy information
#Note this is a new file generated in excel from the output of the command of lines 100-103
#it is a tab-delimited file with 8 columns, the header names are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",
#and the empty cells are filled with the term 'Unassigned'
Ridhdhi_taxa_ordered <- read.delim ("Ridhdhi_dat_tax_noPlants_ordered2.txt", sep = "\t", row.names=1, header=T, blank.lines.skip = FALSE)
Ridhdhi_taxa <- tax_table(as.matrix(Ridhdhi_taxa_ordered))
dim(Ridhdhi_taxa)

#c) The mapping file 
Ridhdhi_map <- sample_data(design)

#d) The phylogenetic tree 
Ridhdhi_tree <- read_tree("Ridhdhi_OTU_table_noPlants_muscle_tree.tre")
#check whether the tree is rooted
is.rooted(Ridhdhi_tree)

#Root the tree
#Identify unique taxa
unique_OTUs <- unique(Ridhdhi_taxa[,2])
unique_OTUs

#OK we could use a Planctomycetes as an outgroup 
tax <- as.data.frame(tax_table(Ridhdhi_taxa))
outgroup <- tax[grepl("p__Planctomycetes", tax$Phylum), ]
dim(outgroup)
outgroup 
#OK now we can identify the top abundant OTU of this group
sort(rowSums(dat_count_noplants[rownames(outgroup), ]))

#With 701 reads, OTUs281 is relatively abundant, let's pick this OTUs as an outgroup
newRoot = c("OTU281")
Ridhdhi_tree <- root(Ridhdhi_tree,newRoot,resolve.root = TRUE)

#let's check it now
is.rooted(Ridhdhi_tree)

#merge the files and create the phyloseq object
Ridhdhi_data_phyloseq <- merge_phyloseq(Ridhdhi_OTU, Ridhdhi_taxa, Ridhdhi_map,  Ridhdhi_tree)

#inspect the generated data
Ridhdhi_data_phyloseq
sum(colSums(otu_table(Ridhdhi_data_phyloseq)))
dim(dat_count_noplants)
sum(colSums(dat_count_noplants))

#############################################################

#############################################################

##########
#Figure 2: alphadiversity calculation (rarefied dataset)
##########
Ridhdhi_data_phyloseq_Soil <- subset_samples(Ridhdhi_data_phyloseq, Habitat=="Soil" )                                             
Ridhdhi_data_phyloseq_Plant <- subset_samples(Ridhdhi_data_phyloseq, Habitat=="Plant")                                                                                  

#generate a rerefied dataset
#Ridhdhi_data_phyloseq_Soil_rare <- rarefy_even_depth(Ridhdhi_data_phyloseq_Soil)
#Ridhdhi_data_phyloseq_Plant_rare <- rarefy_even_depth(Ridhdhi_data_phyloseq_Plant)

#extract and save the OTU table for reproducibiity of the code
#Ridhdhi_data_phyloseq_Soil_rare_table <- as.data.frame(otu_table(Ridhdhi_data_phyloseq_Soil_rare))
#Ridhdhi_data_phyloseq_Plant_rare_table <- as.data.frame(otu_table(Ridhdhi_data_phyloseq_Plant_rare))

#inspect the generated file
#class(Ridhdhi_data_phyloseq_Soil_rare_table)
#class(Ridhdhi_data_phyloseq_Plant_rare_table)

#dim(Ridhdhi_data_phyloseq_Soil_rare_table)
#dim(Ridhdhi_data_phyloseq_Plant_rare_table)

#save the file for the reproducibility of the code

#write.table(Ridhdhi_data_phyloseq_Soil_rare_table, file="Ridhdhi_data_phyloseq_Soil_rare_table_counts.txt", sep="\t")
#write.table(Ridhdhi_data_phyloseq_Plant_rare_table, file="Ridhdhi_data_phyloseq_Plant_rare_table_counts.txt", sep="\t")

#import the rarefied OTU counts (note file name counts2.txt this file has been generated in excel and includes the #OTU ID as header of the first column)
dat_count_Soil_rare <- read.delim("Ridhdhi_data_phyloseq_Soil_rare_table.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
dat_count_Plant_rare <- read.delim("Ridhdhi_data_phyloseq_Plant_rare_table.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#inspect the generated file
dim(dat_count_Soil_rare)
dim(dat_count_Plant_rare)

#total count
colSums(dat_count_Soil_rare)
colSums(dat_count_Plant_rare)

#split the design files
Ridhdhi_map_soil <- Ridhdhi_map[colnames(dat_count_Soil_rare), ]
Ridhdhi_map_plant <- Ridhdhi_map[colnames(dat_count_Plant_rare), ]

#generate a new phyloseq object wich will contain only the rarefied counts and the design file (only these two pieces of information are required for alphadiversity calculation)
Ridhdhi_OTU_Soil_rare <- otu_table(dat_count_Soil_rare, taxa_are_rows=TRUE)
Ridhdhi_data_Soil_rare_phyloseq <- merge_phyloseq(Ridhdhi_OTU_Soil_rare , Ridhdhi_map_soil)
Ridhdhi_OTU_Plant_rare <- otu_table(dat_count_Plant_rare, taxa_are_rows=TRUE)
Ridhdhi_data_Plant_rare_phyloseq <- merge_phyloseq(Ridhdhi_OTU_Plant_rare , Ridhdhi_map_plant)

#Inspect the generated file
Ridhdhi_data_Soil_rare_phyloseq 
Ridhdhi_data_Plant_rare_phyloseq 

#number of reads per sample: original dataset different reads per sample
sample_sums(Ridhdhi_data_phyloseq_Soil)
sample_sums(Ridhdhi_data_phyloseq_Plant)

#number of reads per sample: rarefied dataset, they are all the same
sample_sums(Ridhdhi_data_Soil_rare_phyloseq)
sample_sums(Ridhdhi_data_Plant_rare_phyloseq)

#Index calculation
Ridhdhi_alpha_Soil_rare <-  estimate_richness(Ridhdhi_data_Soil_rare_phyloseq, measures = c("Observed", "Chao1", "Shannon"))
Ridhdhi_alpha_Plant_rare <-  estimate_richness(Ridhdhi_data_Plant_rare_phyloseq, measures = c("Observed", "Chao1", "Shannon"))

#generate a new dataframe for data visualisation and statistical analysis
design_soil <- design[rownames(Ridhdhi_alpha_Soil_rare), ]
Ridhdhi_alpha_Soil_rare_info <- cbind(design_soil, Ridhdhi_alpha_Soil_rare)

design_plant <- design[rownames(Ridhdhi_alpha_Plant_rare), ]
Ridhdhi_alpha_Plant_rare_info <- cbind(design_plant, Ridhdhi_alpha_Plant_rare)

#check the new dataset: it contains both description of the samples and alpha diversity indices 
Ridhdhi_alpha_Soil_rare_info 
Ridhdhi_alpha_Plant_rare_info

#generate a box plot for data visualisation
#re-order the factors
Ridhdhi_alpha_Soil_rare_info$Description <- ordered(Ridhdhi_alpha_Soil_rare_info$Description, levels=c("BkConventional", "RzConventional", "BkConservation", "RzConservation"))
Ridhdhi_alpha_Plant_rare_info$Description <- ordered(Ridhdhi_alpha_Plant_rare_info$Description, levels=c("ShConventional", "RtConventional", "ShConservation", "RtConservation"))

#observed
dev.off()
par(mfrow=c(2,1))
with(Ridhdhi_alpha_Soil_rare_info, boxplot(Observed  ~ Description, xlab = "Compartments", ylab = "Observed OTUs",   main = "Observed OTUs rarefied dataset"))
with(Ridhdhi_alpha_Plant_rare_info, boxplot(Observed  ~ Description, xlab = "Compartments", ylab = "Observed OTUs",   main = "Observed OTUs rarefied dataset"))
#chao 1
dev.off()
par(mfrow=c(2,1))
with(Ridhdhi_alpha_Soil_rare_info, boxplot(Chao1  ~ Description, xlab = "Compartments", ylab = "Chao1",   main = "Chao1 rarefied dataset"))
with(Ridhdhi_alpha_Plant_rare_info, boxplot(Chao1  ~ Description, xlab = "Compartments", ylab = "Chao1",   main = "Chao1 rarefied dataset"))
#shannon
dev.off()
par(mfrow=c(2,1))
with(Ridhdhi_alpha_Soil_rare_info, boxplot(Shannon  ~ Description, xlab = "Compartments", ylab = "Shannon",   main = "Shannon rarefied dataset"))
with(Ridhdhi_alpha_Plant_rare_info, boxplot(Shannon  ~ Description, xlab = "Compartments", ylab = "Shannon",   main = "Shannon rarefied dataset"))

#inspect the normal distribution of the datasets
shapiro.test(Ridhdhi_alpha_Soil_rare_info$Observed)
shapiro.test(Ridhdhi_alpha_Soil_rare_info$Chao1)
shapiro.test(Ridhdhi_alpha_Soil_rare_info$Shannon)

shapiro.test(Ridhdhi_alpha_Plant_rare_info$Observed)
shapiro.test(Ridhdhi_alpha_Plant_rare_info$Chao1)
shapiro.test(Ridhdhi_alpha_Plant_rare_info$Shannon)

#Define the tillage effect in soil habitat
#Observed
t.test(Observed ~ Tillage, data = Ridhdhi_alpha_Soil_rare_info)
#Chao1
t.test(Chao1 ~ Tillage, data = Ridhdhi_alpha_Soil_rare_info)
#Shannon
wilcox.test(Shannon ~ Tillage, data = Ridhdhi_alpha_Soil_rare_info)

#Define the Compartment effect in soil habitat
#Observed
t.test(Observed ~ Compartments, data = Ridhdhi_alpha_Soil_rare_info)
#Chao1
t.test(Chao1 ~ Compartments, data = Ridhdhi_alpha_Soil_rare_info)
#Shannon
wilcox.test(Shannon ~ Compartments, data = Ridhdhi_alpha_Soil_rare_info)

#Define the Tillage effect in plant habitat
#observed
wilcox.test(Observed ~ Tillage, data = Ridhdhi_alpha_Plant_rare_info)
#Chao1
wilcox.test(Chao1 ~ Tillage, data = Ridhdhi_alpha_Plant_rare_info)
#Shannon
wilcox.test(Shannon ~ Tillage, data = Ridhdhi_alpha_Plant_rare_info)

#Define the compartment effect in plant habitat
#Observed
wilcox.test(Observed ~ Compartments, data = Ridhdhi_alpha_Plant_rare_info)
#Chao1
wilcox.test(Chao1 ~ Compartments, data = Ridhdhi_alpha_Plant_rare_info)
#Shannon
wilcox.test(Shannon ~ Compartments, data = Ridhdhi_alpha_Plant_rare_info)

########################################################################################

########################################################################################

##########
#Figure 3: betadiversity calculation 
##########

#Transform the count in relative abundance
Ridhdhi_data_phyloseq_Soil_prop <- transform_sample_counts(Ridhdhi_data_phyloseq_Soil,  function(x) 1e+06 * x/sum(x))
Ridhdhi_data_phyloseq_Plant_prop <- transform_sample_counts(Ridhdhi_data_phyloseq_Plant,  function(x) 1e+06 * x/sum(x))

#PCoA bray distance plant
Ridhdhi_data_phyloseq_Plant_prop_bray <- ordinate(Ridhdhi_data_phyloseq_Plant_prop, "PCoA", "bray")
plot_ordination(Ridhdhi_data_phyloseq_Plant_prop, Ridhdhi_data_phyloseq_Plant_prop_bray, color = "Tillage")
#assign shapes to plant Compartments
p=plot_ordination(Ridhdhi_data_phyloseq_Plant_prop, Ridhdhi_data_phyloseq_Plant_prop_bray , color = "Tillage", shape="Compartments")
p = p + geom_point(size = 5, alpha = 0.75)
p + ggtitle("PCoA 16S data, Bray distance")

#PCoA bray distance Soil
Ridhdhi_data_phyloseq_Soil_prop_bray <- ordinate(Ridhdhi_data_phyloseq_Soil_prop, "PCoA", "bray")
plot_ordination(Ridhdhi_data_phyloseq_Soil_prop, Ridhdhi_data_phyloseq_Soil_prop_bray, color = "Tillage")
#assign shapes to Soil Compartments
p=plot_ordination(Ridhdhi_data_phyloseq_Soil_prop, Ridhdhi_data_phyloseq_Soil_prop_bray , color = "Tillage", shape="Compartments")
p = p + geom_point(size = 5, alpha = 0.75)
p + ggtitle("PCoA 16S data, Bray distance")

#PCoA wunifrac distance plant
Ridhdhi_data_phyloseq_Plant_prop_wunifrac <- ordinate(Ridhdhi_data_phyloseq_Plant_prop, "PCoA", "unifrac", weighted = TRUE)
plot_ordination(Ridhdhi_data_phyloseq_Plant_prop, Ridhdhi_data_phyloseq_Plant_prop_wunifrac, color = "Tillage")
#assign shapes to plant Compartments
p=plot_ordination(Ridhdhi_data_phyloseq_Plant_prop, Ridhdhi_data_phyloseq_Plant_prop_wunifrac , color = "Tillage", shape="Compartments")
p = p + geom_point(size = 5, alpha = 0.75)
p + ggtitle("PCoA 16S data, wunifrac distance")

#PCoA wunifrac distance Soil
Ridhdhi_data_phyloseq_Soil_prop_wunifrac <- ordinate(Ridhdhi_data_phyloseq_Soil_prop, "PCoA", "unifrac", weighted = TRUE)
plot_ordination(Ridhdhi_data_phyloseq_Soil_prop, Ridhdhi_data_phyloseq_Soil_prop_wunifrac, color = "Tillage")
#assign shapes to Soil Compartments
p=plot_ordination(Ridhdhi_data_phyloseq_Soil_prop, Ridhdhi_data_phyloseq_Soil_prop_wunifrac , color = "Tillage", shape="Compartments")
p = p + geom_point(size = 5, alpha = 0.75)
p + ggtitle("PCoA 16S data, wunifrac distance")

#permutational analysis of variance on dissimilarity matrices

#extract the dissimilarity matrices
BC_Soil <- phyloseq::distance(Ridhdhi_data_phyloseq_Soil_prop, "bray")
BC_Plant <- phyloseq::distance(Ridhdhi_data_phyloseq_Plant_prop, "bray")
WU_Soil <- phyloseq::distance(Ridhdhi_data_phyloseq_Soil_prop, "unifrac", weighted= TRUE)
WU_Plant <- phyloseq::distance(Ridhdhi_data_phyloseq_Plant_prop, "unifrac", weighted= TRUE)

#BC distance
adonis(BC_Soil ~ Tillage * Compartments, data = design_soil, permutations = 5000)
adonis(BC_Plant ~ Tillage * Compartments, data = design_plant, permutations = 5000)
#WU distance
adonis(WU_Soil ~ Tillage * Compartments, data= design_soil, permutations = 5000)
adonis(WU_Plant ~ Tillage * Compartments, data= design_plant, permutations = 5000)
############################################################### 

############################################################### 
##########
#Figure4
##########

#create a deseq object

#extract count data and 
Ridhdhi_OTU_counts_integer <- otu_table(Ridhdhi_data_phyloseq)
countData = as.data.frame(Ridhdhi_OTU_counts_integer)

#the design file containing sample information
colData = design

#construct a DESeq dataset combining count data and sample information
Ridhdhi_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData, design= ~ Description)

#execute the differential count analysis with the function DESeq 
Ridhdhi_cds_test <- DESeq(Ridhdhi_cds, fitType="local", betaPrior=FALSE)

#define the OTUs differentially enriched between compartments and tillage enriched in rhizosphere, root and shoot corrected for tillage

#conventional
Bulk_Rhizosphere_conventional <- results(Ridhdhi_cds_test , contrast = c("Description",  "BkConventional", "RzConventional"))
Root_Shoot_conventional <- results(Ridhdhi_cds_test , contrast = c("Description", "RtConventional", "ShConventional"))

#inspect a result file
Bulk_Rhizosphere_conventional 
mcols(Bulk_Rhizosphere_conventional, use.names=TRUE)
mcols(Root_Shoot_conventional, use.names = TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
Bulk_Rhizosphere_conventional_FDR_005 <- Bulk_Rhizosphere_conventional[(rownames(Bulk_Rhizosphere_conventional)[which(Bulk_Rhizosphere_conventional$padj <0.05)]), ]
Root_Shoot_conventional_FDR_005 <- Root_Shoot_conventional[(rownames(Root_Shoot_conventional)[which(Root_Shoot_conventional$padj <0.05)]), ]

#Identify OTUs enriched in the rhizosphere vs soil, and in shoot vs. roots (second term of the comparison)
Rhizosphere_conventional_enriched <-  Bulk_Rhizosphere_conventional[(rownames(Bulk_Rhizosphere_conventional)[which(Bulk_Rhizosphere_conventional$log2FoldChange < 0)]), ]
Shoot_conventional_enriched <-  Root_Shoot_conventional[(rownames(Root_Shoot_conventional)[which(Root_Shoot_conventional$log2FoldChange < 0)]), ]

#Commands on lines 487/493  provides lists of OTUs fulfilling the imposed criteria. To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
Rhizosphere_conventional_enriched_FDR005 <- intersect(rownames(Bulk_Rhizosphere_conventional_FDR_005), rownames(Rhizosphere_conventional_enriched))
Shoot_conventional_enriched_FDR005 <- intersect(rownames(Root_Shoot_conventional_FDR_005), rownames(Shoot_conventional_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(Rhizosphere_conventional_enriched_FDR005)
length(Shoot_conventional_enriched_FDR005)

#Identify OTUs enriched in the soil vs. rhizosphere, and in roots vs. shoot (first term of the comparison)
Soil_conventional_enriched <-  Bulk_Rhizosphere_conventional[(rownames(Bulk_Rhizosphere_conventional)[which(Bulk_Rhizosphere_conventional$log2FoldChange > 0)]), ]
Root_conventional_enriched <-  Root_Shoot_conventional[(rownames(Root_Shoot_conventional)[which(Root_Shoot_conventional$log2FoldChange > 0)]), ]

#Commands on lines 487-489/503-505  provides lists of OTUs fulfilling the imposed criteria. To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
Soil_conventional_enriched_FDR005 <- intersect(rownames(Bulk_Rhizosphere_conventional_FDR_005), rownames(Soil_conventional_enriched))
Root_conventional_enriched_FDR005 <- intersect(rownames(Root_Shoot_conventional_FDR_005), rownames(Root_conventional_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(Soil_conventional_enriched_FDR005)
length(Root_conventional_enriched_FDR005)


##########################################################################################################

#conservation
Bulk_Rhizosphere_conservation <- results(Ridhdhi_cds_test , contrast = c("Description",  "BkConservation", "RzConservation"))
Root_Shoot_conservation <- results(Ridhdhi_cds_test , contrast = c("Description", "RtConservation", "ShConservation"))

#inspect a result file
Bulk_Rhizosphere_conservation 
mcols(Bulk_Rhizosphere_conservation, use.names=TRUE)
mcols(Root_Shoot_conservation, use.names = TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
Bulk_Rhizosphere_conservation_FDR_005 <- Bulk_Rhizosphere_conservation[(rownames(Bulk_Rhizosphere_conservation)[which(Bulk_Rhizosphere_conservation$padj <0.05)]), ]
Root_Shoot_conservation_FDR_005 <- Root_Shoot_conservation[(rownames(Root_Shoot_conservation)[which(Root_Shoot_conservation$padj <0.05)]), ]

#Identify OTUs enriched in the rhizosphere, root and shoot compartments (second term of the comparison)
Rhizosphere_conservation_enriched <-  Bulk_Rhizosphere_conservation[(rownames(Bulk_Rhizosphere_conservation)[which(Bulk_Rhizosphere_conservation$log2FoldChange < 0)]), ]
Shoot_conservation_enriched <-  Root_Shoot_conservation[(rownames(Root_Shoot_conservation)[which(Root_Shoot_conservation$log2FoldChange < 0)]), ]

#Commands on lines 527/533 provides lists of OTUs fulfilling the imposed criteria. To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
Rhizosphere_conservation_enriched_FDR005 <- intersect(rownames(Bulk_Rhizosphere_conservation_FDR_005), rownames(Rhizosphere_conservation_enriched))
Shoot_conservation_enriched_FDR005 <- intersect(rownames(Root_Shoot_conservation_FDR_005), rownames(Shoot_conservation_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(Rhizosphere_conservation_enriched_FDR005)
length(Shoot_conservation_enriched_FDR005)

#Identify OTUs enriched in the soil vs. rhizosphere, and in roots vs. shoot (first term of the comparison)
Soil_conservation_enriched <-  Bulk_Rhizosphere_conservation[(rownames(Bulk_Rhizosphere_conservation)[which(Bulk_Rhizosphere_conservation$log2FoldChange > 0)]), ]
Root_conservation_enriched <-  Root_Shoot_conservation[(rownames(Root_Shoot_conservation)[which(Root_Shoot_conservation$log2FoldChange > 0)]), ]

#Commands on lines 527-529/543-545  provides lists of OTUs fulfilling the imposed criteria. To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
Soil_conservation_enriched_FDR005 <- intersect(rownames(Bulk_Rhizosphere_conservation_FDR_005), rownames(Soil_conservation_enriched))
Root_conservation_enriched_FDR005 <- intersect(rownames(Root_Shoot_conservation_FDR_005), rownames(Root_conservation_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(Soil_conservation_enriched_FDR005)
length(Root_conservation_enriched_FDR005)

#Figure 4
#set up the screen
dev.off()
par(mfrow=c(2,1))
#Figure 4A
plotMA(Bulk_Rhizosphere_conventional, alpha = 0.05, main="Bulk vs. Rhizosphere CT")
#Figure 4C
plotMA(Bulk_Rhizosphere_conservation, alpha = 0.05, main="Bulk vs. Rhizosphere ST")

dev.off()
par(mfrow=c(2,1))
#Figure DB
plotMA(Root_Shoot_conventional, alpha = 0.05, main= "Root vs shoot CT")
#Figure 4D
plotMA(Root_Shoot_conservation, alpha = 0.05, main="Root vs shoot ST")

##########
#Figure5
##########

#note: this diagram now represents the OTU enriched in root verus shoot and shoot versus root in both treatments (and not anymore OTUs enriched versus unplanted soil)

#define the for main areas of the diagram
#Area 1
length(Root_conventional_enriched_FDR005)
#Area 2
length(Root_conservation_enriched_FDR005)
#Area 3
length(Shoot_conventional_enriched_FDR005)
#Area 4
length(Shoot_conservation_enriched_FDR005)
#intersection
#Area 12
Area12 <- intersect(Root_conventional_enriched_FDR005,Root_conservation_enriched_FDR005)
length(Area12)
#Area 34
Area34 <- intersect(Shoot_conventional_enriched_FDR005,Shoot_conservation_enriched_FDR005)
length(Area34)
#Area 14
Area14 <- intersect(Root_conventional_enriched_FDR005,Shoot_conservation_enriched_FDR005)
length(Area14)
#Area 23
Area23 <- intersect(Root_conservation_enriched_FDR005,Shoot_conventional_enriched_FDR005)
length(Area23)
#Area 24
Area24 <- intersect(Root_conservation_enriched_FDR005,Shoot_conservation_enriched_FDR005)
length(Area24)
#Area 13
Area13 <- intersect(Root_conventional_enriched_FDR005,Shoot_conventional_enriched_FDR005)
length(Area13)
#Area 123
Area123 <- intersect(Area12,Shoot_conventional_enriched_FDR005)
length(Area123)
#Area 124
Area124 <- intersect(Area12,Shoot_conservation_enriched_FDR005)
length(Area124)
#Area 134
Area134 <- intersect(Area13,Shoot_conventional_enriched_FDR005)
length(Area134)
#Area 234
Area234 <- intersect(Area23,Shoot_conservation_enriched_FDR005)
length(Area234)
#Area 1234
Area1234 <- intersect(Area123,Shoot_conservation_enriched_FDR005)
length(Area1234)

#drawing
dev.off()
draw.quad.venn(area1 = 368, area2 = 174, area3 = 39, area4 = 51,
               n12 = 163, n13 = 0, n14 = 0, n23 = 1, n24 = 0 , n34 = 25, n123 = 0, n124 = 0, n134 = 0, n234 = 0, n1234 = 0,
               category = c("Root Conventional", "Root Conservation", "Shoot Conventional", "Shoot Conservation"), lty = "blank", 
               fill = c("skyblue", "pink1", "darkblue", "red"))


#end#


################################################################################################################################

#Figure s6
#create a deseq object
#CT vs ST
#define the OTUs significantly enriched in Conservation and conventional tillage corrected for compartments (bulk soil, rhizosphere, root and shoot)
Bulk_conservation_conventional <- results(Ridhdhi_cds_test , contrast = c("Description", "BkConservation", "BkConventional"))
Rhizosphere_conservation_conventional <- results(Ridhdhi_cds_test , contrast = c("Description", "RzConservation", "RzConventional"))
Root_conservation_conventional <- results(Ridhdhi_cds_test , contrast = c("Description", "RtConservation", "RtConventional"))
Shoot_conservation_conventional <- results(Ridhdhi_cds_test , contrast = c("Description", "ShConservation", "ShConventional"))

#inspect a result file
Bulk_conservation_conventional 
mcols(Bulk_conservation_conventional, use.names=TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
Bulk_conservation_conventional_FDR_005 <- Bulk_conservation_conventional[(rownames(Bulk_conservation_conventional)[which(Bulk_conservation_conventional$padj <0.05)]), ]
Rhizosphere_conservation_conventional_FDR_005 <- Rhizosphere_conservation_conventional[(rownames(Rhizosphere_conservation_conventional)[which(Rhizosphere_conservation_conventional$padj <0.05)]), ]
Root_conservation_conventional_FDR_005 <- Root_conservation_conventional[(rownames(Root_conservation_conventional)[which(Root_conservation_conventional$padj <0.05)]), ]
Shoot_conservation_conventional_FDR_005 <- Shoot_conservation_conventional[(rownames(Shoot_conservation_conventional)[which(Shoot_conservation_conventional$padj <0.05)]), ]

#Identify OTUs enriched in the rhizosphere, root and shoot compartments (second term of the comparison)
Bulk_conservation_ST_enriched <- Bulk_conservation_conventional[(rownames(Bulk_conservation_conventional)[which(Bulk_conservation_conventional$log2FoldChange < 0)]), ]
Rhizosphere_conservation_ST_enriched <- Rhizosphere_conservation_conventional[(rownames(Rhizosphere_conservation_conventional)[which(Rhizosphere_conservation_conventional$log2FoldChange < 0)]), ]
Root_conservation_ST_enriched <- Root_conservation_conventional[(rownames(Root_conservation_conventional)[which(Root_conservation_conventional$log2FoldChange < 0)]), ]
Shoot_conservation_ST_enriched <- Shoot_conservation_conventional[(rownames(Shoot_conservation_conventional)[which(Shoot_conservation_conventional$log2FoldChange < 0)]), ]

#Intersection
Bulk_conservation_ST_enriched_FDR005 <- intersect(rownames(Bulk_conservation_conventional_FDR_005), rownames(Bulk_conservation_ST_enriched))
Rhizosphere_conservation_ST_enriched_FDR005 <- intersect(rownames(Rhizosphere_conservation_conventional_FDR_005), rownames(Rhizosphere_conservation_ST_enriched))
Root_conservation_ST_enriched_FDR005 <- intersect(rownames(Root_conservation_conventional_FDR_005), rownames(Root_conservation_ST_enriched))
Shoot_conservation_ST_enriched_FDR005 <- intersect(rownames(Shoot_conservation_conventional_FDR_005), rownames(Shoot_conservation_ST_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(Bulk_conservation_ST_enriched_FDR005)
length(Rhizosphere_conservation_ST_enriched_FDR005)
length(Root_conservation_ST_enriched_FDR005)
length(Shoot_conservation_ST_enriched_FDR005)

#create supplementary database combining statisical and tanonomy information for the enriched OTUs
Rhizosphere_conservation_ST_enriched_taxa <- Ridhdhi_taxa_ordered[Rhizosphere_conservation_ST_enriched_FDR005, ]
Rhizosphere_conservation_ST_enriched_FDR005_taxa <- cbind(as.data.frame(Rhizosphere_conservation_conventional[Rhizosphere_conservation_ST_enriched_FDR005, ]), Rhizosphere_conservation_ST_enriched_taxa)

Root_conservation_ST_enriched_taxa <- Ridhdhi_taxa_ordered[Root_conservation_ST_enriched_FDR005, ]
Root_conservation_ST_enriched_FDR005_taxa <- cbind(as.data.frame(Root_conservation_conventional[Root_conservation_ST_enriched_FDR005, ]), Root_conservation_ST_enriched_taxa)

Shoot_conservation_ST_enriched_taxa <- Ridhdhi_taxa_ordered[Shoot_conservation_ST_enriched_FDR005, ]
Shoot_conservation_ST_enriched_FDR005_taxa <- cbind(as.data.frame(Shoot_conservation_conventional[Shoot_conservation_ST_enriched_FDR005, ]), Shoot_conservation_ST_enriched_taxa)

#save the files
#write.table(Bulk_conservation_ST_enriched_FDR005_taxa, file="Bulk_conservation_ST_enriched_FDR005_taxa_table.txt", sep="\t")
#write.table(Rhizosphere_conservation_ST_enriched_FDR005_taxa, file="Rhizosphere_conservation_ST_enriched_FDR005_taxa_table.txt", sep="\t")
#write.table(Root_conservation_ST_enriched_FDR005_taxa, file="Root_conservation_ST_enriched_FDR005_taxa_table.txt", sep="\t")
#write.table(Shoot_conservation_ST_enriched_FDR005_taxa, file="Shoot_conservation_ST_enriched_FDR005_taxa_table.txt", sep="\t")

#visualise the data
#set up the screen
dev.off()
par(mfrow=c(2,1))
#Figure S6A
plotMA(Bulk_conservation_conventional, alpha = 0.05, main="ST vs CT, Bulk soil")
#Figure s6B
plotMA(Rhizosphere_conservation_conventional, alpha = 0.05, main="ST vs CT, Rhizosphere")

dev.off()
par(mfrow=c(2,1))
#Figure S6C
plotMA(Root_conservation_conventional, alpha = 0.05, main="ST vs CT, Root")
#Figure S6D
plotMA(Shoot_conservation_conventional, alpha = 0.05, main="St vs CT, Shoot")

################################################################################################################################################

#Define the number of OTUs excluded from the plant-associated compartments
Bulk_conservation_ST_depleted <- Bulk_conservation_conventional[(rownames(Bulk_conservation_conventional)[which(Bulk_conservation_conventional$log2FoldChange > 0)]), ]
Bulk_conservation_ST_depleted_FDR005 <- intersect(rownames(Bulk_conservation_conventional_FDR_005), rownames(Bulk_conservation_ST_depleted))
length(Bulk_conservation_ST_depleted_FDR005)

Rhizosphere_conservation_ST_depleted <- Rhizosphere_conservation_conventional[(rownames(Rhizosphere_conservation_conventional)[which(Rhizosphere_conservation_conventional$log2FoldChange > 0)]), ]
Rhizosphere_conservation_ST_depleted_FDR005 <- intersect(rownames(Rhizosphere_conservation_conventional_FDR_005), rownames(Rhizosphere_conservation_ST_depleted))
length(Rhizosphere_conservation_ST_depleted_FDR005)

Root_conservation_ST_depleted <- Root_conservation_conventional[(rownames(Root_conservation_conventional)[which(Root_conservation_conventional$log2FoldChange > 0)]), ]
Root_conservation_ST_depleted_FDR005 <- intersect(rownames(Root_conservation_conventional_FDR_005), rownames(Root_conservation_ST_depleted))
length(Root_conservation_ST_depleted_FDR005)

Shoot_conservation_ST_depleted <- Shoot_conservation_conventional[(rownames(Shoot_conservation_conventional)[which(Shoot_conservation_conventional$log2FoldChange > 0)]), ]
Shoot_conservation_ST_depleted_FDR005 <- intersect(rownames(Shoot_conservation_conventional_FDR_005), rownames(Shoot_conservation_ST_depleted))
length(Shoot_conservation_ST_depleted_FDR005)

#create supplementary database combining statisical and tanonomy information for the depleted OTUs
Rhizosphere_conservation_ST_depleted_taxa <- Ridhdhi_taxa_ordered[Rhizosphere_conservation_ST_depleted_FDR005, ]
Rhizosphere_conservation_ST_depleted_FDR005_taxa <- cbind(as.data.frame(Rhizosphere_conservation_conventional[Rhizosphere_conservation_ST_depleted_FDR005, ]), Rhizosphere_conservation_ST_depleted_taxa)

Root_conservation_ST_depleted_taxa <- Ridhdhi_taxa_ordered[Root_conservation_ST_depleted_FDR005, ]
Root_conservation_ST_depleted_FDR005_taxa <- cbind(as.data.frame(Root_conservation_conventional[Root_conservation_ST_depleted_FDR005, ]), Root_conservation_ST_depleted_taxa)

Shoot_conservation_ST_depleted_taxa <- Ridhdhi_taxa_ordered[Shoot_conservation_ST_depleted_FDR005, ]
Shoot_conservation_ST_depleted_FDR005_taxa <- cbind(as.data.frame(Shoot_conservation_conventional[Shoot_conservation_ST_depleted_FDR005, ]), Shoot_conservation_ST_depleted_taxa)

#save these files
#write.table(Bulk_conservation_ST_depleted_FDR005_taxa, file="Bulk_conservation_ST_depleted_FDR005_taxa_table.txt", sep="\t")
#write.table(Rhizosphere_conservation_ST_depleted_FDR005_taxa, file="Rhizosphere_ST_conservation_depleted_FDR005_taxa_table.txt", sep="\t")
#write.table(Root_conservation_ST_depleted_FDR005_taxa, file="Root_conservation_ST_depleted_FDR005_taxa_table.txt", sep="\t")
#write.table(Shoot_conservation_ST_depleted_FDR005_taxa, file="Shoot_conservation_ST_depleted_FDR005_taxa_table.txt", sep="\t")

########################################################################################################

