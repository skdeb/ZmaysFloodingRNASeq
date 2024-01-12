###############################################################
# Differential gene expression analysis using limma voom
# To build this script I followed tutorials from here (https://diytranscriptomics.com/)
#The script was run seperately for each species in their respective directories. 
###############################################################

# Load libraries

library(rhdf5)
library(tidyverse)
library(tximport)
library(edgeR)
library(ggplot2)
library(viridis)
library(pheatmap)
library(cowplot)
library(ggVennDiagram)
library(here)

# Load metadata file
metadata <- read_tsv("Flooding_meta_data.txt")

# Set variables with different levels as factors

Condition <- factor(metadata$Condition,levels = c("CONTROL", "SHORT", "LONG"))
Tissue <- factor(metadata$Tissue) 
RepGroup <- factor(metadata$RepGroup) 
Group <- factor(paste(Condition, Tissue, sep = "."))

# set the design matrix

design <- model.matrix(~0 + Group)
design #check
colnames(design) <- levels(Group) # fixing column name of the design matrix
design #check again


# Load kallisto output for dge analysis
# set file path
path <- file.path(metadata$Sample, "abundance.tsv")
all(file.exists(path)) # check if all the files exist

# Load gene to transcripts map for gene level quantification

##Zea transcripts to gene map
library(biomaRt)
zea.anno <- useMart("plants_mart",dataset="zmays_eg_gene", host="https://plants.ensembl.org")
zea.attributes <- listAttributes(zea.anno)
Tx.zea <- getBM(attributes=c('ensembl_transcript_id',
                             'ensembl_gene_id'),
                mart = zea.anno)

Tx.zea <- as_tibble(Tx.zea)
Tx.gene <- dplyr::rename(Tx.zea, target_id = ensembl_transcript_id, 
                        gene_name = ensembl_gene_id)

Tx.gene <- Tx.gene[, c(2, 1)] # swap to bring transcirpt id in 1st column

Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx.gene, 
                     txOut = FALSE, 
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)


# Others: Znic, Trip, and Vossia transcripts to gene mapping
# Trinity gene trans map file from Trinity transcriptome assembly would work.
# make sure first column is transcript id and the second column is gene id.
# Trinity output map in other way. genes first and then transcript second.


gene_trans_map <- read_tsv("trinity_out.Trinity.fasta.gene_trans_map.tsv")
Tx.gene <- as_tibble(gene_trans_map)
Tx.gene <- dplyr::rename(Tx.gene, target_id = transcript_id, 
                        gene_name = gene_id)

head(Tx.gene) # check everything ok

Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx.gene, 
                     txOut = FALSE, 
                     countsFromAbundance = "lengthScaledTPM")
                     
###########################################################
sample.labels <- metadata$Sample # get sample names to change the column names later.

# Examine data up to this point

data.tpm <- Txi_gene$abundance
data.counts <- Txi_gene$counts
colSums(data.tpm)
colSums(data.counts)                    

# save myTPM and myCounts data for using them later if required

data.tpm.df <- as_tibble(data.tpm, rownames = "geneID")
colnames(data.tpm.df) <- c("geneID", sample.labels)
write_csv(data.tpm.df,"tpm_data.csv")
dim(data.tpm.df)

data.counts.df <- as_tibble(data.counts, rownames = "geneID")
colnames(data.counts.df) <- c("geneID", sample.labels)
write_csv(data.counts.df,"counts_data.csv")
head(data.counts) #check
dim(data.counts.df) #check

# Create a DGEList object in R using the count data stored in data.counts.

dge.list <- DGEList(data.counts)
save(dge.list, file = "DgeList")
# normalize gene expression count using Counts Per Million (CPM)
# filter data and visualize in each step

cpm <- cpm(dge.list) 
colSums(cpm)
dge.log2.cpm <- cpm(dge.list, log=TRUE)
dge.log2.cpm.df <- as_tibble(dge.log2.cpm, rownames = "geneID")
colnames(dge.log2.cpm.df) <- c("geneID", sample.labels)
dge.log2.cpm.df.pivot <- pivot_longer(dge.log2.cpm.df,
                                  cols = 2:19,
                                  names_to = "samples",
                                  values_to = "expression")
#########
# plot the unfiltered cpm data

p1 <- ggplot(dge.log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 124, 
               size = 6, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Unfiltered, non-normalized") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  coord_flip()
p1

##########
# Filter data
# Remove genes with 0 counts in all samples
# keep genes with more than 0.5 cpm in at least 3 samples

table(rowSums(dge.list$counts==0)==18)
keep <- rowSums(cpm>0.5)>=3
dge.list.filtered <- dge.list[keep,]
dim(dge.list.filtered)

###########
# Hierarchical clustering
dge.list.log2.cpm.filtered <- cpm(dge.list.filtered, log=TRUE)
distance <- dist(t(dge.list.log2.cpm.filtered), method = "euclidean")
clusters <- hclust(distance, method = "average")
plot(clusters, labels=sample.labels)

# PCA analysis using the filtered data
pca.res <- prcomp(t(dge.list.log2.cpm.filtered), scale.=F, retx=T)
ls(pca.res)
summary(pca.res)  
pc.var<-pca.res$sdev^2 
pc.per<-round(pc.var/sum(pc.var)*100, 1)

# converting PCA results as tiblle for plotting
pca.res.df <- as_tibble(pca.res$x)

p2 <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2,color = Condition, shape = Tissue) +
  geom_point(size = 4, stroke = 1) +
  scale_shape_manual(values=c(1, 5)) +
  scale_color_manual(values=c("#ba4343", "#22A884FF", "#440154FF"),
                     labels = c("CONTROL (0h)", "SHORT-TERM (4h)", 
                                "LONG-TERM (72h)")) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) +
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  #labs(title="PCA plot")+
  coord_fixed() +
  theme_bw()
  #theme(plot.title = element_text(hjust = 0.5))

p2
# Use quantile normalize from voom function for normalization of data
voom.dge.list.filtered.norm <- voom(dge.list.filtered, design, plot = TRUE,normalize.method ="quantile")

head (voom.dge.list.filtered.norm$E)
## Visulalize
data.voom.qt.norm <- voom.dge.list.filtered.norm$E
data.voom.qt.norm.df <- as_tibble(data.voom.qt.norm, rownames = "geneID")
colnames(data.voom.qt.norm.df) <- c("geneID", sample.labels)
write_csv(data.voom.qt.norm.df,"data.voom.quantile.norm.csv")
data.voom.qt.norm.df.pivot <- pivot_longer(data.voom.qt.norm.df,
                                  cols = 2:19,
                                  names_to = "samples",
                                  values_to = "expression")
#########

p3 <- ggplot(data.voom.qt.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 124, 
               size = 6, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Filtered, quantile-normalized") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  coord_flip()
p3
plot_grid(p1, p3)


# Get correlation between replicates of samples belonging to same group

corr.data.voom.qt.norm <- voom.dge.list.filtered.norm$E
colnames(corr.data.voom.qt.norm) <- sample.labels
corrSamples <- cor(corr.data.voom.qt.norm)
head(corrSamples)

## Plot
pheatmap(corrSamples, color=viridis::viridis(50),
         display_numbers = TRUE,
         number_color = "black", 
         fontsize_number = 4)



########################################################################
# Differential gene expression analysis

# Do contrasts

contrast.matrix <- makeContrasts(L4h_L0h = SHORT.LEAF - CONTROL.LEAF,
                                 L72h_L0h = LONG.LEAF - CONTROL.LEAF,
                                 R4h_R0h = SHORT.ROOT - CONTROL.ROOT,
                                 R72h_R0h = LONG.ROOT - CONTROL.ROOT,
                                 levels=design)


#### Fit and do Diff Expression

vfit <- lmFit(voom.dge.list.filtered.norm, design) 
vfit2 <- contrasts.fit(vfit, contrast.matrix)
efit <- eBayes(vfit2)

################################
## -- Summary and Venn diagrams
results <- decideTests(efit, adjust.method="BH", p.value=0.05, lfc=1)
summary(results)
DESummary <- t(summary(results))[,-2]
write.csv(x=DESummary,"DESummary.csv",row.names = T)

################################
# Save differentially expressed genes along with their normalized expression values

diffGenes.all <- voom.dge.list.filtered.norm$E[results[,1]!=0 |results[,2]!=0 | results[,3]!=0 | 
                                             results[,4]!=0,]

colnames(diffGenes.all) <- sample.labels
head(diffGenes.all) #check
dim(diffGenes.all) #check
diffGenes.all.df <- as_tibble(diffGenes.all, rownames = "geneID")
write_csv(diffGenes.all.df,"diffGenes_all.csv")

# Get average normalized gene expression values between the replicates

colnames(diffGenes.all) <- metadata$RepGroup # Giving each replicates group ids
diffGenes.AVG <- avearrays(diffGenes.all) # Take the average of replicates in each group
head(diffGenes.AVG) #check
dim(diffGenes.AVG) # check
diffGenes.AVG.df <- as_tibble(diffGenes.AVG, rownames = "geneID")
write_csv(diffGenes.AVG.df,"diffGenes.AVG.csv")

######################################################################
# Saving files of significant DGE list for separate contrasts
# The following chunks of code adapted and modified from https://github.com/plant-plasticity/tomato-root-atlas-2020/tree/master

# directory to save files
DE_list_dir <- here("DE_list")

# save files
DEList <- list()
for (contrast in colnames(contrast.matrix)){
  print(contrast)
  ## Sorting by none ensures all contrasts will be in the same order
  tmp <- topTable(efit,coef=contrast,number = Inf,sort.by = "none")

  ## Write genes that are up or downregulated (logFC > 1; logFC < -1)
  upGenes <- as.data.frame(rownames(tmp[tmp$adj.P.Val < 0.05 & tmp$logFC > 1,]))
  tmpSave <- paste(DE_list_dir, contrast,"_up",".csv",sep="")
  write.csv(x=upGenes,tmpSave,quote = F,row.names = T)
  #
  downGenes <- as.data.frame(rownames(tmp[tmp$adj.P.Val < 0.05 & tmp$logFC < -1,]))
  tmpSave <- paste(DE_list_dir, contrast,"_down",".csv",sep="")
  write.csv(x=downGenes,tmpSave,quote = F,row.names = T)
  #####
  ## Add contrast name to the column names, in case of multiple contrasts.
  colnames(tmp) <- paste(colnames(tmp),contrast,sep = ".")
  
  # Write each contrast to file
  tmpSave <- paste(DE_list_dir,contrast, ".csv",sep="")
  write.csv(x=tmp,tmpSave,quote = F,row.names = T)
  
  # Save result to list
  DEList[[contrast]] <- tmp 
}
tmpSave <- paste(DE_list_dir, "DEList.RData",sep="")
save(DEList,file = tmpSave)

##########
# Save log2FC values in each contrast
# Write a function to condense the results from all the contrasts
condenseListTables <- function(listDFs) {
  cat ("-- make condenseListTables function called \n")
  # First make an empty table
  uniqRows <- Reduce(union,lapply(listDFs,rownames))
  uniqCols <- unlist(sapply(listDFs, colnames))
  zeroTable <- as.data.frame(matrix(0,
                                    nrow = length(uniqRows),
                                    length(uniqCols)),
                             row.names =uniqRows)
  colnames(zeroTable) = uniqCols
  # Then Fill it
  for (each in names(listDFs)){
    cat ("Filling binary table:", each,"\n")
    zeroTable[rownames(listDFs[[each]]),
              colnames(listDFs[[each]])] <- listDFs[[each]]
  }

  return(zeroTable)
  cat ("-- binary table done \n")
  cat ("\n")
}
#################################################

### Use the function to Condense all DE files into a single list

DE_All <- condenseListTables(DEList) 
DE_All <- DE_All[,-grep("t.|B.|P.Value|AveExpr",colnames(DE_All))] #Remove unwanted columns

write.csv(x = DE_All,"DE_All_logFC.csv",quote = F,row.names = T)


#############################################################
#############################################################
###############################
# Venn diagram
##############################
## VennDiagram for both up and down regulated genes in leaf and root tissues
## Taking only contrast between 72h and control condition in both leaf and root.

venn.data <- list (rownames(results)[results[,"R72h_R0h"] %in% c(1, -1)],
                   rownames(results)[results[,"L72h_L0h"] %in% c(1, -1)])
                  
names(venn.data) <- c("Root72h_vs_0h","Leaf72h_vs_0h")
venn <- Venn(venn.data)
dat <- process_data(venn, shape_id == "201f")

ggplot() +
  geom_sf(aes(fill = name), data = venn_region(dat)) +
  geom_sf_text(aes(label = count),size = 6, data = venn_region(dat)) +
  scale_fill_manual(values=c("#99d8c9","#756bb1", "#d9f0a3"))+
  guides(fill = "none") +
  theme_void()












