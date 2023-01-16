# From ming tommy tang r blog https://rpubs.com/crazyhottommy/heatmap_demystified

library(dplyr)
library(tidyr)
library(ggplot2)
set.seed(1)

# repeat the sampling 
# Creates a df with 20 columns, and 10 rows. Randomly assigns a value of 0 or 1 for the rows
mut<- replicate(20, sample(c(0,1), 10, replace=TRUE))
mut<- as.data.frame(mut)
print(mut)

# rename columns "samplex" rename rows "genex"
colnames(mut)<- paste0("sample", 1:20)
mut<- mut %>% mutate(gene=paste0("gene", 1:10))
head(mut)

# transposes the df into long form. Samples in ascending order down (ie. sample1 gene1, sample1 gene2...)
mut.tidy<- mut %>% tidyr::gather(sample, mutated, 1:20)


## change the levels for gene names and sample names so it goes 1,2,3,4... rather than 1, 10...
mut.tidy$gene<- factor(mut.tidy$gene, levels = paste0("gene", 1:10))
mut.tidy$sample<- factor(mut.tidy$sample, levels = paste0("sample", 1:20))

# mutated is a binary value, not continous. Change to as factor so R reads it as 'has' or 'doesn't have'.
mut.tidy$mutated<- factor(mut.tidy$mutated)

## use a white border of size 0.5 unit to separate the tiles
gg<- ggplot(mut.tidy, aes(x=sample, y=gene, fill=mutated)) + geom_tile(color="white", linewidth=0.5)

library(RColorBrewer) ## better color schema

## check all the color pallete and choose one, here we use 'set1'
display.brewer.all()
gg<- gg + scale_fill_brewer(palette = "Set1", direction = -1)

# geom_tile() draws rectangles, add coord_equal to draw squres.
gg<- gg + coord_equal()
## add title
gg<- gg + labs(x=NULL, y=NULL, title="mutation spectrum of 20 breast cancers")

library(ggthemes)
##starting with a base theme of theme_tufte() from the ggthemes package. It removes a lot of chart junk without having to do it manually.
gg <- gg + theme_tufte(base_family="Helvetica")

#We don’t want any tick marks on the axes and I want the text to be slightly smaller than the default.
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + theme(axis.text.x=element_text(angle = 45, hjust = 1))
gg


# If you want to mannually fill the color, you can use scale_fill_manual, and check http://colorbrewer2.org/ to get the HEX representation of the color.

ggplot(mut.tidy, aes(x=sample, y=gene, fill=mutated)) + geom_tile(color="white", size=0.5) +
  coord_equal() +
  labs(x=NULL, y=NULL, title="mutation spectrum of 20 breast cancers") +
  theme_tufte(base_family="Helvetica") +
  scale_fill_manual(values = c("#7570b3", "#1b9e77")) +
  theme(axis.ticks=element_blank()) + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1))



# Using Heatmap to represent continuous variables
# In my last blog post, I showed an example to use heatmap to represent discrete values. 
# I am going to continue the theme to introduce heatmap to represent continuous values.
# As I mentioned before, I will stress on three main points:
#     1.  Whether or not scale your data (center/standardize your data or not).
#     2.  Make sure you know the range of the data and do reasonable color mapping.
#     3.  How to perform the clustering. (which distance measure and linkage method to use).
# In order to make the example reproducible and interesting to biologist, 
# I am going to use an RNAseq data set downloaded by the recount package. 
# The data contains 347 single cell RNA-seq data of 11 different distinct populations. 
# Please read this Nature biotechnology paper:
# Low-coverage single-cell mRNA sequencing reveals cellular heterogeneity and activated signaling pathways in developing cerebral cortex. 
# I will try to reproduce the figures in the paper.



# BiocManager::install("recount")
library(recount)
library(dplyr)

## if you go to the paper you can find the project name
project = 'SRP041736'
download_study(project)

## load the data
load(file.path(project, 'rse_gene.Rdata'))
#We can explore a bit this RangedSummarizedExperiment 
rse_gene

colData(rse_gene) %>% head()



## At the gene level, the row data includes the gene ENTREZ ids, the gene
## symbols and the sum of the reduced exons widths, which can be used for 
## taking into account the gene length.
rowData(rse_gene)


# now to download the metadata needed 
# BiocManager::install("SRAdb")
# cant seem to get this to work, this is a problem for future me
library(SRAdb)
project 
# magrittr pipe
'%>%' <- dplyr:: '%>%'
# Declare SRA directory 
# Here ~ refers to your home directory
sra_db_dir <- file.path("~", "shared-data", "SRAdb")
# Declare database file path
sqlite_file <- file.path(sra_db_dir, "SRAmetadb.sqlite")
# Declare output directory 
output_dir <- file.path("SRA_metadata")
# Declare output file path using the project ID
output_metadata_file <- file.path(output_dir, paste0(project, "_metadata.tsv"))
# Create the directory if it doesn't exist.
if (!dir.exists(output_dir)) {
  dir.create(output_dir,  recursive = TRUE)
}
# Make connection to sqlite database file
sra_con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_file)
study_table <- dplyr::tbl(sra_con, "study")

# Use the sqlite connection table to apply functions to it
sra_df <- sra_table %>% 
  # Filter to the study_id we are looking for
  dplyr::filter(study_accession == project) %>% 
  # Inner join to the study table that has more specific info about the study
  dplyr::inner_join(study_table, by = "study_accession") %>% 
  # We need to do this so the dbplyr queries are transformed to a data frame
  as.data.frame() 





# or download the old fashioned way 
library(readr)
sra.runs <- read_delim("C:/Users/A02342347/Downloads/SraRunTable.txt", delim= ",")

## much more info I need!
head(sra.runs)

## Scale counts by taking into account the total coverage per sample (ie. bp per sample?)
rse <- scale_counts(rse_gene)

## check the count table
assay(rse)[1:6, 1:20]


sra.runs$Run %in% colData(rse)$run %>% table()
colData(rse)$run  %in% sra.runs$Run %>% table()
# Which runs are not available in the count matrix?
sra.runs$Run[! sra.runs$Run %in% colData(rse)$run]

## there should be two same samples, one is deep sequenced, the other is shallow sequenced (downsampled)
colData(rse) %>% as.data.frame() %>% group_by(sample) %>% summarise(n=n()) %>% arrange(n) %>% head()

# How I would find sample names to filter... If I had them! No column for this in the txt that I downloaded
sra.runs %>% dplyr::filter(Run == "SRR1274250") %>% select(SRA_sample)
sra.runs %>% dplyr::filter(Run == "SRR1274344") %>% select(SRA_sample)

# filter these samples from the df
colData(rse) %>% as.data.frame() %>% dplyr::filter(sample == "SRS603223")
colData(rse) %>% as.data.frame() %>% dplyr::filter(sample == "SRS603270")


shallow.run<- colData(rse) %>% as.data.frame() %>% group_by(sample) %>% dplyr::slice(which.min(mapped_read_count)) %>% ungroup()

boxplot(shallow.run$mapped_read_count)
title(main = "number of mapped reads for the shallow sequenced samples")


## merge the meta data
shallow.run <- shallow.run %>% left_join(sra.runs, by= c("run" = "Run"))
## summary of mapped_read_count in shallow ones. median of  ~0.2 million reads
summary(shallow.run$mapped_read_count)


## exclude the shallow runs to find the deep runs
shallow.runs<- shallow.run$run
deep.run<- colData(rse) %>% as.data.frame() %>% dplyr::filter(! run %in% shallow.runs) %>% ungroup()
boxplot(deep.run$mapped_read_count)
title(main = "number of mapped reads for the deep sequenced samples")

## merge the meta data
deep.run <- deep.run %>% left_join(sra.runs, by = c("run" = "Run"))

shallow.run %>% mutate(cell_line = gsub("([^ ]+).+", "\\1", Population)) %>% select(cell_line) %>% table()


#It is 301, the same number of samples used in the paper,
# so I guess they excluded those k562 cells
# exclude those k562 cells to get only 301 single cells.
shallow.run <- shallow.run %>% mutate(cell_line = gsub("([^ ]+).+", "\\1", Population)) 
shallow.run<- shallow.run %>% filter(! cell_line %in% c("chip", "plate"))


## a peek into the count matrix
assay(rse)[1:6, 1:6]

head(rse)


## select out only the counts for shallow ones
shallow.counts<- assay(rse)[, colnames(rse) %in% shallow.run$run]
## I need the cell type, tissue type info for each run in the order of runs (each column is a run) in the matrix
shallow.run.cell<- left_join(data.frame(run = colnames(shallow.counts)), shallow.run) %>% select(run, Cell_Line, Tissue)


## tidy the dataframe a bit, set the NPC tissue from N/A to brain for comparisons with GW lines, and iPS tissue from N/A to pluripotent
shallow.run.cell<- shallow.run.cell %>% 
  mutate(Tissue = gsub("(mammary gland);.+", "\\1", Tissue)) %>% 
  mutate(Tissue = replace(Tissue, Tissue== "N/A" & Cell_Line == "iPS", "pluripotent")) %>% 
  mutate(Tissue = replace(Tissue, Tissue== "N/A" & Cell_Line == "NPC", "brain"))


## use log2 count for PCA and center the data for each sample.
X<- t(scale(t(log2(shallow.counts+1)),center=TRUE,scale=FALSE))
sv<- svd(t(X))
U<- sv$u
V<- sv$v
D<- sv$d
Z<- t(X)%*%V

## put PC1 and PC2 into a dataframe and add the cell line info side by side
pc_dat<- data.frame(cell.line = shallow.run.cell$Cell_Line, PC1 = Z[,1], PC2 = Z[,2])

## make figure with ggplot2
## NOT WORKING
library(ggplot2)
library(ggthemes)
ggplot(pc_dat, aes(x=PC1, y=PC2, col=cell.line)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(color="black", linewidth = 0.6),
        axis.line.y = element_line(color="black", linewidth = 0.6))


p + guides(col = guide_legend(ncol=2)) + theme(legend.position = c(.9, .75))



library(genefilter)
rv<- rowVars(X)
## select the top 250 most variable genes for clustering
idx<- order(-rv)[1:250]

library(ComplexHeatmap) ## a very nice package from Zuguang Gu, he is very responsive.
Heatmap(X[idx,], name = "log2 RNAseq\ncounts scaled", 
        show_row_names = FALSE, show_column_names = FALSE, 
        row_dend_reorder = TRUE, column_dend_reorder = TRUE)



# Changing clustering methods from single/average to complete linkage, as this outperforms the other two
# Changing gene clustering distance to pearson correlation, and sample distance to euclidian 
# Default is euclidian and complete in complexheatmap. NEVER USE SINGLE LINKAGE
Heatmap(X[idx,], name = "log2 RNAseq\ncounts scaled", 
        show_row_names = FALSE, show_column_names = FALSE, 
        row_dend_reorder = TRUE, column_dend_reorder = TRUE, 
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete")


# Now to add cell line and tissue type to the graphic 
df<- data.frame(cell.line = shallow.run.cell$Cell_Line, tissue = shallow.run.cell$Tissue)
## choose some colors to repesent cell lines and tissues
library(RColorBrewer)
## how many we need?
table(df$cell.line)
table(df$tissue)

##  6 for tissues and 11 for cell lines 
tissue.cols<- brewer.pal(6, "Dark2")
cell.cols<- brewer.pal(11, "Set3")

## make a named vector from two vectors
cell.cols.assigned<- setNames(cell.cols, unique(as.character(df$cell.line)))
tissue.cols.assigned<- setNames(tissue.cols, unique(as.character(df$tissue)))

## Heatmap annotation
ha<- HeatmapAnnotation(df = df, 
                       col = list(cell.line = cell.cols.assigned, 
                                  tissue = tissue.cols.assigned))

Heatmap(X[idx,], name = "log2 RNAseq\ncounts scaled", 
        show_row_names = FALSE, show_column_names = FALSE, 
        row_dend_reorder = TRUE, column_dend_reorder = TRUE, 
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        top_annotation = ha)


# Now to focus on brain cells 
brain.run<- shallow.run %>% filter(cell_line %in% c("NPC", "GW16", "GW21", "GW21+3"))
brain.counts<- assay(rse)[, colnames(rse) %in% brain.run$run]

Y<- t(scale(t(log2(brain.counts+1)),center=TRUE,scale=FALSE))

##SVD to get PCs
sv.brain<- svd(t(Y))
U.brain<- sv.brain$u
V.brain<- sv.brain$v
D.brain<- sv.brain$d
Z.brain<- t(Y)%*%V.brain


## x is the matrix D from SVD
variance_explained_each_PC<- function(x){
  var.list= list()
  varex = 0
  cumvar = 0
  denom = sum(x^2)
  for(i in 1:length(x)){
    varex[i] = x[i]^2/denom
    cumvar[i] = sum(x[1:i]^2)/denom
  }
  var.list$varex<- varex
  var.list$cumvar<- cumvar
  var.list
}

### screen plot 
screen.plot<- function(var.list){
  par(mfrow=c(1,2))
  plot(1:length(var.list$varex), var.list$varex *100,type="h",
       lwd=2,xlab="PC",ylab="% Variance Explained")
  plot(1:length(var.list$cumvar),var.list$cumvar,type="h",
       lwd=2,xlab="PC",ylab="Cummulative Variance Explained")
  
}

screen.plot(variance_explained_each_PC(D.brain))








## for PC1, PC2, and PC3, choose the maximum absulute loadings 
variations.3PC<- apply(abs(V.brain[,1:3]), 1, max)

## order according to the loadings and choose the top 500 genes
genes.3PC<- order(-variations.3PC)[1:500]

## cell type/tissue information
brain.run.cell<- shallow.run.cell %>% filter(cell_line %in% c("NPC", "GW16", "GW21", "GW21+3"))

df.brain<- data.frame(cell.line = brain.run.cell$cell_line)
## choose some colors to repesent cell lines
## how many we need?
table(df.brain$cell.line)




##  4 for cell lines, I select the color mannually http://www.color-hex.com/ to match the 
## color in the original figure

## olor space is important for interpolating colors. By default (in ComplexHeatmap pacakge), colors are linearly interpolated in LAB color space, but you can select the color space in colorRamp2() function.

brain.cell.cols<- c("#A6FFB6", "#33a1f6", "#d16303", "#322f64")
#brain.cell.cols<- brewer.pal(4, "Dark2")

## make a named vector from two vectors
brain.cell.cols.assigned<- setNames(brain.cell.cols, unique(as.character(df.brain$cell.line)))

## Heatmap annotation
ha.brain<- HeatmapAnnotation(df = df.brain, 
                             col = list(cell.line = brain.cell.cols.assigned))


library(circlize) ## for colorRamp2 color mapping
Heatmap(Y[genes.3PC, ], name = "log2 RNAseq\ncounts scaled", 
        col = colorRamp2(c(-10, 0, 10), c("Darkblue", "white", "red")),
        show_row_names = FALSE,
        show_column_names = FALSE, 
        row_dend_reorder = TRUE, column_dend_reorder = FALSE, 
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        top_annotation = ha.brain,
        gap = unit(0.5, "mm"))


library("PMA")
## we also look at the first 3 PCs, we centered the data before hand
## sumabsv  
## How sparse do you want v to be? This is the sum of absolute values of elements of v. 
# It must be between 1 and square root of number of columns (No. of genes in this case) of data. 
# The smaller it is, the sparser v will be.
spc<- SPC(t(Y),sumabsv=10, K=3, center = FALSE)


# how many genes selected? we look at the V matrix again, if the weights are zeros, 
# they are not important features. sparse PCA zero them out. For each of the four PCs, 
# how many features are retained?
apply(spc$v!=0, 2, sum)




## genes retained for the PC1, PC2 and PC3
PC1.genes<- which(spc$v[,1] !=0)
PC2.genes<- which(spc$v[,2] !=0)
PC3.genes<- which(spc$v[,3] !=0)

selected.genes<- unique(c(PC1.genes, PC2.genes, PC3.genes))
## total 493 genes for PC1, PC2 and PC3 were selected

## heatmap using these 493 genes, kmeans for genes (k=4), I visually selected for k=4 after testing around.

Heatmap(Y[selected.genes, ], name = "log2 RNAseq\ncounts scaled", 
        col = colorRamp2(c(-10, 0, 10), c("Darkblue", "white", "red")),
        show_row_names = FALSE,
        show_column_names = FALSE, 
        row_dend_reorder = TRUE, column_dend_reorder = TRUE, 
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        top_annotation = ha.brain,
        gap = unit(0.5, "mm"),
        km = 4)





## select genes based on variance across samples
rv.brain<- rowVars(Y)
## select the top 500 most variable genes for clustering
idx.brain<- order(-rv.brain)[1:500]
library(circlize)  ## for colorRamp2 function to map values to colors  

## kmeans for genes k=3
Heatmap(Y[idx.brain, ], name = "log2 RNAseq\ncounts scaled", 
        col = colorRamp2(c(-10, 0, 10), c("Darkblue", "white", "red")),
        show_row_names = FALSE,
        show_column_names = FALSE, 
        row_dend_reorder = TRUE, column_dend_reorder = TRUE, 
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D",
        top_annotation = ha.brain,
        gap = unit(0.5, "mm"),
        km = 3)
# It is not the same as the figure in the paper.
# Now, you probably understand how important the feature selection is
# (which genes we choose to generate heatmaps); how important to center the data or not; 
#and how important to choose clustering distance measures and linkage methods.


# Better color scheme with Viridis package 
library(viridis)

Heatmap(Y[genes.3PC, ], name = "log2 RNAseq\ncounts scaled", 
        col = colorRamp2(c(-10, 0, 10), viridis(3)),
        show_row_names = FALSE,
        show_column_names = FALSE, 
        row_dend_reorder = TRUE, column_dend_reorder = FALSE, 
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        top_annotation = ha.brain,
        gap = unit(0.5, "mm"))

