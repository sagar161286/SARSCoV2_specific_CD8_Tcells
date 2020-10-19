library(Matrix)
library(readr)
library(data.table)
library(dplyr)
library(plyr)
library(HDCytoData)
library(flowCore)
library(CATALYST)

# Specify directories
file_paths <- list(paste0(getwd(),"/fcs_files/"))
# Name file_paths
names(file_paths) <- file_paths

# Search for all files in the specified directories and extract files by a given extension
files_list            <- lapply(file_paths, list.files, recursive=T)
files_by_ext          <- lapply(files_list, function(x){x[endsWith(x, suffix=".fcs")]} )

# Get complete paths to all files
all_file_paths        <- unlist(lapply(seq_along(files_by_ext), function(x) {  paste(names(files_by_ext[x]), files_by_ext[[x]], sep="") } ))
names(all_file_paths) <- lapply(strsplit(all_file_paths,split="/"), function(x) { sub(".fcs","",x[length(x)]) } )
file_names <- unname(unlist(lapply(strsplit(unlist(files_by_ext),split = "/"),tail,1)))

# Read all the FCS files present in the folder fcs/
data_list <- list()
fcs.par <- list()
for (i in c(1:length(file_names))) {
  data_list[[i]] <- read.FCS(as.character(all_file_paths[i]), alter.names = T)
  fcs.par[[i]] <- as.character(data_list[[i]]@parameters@data$name)
}
common.par <- Reduce(intersect, fcs.par)
common.par <- common.par[7:21]

dir.create(paste0(getwd(),"/subset/"))
sub_data_list <- list()
subset_file_names <- list()

for (i in c(1:length(file_names))) {
  sub_data_list[[i]] <- data_list[[i]][,common.par]
  write.FCS(sub_data_list[[i]] , paste(getwd(),"/subset/",paste(sub(".fcs","",file_names[[i]]),"_subset",".fcs", sep = ""), sep = ""))
  subset_file_names[[i]] <- paste(sub(".fcs","",file_names[[i]]),"_subset",".fcs", sep = "")
}

fs <- read.flowSet(path =paste0(getwd(),"/subset/"))

md <- as.data.frame(as.character(unlist(subset_file_names)))
colnames(md) <- c("file_name")
md$patient_id <- read.table(text = md$file_name, sep = "_", as.is = TRUE)$V2
md$condition <-  sub("(\\_).+","", md$file_name)
md$sample_id <- md$condition

md$condition <- factor(md$condition, levels = unique(md$condition))
md$patient_id <- factor(md$patient_id, levels = unique(md$patient_id))

panel <- data.frame(as.character(sub_data_list[[1]]@parameters@data$name),as.character(sub_data_list[[1]]@parameters@data$desc))
rownames(panel) <- c(1:dim(panel)[1])
panel$marker_class <- c("type")
colnames(panel) <- c("fcs_colname",	"antigen", "marker_class")

panel[grep("CD8",panel$antigen),]$marker_class <- c("state")
panel[grep("Via",panel$antigen),]$marker_class <- c("state")
panel[grep("CCR7",panel$antigen),]$marker_class <- c("state")
panel[grep("CD27",panel$antigen),]$marker_class <- c("state")
panel[grep("CD45RA",panel$antigen),]$marker_class <- c("state")
panel[grep("TET",panel$antigen),]$marker_class <- c("state")


sce <- prepData(
  fs,
  panel = panel,
  md = md,
  features = panel$fcs_colname,
  transform = TRUE,
  cofactor = 150,
  by_time = F, FACS = T
)

p <- plotExprs(sce, color_by = "sample_id")
p$facet$params$ncol <- 3
p

plotCounts(sce, group_by = "sample_id", color_by = "condition")
pbMDS(sce, color_by = "condition", label_by = "sample_id")

plotExprHeatmap(sce, scale = "last", hm_pal = rev(hcl.colors(10, "YlGnBu")))

plotNRS(sce, features = "type", color_by = "condition")

# run t-SNE/UMAP on at most 500/1000 cells per sample
set.seed(1234)
sce <- runDR(sce, "TSNE", features = "type", cells = 200)
sce <- runDR(sce, "UMAP", features = "type", cells = 200)

library(ggplot2)
plotDR(sce, "TSNE") + geom_point(size=1)
plotDR(sce, "UMAP")+ geom_point(size=1)

plotDR(sce, "TSNE", color_by = "PD-1")
plotDR(sce, "TSNE", color_by = "KLRG1")
plotDR(sce, "TSNE", color_by = "T-bet")
plotDR(sce, "TSNE", color_by = "TCF-1")
plotDR(sce, "TSNE", color_by = "CD127")
plotDR(sce, "TSNE", color_by = "TOX")
plotDR(sce, "TSNE", color_by = "Bcl-2")
plotDR(sce, "TSNE", color_by = "Eomes")
plotDR(sce, "TSNE", color_by = "IRF4")

cdx <- rownames(sce)[c(1,3:10,12:14)]
#cdx <- rownames(sce)[c(2:12)]
plotDR(sce, scale = T,color_by = cdx, ncol = 4, a_pal = rev(hcl.colors(10, "Spectral"))) + geom_point(size=0.5) 
plotDR(sce, "UMAP",scale = T,color_by = cdx, ncol = 3, a_pal = rev(hcl.colors(10, "Spectral"))) + geom_point(size=0.1) 

plotDR(sce, scale = F,color_by = "PD-1", ncol = 3, a_pal = rev(hcl.colors(10, "Spectral"))) + geom_point(size=2) 
plotDR(sce, scale = T,color_by = "TCF-1", ncol = 3, a_pal = rev(hcl.colors(10, "Spectral"))) + geom_point(size=2) 

plotDR(sce, scale = F,color_by = "TOX", ncol = 3, a_pal = rev(hcl.colors(10, "Spectral"))) + geom_point(size=2) 
plotDR(sce, scale = F,color_by = "TOX", ncol = 3, a_pal = rev(hcl.colors(10, "Spectral"))) + geom_point(size=2) 
