library(Matrix)
library(readr)
library(data.table)
library(dplyr)
library(plyr)
library(HDCytoData)
library(flowCore)
library(CATALYST)
library(ggplot2)


# Specify directories
file_paths <- list(paste0(getwd(),"/fcs/"))

# Name file_paths
names(file_paths) <- file_paths

# Search for all files in the specified directories and extract files by a given extension
files_list            <- lapply(file_paths, list.files, recursive=T)
files_by_ext          <- lapply(files_list, function(x){x[endsWith(x, suffix=".fcs")]} )

# Get complete paths to all files
all_file_paths        <- unlist(lapply(seq_along(files_by_ext), function(x) {  paste(names(files_by_ext[x]), files_by_ext[[x]], sep="") } ))
names(all_file_paths) <- lapply(strsplit(all_file_paths,split="/"), function(x) { sub(".fcs","",x[length(x)]) } )
file_names <- unname(unlist(lapply(strsplit(unlist(files_by_ext),split = "/"),tail,1)))

#read the FCS files in the directory fcs/
data_list <- list()
fcs.par <- list()
for (i in c(1:length(file_names))) {
  data_list[[i]] <- read.FCS(as.character(all_file_paths[i]), alter.names = T)
  fcs.par[[i]] <- as.character(data_list[[i]]@parameters@data$name)
}

# find all the measured FACS parameters
common.par <- Reduce(intersect, fcs.par)

# get all the fluorescent markers measured during FACS analysis
common.par <- common.par[7:21]

#subset the FCS files for fluorescent markers measured during FACS analysis
dir.create(paste0(getwd(),"/subset/"))
sub_data_list <- list()
subset_file_names <- list()

for (i in c(1:length(file_names))) {
  sub_data_list[[i]] <- data_list[[i]][,common.par]
  write.FCS(sub_data_list[[i]] , paste(getwd(),"/subset/",paste(sub(".fcs","",file_names[[i]]),"_subset",".fcs", sep = ""), sep = ""))
  subset_file_names[[i]] <- paste(sub(".fcs","",file_names[[i]]),"_subset",".fcs", sep = "")
}

#read the subset FCS files
fs <- read.flowSet(path =paste0(getwd(),"/subset/"))

#create patient and sample annotations
md <- as.data.frame(as.character(unlist(subset_file_names)))
colnames(md) <- c("file_name")
md$patient_id <- read.table(text = md$file_name, sep = "_", as.is = TRUE)$V2
md$condition <-  sub("(\\_).+","", md$file_name)
md$sample_id <- md$condition

md$condition <- factor(md$condition, levels = unique(md$condition))
md$patient_id <- factor(md$patient_id, levels = unique(md$patient_id))

#create panel annotations
panel <- data.frame(as.character(sub_data_list[[1]]@parameters@data$name),as.character(sub_data_list[[1]]@parameters@data$desc))
rownames(panel) <- c(1:dim(panel)[1])
panel$marker_class <- c("type")
colnames(panel) <- c("fcs_colname",	"antigen", "marker_class")

# exclude the parameters not to be included in the dimensionality reduction
panel[grep("CD8|Via|CCR7|CD27|CD45RA|TET",panel$antigen),]$marker_class <- c("state")

# create single cell experiemnt
sce <- prepData(
  fs,
  panel = panel,
  md = md,
  features = panel$fcs_colname,
  transform = TRUE,
  cofactor = 150,
  by_time = F, FACS = T
)

# plot the normalized expression density for all markers
p <- plotExprs(sce, color_by = "sample_id")
p$facet$params$ncol <- 3
p

# plot the number of cells in each condition
plotCounts(sce, group_by = "sample_id", color_by = "condition")

# MDS plot for all conditions
pbMDS(sce, color_by = "condition", label_by = "sample_id")

# run and plot t-SNE on at most 200 cells per sample
set.seed(1234)
sce <- runDR(sce, "TSNE", features = "type", cells = 200)
plotDR(sce, "TSNE") + geom_point(size=1)

# plot scaled expression of the desired markers
cdx <- rownames(sce)[c(1,3:10,12:14)]
plotDR(sce, scale = T,color_by = cdx, ncol = 4, a_pal = rev(hcl.colors(10, "Spectral"))) + geom_point(size=0.5) 