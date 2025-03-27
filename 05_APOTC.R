gc()
rm(list = ls())
# install.packages("dplyr")
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-1.tar.gz", type = "source", repos = NULL)
new.pkgs <- c("APackOfTheClones", "svglite", "car", "ggpubr", "ggthemes", "dplyr")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }
}
library(ggpubr)
library(ggthemes)
library(APackOfTheClones)
library("gridExtra")

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))
#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25
regressOut_mode <- NULL
features_to_regressOut <- NULL
use.sctransform <- TRUE
vars.to.regress <- c("percent.mt")
cluster.resolution <- 0.5
analysis.round <- "1st"

# input.case <- "with_IGgenes"
input.case <- "without_IGgenes"

PROJECT <- "241219_BSimons_EManiak_v0.1"

outdir <- "/home/hieunguyen/CRC1382/outdir"

path.to.storage <- "/media/hieunguyen/HD01/storage"

path.to.main.input <- file.path(outdir, PROJECT)

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output", input.case)
path.to.04.output <- file.path(path.to.main.output, "04_output", input.case)
path.to.05.output <- file.path(path.to.main.output, "05_output", input.case)
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)

colors <- tableau_color_pal(palette = "Tableau 20")(20)

path.to.s.obj <- file.path(path.to.04.output, "241219_BSimons_EManiak_v0.1.output.s8.addedVDJ.rds")
s.obj <- readRDS(path.to.s.obj)
clonedf <- readxl::read_excel(file.path(path.to.04.output, "clonedf.xlsx"))

top.clones <- subset(clonedf, clonedf$total_count >= 10)
plot.clonedf <- head(top.clones, 10)
plot.clonedf$color <- head(colors, nrow(plot.clonedf))

Idents(s.obj) <- "cca.cluster.0.5"
clone.name <- "CTstrict"
reduction.name <- "cca_UMAP"

s.obj <- RunAPOTC(seurat_obj = s.obj, 
                  reduction_base = reduction.name, 
                  clonecall = clone.name)
tmp.plot <- vizAPOTC(s.obj, clonecall = clone.name,
                     verbose = FALSE,
                     reduction_base = reduction.name,
                     show_labels = TRUE,
                     repulsion_strength = 5,
                     legend_position = "top_right",
                     legend_sizes = 2) %>%
  showCloneHighlight(clonotype =  as.character(plot.clonedf$CloneID),
                     fill_legend = TRUE,
                     color_each = plot.clonedf$color,
                     default_color = "lightgray")
