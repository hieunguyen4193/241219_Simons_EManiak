gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5/processes_src"

source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering_SeuratV5.selectedGenes.R"))

options(future.globals.maxSize = 10000 * 1024^2)

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

PROJECT <- "241219_BSimons_EManiak_v0.1"

outdir <- "/home/hieunguyen/CRC1382/outdir"
path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.02.output <- file.path(path.to.main.output, "02_output")
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

path.to.input.files <- file.path(path.to.main.input, 
                                 sprintf("%s_round", analysis.round), 
                                 sprintf("%s_1st_round", "*"), "s8a_output", 
                                 sprintf("%s.output.s8a.rds", PROJECT))
all.input.files <- Sys.glob(path.to.input.files)
names(all.input.files) <- to_vec(
  for (item in all.input.files){
    str_split(item, "/")[[1]][[8]] %>% str_replace("_1st_round", "")
  }
)

if (file.exists(file.path(path.to.02.output, sprintf("raw_merge_dataset.rds"))) == FALSE){
  data.list <- list()
  for (i in seq(length(all.input.files))){
    sample.id <- names(all.input.files)[[i]]
    print(sprintf("reading in sample %s", sample.id))
    data.list[[i]] <- readRDS(all.input.files[[i]])
  }
  
  s.obj <- merge(data.list[[1]], data.list[2: length(all.input.files)])
  saveRDS(s.obj, file.path(path.to.02.output, sprintf("raw_merge_dataset.rds")))
} else {
  print("merge data object exists, reading in ...")
  s.obj <- readRDS(file.path(path.to.02.output, sprintf("raw_merge_dataset.rds")))
}

vars.to.regress <- c("percent.mt")

DefaultAssay(s.obj) <- "RNA"
s.obj <- JoinLayers(s.obj)

BR_genes_patterns <- c("Ighv", "Ighd", "Ighj", "Ighc", "Igkv",
                       "Igkj", "Igkc", "Iglv", "Iglj", "Iglc") %>% toupper()
BCRgenes.to.exclude <- unlist(lapply(row.names(s.obj), function(x){
  if (substr(x, 1, 4) %in% BR_genes_patterns){
    return(x)
  } else {
    return(NA)
  }
}))
BCRgenes.to.exclude <- subset(BCRgenes.to.exclude, is.na(BCRgenes.to.exclude) == FALSE)

dir.create(file.path(path.to.02.output, "with_IGgenes"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.02.output, "without_IGgenes"), showWarnings = FALSE, recursive = TRUE)

s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = s.obj, 
                                                     save.RDS.s8 = TRUE,
                                                     path.to.output = file.path(path.to.02.output, "with_IGgenes"),
                                                     use.sctransform = TRUE,
                                                     num.PCA = num.PCA,
                                                     num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                     num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                     cluster.resolution = cluster.resolution,
                                                     vars.to.regress = vars.to.regress,
                                                     remove.genes = NULL)

s.obj.integrated.noIG <- s8.integration.and.clustering_V5(s.obj = s.obj, 
                                                     save.RDS.s8 = TRUE,
                                                     path.to.output = file.path(path.to.02.output, "without_IGgenes"),
                                                     use.sctransform = TRUE,
                                                     num.PCA = num.PCA,
                                                     num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                     num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                     cluster.resolution = cluster.resolution,
                                                     vars.to.regress = vars.to.regress,
                                                     remove.genes = BCRgenes.to.exclude)


