##### clean up #####
gc()
rm(list = ls())

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
dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)

path.to.s.obj <- file.path(path.to.02.output, input.case, "s8_output", "241219_BSimons_EManiak_v0.1.output.s8.rds")
s.obj <- readRDS(path.to.s.obj)

#####-------------------------------------------------------------------------------#####
##### ERROR: in previous run of "run_pipeline.R", forgot to set to "B" cells, 
##### correct now, rerun.
##### Remember to set T_or_B = "B" for this data, B cells. 
#####-------------------------------------------------------------------------------#####
path.to.VDJ.input <- file.path(path.to.storage, "241219_BSimons_EManiak", "VDJ")
path.to.VDJ.output <- file.path(path.to.main.output, "VDJ_output_correction")
dir.create(path.to.VDJ.output, showWarnings = FALSE, recursive = TRUE)

path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
path2src <- file.path(path.to.pipeline.src, "processes_src")

source("/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_VDJ_pipeline/main_VDJ_pipeline.R")

source(file.path(path2src, "import_libraries.R"))
source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))

if (file.exists(file.path(path.to.VDJ.output, "corrected.csv")) == FALSE){
  summarize_vdj_data(path.to.input = path.to.VDJ.input, 
                     path.to.output = path.to.VDJ.output, 
                     PROJECT = PROJECT, 
                     removeNA=FALSE, 
                     removeMulti=FALSE, 
                     T_or_B = "B")
  write.csv(data.frame(status = c("finished correcting VDJ output")), file.path(path.to.VDJ.output, "corrected.csv"))
} else {
  print(sprintf("reading in VDJ corrected data from %s", path.to.VDJ.output))
}

#####-------------------------------------------------------------------------------#####
##### get all VDJ files
#####-------------------------------------------------------------------------------#####
all.vdj.files <- Sys.glob(file.path(path.to.VDJ.output, "*annotated_contigs_clonaltype*.csv"))
vdjdf <- data.frame()
for (input.file in all.vdj.files){
  tmpdf <- read.csv(input.file)
  if ("X" %in% colnames(tmpdf)){
    tmpdf <- subset(tmpdf, select = -c(X))
  }
  tmpdf$SampleID <- str_replace(str_replace(basename(input.file), "annotated_contigs_clonaltype_", ""), ".csv", "")
  tmpdf <- tmpdf %>% rowwise() %>%
    mutate(barcode = str_replace(barcode, sprintf("%s_%s", SampleID, SampleID), SampleID))
  vdjdf <- rbind(vdjdf, tmpdf)
}

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
meta.data <- merge(meta.data, vdjdf, by.x = "barcode", by.y = "barcode", all.x = TRUE)

meta.data <- meta.data %>% column_to_rownames("barcode")

meta.data <- meta.data[row.names(s.obj@meta.data),]

old.cols <- colnames(s.obj@meta.data)
new.cols <- setdiff(colnames(meta.data), c(old.cols, "barcode"))

for (c in new.cols){
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data[[c]], col.name = c)
}

all.samples <- unique(s.obj$name)
#####-------------------------------------------------------------------------------#####
##### clonal data analysis
#####-------------------------------------------------------------------------------#####
meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
clonedf <- data.frame(CloneID = unique(meta.data$CTstrict))
for (sample.id in unique(s.obj$name)){
  clonedf[[sprintf("count_%s", sample.id)]] <- unlist(lapply(
    clonedf$CloneID, function(x){
      return(subset(meta.data, meta.data$SampleID == sample.id & meta.data$CTstrict == x) %>% nrow())
    }
  ))
}
clonedf <- clonedf %>% rowwise() %>%
  mutate(total_count = count_K1 + count_K2 + count_P1 + count_PS) %>% 
  subset(is.na(CloneID) == FALSE) %>% arrange(desc(total_count))

writexl::write_xlsx(clonedf, file.path(path.to.04.output, "clonedf.xlsx"))

#####-------------------------------------------------------------------------------#####
##### plot clones on UMAP
#####-------------------------------------------------------------------------------#####
path.to.save.clone.UMAP <- file.path(path.to.04.output, "clones_on_UMAP")
dir.create(path.to.save.clone.UMAP, showWarnings = FALSE, recursive = TRUE)

top.clones <- subset(clonedf, clonedf$total_count >= 10)

for (clone.id in unique(top.clones$CloneID)){
  savename <- str_replace_all(clone.id, ":", "___")
  available.samples <- subset(clonedf, clonedf$CloneID == clone.id)
  plot.barcodes <- list()
  for (sample.id in all.samples){
    if (available.samples[[sprintf("count_%s", sample.id)]] != 0){
      plot.barcodes[[sample.id]] <- subset(meta.data, meta.data$CTstrict == clone.id & meta.data$name == sample.id)$barcode
    } 
  }
  p <- DimPlot(object = s.obj, reduction = "cca_UMAP", label = TRUE, label.box = TRUE, group.by = "cca.cluster.0.5", 
               cells.highlight = plot.barcodes, cols.highlight = hue_pal()(4), sizes.highlight = 2)
  ggsave(plot = p, filename = sprintf("UMAP_%s.svg", savename), path = path.to.save.clone.UMAP, device = "svg", width = 14, height = 10, dpi = 300)
}

# EOF