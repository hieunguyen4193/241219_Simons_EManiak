gc()
rm(list = ls())

my_random_seed <- 42
set.seed(my_random_seed)

# __________VDJ DATA ANYLYSIS PIPELINE__________
PROJECT <- "241219_BSimons_EManiak"
# config.version <- "default"
config.version <- "v0.1"

# install scRepertoire v1
# devtools::install_github("ncborcherding/scRepertoire@v1")

if (config.version == "default"){
  version.name <- "SeuratV5"  
} else {
  version.name <- config.version
}

PROJECT.with.version <- sprintf("%s_%s", PROJECT, version.name)

source("/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_VDJ_pipeline/main_VDJ_pipeline.R")

outdir <- "/home/hieunguyen/CRC1382/outdir"

path.to.storage <- "/media/hieunguyen/HD01/storage"
path.to.main.input <- file.path(path.to.storage, PROJECT)

path.to.main.output <- file.path(outdir, PROJECT.with.version)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.VDJ.input <- file.path(path.to.main.input, "VDJ")
path.to.VDJ.output <- file.path(path.to.main.output, "VDJ_output")
dir.create(path.to.VDJ.output, showWarnings = FALSE, recursive = TRUE)

path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
path2src <- file.path(path.to.pipeline.src, "processes_src")

source(file.path(path2src, "import_libraries.R"))
source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))

summarize_vdj_data(path.to.input = path.to.VDJ.input, 
                   path.to.output = path.to.VDJ.output, 
                   PROJECT = PROJECT, 
                   removeNA=FALSE, 
                   removeMulti=FALSE, T_or_B = "T")

# __________GEX DATA ANALYSIS PIPELINE__________
path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/241219_Simons_EManiak"
source(file.path(path.to.project.src, "project_config.R"))

for (analysis.round in c("1st")){
  
  path2input <- file.path(path.to.main.input, "GEX")
  
  path.to.VDJ.output <- file.path(path.to.main.output, "VDJ_output")
  
  # _____stage lst for single sample.id_____
  stage_lst <- hash()
  
  stage_lst[["K1"]] <- c(K1 = "K1")
  stage_lst[["K2"]] <- c(K2 = "K2")
  stage_lst[["P1"]] <- c(P1 = "P1")
  stage_lst[["PS"]] <- c(PS = "PS")
  
  MINCELLS  <- 50
  MINGENES  <- 5
  
  save.RDS <- list(s1 = TRUE,
                   s2 = TRUE,
                   s3 = TRUE,
                   s4 = TRUE,
                   s5 = TRUE,
                   s6 = TRUE,
                   s7 = TRUE,
                   s8 = TRUE,
                   s8a = TRUE,
                   s9 = FALSE)
  
  sw <- list(s1 = "on",
             s2 = "on",
             s3 = "on",
             s4 = "on",
             s5 = "on",
             s6 = "on",
             s7 = "on",
             s8 = "off",
             s8a = "on",
             s9 = "off")
  
  rerun <- list(s1 = FALSE, 
                s2 = FALSE,
                s3 = FALSE,
                s4 = FALSE,
                s5 = FALSE,
                s6 = FALSE,
                s7 = FALSE,
                s8 = FALSE,
                s8a = FALSE,
                s9 = FALSE)
  
  filter.thresholds <- filter.config.params[[config.version]]
  
  remove_doublet <- FALSE
  path.to.10X.doublet.estimation <- file.path(path.to.storage, "DoubletEstimation10X.csv"
  )  
  #####--------------------------------------------------------------------#####
  ##### IMPORTANT INPUT PARAMS
  #####--------------------------------------------------------------------#####
  # for (sample.id in names(stage_lst)){
  for (sample.id in c("PS")){
    print(sprintf("Working on sample.id %s", sample.id))
    path.to.output <- file.path(path.to.main.output, sprintf("%s_round", analysis.round))
    dir.create(path.to.output, showWarnings = FALSE)
      path.to.anno.contigs <- NULL
      path.to.count.clonaltype <- NULL
      filtered.barcodes <- NULL
      path.to.s3a <- NULL
      input.method <- "normal"
    s.obj <- run_pipeline_GEX(path2src=path2src,
                              path2input=file.path(path2input, sample.id),
                              path.to.logfile.dir=file.path(path.to.output, sprintf("%s_%s_round", sample.id, analysis.round), "logs"),
                              stage_lst=stage_lst[[sample.id]],
                              path.to.10X.doublet.estimation=path.to.10X.doublet.estimation,
                              MINCELLS=MINCELLS,
                              MINGENES=MINGENES,
                              PROJECT=PROJECT.with.version,
                              remove_doublet=remove_doublet,
                              save.RDS=save.RDS,
                              path.to.output=file.path(path.to.output, sprintf("%s_%s_round", sample.id, analysis.round)),
                              rerun=rerun,
                              DE.test="wilcox",
                              num.PCA=num.PCA,
                              num.PC.used.in.UMAP=num.PC.used.in.UMAP,
                              num.PC.used.in.Clustering=num.PC.used.in.Clustering,
                              use.sctransform=use.sctransform,
                              filtered.barcodes=filtered.barcodes,
                              filter.thresholds=filter.thresholds,
                              path.to.anno.contigs=path.to.anno.contigs,
                              path.to.count.clonaltype=path.to.count.clonaltype,
                              input.method = input.method,
                              my_random_seed = my_random_seed,
                              with.VDJ = TRUE, 
                              path.to.s3a.source = path.to.s3a, 
                              regressOut_mode = regressOut_mode,
                              features_to_regressOut = features_to_regressOut,
                              sw = sw,
                              vars.to.regress = vars.to.regress,
                              cluster.resolution = cluster.resolution)
  }
}
#### ALWAYS REMEMBER TO SAVE SESSIONINFO !!!!!!
writeLines(capture.output(sessionInfo()), file.path(path.to.output, sprintf("%s_sessionInfo.txt", PROJECT)))

