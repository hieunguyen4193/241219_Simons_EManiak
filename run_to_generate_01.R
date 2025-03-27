gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/241219_Simons_EManiak"
path.to.rmd <- file.path(path.to.main.src, "01_downstream_analysis_single_sample.Rmd")

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "241219_BSimons_EManiak_v0.1"
output.version <- "20250325"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.save.html <- file.path(path.to.main.output, "html_outputs", output.version)
dir.create(file.path(path.to.save.html, "01_output"), showWarnings = FALSE, recursive = TRUE)

for (sample.id in c("K1", "K2", "P1", "PS")){
  if (file.exists(file.path(path.to.save.html, 
                            "01_output", 
                            sprintf("01_downstream_analysis_single_sample_%s.html", 
                                    sample.id))) == FALSE){
    rmarkdown::render(input = path.to.rmd, 
                      params = list(
                        outdir = outdir,
                        PROJECT = PROJECT,
                        sample.id = sample.id
                      ),
                      output_dir = file.path(path.to.save.html, "01_output"),
                      output_file = sprintf("01_downstream_analysis_single_sample_%s.html", 
                                            sample.id))      
  } else {
    print(sprintf("File %s exists", 
                  file.path(path.to.save.html, 
                            "01_output", 
                            sprintf("01_downstream_analysis_single_sample_%s.html", 
                                    sample.id))))
  }
}
