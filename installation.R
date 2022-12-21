
packages <- c("tidyverse", "plyr", "stringr", "ggfortify", "broom","ggthemes","ggrepel",
              "data.table", "factoextra", "cluster", "EnhancedVolcano", "UpsetR", "ComplexHeatmap",
              "conclust", "matrixStats", "circlize", "cluster", "RColorBrewer", "ggpubr",
              "devtools", "BiocManager", "lemon")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], dependencies = TRUE)
}
devtools::install_github("GoldfarbLab/TMTPurityCorrection")
devtools::install_github("GoldfarbLab/MSstatsHelper")
BiocManager::install("MSstats")
BiocManager::install("MSstatsTMT", version = "2.1")
BiocManager::install("EnhancedVolcano")
BiocManager::install("UpSetR")
BiocManager::install("ComplexHeatmap")
