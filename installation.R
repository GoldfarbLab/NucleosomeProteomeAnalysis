
packages <- c("tidyverse", "plyr", "stringr", "MSstatsTMT", "ggfortify", "broom","ggthemes","ggrepel",
              "data.table", "factoextra", "cluster", "MSstats", "EnhancedVolcano", "UpsetR", "ComplexHeatmap",
              "conclust","ComplexHeatmap", "matrixStats", "circlize", "cluster", "RColorBrewer", "ggpubr")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
library(devtools)
install_github("GoldfarbLab/TMTPurityCorrection")
install_github("GoldfarbLab/MSstatsHelper")
