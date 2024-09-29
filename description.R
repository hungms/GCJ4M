pkgs <- c(
    "tidyverse", "ComplexHeatmap", "cowplot", "viridis",
    "ggpubr", "ggrepel", "readxl", "writexl", "patchwork")

for(x in pkgs){
    usethis::use_package(x)}