library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(readxl)
library(maftools)

# ReadVariant: function to read a variant file(in csv format)
ReadVariant <- function(file_dir) {
    tumor_Sample_Barcode <- file_dir[1] %>% str_extract("P\\d+(-LN-\\d)?")
    total_variant <- read_csv(file_dir[1], show_col_types = F) %>%
        mutate(Tumor_Sample_Barcode = tumor_Sample_Barcode) %>%
        select(-starts_with(c("Otherinfo", "VAF")))
    return(total_variant)
}


# ReadGISTIC: function to read GISTIC copy number results
ReadGISTIC <- function(file_dir, Clinical_Data) {
    gistic_all_lessions <- paste(file_dir[1], "all_lesions.conf_95.txt", sep = "/")
    gistic_amp_genes <- paste(file_dir[1], "amp_genes.conf_95.txt", sep = "/")
    gistic_del_genes <- paste(file_dir[1], "del_genes.conf_95.txt", sep = "/")
    gistic_score_file <- paste(file_dir[1], "scores.gistic", sep = "/")
    gistic_res <- readGistic(
        gistic_all_lessions, gistic_amp_genes,
        gistic_del_genes, gistic_score_file
    )
    total_variant_maf_add_gistic <- read.maf("C:/Users/ydc2020/Desktop/Zhang Lab/total_variant.maf", gisticAllLesionsFile = gistic_all_lessions, gisticAmpGenesFile = gistic_amp_genes, gisticDelGenesFile = gistic_del_genes, gisticScoresFile = gistic_score_file, clinicalData = Clinical_Data)
    return(total_variant_maf_add_gistic)
}

# Read all data ------------------------------------------------------------

setwd("C:\\Users\\ydc2020\\Desktop\\Zhang Lab")
variant_files <- list.files("./Added_VAFs/", pattern = "*.csv")
variant_dir <- paste("./Added_VAFs/", variant_files, sep = "")
total_variant <- map(variant_dir, ReadVariant) %>% purrr::reduce(rbind)

Clinical_Data <- read_xlsx("C:\\Users\\ydc2020\\Clinical.xlsx")

Clinical_Data <- Clinical_Data %>% mutate(Tumor_Sample_Barcode = str_extract(SampleID, "P[0-9]+"))
write_delim(total_variant, file = "total_variant.txt", delim = "\t")
total_variant_maf <- annovarToMaf("total_variant.txt", refBuild = "hg38", tsbCol = "Tumor_Sample_Barcode", basename = "total_variant", MAFobj = T, sampleAnno = Clinical_Data)

TMB <- tmb(total_variant_maf, captureSize = 18.77)
TMB <- TMB %>%
  right_join(Clinical_Data, by = c(Tumor_Sample_Barcode = "Tumor_Sample_Barcode")) %>%
  select(total_perMB, Tumor_Sample_Barcode, total) %>%
  distinct()
  
TMB %>% ggbarplot(x = "Tumor_sample_barcode", y = "total_perMB", sort.val = "desc", palette = "npg", fill = "response", xlab = "Patient", ylab = "TMB", x.text.angle = 90, sort.by.groups = F) +
  geom_hline(yintercept = 10, linetype = "dotted")