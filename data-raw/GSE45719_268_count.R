## code to prepare `GSE45719_268_count` dataset goes here

# setwd("/Users/niexiner/Desktop/")

# BiocManager::install("GEOquery")
library(GEOquery)
# Download raw data.
getGEOSuppFiles(GEO = "GSE45719",
                makeDirectory = TRUE,
                baseDir = ".",
                fetch_files = TRUE,
                filter_regex = NULL)
# info files
data_info <- getGEO(GEO = "GSE45719")
setwd("GSE45719")
# save(data_info, file = "GSE45719_information.RData")
str(data_info)

untar("GSE45719_RAW.tar")
file.remove("GSE45719_RAW.tar")
dir()


# load("GSE45719_information.R")
title <- data_info$GSE45719_series_matrix.txt.gz@phenoData@data$title
title <- as.character(title)
file.dir <- dir() #

# data_new <- read.table(file = file.dir[1], header = FALSE)
# data_agg <- aggregate(data_new[, c(1, 4)], # 4th:count
#                       by=list(Gene_symbol = data_new$V1), mean) # Average the rows of duplicate gene names
# data_tmp <- data.frame("Gene_symbol" = data_agg$Gene_symbol, data_agg$V4)
# colnames(data_tmp)[2] <- title[1]
# data_comb <- data_tmp

# all
for (f in 2:length(title)){
  data_new <- read.table(file = file.dir[f], header = FALSE)

  data_agg <- aggregate(data_new[, c(1, 4)],
                        by = list(Gene_symbol = data_new$V1), mean)

  data_tmp <- data.frame("Gene_symbol" = data_agg$Gene_symbol, data_agg$V4)

  colnames(data_tmp)[2] <- title[f]

  data_comb <- merge(data_comb, data_tmp, by = "Gene_symbol", all = TRUE)
}

dim(data_comb)
rownames(data_comb) <- data_comb[, 1]
data_comb <- data_comb[, -1]

# setwd("..")
# save(data_comb, file = "data_combination_317.RData")

# Deleting
# "liver"
# "twocell"
# "split"
# "fibroblast"
library(dplyr)
GSE45719_268_count <- select(data_comb, -matches(".liver.")) %>%
  select(-matches("twocell")) %>%
  select(-matches("split")) %>%
  select(-matches("fibroblast"))
# colnames(mydata)
# save(mydata, file = "GSE45719_268_count.RData")
GSE45719_268_count <- as.matrix(GSE45719_268_count[,c(1, 55, 67, 94, 118, 145, 154, 184, 199, 256)])

write.table(GSE45719_268_count, file = "data-raw/GSE45719_268_count.txt")
usethis::use_data(GSE45719_268_count, overwrite = TRUE)
usethis::use_build_ignore("data-raw/GSE45719_268_count.R")

