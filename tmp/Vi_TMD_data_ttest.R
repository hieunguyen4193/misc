gc()
rm(list = ls())

library(comprehenr)
library(stringr)
library(dplyr)
library(tidyverse)
library(limma)

path.to.main.input <- "/media/tmp/Vi_TMD_data"
path.to.main.output <- file.path(path.to.main.input, "output")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

feature.name <- "EM"

path.to.save.output <- file.path(path.to.main.output, "output_Vi", feature.name)
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.save.output, "features"), showWarnings = FALSE, recursive =  TRUE)
dir.create(file.path(path.to.save.output, "test_result"), showWarnings = FALSE, recursive =  TRUE)

files <- Sys.glob(file.path(path.to.main.output, "features", sprintf("feature_%s", feature.name), "*"))
for (input.file in files){
  print(basename(input.file))
  region.name <- str_replace(basename(input.file), ".csv", "")
  testdf <- read.csv(input.file)
  group1 <- subset(testdf, testdf$Class == "Benign")$SampleID
  group2 <- subset(testdf, testdf$Class == "BC")$SampleID
  
  testdf <- testdf %>% 
    subset(select = -c(X)) %>% column_to_rownames("SampleID") %>% subset(select = -c(Class)) %>% t()
  
  resdf <- data.frame(feature = row.names(testdf))
  resdf$pvalue <- unlist(lapply(
    resdf$feature, function(x){
      data1 <- testdf[x, group1]
      data2 <- testdf[x, group2]
      data1 <- data1[is.na(data1) == FALSE]
      data2 <- data2[is.na(data2) == FALSE]
      if (length(data1) >= 10 & length(data2) >= 10){
        res <- t.test(data1, data2)
        output <- res$p.value        
      } else {
        output <- NA
      }
      return(output)    
    }
  ))
  
  resdf$pvalue.wilcox <- unlist(lapply(
    resdf$feature, function(x){
      data1 <- testdf[x, group1]
      data2 <- testdf[x, group2]
      data1 <- data1[is.na(data1) == FALSE]
      data2 <- data2[is.na(data2) == FALSE]
      if (length(data1) >= 10 & length(data2) >= 10){
        res <- wilcox.test(data1, data2)
        output <- res$p.value        
      } else {
        output <- NA
      }
      return(output)    
    }
  ))
  
  resdf$adj.pvalue <- p.adjust(resdf$pvalue, method = "BH")
  resdf$adj.pvalue.wilcox <- p.adjust(resdf$pvalue.wilcox, method = "BH")
  write.csv(resdf, file.path(path.to.save.output, "test_result", sprintf("%s.csv", region.name)))
  write.csv(testdf, file.path(path.to.save.output, "features", sprintf("%s.csv", region.name)))
}

