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

inputdf <- read.csv(file.path(path.to.main.input, sprintf("TMfeature_%sx450_1D_full_Class.csv", feature.name)))

all.col.names <- setdiff(colnames(inputdf), c("SampleID", "Class"))

region.names <- to_vec( for(item in all.col.names) str_split(item, "_")[[1]][[2]]) %>% unique()

dir.create(file.path(path.to.main.output, "features", sprintf("feature_%s", feature.name)), showWarnings = FALSE, recursive = TRUE)
split.col.names <- list()
for (region in region.names){
  if (file.exists(file.path(path.to.main.output, sprintf("region_%s.csv", region))) == FALSE){
    print(sprintf(region))
    split.col.names[[region]] <- c(c("SampleID", "Class", to_vec (for (item in all.col.names) if (grepl(region, item) == TRUE) item)))
    write.csv(inputdf[, split.col.names[[region]]], file.path(path.to.main.output, 
                                                              "features", 
                                                              sprintf("feature_%s", feature.name),  
                                                              sprintf("region_%s.csv", region)))    
  }
}

dir.create(file.path(path.to.main.output, "test_res", sprintf("feature_%s", feature.name)), showWarnings = FALSE, recursive = TRUE)

for (region in region.names){
  
  testdf <- read.csv(file.path(path.to.main.output, "features", sprintf("feature_%s", feature.name), sprintf("region_%s.csv", region))) 
  if (dim(testdf[complete.cases(testdf),]) >= 1){
    group1 <- subset(testdf, testdf$Class == "Benign")$SampleID
    group2 <- subset(testdf, testdf$Class == "BC")$SampleID
    testdf <- testdf %>% 
      subset(select = -c(X)) %>% column_to_rownames("SampleID") %>% subset(select = -c(Class)) %>% t()
    
    
    ##### test
    input.metadata <- data.frame(sample = c(group1, group2), 
                                 label = c(
                                   to_vec(for(item in seq(1, length(group1))) "group1"),
                                   to_vec(for(item in seq(1, length(group2))) "group2")
                                 ))
    
    # this is the factor of interest
    g <- factor(input.metadata$label, levels = c("group1", "group2"))
    # use the above to create a design matrix
    design <- model.matrix(~0+label, data=input.metadata)
    colnames(design) <- levels(g)
    fit <- lmFit(testdf[, input.metadata$sample], design)
    # create a contrast matrix for specific comparisons
    contMatrix <- makeContrasts(group1-group2,
                                levels=design)
    fit2 <- contrasts.fit(fit, contMatrix)
    fit2 <- eBayes(fit2)
    
    res <- topTable(fit2, num=Inf, coef=1)
    
    writexl::write_xlsx(as.data.frame(res), file.path(path.to.main.output, 
                                                      "test_res", 
                                                      sprintf("feature_%s", feature.name), 
                                                      sprintf("region_%s.xlsx", region)))
  }
}
