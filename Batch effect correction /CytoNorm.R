setwd("~/Desktop/Flow_Cytometry_Trail/Scripts_and_Data_for_MDS and markerintensity_plots/T cell panel/CytoNorm")
library(CytoNorm)
library(flowCore)
library(FlowSOM)
dir <- system.file("extdata", package = "CytoNorm")
files <- list.files(dir, pattern = "fcs$")
Batch_pre = stringr::str_match(files, "[_][0-9][_]")[, 1]
library(tidyr)
Batch = readr::parse_number(Batch_pre)

data <- data.frame(File = files,
                   Path = file.path(dir, files),
                   Type = stringr::str_match(files, "([Ns])")[,2],
                   Batch,
                   stringsAsFactors = FALSE)
data$Type <- c("N" = "Train", "s" = "Validation")[data$Type]

train_data <- dplyr::filter(data, Type == "Train")
validation_data <- dplyr::filter(data, Type == "Validation")

ff <- flowCore::read.FCS(data$Path[1],truncate_max_range = FALSE)
channels <- flowCore::colnames(ff)
transformList <- flowCore::transformList(channels,arcsinhTransform(transformationId = "defaultArchsinTransform",a =1, b =1,c=0))
transformList.reverse <- flowCore::transformList(channels,cytofTransform.reverse)

fsom <- prepareFlowSOM(train_data$Path,
                       channels,
                       nCells = 6000,
                       FlowSOM.params = list(xdim = 5,
                                             ydim = 5,
                                             nClus = 10,
                                             scale = FALSE),
                       transformList = transformList,
                       seed = 1)

fsom2 <- fsom[1]$FlowSOM
cvs <- testCV(fsom2,cluster_values = c(5, 10, 15))

model <- CytoNorm.train(files = train_data$Path,
                        labels = train_data$Batch,
                        channels = channels,
                        transformList = transformList,
                        FlowSOM.params = list(nCells = 8000, 
                                              xdim = 5,
                                              ydim = 5,
                                              nClus = 10,
                                              scale = FALSE),
                        normMethod.train = QuantileNorm.train ,
                        normParams = list(nQ = 101,
                                          goal = "mean"),
                        seed = 1,
                        verbose = TRUE)

CytoNorm.normalize(model = model,
                   files = validation_data$Path,
                   labels = validation_data$Batch,
                   transformList = transformList,
                   transformList.reverse = transformList.reverse,
                   normMethod.normalize = QuantileNorm.normalize,
                   outputDir = "Normalized",
                   prefix = "Norm_",
                   clean = TRUE,
                   verbose = TRUE)
          