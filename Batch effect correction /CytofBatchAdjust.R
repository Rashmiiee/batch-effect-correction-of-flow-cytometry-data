setwd("~/Desktop/Flow_Cytometry_Trail/Scripts_and_Data_for_MDS and markerintensity_plots/T cell panel/CytofBatchAdjust_try")
library("flowCore")
source("BatchAdjust.R")
BatchAdjust(
  basedir=".",
  outdir="./Results_2",
  channelsFile = "ChannelsToAdjust.txt",
  batchKeyword="T_live_",
  anchorKeyword = "Ns",
  method="95p",
  transformation=FALSE,
  addExt=NULL,
  plotDiagnostics=TRUE)