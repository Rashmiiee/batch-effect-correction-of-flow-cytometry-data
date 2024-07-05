setwd("~/Desktop/Flow_Cytometry_Trail/Scripts_and_Data_for_MDS and markerintensity_plots/T cell panel/Batch_number_in_Stimulated")
library(readxl)
library(flowCore)
md <- read_excel("metadata_file.xlsx")
md$condition <- factor(md$bunch_id,levels = c("1","2","3","4","5","6"))
color_conditions <- c("#6A3D9A", "#FF7F00","#9a403d","#549a3d","#27b5f2","#f2f241")
names(color_conditions) <- levels(md$bunch_id)
fcs_filename <- "fcs_files.zip"
unzip(fcs_filename)
head(md)
fcs_raw <- read.flowSet(md$filename, transformation = FALSE, truncate_max_range = FALSE)
panel_filename <- "panel.xlsx"
panel <- read_excel(panel_filename)
head(data.frame(panel))
panel$Antigen <- gsub("-", "_", panel$Antigen)
panel_fcs <- pData(parameters(fcs_raw[[1]]))
head(panel_fcs)
panel_fcs$desc
parameters(fcs_raw[[1]])
summary(fcs_raw[1])
panel_fcs$desc <- gsub("-", "_", panel_fcs$desc)
(lineage_markers <- panel$Antigen[panel$Lineage == 1])
(functional_markers <- panel$Antigen[panel$Functional == 1])
all(lineage_markers %in% panel_fcs$desc)
all(functional_markers %in% panel_fcs$desc)
fcs <- fsApply(fcs_raw, function(x, cofactor = 150){
  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  expr <- asinh(expr[, c(lineage_markers, functional_markers)] / cofactor)
  exprs(x) <- expr
  x
})
fcs
expr <- fsApply(fcs, exprs)
expr <- expr[, -2]
expr <- expr[, -6]
library(matrixStats)
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1
hist(expr)
View(expr)
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
library(ggplot2)
library(reshape2)

ggdf <- data.frame(sample_id = sample_ids, expr)
ggdf <- melt(ggdf, id.var = "sample_id", 
             value.name = "expression", variable.name = "antigen")
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$bunch_id[mm]

ggplot(ggdf, aes(x = expression, color = as.factor(condition), 
                 group = sample_id)) +
  geom_density() +
  facet_wrap(~ antigen, nrow = 4, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  scale_color_manual(values = color_conditions)
library(dplyr)

expr_median_sample_tbl <- data.frame(sample_id = sample_ids, expr) %>%
  group_by(sample_id) %>% 
  summarize_all(funs(median))

expr_median_sample <- t(expr_median_sample_tbl[, -1])
colnames(expr_median_sample) <- expr_median_sample_tbl$sample_id

library(limma)
mds <- plotMDS(expr_median_sample, plot = FALSE)

library(ggrepel)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y, 
                   sample_id = colnames(expr_median_sample))
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$bunch_id[mm]

ggplot(ggdf, aes(x = MDS1, y = MDS2, color = as.factor(condition))) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = sample_id)) +
  theme_bw() +
  scale_color_manual(values = color_conditions) 
