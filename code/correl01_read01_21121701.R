# Include Libraries -------------------------------------------------------
library(ggplot2)
library(FactoMineR)
library(corrplot)
library(xlsx)

Sys.setlocale(category = "LC_ALL", locale = "Russian, Russia")

source("./code/read01_aa_clin_21121601.R")

rm(list = setdiff(ls(), c("ms_clinical_data", 
                          "ms_clinical_names")))

# Functions ---------------------------------------------------------------

# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = 'spearman', ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

chi_square_fun <- function(chi_square_data, v1, v2){
  result <- chisq.test(unlist(chi_square_data[v1]), 
                       unlist(chi_square_data[v2]), 
                       simulate.p.value = T, correct = T)
  return(result$p.value)
  
}

cv.test = function(chi_square_data, v1, v2) {
  x <- unlist(chi_square_data[v1])
  y <- unlist(chi_square_data[v2])
  CV = sqrt(chisq.test(x, y, correct=T)$statistic /
              (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  return(as.numeric(CV))
}


# Clean Data --------------------------------------------------------------

path_score <- ms_clinical_data %>%
  select(c194, c195, c196, c197,	c198,	c199,	c200,	c201,	c202,
         c203, c204,	c205,	c206,	c207,	c208,	c209,	c210,	c211,	c212,	
         c213, c214,	c215,	c216,	c217,	c218) %>%
  rowSums()

ms_data_aa_pneu <- ms_clinical_data %>%
  mutate(path_score = path_score) %>%
  select(-c183, -c184, -c194, -c195, -c196, -c197, -c198,	-c199,	-c200,	-c201,	-c202,
         -c203, -c204, -c205,	-c206, -c207,	-c208,	-c209,	-c210,	-c211,	-c212,	
         -c213, -c214, -c215,	-c216, -c217,	-c218) %>%
  dplyr::rename(group = path_score) %>%
  mutate(group = factor(group)) %>%
  gather(miR, deltaCT, -group) %>%
  mutate(miR = factor(miR, levels = ms_clinical_names$name2,
                      labels = ms_clinical_names$name1))



np_data_clean <- ms_data_aa_pneu


# np_data_clean ----------------------------------------------------------------
np_data_cor_res <- round(cor(np_data_clean, method = 'spearman'), 2)
rownames(np_data_cor_res) <-  factor(rownames(np_data_cor_res),
                                     levels = ms_clinical_names$name2,
                                     labels = ms_clinical_names$name1)
colnames(np_data_cor_res) <-  factor(colnames(np_data_cor_res),
                                     levels = ms_clinical_names$name2,
                                     labels = ms_clinical_names$name1)
write.xlsx(np_data_cor_res, 
           './writeup/correl01_pathology_score_spear_21121701.xlsx')

np_p.mat_res <- round(cor.mtest(np_data_clean), 4)

rownames(np_p.mat_res) <-  factor(rownames(np_p.mat_res),
                                     levels = ms_clinical_names$name2,
                                     labels = ms_clinical_names$name1)
colnames(np_p.mat_res) <-  factor(colnames(np_p.mat_res),
                                     levels = ms_clinical_names$name2,
                                     labels = ms_clinical_names$name1)

write.xlsx(np_p.mat_res, 
           './writeup/correl01_p-value_pathology_score_spear_21121701.xlsx')

## Make Figures ##
png("./figures/correl01_pathology_score_spear_21121701.png", w=2000, h=2000, pointsize=25)
corrplot(np_data_cor_res, method = 'circle', type = 'upper', 
         order = 'original', p.mat = np_p.mat_res, sig.level = 0.05,
         tl.col="black", tl.srt=45, tl.cex = 1.3, #Text label color and rotation
         cl.cex = 1, number.cex = 1, pch.cex = 2,
         mar = c(0, 0, 0, 0))
dev.off()

