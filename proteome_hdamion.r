library(readxl)
library(tidyverse)
library(dplyr)
# library(readr)
# library(chromoMap)
# library(caroline)
# library(xlsx)
library(limma)
library(fgsea)

##### input files ##### 
new_ass <- read_excel("New assembly annotation.xlsx") 
yusurika <- read_excel("19c15301_Proteome_yusurika_data_2019 working.xlsx", sheet = 2) 

names(yusurika)[1] = "Transcript"
proteome <- left_join(yusurika, new_ass, by = "Transcript")

##### cleaning NA s and duplicates #####

proteome %>% distinct(Transcript, .keep_all = TRUE) -> proteome

# если пропущено одно значение из 4 - заменяем средним значением
# использовать медиану?
for (i in 1:length(proteome$`Abundances Scaled F5 Sample yusurika_T24`)){
  if (is.na(proteome$`Abundances Scaled F5 Sample yusurika_T24`[i])){
    proteome$`Abundances Scaled F5 Sample yusurika_T24`[i] = 
      (proteome$`Abundances Scaled F6 Sample yusurika_T24`[i] +
         proteome$`Abundances Scaled F7 Sample yusurika_T24`[i] +
         proteome$`Abundances Scaled F8 Sample yusurika_T24`[i]) / 3
  }
}

for (i in 1:length(proteome$`Abundances Scaled F6 Sample yusurika_T24`)){
  if (is.na(proteome$`Abundances Scaled F6 Sample yusurika_T24`[i])){
    proteome$`Abundances Scaled F6 Sample yusurika_T24`[i] = 
      (proteome$`Abundances Scaled F5 Sample yusurika_T24`[i] +
         proteome$`Abundances Scaled F7 Sample yusurika_T24`[i] +
         proteome$`Abundances Scaled F8 Sample yusurika_T24`[i]) / 3
  }
}

for (i in 1:length(proteome$`Abundances Scaled F7 Sample yusurika_T24`)){
  if (is.na(proteome$`Abundances Scaled F7 Sample yusurika_T24`[i])){
    proteome$`Abundances Scaled F7 Sample yusurika_T24`[i] = 
      (proteome$`Abundances Scaled F6 Sample yusurika_T24`[i] +
         proteome$`Abundances Scaled F5 Sample yusurika_T24`[i] +
         proteome$`Abundances Scaled F8 Sample yusurika_T24`[i]) / 3
  }
}

for (i in 1:length(proteome$`Abundances Scaled F8 Sample yusurika_T24`)){
  if (is.na(proteome$`Abundances Scaled F8 Sample yusurika_T24`[i])){
    proteome$`Abundances Scaled F8 Sample yusurika_T24`[i] = 
      (proteome$`Abundances Scaled F6 Sample yusurika_T24`[i] +
         proteome$`Abundances Scaled F7 Sample yusurika_T24`[i] +
         proteome$`Abundances Scaled F5 Sample yusurika_T24`[i]) / 3
  }
}

for (i in 1:length(proteome$`Abundances Scaled F1 Sample yusurika_T0`)){
  if (is.na(proteome$`Abundances Scaled F1 Sample yusurika_T0`[i])){
    proteome$`Abundances Scaled F1 Sample yusurika_T0`[i] = 
      (proteome$`Abundances Scaled F2 Sample yusurika_T0`[i] +
         proteome$`Abundances Scaled F3 Sample yusurika_T0`[i] +
         proteome$`Abundances Scaled F4 Sample yusurika_T0`[i]) / 3
  }
}

for (i in 1:length(proteome$`Abundances Scaled F2 Sample yusurika_T0`)){
  if (is.na(proteome$`Abundances Scaled F2 Sample yusurika_T0`[i])){
    proteome$`Abundances Scaled F2 Sample yusurika_T0`[i] = 
      (proteome$`Abundances Scaled F1 Sample yusurika_T0`[i] +
         proteome$`Abundances Scaled F3 Sample yusurika_T0`[i] +
         proteome$`Abundances Scaled F4 Sample yusurika_T0`[i]) / 3
  }
}

for (i in 1:length(proteome$`Abundances Scaled F3 Sample yusurika_T0`)){
  if (is.na(proteome$`Abundances Scaled F3 Sample yusurika_T0`[i])){
    proteome$`Abundances Scaled F3 Sample yusurika_T0`[i] = 
      (proteome$`Abundances Scaled F2 Sample yusurika_T0`[i] +
         proteome$`Abundances Scaled F1 Sample yusurika_T0`[i] +
         proteome$`Abundances Scaled F4 Sample yusurika_T0`[i]) / 3
  }
}

for (i in 1:length(proteome$`Abundances Scaled F4 Sample yusurika_T0`)){
  if (is.na(proteome$`Abundances Scaled F4 Sample yusurika_T0`[i])){
    proteome$`Abundances Scaled F4 Sample yusurika_T0`[i] = 
      (proteome$`Abundances Scaled F2 Sample yusurika_T0`[i] +
         proteome$`Abundances Scaled F3 Sample yusurika_T0`[i] +
         proteome$`Abundances Scaled F1 Sample yusurika_T0`[i]) / 3
  }
}

# больше 1 пропущено - удаляем
proteome %>% drop_na(`Abundances Scaled F1 Sample yusurika_T0`,
                     `Abundances Scaled F2 Sample yusurika_T0`,
                     `Abundances Scaled F3 Sample yusurika_T0`,
                     `Abundances Scaled F4 Sample yusurika_T0`,
                     `Abundances Scaled F5 Sample yusurika_T24`,
                     `Abundances Scaled F6 Sample yusurika_T24`,
                     `Abundances Scaled F7 Sample yusurika_T24`,
                     `Abundances Scaled F8 Sample yusurika_T24`) -> proteome

##### creating fc #####

proteome$t24_mean <- (proteome$`Abundances Scaled F5 Sample yusurika_T24` + proteome$`Abundances Scaled F6 Sample yusurika_T24` +
                        proteome$`Abundances Scaled F7 Sample yusurika_T24` + proteome$`Abundances Scaled F8 Sample yusurika_T24`) / 4
proteome$t0_mean <- (proteome$`Abundances Scaled F1 Sample yusurika_T0` + proteome$`Abundances Scaled F2 Sample yusurika_T0` +
                       proteome$`Abundances Scaled F3 Sample yusurika_T0` + proteome$`Abundances Scaled F4 Sample yusurika_T0`) / 4
proteome$fc <- proteome$t24_mean / proteome$t0_mean

##### limma #####

load("/cloud/project/kegg_definitions_list.RData")

design <- model.matrix(~ 0+factor(c(1,1,1,1,2,2,2,2)))
colnames(design) <- c("T0", "T24")
fit <- lmFit(proteome[3:10], design)
contrast.matrix <- makeContrasts(T24-T0, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2)
# plotSA(fit2)
fit.data <- as.data.frame(fit2)
# topTable(fit2, coef=1, adjust="BH")
glimpse(fit.data)

##### fgsea #####

test <- c(fit2$coefficients)
names(test) <- proteome$Transcript

fgsea(
  pathways = kegg_definitions_grouped_list,
  stats    = test,
  # minSize  = 2,
  # maxSize  = 500,
  eps = 0.0,
  scoreType = "pos"
)

# Warning message:
# In preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :
# There are ties in the preranked stats (13.28% of the list).
# The order of those tied genes will be arbitrary, which may produce unexpected results.

# trying random
rand <- rnorm(nrow(proteome))
names(rand) <- proteome$Transcript

fgsea(
  pathways = kegg_definitions_grouped_list,
  stats    = rand,
  # minSize  = 2,
  # maxSize  = 500,
  eps = 0.0
  # scoreType = "pos"
)

# Empty data.table (0 rows and 8 cols): pathway,pval,padj,log2err,ES,NES...