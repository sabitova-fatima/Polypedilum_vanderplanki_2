# install.packages("chromoMap")
library(readxl)
library(tidyverse)
library(dplyr)
library(readr)
library(chromoMap)

##### input files ##### 
new_ass <- read_excel("New assembly annotation.xlsx") 
yusurika <- read_excel("19c15301_Proteome_yusurika_data_2019 working.xlsx", sheet = 2) 

# names(yusurika)[1] = "Transcript"
proteome <- left_join(yusurika, new_ass, by = "Transcript")

##### visualizing chromosomes ##### 

chr_1 <- filter(proteome, `Scaffold Id` == "chr_1")
chr_2 <- filter(proteome, `Scaffold Id` == "chr_2")
chr_3 <- filter(proteome, `Scaffold Id` == "chr_3")
chr_4 <- filter(proteome, `Scaffold Id` == "chr_4")

chr <- data.frame (
  c("chr_1", "chr_2", "chr_3", "chr_4"), 
  c(1, 1, 1, 1), 
  c(max(chr_1$End), max(chr_2$End), max(chr_3$End), max(chr_4$End)))

anno <- proteome %>% select(Transcript, `Scaffold Id`, Start, End)

write.delim(df, "chr_file.txt", col.names = FALSE, sep = "\t")
write.delim(anno, "anno_file.txt", col.names = FALSE, sep = "\t")

chromoMap("chr_file.txt",  "anno_file.txt")

##### how close proteins are to each other ##### 

proteome <-proteome[order(proteome$Start),]
proteome <-proteome[order(proteome$Start),]

proteome %>% filter(Strand == "+") %>% filter(`Scaffold Id` == "chr_1") -> chr_1_plus
proteome %>% filter(Strand == "-") %>% filter(`Scaffold Id` == "chr_1") -> chr_1_minus

x_plus = c()

for(i in 1:nrow(chr_1_plus)) {
  x_plus[i] = (chr_1_plus$Start[i + 1] - chr_1_plus$End[i])
}

nrow(chr_1_plus)

mean(x_plus, na.rm=TRUE)
# 13241.96
min(x_plus, na.rm=TRUE)
# -51577
max(x_plus, na.rm=TRUE)
# 474558

##
x_minus = c()

for(i in 1:nrow(chr_1_minus)) {
  x_minus[i] = (chr_1_minus$Start[i + 1] - chr_1_minus$End[i])
}


nrow(chr_1_minus)
# 
mean(x_minus, na.rm=TRUE)
# 13539.53
min(x_minus, na.rm=TRUE)
# -58113
max(x_minus, na.rm=TRUE)

