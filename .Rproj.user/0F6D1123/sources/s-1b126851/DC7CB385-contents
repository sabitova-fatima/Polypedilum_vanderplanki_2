# install.packages("chromoMap")
library(readxl)
library(tidyverse)
library(dplyr)
library(readr)
library(chromoMap)

##### input files ##### 
new_ass <- read_excel("New assembly annotation.xlsx") 
yusurika <- read_excel("19c15301_Proteome_yusurika_data_2019 working.xlsx", sheet = 2) 

names(yusurika)[1] = "Transcript"
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

write.delim(chr, "chr_file.txt", col.names = FALSE, sep = "\t")
write.delim(anno, "anno_file.txt", col.names = FALSE, sep = "\t")

chromoMap("chr_file.txt",  "anno_file.txt")

##### how close proteins are to each other ##### 

proteome <-proteome[order(proteome$Start),]
proteome <-proteome[order(proteome$Start),]

##### chromosome 1 #####

proteome %>% filter(Strand == "+") %>% filter(`Scaffold Id` == "chr_1") -> chr_1_plus
proteome %>% filter(Strand == "-") %>% filter(`Scaffold Id` == "chr_1") -> chr_1_minus

x_plus = c()

for(i in 1:nrow(chr_1_plus)) {
  x_plus[i] = (chr_1_plus$Start[i + 1] - chr_1_plus$End[i])
}

summary(x_plus, na.rm=TRUE)

x_minus = c()

for(i in 1:nrow(chr_1_minus)) {
  x_minus[i] = (chr_1_minus$Start[i + 1] - chr_1_minus$End[i])
}

summary(x_minus, na.rm=TRUE)

##### chromosome 2 #####

proteome %>% filter(Strand == "+") %>% filter(`Scaffold Id` == "chr_2") -> chr_2_plus
proteome %>% filter(Strand == "-") %>% filter(`Scaffold Id` == "chr_2") -> chr_2_minus

x_plus = c()

for(i in 1:nrow(chr_2_plus)) {
  x_plus[i] = (chr_2_plus$Start[i + 1] - chr_2_plus$End[i])
}

summary(x_plus, na.rm=TRUE)

x_minus = c()

for(i in 1:nrow(chr_2_minus)) {
  x_minus[i] = (chr_2_minus$Start[i + 1] - chr_2_minus$End[i])
}

summary(x_minus, na.rm=TRUE)

##### chromosome 3 #####

proteome %>% filter(Strand == "+") %>% filter(`Scaffold Id` == "chr_3") -> chr_3_plus
proteome %>% filter(Strand == "-") %>% filter(`Scaffold Id` == "chr_3") -> chr_3_minus

x_plus = c()

for(i in 1:nrow(chr_3_plus)) {
  x_plus[i] = (chr_3_plus$Start[i + 1] - chr_3_plus$End[i])
}


for(i in 1:nrow(chr_4_plus)) {
  x_plus[i] = (chr_4_plus$Start[i + 1] - chr_4_plus$End[i])
}

for (i in 1:nrow(chr_4_plus)){
  n <- chr_4_plus$Start[i + 1] - chr_4_plus$End[i]
  if (n < quantile(x_minus, na.rm=TRUE)[2]){
    print(o)
    o = o + 1
  }
  else{
    o = 1
  }
}

summary(x_plus, na.rm=TRUE)

x_minus = c()

for(i in 1:nrow(chr_3_minus)) {
  x_minus[i] = (chr_3_minus$Start[i + 1] - chr_3_minus$End[i])
}

summary(x_minus, na.rm=TRUE)

##### chromosome 4 #####

proteome %>% filter(Strand == "+") %>% filter(`Scaffold Id` == "chr_4") -> chr_4_plus
proteome %>% filter(Strand == "-") %>% filter(`Scaffold Id` == "chr_4") -> chr_4_minus

x_plus = c()

for(i in 1:nrow(chr_4_plus)) {
  x_plus[i] = (chr_4_plus$Start[i + 1] - chr_4_plus$End[i])
}

for (i in 1:nrow(chr_4_plus)){
  n <- chr_4_plus$Start[i + 1] - chr_4_plus$End[i]
  if (n < quantile(x_minus, na.rm=TRUE)[2]){
    print(o)
    o = o + 1
  }
  else{
    o = 1
  }
}

summary(x_plus, na.rm=TRUE)

x_minus <- c()

for(i in 1:nrow(chr_4_minus)) {
  x_minus[i] <- (chr_4_minus$Start[i + 1] - chr_4_minus$End[i])
}

summary(x_minus, na.rm=TRUE)

# 25% of proteins are closer
quantile(x_minus, na.rm=TRUE)[2]
sum(x_minus < quantile(x_minus, na.rm=TRUE)[2], na.rm=TRUE)

o = 1;

for (i in 1:nrow(chr_4_minus)){
  n <- chr_4_minus$Start[i + 1] - chr_4_minus$End[i]
  if (n < quantile(x_minus, na.rm=TRUE)[2]){
    print(o)
    o = o + 1
  }
  else{
    o = 1
  }
}


