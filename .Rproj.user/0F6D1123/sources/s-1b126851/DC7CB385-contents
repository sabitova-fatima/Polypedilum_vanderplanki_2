# install.packages("chromoMap")
library(readxl)
library(tidyverse)
library(dplyr)
library(readr)
library(chromoMap)

# input files
new_ass <- read_excel("New assembly annotation.xlsx") 
yusurika <- read_excel("19c15301_Proteome_yusurika_data_2019 working.xlsx", sheet = 2) 

# names(yusurika)[1] = "Transcript"
proteome <- left_join(yusurika, new_ass, by = "Transcript")

chr_1 <- filter(proteome, `Scaffold Id` == "chr_1")
chr_2 <- filter(proteome, `Scaffold Id` == "chr_2")
chr_3 <- filter(proteome, `Scaffold Id` == "chr_3")
chr_4 <- filter(proteome, `Scaffold Id` == "chr_4")

first_column <- c("chr_1", "chr_2", "chr_3", "chr_4")
second_column <- c(1, 1, 1, 1)
third_column <- c(max(chr_1$End), max(chr_2$End), max(chr_3$End), max(chr_4$End))
chr <- data.frame(first_column, second_column, third_column)

anno <- proteome %>% select(Transcript, `Scaffold Id`, Start, End)

write.delim(df, "chr_file.txt", col.names = FALSE, sep = "\t")
write.delim(anno, "anno_file.txt", col.names = FALSE, sep = "\t")

head(read.table("chr_file.txt", sep = "\t"))
head(read.table("anno_file.txt", sep = "\t"))

chromoMap("chr_file.txt",  "anno_file.txt")