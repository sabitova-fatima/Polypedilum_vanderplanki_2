


install.packages("chromoMap")
# install.packages("xlsx")
# install.packages("RecordLinkage")
# install.packages("Bioconductor")
# library(RecordLinkage)
# setwd("/cloud/project/Polypedilum_vanderplanki_2")

# library(RecordLinkage)
library(readxl)
library(tidyverse)
library(dplyr)
library(readr)
library(chromoMap)
library(caroline)
library(xlsx)

##### input files ##### 
new_ass <- read_excel("New assembly annotation.xlsx") 
yusurika <- read_excel("19c15301_Proteome_yusurika_data_2019 working.xlsx", sheet = 2) 

names(yusurika)[1] = "Transcript"
proteome <- left_join(yusurika, new_ass, by = "Transcript")

##### cleaning NA s #####

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


##### creating fc #####

proteome$t24_mean <- (proteome$`Abundances Scaled F5 Sample yusurika_T24` + proteome$`Abundances Scaled F6 Sample yusurika_T24` +
proteome$`Abundances Scaled F7 Sample yusurika_T24` + proteome$`Abundances Scaled F8 Sample yusurika_T24`) / 4
proteome$t0_mean <- (proteome$`Abundances Scaled F1 Sample yusurika_T0` + proteome$`Abundances Scaled F2 Sample yusurika_T0` +
                        proteome$`Abundances Scaled F3 Sample yusurika_T0` + proteome$`Abundances Scaled F4 Sample yusurika_T0`) / 4
proteome$fc <- proteome$t24_mean / proteome$t0_mean

# 
# ##### visualizing chromosomes ##### 
# 
# chr_1 <- filter(proteome, `Scaffold Id` == "chr_1")
# chr_2 <- filter(proteome, `Scaffold Id` == "chr_2")
# chr_3 <- filter(proteome, `Scaffold Id` == "chr_3")
# chr_4 <- filter(proteome, `Scaffold Id` == "chr_4")
# 
# chr <- data.frame (
#   c("chr_1", "chr_2", "chr_3", "chr_4"), 
#   c(1, 1, 1, 1), 
#   c(max(chr_1$End), max(chr_2$End), max(chr_3$End), max(chr_4$End)))
# 
# # все белки
# # anno <- proteome_meaningful %>% select(Transcript, `Scaffold Id`, Start, End, fc)
# 
# # только те белки, которые группами
# anno <- proteome %>% select(Transcript, `Scaffold Id`, Start, End, fc)
# write.delim(chr, "chr_file.txt", col.names = FALSE, sep = "\t")
# write.delim(anno, "anno_file.txt", col.names = FALSE, sep = "\t")
# 
# chromoMap("chr_file.txt",  "anno_file.txt")
# chromoMap("chr_file.txt",  "anno_file.txt", segment_annotation = T)
# 
# 
# chromoMap("chr_file.txt",  "anno_file.txt",
#           labels=T,
#           label_angle = -60,
#           chr_length = 8,
#           chr_width = 10,
#           canvas_width = 1000)
# 
# chromoMap("chr_file.txt",  "anno_file.txt",
#           data_type = "categorical",
#           plots = "bar")
# 
# chromoMap("chr_file.txt",  "anno_file.txt",
#           data_based_color_map = T,
#           data_type = "categorical",
#           plots = "scatter")
# 
# # намана
# chromoMap("chr_file.txt",  "anno_file.txt",
#           data_based_color_map = T,
#           segment_annotation = T,
#           chr_color="green",
#           data_colors = list("green"),
#           # anno_col = list("red"),
#           plot_color = "light blue",
#           data_type = "numeric",
#           plots = "bar",
#           ref_line = T,
#           refl_pos = 15,
#           heat_map = F)


##### how close proteins are to each other ##### 

proteome <-proteome[order(proteome$Start),]
proteome <-proteome[order(proteome$Start),]

##### function definitions ##### 

# shows how many proteins in a row come together
# if (o >= 4) - show only group of 4
# quantile(x, na.rm=TRUE)[2] - change 2 to 3 to use median

belki_sosedi <- function(chr1, x) {
  o = 1;
  chr$close <- 0
  chr <- as.data.frame(chr1)
  for (i in 1:nrow(chr)-1){
    n <- chr$Start[i + 1] - chr$End[i]
    # n - расстояние от одного белка до другого
    if (!is.na(n) && n <= quantile(x, na.rm=TRUE)[3]){
      if (o >= 8){
        print(o)
      }
      o = o + 1
    }
    else
      o = 1
  }
  return (chr)
}

##### chromosome 1 #####

proteome %>% filter(Strand == "+") %>% filter(`Scaffold Id` == "chr_1") -> chr_1_plus
proteome %>% filter(Strand == "-") %>% filter(`Scaffold Id` == "chr_1") -> chr_1_minus

x_plus = c()

for(i in 1:nrow(chr_1_plus)) {
  x_plus[i] = (chr_1_plus$Start[i + 1] - chr_1_plus$End[i])
}

chr <- belki_sosedi(chr_1_plus, x_plus)
# summary(x_plus, na.rm=TRUE)

x_minus = c()

for(i in 1:nrow(chr_1_minus)) {
  x_minus[i] = (chr_1_minus$Start[i + 1] - chr_1_minus$End[i])
}

belki_sosedi(chr_1_minus, x_minus)
# summary(x_minus, na.rm=TRUE)

##### chromosome 2 #####

proteome %>% filter(Strand == "+") %>% filter(`Scaffold Id` == "chr_2") -> chr_2_plus
proteome %>% filter(Strand == "-") %>% filter(`Scaffold Id` == "chr_2") -> chr_2_minus

x_plus = c()

for(i in 1:nrow(chr_2_plus)) {
  x_plus[i] = (chr_2_plus$Start[i + 1] - chr_2_plus$End[i])
}

chr <- belki_sosedi(chr_2_plus, x_plus)
# summary(x_plus, na.rm=TRUE)

x_minus = c()

for(i in 1:nrow(chr_2_minus)) {
  x_minus[i] = (chr_2_minus$Start[i + 1] - chr_2_minus$End[i])
}

belki_sosedi(chr_2_plus, x_minus)
# summary(x_minus, na.rm=TRUE)

##### chromosome 3 #####

proteome %>% filter(Strand == "+") %>% filter(`Scaffold Id` == "chr_3") -> chr_3_plus
proteome %>% filter(Strand == "-") %>% filter(`Scaffold Id` == "chr_3") -> chr_3_minus

x_plus = c()

for(i in 1:nrow(chr_3_plus)) {
  x_plus[i] = (chr_3_plus$Start[i + 1] - chr_3_plus$End[i])
}

belki_sosedi(chr_3_plus, x_plus)
# summary(x_plus, na.rm=TRUE)

x_minus = c()

for(i in 1:nrow(chr_3_minus)) {
  x_minus[i] = (chr_3_minus$Start[i + 1] - chr_3_minus$End[i])
}

belki_sosedi(chr_3_minus, x_minus)
# summary(x_minus, na.rm=TRUE)

##### chromosome 4 #####

proteome %>% filter(Strand == "+") %>% filter(`Scaffold Id` == "chr_4") -> chr_4_plus
proteome %>% filter(Strand == "-") %>% filter(`Scaffold Id` == "chr_4") -> chr_4_minus

x_plus = c()

for(i in 1:nrow(chr_4_plus)) {
  x_plus[i] = (chr_4_plus$Start[i + 1] - chr_4_plus$End[i])
}

# summary(x_plus, na.rm=TRUE)
belki_sosedi(chr_4_plus, x_plus)

x_minus <- c()
for(i in 1:nrow(chr_4_minus)) {
  x_minus[i] <- (chr_4_minus$Start[i + 1] - chr_4_minus$End[i])
}

belki_sosedi(chr_4_minus, x_minus)
# summary(x_minus, na.rm=TRUE)



#####  plots   #####


quantile(proteome$fc, na.rm = T)

ggplot(proteome, aes(x = Transcript, y = fc)) +
  geom_point()

proteome_copy <-data.frame(proteome)

median(proteome$fc, na.rm = T)
max(proteome$fc, na.rm = T)

proteome_copy$big_fc = 0

##### общий fc #####
res = 0 
proteome_copy$big_fc = 0
for (i in 1:length(proteome_copy$big_fc)){
  if (!(is.na(proteome_copy$fc[i])) && proteome_copy$fc[i] >= 1.5){
    res <- res + 1
    proteome_copy$big_fc[i] <- 1
  }
}

# 10 - 15
# 9 - 16
# 8 - 20
# 7 - 21
# 6 - 28
# 5 - 33
# 4 - 44
# 3 - 57
# 2 - 107
# 1.5 - 227
# 1 - 1705
# all - 4329

proteome_only_big <- filter(proteome_copy, big_fc == 1)

ggplot(proteome_only_big, aes(x = Transcript, y = fc)) +
  geom_point()

##### кластеры близких генов #####

proteome %>% filter(Strand == "+") %>% filter(`Scaffold Id` == "chr_1") -> chr_1_plus
proteome %>% filter(Strand == "-") %>% filter(`Scaffold Id` == "chr_1") -> chr_1_minus
proteome %>% filter(Strand == "+") %>% filter(`Scaffold Id` == "chr_2") -> chr_2_plus
proteome %>% filter(Strand == "-") %>% filter(`Scaffold Id` == "chr_2") -> chr_2_minus
proteome %>% filter(Strand == "+") %>% filter(`Scaffold Id` == "chr_3") -> chr_3_plus
proteome %>% filter(Strand == "-") %>% filter(`Scaffold Id` == "chr_3") -> chr_3_minus
proteome %>% filter(Strand == "+") %>% filter(`Scaffold Id` == "chr_4") -> chr_4_plus
proteome %>% filter(Strand == "-") %>% filter(`Scaffold Id` == "chr_4") -> chr_4_minus

count_dist <- function(str){
  distance <- c()
  for (i in 1:nrow(str)){
    distance[i] = str$Start[i + 1] - str$End[i]
  }
  return(distance)
}

# 25% белков ближе, чем на q_25 нуклеотидов
quant <- quantile(c(count_dist(chr_1_plus),
                   count_dist(chr_1_minus),
                   count_dist(chr_2_plus),
                   count_dist(chr_2_minus),
                   count_dist(chr_3_plus),
                   count_dist(chr_3_minus),
                   count_dist(chr_4_plus),
                   count_dist(chr_4_minus)), na.rm = T)

# менять, чтобы изменить условия отбора групп белков

# quant[3] - 50%, quant[2] - 25%
stat = quant[3]

# number of proteins in a row
treshold = 5

write_friends <- function(chr, i, count, value){
  n = i - count
  while (n < i)
  {
    chr$friends[n] <- value
    n = n + 1
  }
  return (chr)
}

find_friends <- function(chr, stat, treshold){
  count = 1;
  old_count = 1
  value = 1
  for (i in 1:(nrow(chr)-1)){
    if (chr$Start[i + 1] - chr$End[i] <= stat)
    {
      count = count + 1
    }
    else
    {
      count = 1
    }
    if (count >= treshold)
    {
      old_count = count
    }
    if (old_count >= treshold && count == 1)
    {
      group = paste("group", value, sep = "_")
      chr <- write_friends(chr, i, old_count, group)
      value = value + 1
      old_count = 1
    }
  }
  return (chr)
}

calc_mean_fc <- function(chr, start, stop)
{
  for (i in start:stop)
  {
    temp = filter(chr, friends == paste("group", i, sep = "_"))
    #  print(mean(temp$fc,  na.rm = T) >= mean(proteome$fc, na.rm = T))
    chr$mean_group_fc[chr$friends == paste("group", i, sep = "_")] <- mean(temp$fc,  na.rm = T)
  }
  return (chr)
}

# создаем и зануляем колонки

chr_1_plus$friends <- "none"
chr_1_minus$friends <- "none"
chr_2_plus$friends <- "none"
chr_2_minus$friends <- "none"
chr_3_plus$friends <- "none"
chr_3_minus$friends <- "none"
chr_4_plus$friends <- "none"
chr_4_minus$friends <- "none"

chr_1_plus$mean_group_fc <- 0
chr_1_minus$mean_group_fc <- 0
chr_2_plus$mean_group_fc <- 0
chr_2_minus$mean_group_fc <- 0
chr_3_plus$mean_group_fc <- 0
chr_3_minus$mean_group_fc <- 0
chr_4_plus$mean_group_fc <- 0
chr_4_minus$mean_group_fc <- 0

# пишет "group_1" итд, если n белков расположены ближе, чем 50% других
chr_1_plus <- find_friends(chr_1_plus, stat, treshold)
chr_1_minus <- find_friends(chr_1_minus, stat, treshold)
chr_2_plus <- find_friends(chr_2_plus, stat, treshold)
chr_2_minus <- find_friends(chr_2_minus, stat, treshold)
chr_3_plus <- find_friends(chr_3_plus, stat, treshold)
chr_3_minus <- find_friends(chr_3_minus, stat, treshold)
chr_4_plus <- find_friends(chr_4_plus, stat, treshold)
chr_4_minus <- find_friends(chr_4_minus, stat, treshold)

# добавляем mean изменение конц белка для каждой группы
chr_1_plus <- calc_mean_fc(chr_1_plus, 1, 21) # числа взяты по количеству групп ???
chr_1_minus <- calc_mean_fc(chr_1_minus, 1, 24)
chr_2_plus <- calc_mean_fc(chr_2_plus, 1, 15)
chr_2_minus <- calc_mean_fc(chr_2_minus, 1, 20)
chr_3_plus <- calc_mean_fc(chr_3_plus, 1, 27)
chr_3_minus <- calc_mean_fc(chr_3_minus, 1, 28)
chr_4_plus <- calc_mean_fc(chr_4_plus, 1, 4)
chr_4_minus <- calc_mean_fc(chr_4_minus, 1, 5)

7 + 7 + 3 + 4 + 11 + 7 + 1 + 2

proteome_bind <- bind_rows(
          chr_1_plus, 
          chr_1_minus,
          chr_2_plus,
          chr_2_minus,
          chr_3_plus,
          chr_3_minus,
          chr_4_plus,
          chr_4_minus)

# df в котором только белки, разделенные на группы (= стоящие близкими группами) 
proteome_groups_only <- filter(proteome_bind, mean_group_fc != 0)

# df в котором только белки, разделенные на группы (= стоящие близкими группами) 
# и изменение которых больше среднего изменения по протеому
proteome_meaningful<- filter(proteome_groups_only, mean_group_fc > 
                               mean(proteome$fc, na.rm = T))

# write.xlsx(proteome_meaningful, "proteome_meaningful.xlsx", sheetName="Sheet1", 
#            col.names = TRUE, row.names = TRUE, append = FALSE, showNA = TRUE)

##### сравнение с генами #####

ncomms <- read_excel("ncomms (1).xlsx") 

names(new_ass)[8] <- "Genes"
names(ncomms)[1] <- "Genes"
genome <- left_join(ncomms, new_ass, by = "Genes")

proteome_genome <- inner_join(proteome, genome, by = "Transcript")
proteome_genome$fc_gene <- proteome_genome$`RPKM (PvD0` / proteome_genome$`RPKM (PvD24)`

proteome_genome  %>% filter(fc < 10) -> proteome_genome_2

p <- ggplot(proteome_genome_2, aes(x=fc, y=fc_gene))

p + geom_jitter(size=0.5, color = "dark blue") +
labs( x = "Изменение концентрации белка", 
      y = "Изменение экспрессии гена",
      title ="Зависимость концентрации белка от экспрессии его гена",
      subtitle ="при индукции ангидробиоза") +
      theme_bw()


proteome_genome_2 <- select(proteome_genome_2, fc, fc_gene)

boxplot(proteome_genome_2,
        xlab = "Кратность изменения",
        main = "Изменение концентрации белков и экспрессии 
        их генов при индукции ангидробиоза",
        col = c("steelblue4", "violetred3"), horizontal = TRUE,
        las = 2, 
        log10="y",
        names = c("Белки",
                  "Гены"))

# удалено несколько крайних значений


mean(abs(proteome_genome$fc_gene - proteome_genome$fc), na.rm = T)
mean((proteome_genome$fc_gene - proteome_genome$fc), na.rm = T)
# в среднем изменение гена отличается от изменения его белка на 50%

temp <- filter(proteome_genome, abs(fc_gene - fc) > 5)

write.xlsx(temp, "5_times_changed_genes.xlsx", sheetName="Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE, showNA = TRUE)

write.xlsx(proteome, "proteome.xlsx", sheetName="Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE, showNA = TRUE)

proteome$Protein[1]
proteome$Protein[2]

##### поиск схожих белков #####

levenshteinSim_1 <- 0 * nrow(proteome) * nrow(proteome)

df <- data.frame(first_column, second_column)


proteome_filtered <- filter(proteome, fc > 2)

# 4490
# 20160100
max(levenshteinSim_1)
which(levenshteinSim_1 == max(levenshteinSim_1))

for (j in 2:nrow(proteome_filtered))
{
  for(i in 2:nrow(proteome_filtered))
  {
    levenshteinSim_1[j * i] <- levenshteinSim(proteome_filtered$Protein[j], proteome_filtered$Protein[i])
  }
}

levenshteinSim_1[levenshteinSim_1 >= 0.3 && levenshteinSim_1 < 1]
# 
# 0%       25%       50%       75%      100% 
# 0.0000000 0.1700137 0.2000000 0.2189781 1.0000000


which(str_detect(proteome$Protein, "X", negate = FALSE))

proteome_2 <- filter(proteome, !str_detect(proteome$Protein, "X", negate = FALSE))
# proteome_3 <- filter(proteome_2, is.string(proteome$Protein))
write.csv(proteome_2, "prot.csv")

for (i in 1:nrow(proteome_2))
{
  if (!is.string(proteome$Protein[i]))
  {
    print(i)
  }
}

##### plots from python #####

from_python <- read.csv("proteome_from_python.csv")
view(from_python)

from_python <- filter(from_python, molecular_weight < 750000)
# выкинули два крайниз значения 

p <- ggplot(from_python, aes(x=gravy, y=fc))
p + geom_jitter(size=0.5, color = "dark blue") +
  labs( x = "значение GRAVY", 
        y = "Изменение концентрации белка",
        title ="Зависимость концентрации белка от значения GRAVY",
        subtitle ="при индукции ангидробиоза") +
  theme_bw() + 
  scale_y_log10() 

p <- ggplot(from_python, aes(x=molecular_weight, y=fc))
p + geom_jitter(size=0.5, color = "dark blue") +
  labs( x = "Молекулярная масса белка", 
        y = "Изменение концентрации белка",
        title ="Зависимость изменения концентрации белка\nот его молекулярной массы",
        subtitle ="при индукции ангидробиоза") +
  theme_bw() +
  scale_y_log10() +
  scale_x_log10()

cor(data$log_income, data$wfood, method = "pearson")


# cor(data$log_income, data$wfood, method = "pearson")

from_python <- filter(from_python, is.na(from_python$fc) != T)
from_python <- filter(from_python, is.na(from_python$gravy) != T)
from_python <- filter(from_python, is.na(from_python$molecular_weight) != T)
proteome_genome_2 <- filter (proteome_genome_2, is.na(proteome_genome_2$fc_gene) != T)
proteome_genome_2 <- filter (proteome_genome_2, is.na(proteome_genome_2$fc) != T)

cor(from_python$fc, from_python$molecular_weight, method = "spearman")
cor(proteome_genome_2$fc, proteome_genome_2$fc_gene, method = "pearson")

from_python <- filter(from_python, from_python$fc > 10)
cor(from_python$fc, from_python$gravy, method = "pearson")

p8 <- ggplot(proteome_genome_2) + 
  geom_density(aes(x = fc), fill = "steelblue4", alpha = 0.6) +
  geom_density(aes(x = fc_gene),  fill = "violetred3", alpha = 0.6) +
  scale_x_log10() +
  theme_bw() +
  labs( x = "Изменение концентрации. Синее - белка, красное - гена", 
        y = "плотность",
        title ="Распределение изменеия концентраций белков и генов",
        subtitle ="при индукции ангидробиоза") +
  theme_bw()

p8


