
library(caroline)

names(yusurika)[1] = "Transcript"

proteome <- left_join(yusurika, new_ass, by = "Transcript")
names(proteome)

chr_1 <- filter(proteome, `Scaffold Id` == "chr_1")
chr_2 <- filter(proteome, `Scaffold Id` == "chr_2")
chr_3 <- filter(proteome, `Scaffold Id` == "chr_3")
chr_4 <- filter(proteome, `Scaffold Id` == "chr_4")

max(chr_4$End)

#chromosome file
# chr name, start (1), end

# annotation file
# 1 - Transcript
# 2 - chr_4
# 3 - start
# 4 - end

# new_chr <- chr_1 %>% select(Transcript, `Scaffold Id`, Start, End)


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
