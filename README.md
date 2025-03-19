# R-assignment

## Data Inspection

### assigning txt files to fang or snp
fang <- fang_et_al_genotypes
snp <- snp_position

dim(fang) 
dim(snp) 
#### dim 2783, 986. this lists colms and rows
#### dim 983, 15. this lists colms and rows

ncol(fang) 
ncol(snp) 
#### number of columns is 986
#### number of columns is 15

nrow(fang) 
nrow(snp) 
#### 2783 rows
#### 983 rows

names(snp) 
#### tells the names of column labels

head(fang, 10) 
head(snp, 10)
#### shows first 10 rows of data for fang and snp files

tail(fang, 10) 
tail(snp, 10)
#### shows last 10 rows of data for fang and snp files

### this forloop displays the contents of each variable in the snp file
for(i in 1:ncol(snp)) {
  print(names(snp)[i])
  print(table(snp[, i], useNA = "ifany"))
  cat("\n") # Adds a line before the next iteration
} 

### same thing but for fang file
for(i in 1:ncol(fang)) {
  print(names(fang)[i])
  print(table(fang[, i], useNA = "ifany"))
  cat("\n") # Adds a line before the next iteration
} 

summary(snp) 
#### this shows the length, class, mode, min's/maxes, quantiles of the file

## Data Processing

### read the file 
genotype_data <- read.table(file = "fang_et_al_genotypes.txt", header = TRUE)

### For Maize: filter these specified groups, ZMMIL, ZMMLR, ZMMR and assign to genotypes_maize. Essentially, this is extracting the rows when the group is one of the three indicated groups
genotypes_maize <- filter(genotype_data, Group == "ZMMIL" | Group == "ZMMLR" | Group == "ZMMR")

print(genotypes_maize)

### Same thing but for teosinte data, extracting the specified groups ZMPBA, ZMPIL, and ZMPJA
genotypes_teosinte <- filter(genotype_data, Group == "ZMPBA" | Group == "ZMPIL" | Group == "ZMPJA")

print(genotypes_teosinte)

### Transpose the data for maize and teosinte
maize_transpose <- t(genotypes_maize)
teosinte_transpose <- t(genotypes_teosinte)

print(maize_transpose)
print(teosinte_transpose)

### remove the unwanted columns from the snp data as indicated by the snp[-c(2,5:15)]. This means to remove column 2 AND 5 through 15.
SNP_clean <- snp[-c(2,5:15)]
SNP_clean 
#### we can see we have the columns we want: SNP_ID, Chromosome, and Position

### label rownames as SNP_ID so we can then merge SNP_clean and the teosinte and maize files by row.names
rownames(SNP_clean) <- SNP_clean$SNP_ID
rownames(SNP_clean)

### merge maize and snp_clean
joined_maize <- merge(SNP_clean, maize_transpose, by = "row.names")

### merge teosinte and snp_clean
joined_teosinte <- merge(SNP_clean, teosinte_transpose, by = "row.names")

### create a for-loop that cycles through chromosomes 1-10 in the joined_maize data file in my global env
for (i in 1:10) {chr_data <- subset(joined_maize, joined_maize[[3]] == i)  # Filters rows where the third column matches the chromosome number
  chr_data[is.na(chr_data)] <- "?" # Replaces NA values with "?"
  chr_data <- chr_data[order(chr_data[[2]]), ]  # Sorts the data by the second column in increasing order
  assign(paste0("chr", i, "_increasing_maize"), chr_data)
}  # Saves the result as a new dataframe in the environment

### same thing but for teosinte data
for (i in 1:10) {chr_teosinte <- subset(joined_teosinte, joined_teosinte[[3]] == i)
  chr_teosinte[is.na(chr_teosinte)] <- "?"  
  chr_teosinte <- chr_teosinte[order(chr_teosinte[[2]]), ] 
  assign(paste0("chr", i, "_increasing_teosinte"), chr_teosinte)
}

### Maize in decreasing order replacing "?" with "-"

for (i in 1:10) {
  chr_maize <- subset(joined_maize, joined_maize[[3]] == i) # Creates a new variable for each chromosome dataset
  chr_maize[chr_maize == "?"] <- "-" # Replaces "?" with "-"
  chr_maize <- chr_maize[order(chr_maize[[2]], decreasing = TRUE), ] # Sorts by the second column in decreasing order
  assign(paste0("chr", i, "_decreasing_maize"), chr_maize) # Saves as a separate dataframe in the environment
  write.table(chr_maize, file = paste0("chr", i, "_decreasing_maize.txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
} # Writes to a text file to save new data created

### same thing but for teosinte data
for (i in 1:10) {
  dec_teosinte <- subset(joined_teosinte, joined_teosinte[[3]] == i)
  dec_teosinte[dec_teosinte == "?"] <- "-"
  dec_teosinte <- dec_teosinte[order(dec_teosinte[[2]], decreasing = TRUE), ]
  assign(paste0("chr", i, "_decreasing_teosinte"), dec_teosinte)
  write.table(dec_teosinte, file = paste0("chr", i, "_decreasing_teosinte.txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

### unknown positions for maize, this is looking at rows where the third column has a "?", had an issue with this in UNIX assingment not creating a file with content and it did the same thing for R so I must be doing something incorrect
unknown_positions_maize <- subset(joined_maize, joined_maize[[3]] == "?")

write.table(unknown_positions_maize, file = "unknown_positions_maize.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

### unknown positions for teosinte, this is also looking at rows where the third column has a "?", same issue here as above
unknown_positions_teosinte <- subset(joined_teosinte, joined_teosinte[[3]] == "?")

write.table(unknown_positions_teosinte, file = "unknown_positions_teosinte.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## Data Visualization: Part II

library(ggplot2)
library(dplyr)

### Combine maize and teosinte data with a new "type" column, which is either maize or teosinte
maize_type_data <- joined_maize %>% mutate(Type = "Maize")
teosinte_type_data <- joined_teosinte %>% mutate(Type = "Teosinte")

### Merge teosinte and maize datasets, naming this new data combined_data
combined_data <- bind_rows(maize_type_data, teosinte_type_data)

### This creates the plot that shows the distribution of SNPs across chromosomes in both species, the first line is indicating to use combined_data as the data for this graph, and the aes() command is indicating to use the Chromosome data as a factor and fill is the "type" which will be maize or teosinte. The geom_bar() line is telling the position of the bars to be side by side rather than a different arrangement. The lab() is simply labeling the title, x and y axes, and the legend for the "fill" which is species. The last two lines are indicating the theme of the graph produced and the colors of the bars for each species.
ggplot(combined_data, aes(x = as.factor(combined_data$Chromosome), fill = Type)) +
  geom_bar(position = "dodge") +
  labs(title = "SNPs across Chromosomes",
       x = "Chromosome",
       y = "SNP #",
       fill = "Species") +
  theme_minimal() +
  scale_fill_manual(values = c("Maize" = "blue", "Teosinte" = "red"))
# From what I can tell based on the graph it seems as though teosinte and maize have the same number of SNPs in each chromosome





         
