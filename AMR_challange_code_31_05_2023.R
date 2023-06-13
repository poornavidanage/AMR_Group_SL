setwd("C:\\Users\\poorn\\OneDrive - University of Moratuwa\\Research_Projects\\AMR Challenge 2023\\Data Challenge - SL Team\\ATLAS_Antibiotics\\AMR_submission")
list.files()
data <- read.csv("2022_04_28_Vivli_Atlas_Data_Antibiotics.csv")

#.................................................................
#Count the most abundant 10 Source (Blood/Urine etc)
# Count the frequency of each unique species
df2 <- table(data$Source)

df2 <- data %>%
  group_by(Source) %>%
  summarise(count = n())
#make the counts in descending order
df2 <- arrange(df2, desc(df2$count))

# Calculate the cumulative percentages
df2$CumulativePercentage <- cumsum(df2$count) / sum(df2$count) * 100

# Save the first 10 entries of column 1 in a data frame
top_10_Source <- data.frame(First10Entries = df2$Source[1:10])
write.csv(top_10_Source,"top_10_Source.csv")

# ...........................................................
#Count the most abundant 10 Species 
# Count the frequency of each unique species
df2 <- table(data$Species)

df2 <- data %>%
  group_by(Species) %>%
  summarise(count = n())
#make the counts in decending order
df2 <- arrange(df2, desc(df2$count))

# Calculate the cumulative percentages
df2$CumulativePercentage <- cumsum(df2$count) / sum(df2$count) * 100

# Save the first 10 entries of column 1 in a data frame
top_10_Species <- data.frame(First10Entries = df2$Species[1:10])
write.csv(top_10_Species,"top_10_Species.csv")

# ...........................................................

# Screen the entries in "data" based on matching values in Species and Source columns with df2 and df3
df4 <- data[data$Species %in% top_10_Species$First10Entries & data$Source %in% top_10_Source$First10Entries, ]


# ...........................................................
#Getting all the columns with AMR resistance data
df1 <- df4[ , grepl( "_I" , names( df4 ) ) ]
write_xlsx(df1,"AMRlist.xlsx")
# Replace "Intermediate" with "Resistant" in all columns

df2 <- df1 %>% mutate_all(~ifelse(. == "Intermediate", "Resistant", .))

df2 <- cbind(df4$Year, df2)
df2 <- cbind(df4$Country, df2)

# Rename the column header

colnames(df2)[1] <- "Country"
colnames(df2)[2] <- "Year"
library(dplyr)



df2 <- df2[, -2]  # Remove the second column : years

# Calculate the resistant percentage of each 
result_percentage <- df2 %>%  
  group_by(Country) %>%
  mutate_at(vars(starts_with("Amikacin_I"):starts_with("Meropenem.vaborbactam_I")),
            ~ sum(. == "Resistant") / sum(. %in% c("Resistant", "Susceptible"))) %>%
  rename_with(~paste0("Resistant_%", .), starts_with("Amikacin_I"):starts_with("Meropenem.vaborbactam_I")) %>%
  distinct()
# get the summary for all the Country to get the AMR profile 
df5 <- result_percentage #load result to df5 
write_xlsx(df5,"result_percentage.xlsx")
# replace "Resistant_%" in the column headers
colnames(df5) <- gsub("Resistant_%", "", colnames(df5))

# Remove columns with all NAs
df_filtered <- df5[, !apply(is.na(df5), 2, all)]

write_xlsx(df_filtered,"country_profile_AMR.xlsx")

#screening the antibiotics based on the selected_list
# ths file is manually curated based on the discussion with Medical Microbiologist
list.files()
screened_antibiotics <- read.csv("AMRlist_selected.csv", header = TRUE)

# Step 1: Extract column names from df2

column_names <- screened_antibiotics$antibiotic

# Step 2: Find common column names between df1 and df2
common_columns <- intersect(column_names, colnames(df_filtered))

# Step 3: Subset columns in df1 based on common column names
df1_subset <- df_filtered[, common_columns, drop = FALSE]
# add the country names
country_profile <- cbind(df_filtered$Country, df1_subset)
colnames(new_country_profile)[1] <- "Country"
# save the data set with antibiotics filterd 

write_xlsx(country_profile,"country_profile_AMR_23.xlsx")
colnames(country_profile)[1] <- "Country"


# Function to count values between 0 and 1 in a row
count_between_0_1 <- function(row) {
  sum(row > 0 & row < 1, na.rm = TRUE)
}

# Count cells with values between 0 and 1 for each row
num_cells_per_row <- apply(country_profile, 1, count_between_0_1)

# Subset rows with more than 17 cells between 0 and 1
new_country_profile <- country_profile[num_cells_per_row > 17, ] # new_country_profile generated here 

#rearrange the column names
colnames(new_country_profile) <- gsub("_I", "", colnames(new_country_profile))
colnames(new_country_profile)[1] <- "Country"

write_xlsx(new_country_profile,"country_61_profile_AMR_23_1.xlsx")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Remove new_country_profile : Moxifloxacin , and test_1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_1 <- new_country_profile[, -which(names(new_country_profile) == "Moxifloxacin")]
write_xlsx(test_1,"country_61_profile_AMR_23_1.xlsx")

new_country_profile <- test_1
# Set the first column as row headers
rownames(new_country_profile) <- new_country_profile$Country

# Remove the first column
new_country_profile <- new_country_profile[, -1]
#Heat Map with dendrogram

#______________________________________________________________
#get the 56 countries refined manually  to plot the heat map


heat_map_data <- read.csv("heatmap_56_countries_AMR.csv", header = TRUE)
# Set the first column as row headers
colnames(heat_map_data)[1] <- "Country"
rownames(heat_map_data) <- heat_map_data$Country

# Remove the first column
heat_map_data <- heat_map_data[, -1]
#Heat Map with dendrogram
# Get the maximum value of the data frame
max_value <- max(unlist(apply(heat_map_data, 2, max, na.rm = TRUE)))

# Print the maximum value
print(max_value)
list.files()
library(circlize)
mycols <- colorRamp2(breaks = c(0, 0.5, 0.906), 
                     colors = c("blue", "green", "red"))
#Transforming the data frame into a matrix
hm_matrx <- data.matrix(heat_map_data)

library(dendextend)

row_dend = hclust(dist(hm_matrx)) # row clustering
col_dend = hclust(dist(t(hm_matrx))) # column clustering
Heatmap(hm_matrx, name = "AMR Percentage (%)", 
        column_title = "Antibiotic", row_title = "Countries",
        row_names_gp = gpar(fontsize = 10),col = mycols,
        #cluster_rows = color_branches(row_dend, k = 3),
        #cluster_columns = color_branches(col_dend, k = 4)
        )

#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Network Analysis
#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

library(dplyr)
data <- read.csv("2022_04_28_Vivli_Atlas_Data_Antibiotics.csv")
#_________________________________________________________________________________
#1
#________________________________________________________________________________________________
#Species == "Staphylococcus aureus", Phenotype == "MRSA"

# Filter the data for the years 2010 to 2020 and the given species
filtered_df <- data %>%
  filter(Year >= 2004 & Year <= 2020, Species == "Staphylococcus aureus", Phenotype == "MRSA")

# Get the count of the given species for each year and country separately
species_count <- filtered_df %>%
  group_by(Year, Country) %>%
  summarise(count = n())
# Pivot the pass rate data frame back to wide format
Str_AU_MRSA <- species_count %>%
  pivot_wider(names_from = Country, values_from = count)
# normalizing with the Population 
popl_78 <- read.csv("population_data_78_countries.csv" )

df_pop <- Str_AU_MRSA
for (country in names(Str_AU_MRSA)[-1]) {
  df_pop[[country]] <- Str_AU_MRSA[[country]]*1E6 / popl_78$Population[popl_78$ï..Country == country]
}

#remove the columns
# Calculate the total number of cells in each column
total_cells <- colSums(!is.na(df_pop))

# Calculate the total number of NA values in each column
na_counts <- colSums(is.na(df_pop))

# Calculate the percentage of NA values in each column
#17 = observation of thr df_pop
na_percentage <- na_counts*100 / 17

# Set the threshold as 40%
threshold <- 40

# Identify columns exceeding the threshold
columns_to_remove <- names(na_percentage[na_percentage > threshold])

# Remove the identified columns from the data frame
df_to_NW <- df_pop[, !(names(df_pop) %in% columns_to_remove)]

# get the average of df_to_NW to load to cytoescape 

# Calculate the average of each column excluding NA values, omitting column 1
averages_Str_AU_MRSA <- colMeans(df_to_NW[, -1], na.rm = TRUE)
averages_Str_AU_MRSA <- data.frame(averages_Str_AU_MRSA)

write.csv(averages_Str_AU_MRSA, file = 'averages_Str_AU_MRSA.csv')

#input for the network analysis

data2 <- df_to_NW

my_data_2 <- data2[-1]
row.names(my_data_2) <- data2[,1]

my_data_2

dim(my_data_2)

df <- data.frame(matrix(ncol = 17, nrow = 40)) #change the number of columns and rows
samples = data2[,1]


res <- cor(my_data_2)
round(res, 2)

cor(my_data_2, use = "complete.obs")

library("Hmisc")
res2 <- rcorr(as.matrix(my_data_2))
# Extract the correlation coefficients
res2$r
# Extract p-values
res2$P

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

library(Hmisc)
res2<-rcorr(as.matrix(my_data_2[,1:26])) # set the column numbers 1:N of my_data
flattenCorrMatrix(res2$r, res2$P)
sq_cor_p <-  flattenCorrMatrix(res2$r, res2$P)

head(sq_cor_p)

data_P_fil <- sq_cor_p[sq_cor_p$p<=0.05,]

max(data_P_fil$cor)
min(data_P_fil$cor)


data_R_min <- subset(data_P_fil, data_P_fil$cor >= -1 & data_P_fil$cor <= -0.5)
data_R_max <- subset(data_P_fil, data_P_fil$cor >= 0.5 & data_P_fil$cor <= 1)

df_to_ctoescape <- rbind(data_R_min, data_R_max)

dim(df_to_ctoescape) ### input for the cytoescape 

write.csv(df_to_ctoescape, file = 'Str_AU_MRSA_input_to_cytoescape.csv')
#_________________________________________________________________________________
#2

#________________________________________________________________________________________________
#Species == "Klebsiella pneumoniae", Phenotype == "ESBL"

# Filter the data for the years 2010 to 2020 and the given species
filtered_df <- data %>%
  filter(Year >= 2004 & Year <= 2020, Species == "Klebsiella pneumoniae", Phenotype == "ESBL")

# Get the count of the given species for each year and country separately
species_count <- filtered_df %>%
  group_by(Year, Country) %>%
  summarise(count = n())
# Pivot the pass rate data frame back to wide format
K_pnu_ESBL <- species_count %>%
  pivot_wider(names_from = Country, values_from = count)
# normalizing with the Population 
popl_78 <- read.csv("population_data_78_countries.csv" )

df_pop <- K_pnu_ESBL
for (country in names(K_pnu_ESBL)[-1]) {
  df_pop[[country]] <- K_pnu_ESBL[[country]]*1E6 / popl_78$Population[popl_78$ï..Country == country]
}

#remove the columns
# Calculate the total number of cells in each column
total_cells <- colSums(!is.na(df_pop))

# Calculate the total number of NA values in each column
na_counts <- colSums(is.na(df_pop))

# Calculate the percentage of NA values in each column
#17 = observation of thr df_pop
na_percentage <- na_counts*100 / 17

# Set the threshold as 10%
threshold <- 40

# Identify columns exceeding the threshold
columns_to_remove <- names(na_percentage[na_percentage > threshold])

# Remove the identified columns from the data frame
df_to_NW <- df_pop[, !(names(df_pop) %in% columns_to_remove)]

# get the average of df_to_NW to load to cytoescape 

# Calculate the average of each column excluding NA values, omitting column 1
averages_K_pnu_ESBL <- colMeans(df_to_NW[, -1], na.rm = TRUE)
averages_K_pnu_ESBL <- data.frame(averages_K_pnu_ESBL)

write.csv(averages_K_pnu_ESBL, file = 'averages_K_pnu_ESBL.csv')

#input for the network analysis

data2 <- df_to_NW

my_data_2 <- data2[-1]

my_data_2

dim(my_data_2)

df <- data.frame(matrix(ncol = 17, nrow = 38)) #change the number of columns and rows
samples = data2[,1]


res <- cor(my_data_2)
round(res, 2)

cor(my_data_2, use = "complete.obs")

library("Hmisc")
res2 <- rcorr(as.matrix(my_data_2))
# Extract the correlation coefficients
res2$r
# Extract p-values
res2$P

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

library(Hmisc)
res2<-rcorr(as.matrix(my_data_2[,1:38])) # set the column numbers 1:N of my_data
flattenCorrMatrix(res2$r, res2$P)
sq_cor_p <-  flattenCorrMatrix(res2$r, res2$P)

head(sq_cor_p)

data_P_fil <- sq_cor_p[sq_cor_p$p<=0.05,]

max(data_P_fil$cor)
min(data_P_fil$cor)


data_R_min <- subset(data_P_fil, data_P_fil$cor >= -1 & data_P_fil$cor <= -0.5)
data_R_max <- subset(data_P_fil, data_P_fil$cor >= 0.5 & data_P_fil$cor <= 1)

df_to_ctoescape <- rbind(data_R_min, data_R_max)

dim(df_to_ctoescape) ### input for the cytoescape averages_K_pnu_ESBL

write.csv(df_to_ctoescape, file = 'K_pnu_ESBL_input_to_cytoescape.csv')
#_________________________________________________________________________________
#3

#________________________________________________________________________________________________
#Species == "Escherichia coli", Phenotype == "ESBL"

# Filter the data for the years 2010 to 2020 and the given species
filtered_df <- data %>%
  filter(Year >= 2004 & Year <= 2020, Species == "Escherichia coli", Phenotype == "ESBL")

# Get the count of the given species for each year and country separately
species_count <- filtered_df %>%
  group_by(Year, Country) %>%
  summarise(count = n())
# Pivot the pass rate data frame back to wide format
E_col_ESBL <- species_count %>%
  pivot_wider(names_from = Country, values_from = count)
# normalizing with the Population 
popl_78 <- read.csv("population_data_78_countries.csv" )

df_pop <- E_col_ESBL
for (country in names(E_col_ESBL)[-1]) {
  df_pop[[country]] <- E_col_ESBL[[country]]*1E6 / popl_78$Population[popl_78$ï..Country == country]
}

#remove the columns
# Calculate the total number of cells in each column
total_cells <- colSums(!is.na(df_pop))

# Calculate the total number of NA values in each column
na_counts <- colSums(is.na(df_pop))

# Calculate the percentage of NA values in each column
#17 = observation of thr df_pop
na_percentage <- na_counts*100 / 17

# Set the threshold as 10%
threshold <- 40

# Identify columns exceeding the threshold
columns_to_remove <- names(na_percentage[na_percentage > threshold])

# Remove the identified columns from the data frame
df_to_NW <- df_pop[, !(names(df_pop) %in% columns_to_remove)]

# get the average of df_to_NW to load to cytoescape 

# Calculate the average of each column excluding NA values, omitting column 1
averages_E_col_ESBL <- colMeans(df_to_NW[, -1], na.rm = TRUE)
averages_E_col_ESBL <- data.frame(averages_E_col_ESBL)

write.csv(averages_E_col_ESBL, file = 'averages_E_col_ESBL.csv')

#input for the network analysis

data2 <- df_to_NW

my_data_2 <- data2[-1]

my_data_2

dim(my_data_2)

df <- data.frame(matrix(ncol = 17, nrow = 39)) #change the number of columns and rows
samples = data2[,1]


res <- cor(my_data_2)
round(res, 2)

cor(my_data_2, use = "complete.obs") # remove - error

library("Hmisc")
res2 <- rcorr(as.matrix(my_data_2))
# Extract the correlation coefficients
res2$r
# Extract p-values
res2$P

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

library(Hmisc)
res2<-rcorr(as.matrix(my_data_2[,1:38])) # set the column numbers 1:N of my_data
flattenCorrMatrix(res2$r, res2$P)
sq_cor_p <-  flattenCorrMatrix(res2$r, res2$P)

head(sq_cor_p)

data_P_fil <- sq_cor_p[sq_cor_p$p<=0.05,]

max(data_P_fil$cor)
min(data_P_fil$cor)


data_R_min <- subset(data_P_fil, data_P_fil$cor >= -1 & data_P_fil$cor <= -0.5)
data_R_max <- subset(data_P_fil, data_P_fil$cor >= 0.5 & data_P_fil$cor <= 1)

df_to_ctoescape <- rbind(data_R_min, data_R_max)

dim(df_to_ctoescape) ### input for the cytoescape averages_K_pnu_ESBL

write.csv(df_to_ctoescape, file = 'E_col_ESBL_input_to_cytoescape.csv')
#_________________________________________________________________________________

#4

#________________________________________________________________________________________________
#Species == "Pseudomonas aeruginosa", Phenotype == "ESBL"

# Filter the data for the years 2010 to 2020 and the given species
filtered_df <- data %>%
  filter(Year >= 2004 & Year <= 2020, Species == "Pseudomonas aeruginosa", Phenotype == "ESBL")

# Get the count of the given species for each year and country separately
species_count <- filtered_df %>%
  group_by(Year, Country) %>%
  summarise(count = n())
# Pivot the pass rate data frame back to wide format
Ps_aer_ESBL <- species_count %>%
  pivot_wider(names_from = Country, values_from = count)
# normalizing with the Population 
popl_78 <- read.csv("population_data_78_countries.csv" )

df_pop <- Ps_aer_ESBL
for (country in names(Ps_aer_ESBL)[-1]) {
  df_pop[[country]] <- Ps_aer_ESBL[[country]]*1E6 / popl_78$Population[popl_78$ï..Country == country]
}

#remove the columns
# Calculate the total number of cells in each column
total_cells <- colSums(!is.na(df_pop))

# Calculate the total number of NA values in each column
na_counts <- colSums(is.na(df_pop))

# Calculate the percentage of NA values in each column
#17 = observation of thr df_pop
na_percentage <- na_counts*100 / 17

# Set the threshold as 10%
threshold <- 40

# Identify columns exceeding the threshold
columns_to_remove <- names(na_percentage[na_percentage > threshold])

# Remove the identified columns from the data frame
df_to_NW <- df_pop[, !(names(df_pop) %in% columns_to_remove)]

# get the average of df_to_NW to load to cytoescape 

# Calculate the average of each column excluding NA values, omitting column 1
averages_Ps_aer_ESBL <- colMeans(df_to_NW[, -1], na.rm = TRUE)
averages_Ps_aer_ESBL <- data.frame(averages_Ps_aer_ESBL)

write.csv(averages_E_col_ESBL, file = 'averages_Ps_aer_ESBL.csv')

#input for the network analysis

data2 <- df_to_NW

my_data_2 <- data2[-1]

my_data_2

dim(my_data_2)

df <- data.frame(matrix(ncol = 9, nrow = 22)) #change the number of columns and rows
samples = data2[,1]


res <- cor(my_data_2)
round(res, 2)

cor(my_data_2, use = "complete.obs")

library("Hmisc")
res2 <- rcorr(as.matrix(my_data_2))
# Extract the correlation coefficients
res2$r
# Extract p-values
res2$P

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

library(Hmisc)
res2<-rcorr(as.matrix(my_data_2[,1:22])) # set the column numbers 1:N of my_data
flattenCorrMatrix(res2$r, res2$P)
sq_cor_p <-  flattenCorrMatrix(res2$r, res2$P)

head(sq_cor_p)

data_P_fil <- sq_cor_p[sq_cor_p$p<=0.05,]

max(data_P_fil$cor)
min(data_P_fil$cor)


data_R_min <- subset(data_P_fil, data_P_fil$cor >= -1 & data_P_fil$cor <= -0.5)
data_R_max <- subset(data_P_fil, data_P_fil$cor >= 0.5 & data_P_fil$cor <= 1)

df_to_ctoescape <- rbind(data_R_min, data_R_max)

dim(df_to_ctoescape) ### input for the cytoescape averages_K_pnu_ESBL

write.csv(df_to_ctoescape, file = 'Ps_aer_ESBL_input_to_cytoescape.csv')
#_________________________________________________________________________________

