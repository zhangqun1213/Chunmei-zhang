setwd("D:/R/TCGA_Clinical_Data")  # Set working directory

# Install required R packages ------------------------------------------------
# Check if 'compareGroups' is installed
if (!requireNamespace("compareGroups", quietly = TRUE)) {
  install.packages("compareGroups")
}

# Load required packages -----------------------------------------------------
library(compareGroups)

# Load example data ----------------------------------------------------------
load("TCGA-ESCA_clinical.rdata")

identical(rownames(clin_info), colnames(mrna_expr_counts))  # Check if row and column names match

# Extract the columns of interest --------------------------------------------
colnames(clin_info)  # View column names to decide which data to extract
# Example: xx = clin_info$xx (vector_name = your_clinical_dataframe$column_name)

# 1. Extract patient ID
ID <- clin_info$barcode
# 2. Extract patient survival status
status <- clin_info$vital_status
# 3. Extract patient survival time
time1 <- clin_info$days_to_last_follow_up
time2 <- clin_info$days_to_death
# Use 'days_to_death' for survival analysis; if NA, use 'days_to_last_follow_up'
# 4. Extract patient gender
gender <- clin_info$gender
# 5. Extract patient age
age <- clin_info$age_at_index
# 6. Extract patient staging information
stage <- clin_info$ajcc_pathologic_stage

stage_T <- clin_info$ajcc_pathologic_t
stage_N <- clin_info$ajcc_pathologic_n
stage_M <- clin_info$ajcc_pathologic_m

# Merge survival information -------------------------------------------------
mydata <- data.frame(ID,
                     status,
                     time1,
                     time2,
                     gender,
                     age,
                     stage,
                     stage_T,
                     stage_N,
                     stage_M)
data <- mydata  # Create a backup copy

# Process and filter clinical data -------------------------------------------
# 1. Set row names as sample names
rownames(data) <- data$ID
data <- subset(data, select = -ID)  # Remove the ID column

# 2. Process survival status (convert to factor/numeric variable)
table(data$status)
data$status <- as.factor(data$status)
data$status <- ifelse(data$status == "Alive", 0, 1)  # Convert to numeric

# 3. Process survival time (merge two time variables)
data$time <- ifelse(is.na(data$time2), data$time1, data$time2)
data <- subset(data, select = -c(time1, time2))  # Remove unnecessary columns

# 4. Process gender (binary variable)
table(data$gender)
data$gender <- as.factor(data$gender)
data$gender <- ifelse(data$gender == "female", 0, 1)  # Convert to numeric

# 5. Process age (continuous variable)
summary(data$age)
data$age <- as.numeric(data$age)
# Handle missing values, e.g., fill with mean/median
data[1, 5] <- NA  # Example: Set a missing value
data$age[is.na(data$age)] <- 61  # Fill missing values with 61

# 6. Process staging (categorical variables)
table(data$stage)
# Check for missing values (if T/N/M stages exist but overall stage is missing)
data[8, 4] <- "example"  # Example: Set a specific value

# The simplest way to handle missing values is to remove those samples
newdata <- data[!is.na(data$stage), ]  # Create a new dataframe excluding missing samples

# If you want to fill a specific value, e.g., replace unknown with 'Mx' in M stage
data$stage_M[is.na(data$stage_M)] <- "Mx"
data$stage[is.na(data$stage)] <- "stage x"

# Group the stages into fewer categories
table(data$stage)
data$stage_copy <- ifelse(data$stage %in% c("Stage I", "Stage IA", "Stage IB"), "Stage I",
                          ifelse(data$stage %in% c("Stage II", "Stage IIA", "Stage IIB"), "Stage II",
                                 ifelse(data$stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"), "Stage III",
                                        ifelse(data$stage %in% c("Stage IV", "Stage IVA"), "Stage IV",
                                               ifelse(data$stage == "stage x", "stage x", NA)))))
table(data$stage_copy)

# Save data and generate summary tables --------------------------------------
colnames(data)
data <- data[, c("status", "time", "gender", "age", "stage_copy", 
                 "stage_T", "stage_N", "stage_M")]
write.csv(data, "Clinical_Data.csv")

# Generate a simple summary table
tab_1 <- descrTable(status ~ .,  # Describe variables grouped by 'status'
                    data = data)
print(tab_1)

# Export the table to Word format
export2word(tab_1, file = "Clinical_Data_Table.docx")
