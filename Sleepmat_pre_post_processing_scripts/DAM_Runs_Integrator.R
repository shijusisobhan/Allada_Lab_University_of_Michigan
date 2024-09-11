rm(list=ls())


packages = c("readxl", "dplyr", "openxlsx")

## Load or install the required packges
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


file_path<-file.choose()
mainDir <- dirname(file_path)


read_with_default <- function(prompt_text, default_value) {
  input <- readline(paste0(prompt_text, " [", default_value, "]: "))
  if (input == "") {
    return(default_value)
  } else {
    return(input)
  }
}

# Usage example
GT_col <- read_with_default("Enter the number of columns to select in the genotype specification file", "6")

GT_col<-as.numeric(GT_col)



subDir<-readline(prompt = "Enter a folder name to save the results: "); 
dir.create(file.path(mainDir, subDir))
New_file_location<-file.path(mainDir, subDir)


file_neme<-readline(prompt = "Enter a file name to save the combined geneotype: "); 
file_neme_to_save<-paste0(file_neme,'.xlsx')

#****************************************************************************************




data <- read_excel(file_path)

new_column_names <- c("Run_number", "Monitor", "GT", "start_ch", "end_ch", "Analysis_date", "Triage date", "Control_GT", "Environment_monitor")

# new_column_names <- c("Run_number", "Monitor", "GT", "start_ch", "end_ch", "Analysis_date")



dataframes <- list()

# Function to rename files in a folder
rename_files <- function(run_number, folder_location) {
  # List all files in the folder
  files <- list.files(folder_location, full.names = TRUE)
  
  merged_data <- data.frame()  # Initialize an empty dataframe for merging
  
  # Loop through each file and rename if it contains "Monitor"
  for (file in files) {
    # Extract file extension
    file_extension <- tools::file_ext(file)
    
    if (grepl("Monitor", basename(file), ignore.case = TRUE)) {
      
      
      # Generate new file name
      new_file_name <- gsub("Monitor", paste0("Monitor", run_number), basename(file), ignore.case = TRUE)
      new_file_path <- file.path(New_file_location, new_file_name)
      
      # Rename the file and copy to another location
      # file.rename(file, new_file_path)
      file.copy(file, new_file_path)
      
      
    }
    
    # If the file is an Excel file, read and merge it
    if (file_extension == "xlsx" || file_extension == "xls") {
      df <- read_excel(file)  # Read file and exclude the header
      # Select the first six columns
      df_selected <- df %>% select(1:GT_col)
      # Rename the columns
      colnames(df_selected) <- new_column_names[1:GT_col]
      
      df_selected$Monitor<- as.numeric(paste0(run_number,df_selected$Monitor))
      
      
    }
    
  }
  
  return(df_selected)
}



# Iterate over each row in the data
for (i in 1:nrow(data)) {
  run_number <- data[[1]][i]  # Assuming the first column is run_number
  folder_location <- data[[2]][i]  # Assuming the second column is folder location
  
  # Rename, copy, and merge files in the folder
  df_selected <- rename_files(run_number, folder_location)
  
  # Append the dataframe to the list
  dataframes[[length(dataframes) + 1]] <- df_selected
  
}


dataframes <- lapply(dataframes, function(df) {
  df$Run_number <- as.character(df$Run_number)
  return(df)
})


# Combine all dataframes by row
combined_df <- bind_rows(dataframes)







xl_output<-file.path(New_file_location, file_neme_to_save)
# Save the combined dataframe to a new Excel file
write.xlsx(combined_df, xl_output, rowNames = FALSE)



print("Files integration has been sucessfully completed")



