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
subDir<-'Combined_Run_files'


dir.create(file.path(mainDir, subDir))
New_file_location<-file.path(mainDir, subDir)

data <- read_excel(file_path)

new_column_names <- c("Run_number", "Monitor", "GT", "start_ch", "end_ch", "Analysis_date")

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
      df_selected <- df %>% select(1:6)
      # Rename the columns
      colnames(df_selected) <- new_column_names
      
      df_selected$Monitor<- as.numeric(paste0(run_number,df_selected$Monitor))
      
      
    }
    
  }
  
  return(df_selected)
}





# Iterate over each row in the data
for (i in 1:nrow(data)) {
  run_number <- data[[1]][i]  # Assuming the first column is "A"
  folder_location <- data[[2]][i]  # Assuming the second column is "B"
  
  # Rename, copy, and merge files in the folder
  df_selected <- rename_files(run_number, folder_location)
  
  # Append the dataframe to the list
  dataframes[[length(dataframes) + 1]] <- df_selected
  
  # Combine all dataframes by row
  combined_df <- bind_rows(dataframes)
  
  #combined_df$Monitor <- as.numeric(paste0(combined_df$Run_number, combined_df$Monitor))
  xl_output<-file.path(New_file_location, "combined_Genotype.xlsx")
  # Save the combined dataframe to a new Excel file
  write.xlsx(combined_df, xl_output, rowNames = FALSE)
  
}

print("Files integration has been sucessfully completed and rsults are in the folder combined_Run_files")



