
#*******************************************************
# created by shiju Sisobhan
# 9/25/2024
# Allada Lab, University of Michigan
# shijusis@umich.edu
#***********************************************************

#************************************************************************************************************************
# This is the program for combining the screening results for different RUN (SleepMat output)
# Requirement: Sleepmat results folders which have screening results for each run
# Folder location specification file (.xlsx) - information regarding the location of folders where screening results are
# It has two column- column-1 : Run number, column-2: folder location (path) to the sleepmat screening results
#*****************************************************************************************************************************

rm(list = ls())


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



main_excel_file<-file.choose() # Browse for folder location file


# Main function to combine data from all folders
combine_all_folders <- function(main_excel_file) {
  # Read the main Excel file to get the folder paths
  main_data <- read_excel(main_excel_file)
  
  # Initialize an empty data frame to store the combined results
  final_result <- data.frame()
  
  # Loop through each row in the main Excel file
  for (i in 1:nrow(main_data)) {
    folder_path <- main_data[[2]][i]
    run <- main_data[[1]][i]
    
    # List all Excel files in the folder
    files <- list.files(folder_path, pattern = "*.xls", full.names = TRUE)
    
    # Read the data from both files
    folder_data <- read_excel(files)
    
    # Add the Run column
    folder_data$Run <- run
    
    # Append the combined data to the final result
    final_result <- bind_rows(final_result, folder_data)
  }
  
  return(final_result)
}


# Combine data from all folders
final_result <- combine_all_folders(main_excel_file)

# Write the results
mainDir <- dirname(main_excel_file)

file_Name<-readline(prompt = "Enetr a file name to save the results: "); 

xlsx_fileName<-paste0(file_Name,'.xlsx')

xl_output<-file.path(mainDir, xlsx_fileName)

# Save the combined dataframe to a new Excel file
write.xlsx(final_result, xl_output, rowNames = FALSE)