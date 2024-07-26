# Load necessary libraries
rm(list = ls())

library(readxl)
library(dplyr)
library(openxlsx)


# Specify the path to the main Excel file
# main_excel_file <- "C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/GitHub_folder/Allada_Lab_University_of_Michigan/sleepmat_Windows/Example/folder_location.xlsx"

main_excel_file<-file.choose()

# Function to combine data from a single folder
combine_folder_data <- function(folder_path) {
  # List all Excel files in the folder
  files <- list.files(folder_path, pattern = "*.xls", full.names = TRUE)
  
  # Read the data from both files
  data1 <- read_excel(files[1])
  data2 <- read_excel(files[2])
  
  # Merge the data
  combined_data <- inner_join(data1, data2, by = "Genotype")
  
  return(combined_data)
}

# Main function to combine data from all folders
combine_all_folders <- function(main_excel_file) {
  # Read the main Excel file to get the folder paths
  main_data <- read_excel(main_excel_file)
  
  # Initialize an empty data frame to store the combined results
  final_result <- data.frame()
  
  # Loop through each row in the main Excel file
  for (i in 1:nrow(main_data)) {
    folder <- main_data$folder[i]
    run <- main_data$Run[i]
    # Combine data from the current folder
    folder_data <- combine_folder_data(folder)
    
    # Add the Run column
    folder_data$Run <- run
    
    # Append the combined data to the final result
    final_result <- bind_rows(final_result, folder_data)
  }
  
  return(final_result)
}






# Combine data from all folders
final_result <- combine_all_folders(main_excel_file)


final_result_arranged<-tryCatch(
  
  {
    
    # final_result[,c("Run","Genotype","Total_sleep","Sleep_Pval","BoutNumber","BoutNumber_Pval",
    #                             "BoutLength","BoutLength_Pval",
    #                             "Total_activity","TotalActivity_Pval",
    #                             "Activity/waking min","Activity_wakingMin_Pval","Latency","Latency_Pval",
    #                             "Sleep_change_during_activation", "Sleep_change_during_activataion_Pval",
    #                             "Sleep_change_after activation","Sleep_change_after_activation_Pval",
    #                             "Latency_after_activation","Latency_after_activation_Pval")]
    
    
    final_result[,c("Run","Genotype","Total_sleep","Sleep_Pval",
                    "Sleep_change_during_activation", "Sleep_change_during_activataion_Pval",
                    "Sleep_change_after activation","Sleep_change_after_activation_Pval")]
    
    
  },
  
  error=function(e) {
    
  # final_result[,c("Run","Genotype","Total_sleep","Sleep_Pval","BoutNumber","BoutNumber_Pval",
  #                               "BoutLength","BoutLength_Pval","Total_activity","TotalActivity_Pval",
  #                               "Activity/waking min","Activity_wakingMin_Pval","Latency","Latency_Pval")]
    
    final_result[,c("Run","Genotype","Total_sleep","Sleep_Pval")]
    
  }
  
)

# Writr the results
mainDir <- dirname(main_excel_file)

file_Name<-readline(prompt = "Enetr a file name to save the results: "); 

xlsx_fileName<-paste0(file_Name,'.xlsx')

xl_output<-file.path(mainDir, xlsx_fileName)

# Save the combined dataframe to a new Excel file
write.xlsx(final_result_arranged, xl_output, rowNames = FALSE)