
#********************************************************************************************
# This is the program for R-GUI which perform nparACT and carv analysis for Actigraphy files
# (c) Shiju Sisobhan (shiju.sisobhan@northwestern.edu)
# This version updated on 8/6/2020 
#********************************************************************************************





rm(list=ls()) 

## List of required packages

packages = c("shiny", "ggplot2", "na.tools", "shinyjs", "gridExtra", "VIM",
             "mice", "nparACT", "lattice", "stringr", "vctrs")


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


# Load required functions

source("Gui_Run_npar.R")
source("readActiwareCSVs_action.R")
source("subject_finder.R")
source("Run_carv_individual.R")
source("nparACT_carv.R")
source("nparACT_plot_minaverage_mod.R")
source("nparACT_plot_hraverage_mod.R")
source("Sujcet_info_finder.R")

## **********************************************************************##

# GUI Front Pannel strts

shinyUI(fluidPage(
  
  titlePanel("nparACT & carv Analysis"),
  
  sidebarLayout(
    
    sidebarPanel(("Enter the input"), 
                 
                

                 fileInput("file","Upload the files", multiple = TRUE), # Load actware files
                 options(shiny.maxRequestSize=100*1024^2),
                
                 textInput("ind_st_T", "Enter the start time for individual file run ", "12:00:00 PM"),
                 actionButton("I_Run", "Click for individual file run"),
                 
                 fileInput("spec.file","Upload the run specification file"), 
                 
                 
                 useShinyjs(),
                 hr(),
                 
                 
                  textInput("subject", "Enter the subject ID: "),
                 
                 
                 
                 
                 
                 actionButton("action", "Click for file information"),
                 
                 
                 
                 textInput("start.day", "Enter the start date of analysis(M/D/YYYY): "),
                 textInput("end.day", "Enter the end day of analysis(M/D/YYYY): "),
                 textInput("start.time", "Enter the start time of analysis(H:M:S PM/AM): ", "12:00:00 PM"),
                 textInput("end.time", "Enter the end time of analysis(H:M:S PM/AM): ", "12:00:00 PM"),
                 textInput("sampling.time", "Enter the sampling time (Minute): ")
                
                
                 
    ),
    
    
    mainPanel(("Results"),
              
              
              
              tabsetPanel(type="tab",
                          tabPanel("subject ID", tableOutput("subject.N"))),
              
              
              tabsetPanel(type="tab",
                          tabPanel("Run specification", tableOutput("Run.spec"))),
              
              

              tabsetPanel(type="tab",
                          tabPanel("Individual Result carv", tableOutput("Individual_Result_carv"))),
              
              
              tabsetPanel(type="tab",
                          tabPanel("Individual Result nparACT", tableOutput("Individual_Result_napar"))),
              
              downloadButton(outputId ="Result_individual.csv", label = "Download Individual file Results"),

              
              
              tabsetPanel(type="tab",
                          tabPanel("Subject informations", tableOutput("subject.information"))),
              
              
              textOutput("txt_error"),
              
              tags$head(tags$style("#txt_error{color: red;
                                   font-size: 20px;
                                   }"
                         )
              ),
              
              

              tabsetPanel(type="tab",
                          tabPanel("Parameter estimates of curve fit", tableOutput("Result_carv"))),
              
              
              
              textOutput("txt_warning"),
              
              tags$head(tags$style("#txt_warning{color: orange;
                                   font-size: 20px;
                                   }"
                         )
              ), 
              
              
              
              tabsetPanel(type="tab",
                          tabPanel("nparACT variables (Activity)", tableOutput("Run_npar_Activity"))),
              
              
             
              
              tabsetPanel(type="tab",
                          tabPanel("nparACT variables (White Light)", tableOutput("Run_npar_WL"))),
              


    
              
              downloadButton(outputId ="Result.csv", label = "Download Results"),
              


             downloadButton("plot.pdf", "Download the plot")
             
              

    ))))