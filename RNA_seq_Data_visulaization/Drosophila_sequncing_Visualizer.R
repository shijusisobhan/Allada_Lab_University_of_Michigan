rm(list=ls())

packages = c("shiny", "ggplot2", "ggpubr",  "gridExtra")

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



ui <- fluidPage(
  titlePanel("Gene expression plot"),
  sidebarLayout(position = "left",
                sidebarPanel("select the data set",
                             checkboxInput("MB247", "MB 247", value = T),
                             checkboxInput("MB_alphabeta", "MB alphabeta", value = F),
                             checkboxInput("MB_gamma", "MB gamma", value = F),
                             checkboxInput("R21H11", "R21H11", value = F),
                             checkboxInput("R65H02", "R65H02", value = F),
                             checkboxInput("R85C10", "R85C10", value = F),
                             checkboxInput("Brain_HC", "Brain HC", value = F),
                             checkboxInput("Mech_SD_3hr", "Mech SD 3hr", value = F),
                             checkboxInput("Mech_SD_6hr", "Mech SD 6hr", value = F),
                             checkboxInput("Mech_SD_12hr", "Mech SD 12hr", value = F),
                             checkboxInput("LNv", "LNv", value = F),
                             checkboxInput("DN1", "DN1", value = F),
                             checkboxInput("FatBody_New", "FatBody New", value = F),
                             checkboxInput("FatBody_old", "FatBody old", value = F),
                             checkboxInput("Circadian_prof", "Circadian profile", value = F),
                             textInput("gene_1", "Enter the gene name"),
                             radioButtons('Plot_type', 'Select the plot type', c("Bar plot", "Time series"))
                             
                ),
                mainPanel("main panel",
                          column(6,plotOutput(outputId="plotgraph", width="500px",height="400px"))
                )))




server <- function(input, output) {
  
 
  pt1 <- reactive({
    
    if (input$Plot_type=='Bar plot'){
      
    if (!input$MB247) return(NULL)
    Data<-read.csv('DGE_table_MB_247_all.csv')
    gene_1<-input$gene_1
    
    gene_ind<-match(gene_1,Data$gene)
    
    
    
    DF<-Data[gene_ind,]
    DF<-DF[-1]
    DF<-DF[1:((length(DF))-2)]
    
    DF=stack(DF)
    DF$ind=as.character(DF$ind)
    DF$ind=c(rep(DF$ind[1],12),rep(DF$ind[13],4))
    
  
    
    p=ggbarplot(DF, x = "ind", y = "values",
                xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                title = 'MB 247',
                add = c("mean_se", "jitter"),
                fill = "ind", palette = "jco")
    
    p + theme(legend.position = "none")}
    
    else{ 
      
      if (!input$MB247) return(NULL)
      
      Data<-read.csv('DGE_table_MB_247_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      
      p_val<-DF$ADJ.P
      q_val<-DF$BH.Q
      DF<-DF[1:((length(DF))-2)]
      mx_exp<-max(DF)
      
      
      
      DF<-DF[1:12]
      t_MB247<-seq(0,22,2)
      
      TPM=as.numeric(t(DF))
      DF_new<-data.frame(ZT=t_MB247,TPM=TPM)
      
      ggplot(DF_new, aes(x=ZT, y=TPM)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('MB 247') + 
        annotate("text", x=15, y=mx_exp, label= paste("p=",p_val, ', q=', q_val))
      
      
      }
    
  })
  
  pt2 <- reactive({
    
    if(input$Plot_type=='Bar plot') {
      
    if (!input$MB_alphabeta) return(NULL)
   
    Data<-read.csv('DGE_table_MB_alphabeta_all.csv')
    gene_1<-input$gene_1
    
    gene_ind<-match(gene_1,Data$gene)
    
    
    
    DF<-Data[gene_ind,]
    DF<-DF[-1]
    DF<-DF[1:((length(DF))-2)]
    
    DF=stack(DF)
    DF$ind=as.character(DF$ind)
    DF$ind=c(rep(DF$ind[1],12),rep(DF$ind[13],4))
    
    
    
    p=ggbarplot(DF, x = "ind", y = "values",
                xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                title = 'MB alphabeta',
                add = c("mean_se", "jitter"),
                fill = "ind", palette = "jco")
    
    p + theme(legend.position = "none")
    }
    
    else{
      if (!input$MB_alphabeta) return(NULL)
      
      Data<-read.csv('DGE_table_MB_alphabeta_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      p_val<-DF$ADJ.P
      q_val<-DF$BH.Q
      DF<-DF[1:((length(DF))-2)]
      mx_exp<-max(DF)
      
      
      DF<-DF[1:12]
      t_MBalphabeta<-seq(0,22,2)
      
      TPM=as.numeric(t(DF))
      DF_new<-data.frame(ZT=t_MBalphabeta,TPM=TPM)
      
      ggplot(DF_new, aes(x=ZT, y=TPM)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('MB alphabeta')+
        annotate("text", x=15, y=mx_exp, label= paste("p=",p_val, ', q=', q_val))
      
    }
    
  }) 
  
  
  pt3 <- reactive({
    
    if (input$Plot_type=='Bar plot'){
      
      if (!input$MB_gamma) return(NULL)
      Data<-read.csv('DGE_table_MB_gamma_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      DF<-DF[1:((length(DF))-2)]
      
      DF=stack(DF)
      DF$ind=as.character(DF$ind)
      DF$ind=c(rep(DF$ind[1],12),rep(DF$ind[13],4))
      
   
      
      p=ggbarplot(DF, x = "ind", y = "values",
                  xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                  title = 'MB gamma',
                  add = c("mean_se", "jitter"),
                  fill = "ind", palette = "jco")
      
      p + theme(legend.position = "none")}
    
    else{ 
      
      if (!input$MB_gamma) return(NULL)
      
      Data<-read.csv('DGE_table_MB_gamma_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      p_val<-DF$ADJ.P
      q_val<-DF$BH.Q
      DF<-DF[1:((length(DF))-2)]
      mx_exp<-max(DF)
      
      DF<-DF[1:12]
      t_MBgamma<-seq(0,22,2)
      
      TPM=as.numeric(t(DF))
      DF_new<-data.frame(ZT=t_MBgamma,TPM=TPM)
      
      ggplot(DF_new, aes(x=ZT, y=TPM)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('MB gamma')+
        annotate("text", x=15, y=mx_exp, label= paste("p=",p_val, ', q=', q_val))
      
    }
    
  })
  
  
  
  pt4 <- reactive({
    
    if (input$Plot_type=='Bar plot'){
      
      if (!input$R21H11) return(NULL)
      Data<-read.csv('DGE_table_R21H11_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      
      DF=stack(DF)
      DF$ind=as.character(DF$ind)
      DF$ind=c(rep(DF$ind[1],6),rep(DF$ind[7],3))
      
    
      
      p=ggbarplot(DF, x = "ind", y = "values",
                  xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                  title = 'R21H11',
                  add = c("mean_se", "jitter"),
                  fill = "ind", palette = "jco")
      
      p + theme(legend.position = "none")}
    
    else{ 
      
      return(NULL)      
    }
    
  })
  
  pt5 <- reactive({
    
    if (input$Plot_type=='Bar plot'){
      
      if (!input$R65H02) return(NULL)
      Data<-read.csv('DGE_table_R65H02_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      
      DF=stack(DF)
      DF$ind=as.character(DF$ind)
      DF$ind=c(rep(DF$ind[1],7),rep(DF$ind[8],3))
      
     
      
      p=ggbarplot(DF, x = "ind", y = "values",
                  xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                  title = 'R65H02',
                  add = c("mean_se", "jitter"),
                  fill = "ind", palette = "jco")
      
      p + theme(legend.position = "none")}
    
    else{ 
      
      return(NULL)      
    }
    
  })
  
  
  pt6 <- reactive({
    
    if (input$Plot_type=='Bar plot'){
      
      if (!input$R85C10) return(NULL)
      Data<-read.csv('DGE_table_R85C10_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      
      DF=stack(DF)
      DF$ind=as.character(DF$ind)
      DF$ind=c(rep(DF$ind[1],6),rep(DF$ind[7],3))
      
    
      
      p=ggbarplot(DF, x = "ind", y = "values",
                  xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                  title = 'R85C10',
                  add = c("mean_se", "jitter"),
                  fill = "ind", palette = "jco")
      
      p + theme(legend.position = "none")}
    
    else{ 
      
      return(NULL)      
    }
    
  })
  
  pt7 <- reactive({
    
    if (input$Plot_type=='Bar plot'){
      
      if (!input$Brain_HC) return(NULL)
      Data<-read.csv('DGE_table_Brain_HC_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      
      DF=stack(DF)
      DF$ind=as.character(DF$ind)
      DF$ind=c(rep(DF$ind[1],6),rep(DF$ind[7],3))
      
  
      
      p=ggbarplot(DF, x = "ind", y = "values",
                  xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                  title = 'Brain HC',
                  add = c("mean_se", "jitter"),
                  fill = "ind", palette = "jco")
      
      p + theme(legend.position = "none")}
    
    else{ 
      
      return(NULL)      
    }
    
  })
  
  
  pt8 <- reactive({
    
    if (input$Plot_type=='Bar plot'){
      
      if (!input$Mech_SD_3hr) return(NULL)
      Data<-read.csv('DGE_table_Mech_SD_3hr_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      
      DF=stack(DF)
      DF$ind=as.character(DF$ind)
      DF$ind=c(rep(DF$ind[1],3),rep(DF$ind[4],3))
      
  
      
      p=ggbarplot(DF, x = "ind", y = "values",
                  xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                  title = 'Mech SD 3hr',
                  add = c("mean_se", "jitter"),
                  fill = "ind", palette = "jco")
      
      p + theme(legend.position = "none")}
    
    else{ 
      
      return(NULL)      
    }
    
  })
  
  
  pt9 <- reactive({
    
    if (input$Plot_type=='Bar plot'){
      
      if (!input$Mech_SD_6hr) return(NULL)
      Data<-read.csv('DGE_table_Mech_SD_6hr_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      
      DF=stack(DF)
      DF$ind=as.character(DF$ind)
      DF$ind=c(rep(DF$ind[1],3),rep(DF$ind[4],3))
      
    
      
      p=ggbarplot(DF, x = "ind", y = "values",
                  xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                  title = 'Mech SD 6hr',
                  add = c("mean_se", "jitter"),
                  fill = "ind", palette = "jco")
      
      p + theme(legend.position = "none")}
    
    else{ 
      
      return(NULL)      
    }
    
  })
  
  pt10 <- reactive({
    
    if (input$Plot_type=='Bar plot'){
      
      if (!input$Mech_SD_12hr) return(NULL)
      Data<-read.csv('DGE_table_Mech_SD_12hr_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      
      DF=stack(DF)
      DF$ind=as.character(DF$ind)
      DF$ind=c(rep(DF$ind[1],3),rep(DF$ind[4],3))
      
     
      
      p=ggbarplot(DF, x = "ind", y = "values",
                  xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                  title = 'Mech SD 12hr',
                  add = c("mean_se", "jitter"),
                  fill = "ind", palette = "jco")
      
      p + theme(legend.position = "none")}
    
    else{ 
      
      return(NULL)      
    }
    
  })
  
 #************************************************************************
 
  pt11 <- reactive({
    
    if (input$Plot_type=='Bar plot'){
      
      if (!input$LNv) return(NULL)
      Data<-read.csv('DGE_table_LNv_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)

      DF<-Data[gene_ind,]
      DF<-DF[-1]
      DF<-DF[-c(13,14,27,28,41,42)]
      
      DF=stack(DF)
      DF$ind=as.character(DF$ind)
      DF$ind=c(rep(DF$ind[1],12),rep(DF$ind[13],12),rep(DF$ind[25],12))
 
      
      p=ggbarplot(DF, x = "ind", y = "values",
                  xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                  title = 'LNv',
                  add = c("mean_se", "jitter"),
                  fill = "ind", palette = "jco")
      
      p + theme(legend.position = "none")}
    
    else{ 
      
      if (!input$LNv) return(NULL)
      
      Data<-read.csv('DGE_table_LNv_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      
      HC18_p_val<-DF$HC18_ADJ.P
      HC18_q_val<-DF$HC18_BH.Q
      HC25_p_val<-DF$HC25_ADJ.P
      HC25_q_val<-DF$HC25_BH.Q
      LC25_p_val<-DF$LC25_ADJ.P
      LC25_q_val<-DF$LC25_BH.Q
      
      DF<-DF[-c(13,14,27,28,41,42)]
      
      mx_exp<-max(DF)
    
      DF1<-DF[1:12]
      DF2<-DF[13:24]
      DF3<-DF[25:36]
      mx_exp_HC18<-max(DF1)
      mx_exp_HC25<-max(DF2)
      mx_exp_LC25<-max(DF3)
      
      t_LNv<-seq(2,24,2)
      
      HC18=as.numeric(t(DF1))
      DF_new1<-data.frame(ZT=t_LNv,HC18=HC18)
      
      pp1<-ggplot(DF_new1, aes(x=ZT, y=HC18)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('LNv_HC18')+
        annotate("text", x=15, y=mx_exp_HC18, label= paste("p=",HC18_p_val, ', q=', HC18_q_val))
      
      HC25=as.numeric(t(DF2))
      DF_new2<-data.frame(ZT=t_LNv,HC25=HC25)
      
      pp2<-ggplot(DF_new2, aes(x=ZT, y=HC25)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('LNv_HC25')+
        annotate("text", x=15, y=mx_exp_HC25, label= paste("p=",HC25_p_val, ', q=', HC25_q_val))
      
      LC25=as.numeric(t(DF3))
      DF_new3<-data.frame(ZT=t_LNv,LC25=LC25)
      
      pp3<-ggplot(DF_new3, aes(x=ZT, y=LC25)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('LNv_LC25')+
        annotate("text", x=15, y=mx_exp_LC25, label= paste("p=",LC25_p_val, ', q=', LC25_q_val))
      
      grid.arrange(pp1,pp2,pp3,nrow=3)
      
    }
    
  })
  
  
#***************************************************************************
  
  pt12 <- reactive({
    
    if (input$Plot_type=='Bar plot'){
      
      if (!input$DN1) return(NULL)
      Data<-read.csv('DGE_table_DN1_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      DF<-DF[-c(12,13,26,27)]
      
      DF=stack(DF)
      DF$ind=as.character(DF$ind)
      DF$ind=c(rep(DF$ind[1],11),rep(DF$ind[12],12))
      
      
      
      p=ggbarplot(DF, x = "ind", y = "values",
                  xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                  title = 'DN1',
                  add = c("mean_se", "jitter"),
                  fill = "ind", palette = "jco")
      
      p + theme(legend.position = "none")}
    
    else{ 
      
      if (!input$DN1) return(NULL)
      
      Data<-read.csv('DGE_table_DN1_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      HC25_p_val<-DF$HC25_ADJ.P
      HC25_q_val<-DF$HC25_BH.Q
      LC25_p_val<-DF$LC25_ADJ.P
      LC25_q_val<-DF$LC25_BH.Q
      DF<-DF[-c(12,13,26,27)]
      
      mx_exp<-max(DF)
      
      
      DF1<-DF[1:11]
      DF2<-DF[12:23]
      
      HC25_mx_exp<-max(DF1)
      LC25_mx_exp<-max(DF2)
      
      
      t_HC<-seq(4,24,2)
      t_LC<-seq(2,24,2)
      
 
      
      HC25=as.numeric(t(DF1))
      DF_new1<-data.frame(ZT=t_HC,HC25=HC25)
      
      pp1<-ggplot(DF_new1, aes(x=ZT, y=HC25)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('DN1_HC25') +
        annotate("text", x=15, y=HC25_mx_exp, label= paste("p=",HC25_p_val, ', q=', HC25_q_val))
      
      
      
      LC25=as.numeric(t(DF2))
      DF_new2<-data.frame(ZT=t_LC,LC25=LC25)
      
      pp2<-ggplot(DF_new2, aes(x=ZT, y=LC25)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('DN1_LC25')+
        annotate("text", x=15, y=LC25_mx_exp, label= paste("p=",LC25_p_val, ', q=', LC25_q_val))
      
      
      
      grid.arrange(pp1,pp2,nrow=2)
      
    }
    
  })
  
  
  
  pt13 <- reactive({
    
    if (input$Plot_type=='Bar plot'){
      
      if (!input$FatBody_New) return(NULL)
      Data<-read.csv('DGE_table_fatbody_new_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      DF<-DF[-c(19,20,26,27,52,53,80,81)]
      
      DF=stack(DF)
      DF$ind=as.character(DF$ind)
      DF$ind=c(rep(DF$ind[1],18),rep(DF$ind[19],5),rep(DF$ind[24],24),rep(DF$ind[48],26))
      
    
      
      p=ggbarplot(DF, x = "ind", y = "values",
                  xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                  title = 'FB_new',
                  add = c("mean_se", "jitter"),
                  fill = "ind", palette = "jco")
      
      p + theme(legend.position = "none")}
    
    else{ 
      
      if (!input$FatBody_New) return(NULL)
      
      Data<-read.csv('DGE_table_fatbody_new_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      D12T18_p_val<-DF$D12T18_ADJ.P
      D12T18_q_val<-DF$D12T18_BH.Q
      D12T25_p_val<-DF$D12T25_ADJ.P
      D12T25_q_val<-DF$D12T25_BH.Q
      D45T18_p_val<-DF$D45T18_ADJ.P
      D45T18_q_val<-DF$D45T18_BH.Q
      D45T25_p_val<-DF$D45T25_ADJ.P
      D45T25_q_val<-DF$D45T25_BH.Q
      
      DF<-DF[-c(19,20,26,27,52,53,80,81)]
      
      
      DF1<-DF[1:18]
      DF2<-DF[19:23]
      DF3<-DF[24:47]
      DF4<-DF[48:73]
      
      D12T18_mx_exp<-max(DF1)
      D12T25_mx_exp<-max(DF2)
      D45T18_mx_exp<-max(DF3)
      D45T25_mx_exp<-max(DF4)
      
      
      t_D12T18<-c(seq(0,16,2), seq(20,36,2))
      t_D12T25<-c(seq(0,24,6))
      t_D45T18<-c(seq(6,42,2),54,60,66,78,84)
      t_D45T25<-c(seq(6,42,2),seq(54,90,6))
      
      D12T18=as.numeric(t(DF1))
      DF_new1<-data.frame(ZT=t_D12T18,D12T18=D12T18)
      
      pp1<-ggplot(DF_new1, aes(x=ZT, y=D12T18)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('FB_new-D12T18')+
        annotate("text", x=15, y=D12T18_mx_exp, label= paste("p=",D12T18_p_val, ', q=', D12T18_q_val))
      
      
      D12T25=as.numeric(t(DF2))
      DF_new2<-data.frame(ZT=t_D12T25,D12T25=D12T25)
      
      pp2<-ggplot(DF_new2, aes(x=ZT, y=D12T25)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('FB_new-D12T25')+
        annotate("text", x=15, y=D12T25_mx_exp, label= paste("p=",D12T25_p_val, ', q=', D12T25_q_val))
      
      D45T18=as.numeric(t(DF3))
      DF_new3<-data.frame(ZT=t_D45T18,D45T18=D45T18)
      
      
      pp3<-ggplot(DF_new3, aes(x=ZT, y=D45T18)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('FB_new-D45T18')+
        annotate("text", x=45, y=D45T18_mx_exp, label= paste("p=",D45T18_p_val, ', q=', D45T18_q_val))
      
      
      
      D45T25=as.numeric(t(DF4))
      DF_new4<-data.frame(ZT=t_D45T25,D45T25=D45T25)
      
      
      pp4<-ggplot(DF_new4, aes(x=ZT, y=D45T25)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('FB_new-D45T25')+
        annotate("text", x=45, y=D45T25_mx_exp, label= paste("p=",D45T25_p_val, ', q=', D45T25_q_val))
      
      grid.arrange(pp1,pp2,pp3,pp4,nrow=4)
      
    }
    
  })
  
  
#****************************************************************************

  
  pt14 <- reactive({
    
    if (input$Plot_type=='Bar plot'){
      
      if (!input$FatBody_old) return(NULL)
      Data<-read.csv('DGE_table_fatbody_old_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      DF<-DF[-c(13,14,38,39,64,65)]
      
      DF=stack(DF)
      DF$ind=as.character(DF$ind)
      DF$ind=c(rep(DF$ind[1],12),rep(DF$ind[13],23),rep(DF$ind[36],24))
   
      
      p=ggbarplot(DF, x = "ind", y = "values",
                  xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                  title = 'FB_old',
                  add = c("mean_se", "jitter"),
                  fill = "ind", palette = "jco")
      
      p + theme(legend.position = "none")}
    
    else{ 
      
      if (!input$FatBody_old) return(NULL)
      
      Data<-read.csv('DGE_table_fatbody_old_all.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      HC18_p_val<-DF$HC18_ADJ.P
      HC18_q_val<-DF$HC18_BH.Q
      HC25_p_val<-DF$HC25_ADJ.P
      HC25_q_val<-DF$HC25_BH.Q
      LC25_p_val<-DF$LC25_ADJ.P
      LC25_q_val<-DF$LC25_BH.Q
      
      DF<-DF[-c(13,14,38,39,64,65)]
      
      
      DF1<-DF[1:12]
      DF2<-DF[13:35]
      DF3<-DF[36:59]
      
      HC18_mx_exp<-max(DF1)
      HC25_mx_exp<-max(DF2)
      LC25_mx_exp<-max(DF3)
      
      
      
      t_HC18<-seq(2,24,2)
      t_HC25<-c(seq(2,24,2),seq(28,48,2))
      t_LC25<-seq(2,48,2)
      
      HC18=as.numeric(t(DF1))
      DF_new1<-data.frame(ZT=t_HC18,HC18=HC18)
      
      pp1<-ggplot(DF_new1, aes(x=ZT, y=HC18)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('FB_old_HC18')+
        annotate("text", x=15, y=HC18_mx_exp, label= paste("p=",HC18_p_val, ', q=', HC18_q_val))
      
      HC25=as.numeric(t(DF2))
      DF_new2<-data.frame(ZT=t_HC25,HC25=HC25)
      
      pp2<-ggplot(DF_new2, aes(x=ZT, y=HC25)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('FB_old_HC25')+
        annotate("text", x=30, y=HC25_mx_exp, label= paste("p=",HC25_p_val, ', q=', HC25_q_val))
      
      LC25=as.numeric(t(DF3))
      DF_new3<-data.frame(ZT=t_LC25,LC25=LC25)
      
      pp3<-ggplot(DF_new3, aes(x=ZT, y=LC25)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('FB_old_LC25')+
        annotate("text", x=30, y=LC25_mx_exp, label= paste("p=",LC25_p_val, ', q=', LC25_q_val))
      
      grid.arrange(pp1,pp2,pp3,nrow=3)
      
    }
    
  })
  
  
  pt15 <- reactive({
    
    if(input$Plot_type=='Bar plot') {
      
      if (!input$Circadian_prof) return(NULL)
      
      Data<-read.csv('DGE_table_Circadian_profile.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      DF<-DF[1:((length(DF))-2)]
      
      DF=stack(DF)
      DF$ind=as.character(DF$ind)
      DF$ind=c(rep(DF$ind[1],18))
      
      
      
      p=ggbarplot(DF, x = "ind", y = "values",
                  xlab = 'condition', ylab =paste(gene_1,'(TPM)'),
                  title = 'Circadian profile',
                  add = c("mean_se", "jitter"),
                  fill = "ind", palette = "jco")
      
      p + theme(legend.position = "none")
    }
    
    else{
      if (!input$Circadian_prof) return(NULL)
      
      Data<-read.csv('DGE_table_Circadian_profile.csv')
      gene_1<-input$gene_1
      
      gene_ind<-match(gene_1,Data$gene)
      
      
      
      DF<-Data[gene_ind,]
      DF<-DF[-1]
      p_val<-DF$CR_ADJ.P
      q_val<-DF$CR_BH.Q
      DF<-DF[1:((length(DF))-2)]
      mx_exp<-max(DF)
      
      
      DF<-DF[1:18]
      t_CR<-seq(0,68,4)
      
      TPM=as.numeric(t(DF))
      DF_new<-data.frame(ZT=t_CR,TPM=TPM)
      
      ggplot(DF_new, aes(x=ZT, y=TPM)) + 
        geom_line()+xlab('ZT')+ylab(paste(gene_1,'(TPM)'))+ggtitle('Circadian profile')+
        annotate("text", x=40, y=mx_exp, label= paste("p=",p_val, ', q=', q_val))
      
    }
    
  }) 
  
  
  #*************************************************************************************
  
  
  
  
  
  
  
  output$plotgraph = renderPlot({
    ptlist <- list(pt1(),pt2(), pt3(), pt4(),pt5(),pt6(),pt7(), pt8(),pt9(),pt10(),pt11(),
                   pt12(),pt13(),pt14(), pt15())
    #wtlist <- c(input$wt1,input$wt2,input$wt3)
    # remove the null plots from ptlist and wtlist
    to_delete <- !sapply(ptlist,is.null)
    ptlist <- ptlist[to_delete] 
    #wtlist <- wtlist[to_delete]
    n_pt<-length(ptlist)
    
    if (n_pt==0) return(NULL)
    grid.arrange(grobs=ptlist,ncol=3)
    
  }) 
  
  
}
  
 
# Create Shiny app ----
shinyApp(ui, server)
