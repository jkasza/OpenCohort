#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
library(shinyjs)
library(matrixcalc)
library(swCRTdesign)


source("OpenCohortFunctions.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Open cohort longitudinal cluster randomised trials"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        selectInput("design", label = ("Design type:"), 
                    choices = list("Stepped wedge" = 1, "CRXO" = 4, "Parallel" = 2, "Parallel with baseline" =5), selected = 1),
        numericInput("T",
                     "Number of periods:",
                     min = 1,
                     step = 1,
                     value = 4),
        numericInput("kseq",
                     "Number of clusters assigned to each sequence:",
                     min = 1,
                     step = 1,
                     value = 1),
        
        fileInput('file1', 'Upload a design matrix instead:',
                  accept=c('text/plain', '.txt', '.csv')),
        helpText("The file must be a comma separated .csv or .txt file consisting of 0s and 1s, with a column for each time period. Do not include row or column names."),
        actionButton('reset', 'Clear file'),
        
        numericInput("m",
                     "Number of participants/cluster-period:",
                     min = 1,
                     step = 1,
                     value = 10),
        selectInput("scheme", label = ("Sampling scheme:"), 
                    choices = list("Core group" = 1, "Closed population" = 2, "Rotation (in-for-p)" = 3), selected = 1),
        
        #radioButtons("TEhetYN", "Allow for treatment effect heterogeneity?",
        #             choices = list("No" =1, "Yes"=2), selected = 1),
        # conditionalPanel(
        #  condition = "(input.ClusDecay == '1') && (input$SubjDecay == `1')",
              numericInput("rho",
                    label = HTML(paste("Intra cluster correlation, &rho;:", sep="")),
                    min = 0.001,
                    max=1,
                    step = 0.001,
                    value = 0.33),
              numericInput("pi",
                    label = HTML(paste("Cluster autocorrelation, &pi;:", sep="")),
                    min = 0,
                    max=1,
                    step = 0.001,
                    value = 0.9,
                    ),
              numericInput("tau",
                    label = HTML(paste("Individual autocorrelation, &tau;:", sep="")),
                    min = 0,
                    max=1,
                    step = 0.001,
                    value = 0.7),
              numericInput("totvar",
                     label = HTML(paste("Total variance:", sep="")),
                     min = 0,
                     step = 0.001,
                     value =25),
        helpText("When allowing for decaying between-period correlations, 
                the cluster autocorrelation must be less than 1: if the 
                 cluster autocorrelation is 1, an error will occur."),
        radioButtons("ClusDecay", "Allow for decaying between period correlations?",
                     choices = list("No" =0, "Yes"=1), selected = 0),
        radioButtons("SubjDecay", "Allow for autoregressive subject-level errors?",
                     choices = list("No" =0, "Yes"=1), selected = 0),
        #conditionalPanel(
        #  condition = "input.TEhetYN == '2'",
        #  numericInput("icc_cc",
        #               label= HTML(paste("&rho;",tags$sub("CC"),": ICC of subjects in the same control period", sep="")),
        #               min = 0,
        #               max=1,
        #               step = 0.005,
        #               value = 0.05),
        #  numericInput("icc_tt",
        #               label= HTML(paste("&rho;",tags$sub("TT"),": ICC of subjects in the same intervention period", sep="")),
        #               min = 0,
        #               max = 1,
        #               step = 0.005,
        #               value = 0.05),
        #  numericInput("icc_ct",
        #               label= HTML(paste("&rho;",tags$sub("CT"),": ICC of subjects, one in control, the other in intervention", sep="")),
        #               min = 0,
        #               max=1,
        #               step = 0.005,
        #               value = 0.05),
        #  numericInput("tau2",
        #               label = HTML(paste("Individual autocorrelation, &tau;:", sep="")),
        #               min = 0,
        #               max=1,
        #               step = 0.001,
        #               value = 0.7),
        #  numericInput("totvar2",
        #               label = HTML(paste("Total variance in an intervention period:", sep="")),
        #               min = 0,
        #               step = 0.001,
        #               value =25) 
        #),
        numericInput("effsize",
                     label = HTML(paste("Difference to detect:", sep="")),
                     min = 0,
                     step = 0.001,
                     value = 2),
        numericInput("alpha",
                     label = HTML(paste("Type I error rate,  &alpha;:", sep="")),
                     min = 0,
                     max=1,
                     step = 0.001,
                     value = 0.05),
        conditionalPanel(condition="input.conditionedPanels==1",
                         numericInput("power",
                                      label = HTML(paste("Power:", sep="")),
                                      min = 0,
                                      max=1,
                                      step = 0.001,
                                      value = 0.8)
                         
        )
        
        
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(
          tabPanel("Number of clusters per sequence", value=1, 
                   plotlyOutput("samplesizeplot"), textOutput("sampsizedeets")),
           tabPanel("Power", value=2, 
                    plotlyOutput("powerplot"), textOutput("powerdeets")),
          tabPanel("Design matrix", value=2, tableOutput("designmatrix")),
          tabPanel("Contact Details", value=3, textOutput("Details"))
          
          , id = "conditionedPanels"
        )
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  mymatrix <- reactiveValues(data=NULL) 
  
  observe({
    req(input$file1)
    mymatrix$data <- read.csv(input$file1$datapath)
  })
  
  observeEvent(input$reset, {
    mymatrix$data <- NULL
    reset('file1')
  })
  
 
  
  
  output$designmatrix <- renderTable({
    if(input$design == 1){  
      Xmat <- SWdesmat(input$T)
      #replicate the rows kseq times:
      Xmat <- Xmat[rep(1:nrow(Xmat), each = input$kseq),]
    }
    else if (input$design == 2){
      Xmat <- plleldesmat(input$T)
      #replicate the rows kseq times:
      Xmat <- Xmat[rep(1:nrow(Xmat), each = input$kseq),]
    }
    else if (input$design == 4){
      Xmat <- CRXOdesmat(input$T)
      #replicate the rows kseq times:
      Xmat <- Xmat[rep(1:nrow(Xmat), each = input$kseq),]
    }
    else if (input$design == 5){
      Xmat <- pllelBLdesmat(input$T)
      #replicate the rows kseq times:
      Xmat <- Xmat[rep(1:nrow(Xmat), each = input$kseq),]
    }
    
    
    if(!is.null(mymatrix$data)) Xmat <-  mymatrix$data
    
    Xmat
  },digits=0)
  
  
  output$Details <- renderText({
    "This R Shiny app investigates the impact of the degree of openness of open cohort
    longitudinal cluster randomised trials on required sample sizes and power of planned studies.
    It is currently under development. Please contact Jessica Kasza, jessica.kasza@monash.edu
    with comments or questions."
  })
  
  output$powerdeets <- renderText({
    "This figure plots the power of the user-specified design to detect the given difference, for 
    expected retention rates ranging from 0 (each subject provides one measurement only) to 1 
    (closed cohort, with each subject providing measurements in all trial periods). When a rotation
    sampling scheme is selected, the power for all in-for-p designs, with p=1,...,T is plotted."
  })
  
  output$sampsizedeets <- renderText({
    "This figure plots the number of clusters that must be assigned to each sequence to detect
      the given difference with the specified power, for 
    expected retention rates ranging from 0 (each subject provides one measurement only) to 1 
    (closed cohort, with each subject providing measurements in all trial periods). When a rotation
    sampling scheme is selected, the number of clusters for all in-for-p designs, with p=1,...,T is plotted. Note that 
    if the specified design matrix has repeated rows, these are considered separately in calculations."
  })
  
  output$samplesizeplot <- renderPlotly({
    if(input$design == 1){  
      Xmat <- SWdesmat(input$T)
      #replicate the rows kseq times:
      Xmat <- Xmat[rep(1:nrow(Xmat), each = input$kseq),]
      }
    else if (input$design == 2){
      Xmat <- plleldesmat(input$T)
      #replicate the rows kseq times:
      Xmat <- Xmat[rep(1:nrow(Xmat), each = input$kseq),]
    }
    else if (input$design == 4){
      Xmat <- CRXOdesmat(input$T)
      #replicate the rows kseq times:
      Xmat <- Xmat[rep(1:nrow(Xmat), each = input$kseq),]
    }
    else if (input$design == 5){
      Xmat <- pllelBLdesmat(input$T)
      #replicate the rows kseq times:
      Xmat <- Xmat[rep(1:nrow(Xmat), each = input$kseq),]
    }
   
    
    if(!is.null(mymatrix$data)) Xmat <-  mymatrix$data
    
    
    
    
  #Note that sampsize outputs the total number of clusters, so is later divided by the
    #number of rows in Xmat
    if(input$ClusDecay==0 && input$SubjDecay == 0 ) {
      
      if(input$scheme==1){
        if(input$m >= 100)  x<- seq(0,1, 0.01)
        else if(input$m < 100) x <- seq(0,1,1/input$m)
        sampsize <- sapply(x, SampleSize_open, effsize=input$effsize, 
                         totvar=input$totvar, power=input$power, alpha=input$alpha, 
                         m=input$m, desmat = Xmat, rho=input$rho, pi=input$pi, tau=input$tau)
        sampsize <- sampsize/nrow(Xmat)
      }
      
      else if(input$scheme==2){
        x <- seq(0,1,0.01)
        sampsize <- sapply(x, SampleSize_open, effsize=input$effsize, 
                           totvar=input$totvar, power=input$power, alpha=input$alpha, 
                           m=input$m, desmat = Xmat, rho=input$rho, pi=input$pi, tau=input$tau)
        sampsize <- sampsize/nrow(Xmat)
      }
      else if(input$scheme==3){
        nT <- ncol(Xmat)
        x<- seq(1, nT, 1)
        p<-x
        sampsize <- sapply(p, INFORP_SampleSize_open, effsize=input$effsize, 
                           power=input$power, alpha=input$alpha, desmat = Xmat,
                           m=input$m, totvar=input$totvar, rho=input$rho, pi=input$pi, tau=input$tau,
                           clusdecay=input$ClusDecay, subjdecay=input$SubjDecay)
        sampsize <- sampsize/nrow(Xmat)
      }
      
    }
   
     else if(input$ClusDecay==1 || input$SubjDecay == 1){
      
       if(input$scheme==1){
        if(input$m >= 100)  x<- seq(0,1, 0.01)
        else if(input$m < 100) x <- seq(0,1,1/input$m)
        sampsize <- sapply(x, SampleSize_open_decays, effsize=input$effsize, 
                         power=input$power, alpha=input$alpha, desmat = Xmat,
                         m=input$m, totvar=input$totvar, rho=input$rho, pi=input$pi, tau=input$tau,
                         clusdecay=input$ClusDecay, subjdecay=input$SubjDecay)
       sampsize <- sampsize/nrow(Xmat)
      }
      else if(input$scheme==2){
         x<- seq(0,1, 0.01)
        sampsize <- sapply(x, SampleSize_open_decays, effsize=input$effsize, 
                           power=input$power, alpha=input$alpha, desmat = Xmat,
                           m=input$m, totvar=input$totvar, rho=input$rho, pi=input$pi, tau=input$tau,
                           clusdecay=input$ClusDecay, subjdecay=input$SubjDecay)
        sampsize <- sampsize/nrow(Xmat)
      }
      else if(input$scheme==3){
        nT <- ncol(Xmat)
        x<- seq(1, nT, 1)
        p<-x
        sampsize <- sapply(p, INFORP_SampleSize_open, effsize=input$effsize, 
                           power=input$power, alpha=input$alpha, desmat = Xmat,
                           m=input$m, totvar=input$totvar, rho=input$rho, pi=input$pi, tau=input$tau,
                           clusdecay=input$ClusDecay, subjdecay=input$SubjDecay)
  
        sampsize <- sampsize/nrow(Xmat)
      }
       
    }     
#  if(input$TEhetYN==1) {
#   sampsize <- sapply(x, SampleSize_open, effsize=input$effsize, 
#                    totvar=input$totvar, power=input$power, alpha=input$alpha, 
#                    m=input$m, desmat = Xmat, rho=input$rho, pi=input$pi, tau=input$tau)
#   sampsize <- sampsize/nrow(Xmat)
#    }
#    else if(input$TEhetYN==2) {
#    #something else!
#      #Generate the required variances from the ICCs
#      sig2eta <- input$totvar2*(1-input$icc_tt)*input$tau2
#      sige2   <- input$totvar2*(1-input$icc_tt)*(1-input$tau2)
#      sigu2   <- input$totvar2*(1-input$icc_tt)*input$icc_cc/(1-input$icc_cc)
#      siguv   <- input$totvar2*(input$icc_ct*sqrt((1-input$icc_tt)/(1-input$icc_cc))
#                                -input$icc_cc*(1-input$icc_tt)/(1-input$icc_cc))
#      sigv2   <- input$totvar2*(input$icc_tt 
#                                -2*input$icc_ct*sqrt((1-input$icc_tt)/(1-input$icc_cc))
#                                +input$icc_cc*(1-input$icc_tt)/(1-input$icc_cc))
#      
#      myvar <- matrix(c(sigu2, siguv, siguv, sigv2), nrow=2)
#      
#      validate(
#        need(is.positive.semi.definite(myvar), "Selected intra-cluster correlation values do not lead to a positive semi-definite variance matrix. Try a smaller value for rho_CT.")
#      )
#      sampsize <- sapply(x, SampleSize_open_TEhet, effsize=input$effsize, 
#                         power=input$power, alpha=input$alpha, Xmat = Xmat,
#                         m=input$m, sigu2=sigu2, sigv2=sigv2, siguv =siguv, 
#                         sig2E=sige2, sig2eta = sig2eta)
#      sampsize <- sampsize/nrow(Xmat)
#      
#    }     
    mydata <- as.data.frame(cbind(sampsize, x))
    
    if(input$scheme == 1) {
    myplot <- plot_ly(data=mydata, x=~x, y=~sampsize, type = 'scatter', mode = 'lines',
                      hoverinfo = 'text',
                      text = ~paste('</br> No. of clusters: ', round(sampsize,3),
                                    '</br> Exp. retention rate: ', x)) %>%
              layout(xaxis = list(title = "Expected retention rate"), yaxis = list(title = "Number of clusters per sequence")) 
    }
    else if(input$scheme == 2) {
      myplot <- plot_ly(data=mydata, x=~x, y=~sampsize, type = 'scatter', mode = 'lines',
                        hoverinfo = 'text',
                        text = ~paste('</br> No. of clusters: ', round(sampsize,3),
                                      '</br> Exp. proportion sampled: ', x)) %>%
        layout(xaxis = list(title = "Expected proportion sampled"), yaxis = list(title = "Number of clusters per sequence")) 
    }
    else if(input$scheme == 3) {
      myplot <- plot_ly(data=mydata, x=~x, y=~sampsize, type = 'scatter', mode = 'lines',
                        hoverinfo = 'text',
                        text = ~paste('</br> No. of clusters: ', round(sampsize,3),
                                      '</br> p: ', x)) %>%
        layout(xaxis = list(title = "p (of in-for-p)"), yaxis = list(title = "Number of clusters per sequence")) 
    }
    
    
     myplot 
    
  })
   
   
   output$powerplot <- renderPlotly({

     if(input$design == 1){  
       Xmat <- SWdesmat(input$T)
       #replicate the rows kseq times:
       Xmat <- Xmat[rep(1:nrow(Xmat), each = input$kseq),]
     }
     else if (input$design == 2){
       Xmat <- plleldesmat(input$T)
       #replicate the rows kseq times:
       Xmat <- Xmat[rep(1:nrow(Xmat), each = input$kseq),]
     }
     else if (input$design == 4){
       Xmat <- CRXOdesmat(input$T)
       #replicate the rows kseq times:
       Xmat <- Xmat[rep(1:nrow(Xmat), each = input$kseq),]
     }
     else if (input$design == 5){
       Xmat <- pllelBLdesmat(input$T)
       #replicate the rows kseq times:
       Xmat <- Xmat[rep(1:nrow(Xmat), each = input$kseq),]
     }
     
     
     if(!is.null(mymatrix$data)) Xmat <-  mymatrix$data
     
     if(input$scheme == 1) {
       if(input$m >= 100)  x <- seq(0,1, 0.01)
        else if(input$m < 100) x <- seq(0,1,1/input$m)     
     
       power <- sapply(x, Power_open_decays, effsize=input$effsize, 
                       totvar=input$totvar, alpha=input$alpha, 
                       m=input$m, desmat = Xmat, rho=input$rho, pi=input$pi, tau=input$tau,
                       clusdecay = input$ClusDecay, subjdecay=input$SubjDecay)
     }
     else if(input$scheme == 2) {
       x <- seq(0,1, 0.01)
       power <- sapply(x, Power_open_decays, effsize=input$effsize, 
                       totvar=input$totvar, alpha=input$alpha, 
                       m=input$m, desmat = Xmat, rho=input$rho, pi=input$pi, tau=input$tau,
                       clusdecay = input$ClusDecay, subjdecay=input$SubjDecay)
     }
     
     else if(input$scheme == 3) {
       nT <- ncol(Xmat)
       x<- seq(1, nT, 1)
       p<-x
       power <- sapply(p, Power_open_decaysINFORP, effsize=input$effsize, 
                       totvar=input$totvar, alpha=input$alpha, 
                       m=input$m, desmat = Xmat, rho=input$rho, pi=input$pi, tau=input$tau,
                       clusdecay = input$ClusDecay, subjdecay=input$SubjDecay)
     }
    # if(input$TEhetYN==1) {
    #  power <- sapply(x, Power_open, effsize=input$effsize, 
    #                 totvar=input$totvar, alpha=input$alpha, 
    #                 m=input$m, desmat = Xmat, rho=input$rho, pi=input$pi, tau=input$tau)
    # }
    # else if(input$TEhetYN==2){
    #   #Generate the required variances from the ICCs
    #   sig2eta <- input$totvar2*(1-input$icc_tt)*input$tau2
    #   sige2   <- input$totvar2*(1-input$icc_tt)*(1-input$tau2)
    #   sigu2   <- input$totvar2*(1-input$icc_tt)*input$icc_cc/(1-input$icc_cc)
    #   siguv   <- input$totvar2*(input$icc_ct*sqrt((1-input$icc_tt)/(1-input$icc_cc))
    #                             -input$icc_cc*(1-input$icc_tt)/(1-input$icc_cc))
    #   sigv2   <- input$totvar2*(input$icc_tt 
    #                  -2*input$icc_ct*sqrt((1-input$icc_tt)/(1-input$icc_cc))
    #                  +input$icc_cc*(1-input$icc_tt)/(1-input$icc_cc))
    #   
    #   
    #   myvar <- matrix(c(sigu2, siguv, siguv, sigv2), nrow=2)
    #   
    #   validate(
    #     need(is.positive.semi.definite(myvar), "Selected intra-cluster correlation values do not lead to a positive semi-definite variance matrix. Try a smaller value for rho_CT.")
    #   )
    #   power <- sapply(x, Power_open_TEhet, effsize=input$effsize, 
    #                      alpha=input$alpha, Xmat = Xmat,
    #                      m=input$m, sigu2=sigu2, sigv2=sigv2, siguv =siguv, 
    #                      sig2E=sige2, sig2eta = sig2eta)
    # }
     
     mydata <- as.data.frame(cbind(power, x))
     
     if(input$scheme == 1) {
    myplot <- plot_ly(data=mydata, x=~x, y=~power, type = 'scatter', mode = 'lines',
                      hoverinfo = 'text',
                      text = ~paste('</br> Power: ', round(power,4),
                                    '</br> Exp. retention rate: ', x))  %>%
      layout(xaxis = list(title = "Expected retention rate"), yaxis = list(title = "Power")) 
     }
     else if(input$scheme == 2) {
       myplot <- plot_ly(data=mydata, x=~x, y=~power, type = 'scatter', mode = 'lines',
                         hoverinfo = 'text',
                         text = ~paste('</br> Power: ', round(power,4),
                                       '</br> Exp. proportion sampled: ', x))  %>%
         layout(xaxis = list(title = "Expected proportion sampled"), yaxis = list(title = "Power")) 
     }
     else if(input$scheme == 3) {
       myplot <- plot_ly(data=mydata, x=~x, y=~power, type = 'scatter', mode = 'lines',
                         hoverinfo = 'text',
                         text = ~paste('</br> Power: ', round(power,4),
                                       '</br> p: ', x))  %>%
         layout(xaxis = list(title = "p (of in-for-p)"), yaxis = list(title = "Power")) 
     }
    
    
    myplot 
    
   })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

