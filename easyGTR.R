
# install.packages("shiny")
# package shiny is used to build interactive web applications (apps) straight from R

# version control
# R version 4.2.3 (2023-03-15)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 11.7.3


library(shiny)
library('phylosim')
library(ggplot2)


parameter_tabs <- tabsetPanel(
  id = "parameter_tabs",
  type = "hidden",
  tabPanel("GTR",
           value = 'GTR',
           sliderInput("a", "a", min = 1, max = 10, value = 1),
           sliderInput("b", "b", min = 1, max = 10, value = 2),
           sliderInput("c", "c", min = 1, max = 10, value = 3),
           sliderInput("d", "d", min = 1, max = 10, value = 1),
           sliderInput("e", "e", min = 1, max = 10, value = 2),
           sliderInput("f", "f", min = 1, max = 10, value = 3)
  ),
  tabPanel("K80", 
           value =  "K80",
           sliderInput("Alpha",  "Alpha",
                       min = 1, max = 10, value = 6),
           sliderInput("Beta", "Beta",
                       min = 1, max = 10, value = 2)
  ),
  tabPanel("HKY", 
           value =  "HKY",
           sliderInput("Alpha",  "Alpha",
                       min = 1, max = 10, value = 10),
           sliderInput("Beta", "Beta",
                       min = 1, max = 10, value = 2)
  ),
  tabPanel("TN93", 
           value =  "TN93",
           sliderInput("Alpha1",  "Alpha1",
                       min = 1, max = 10, value = 4),
           sliderInput("Alpha2",  "Alpha2",
                       min = 1, max = 10, value = 3),
           sliderInput("Beta", "Beta",
                       min = 1, max = 10, value = 2)
  ),
  tabPanel("F84", 
           value =  "F84",
           sliderInput("Kappa",  "Kappa",
                       min = 1, max = 10, value = 2)
  ),
  tabPanel(title = "Upload parameters file", 
           value = "upload",
           fileInput("file1", "Choose txt File",
                     multiple = FALSE,
                     accept = c("text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv"))
  )
)

ins_del_tabs <- tabsetPanel(
  id = "ins_del_tabs",
  tabPanel("Root sequence", 
           value =  "root",
           sliderInput("length",  "Length of the root sequence",
                       min = 1, max = 100, value = 20),
           
           sliderInput("trees",  "Number of generated random trees",
                       min = 1, max = 10, value = 3)
  ),
  tabPanel("Insertion and Deletion", 
           value =  "ins_del",
           sliderInput("deletion_rate","Rate of deletion",
                       min = 0, max = 1, value = 1),
           sliderInput("insertion_rate","Rate of insertion",
                       min = 0, max = 1, value = 1)
  ),
  
)

# ---- User interface ----
ui <- fluidPage(
  titlePanel('EasyGTR: A web application to simulate DNA by selected 
             type of GTR model'),
  sidebarLayout(
    sidebarPanel(
      helpText("The meaning for parameters of each model
               can be found in example_input.txt file on EasyGTR
               Github page."),
      selectInput("model_type",
                  label = "Choose a model",
                  choices = c("GTR", "K80","HKY","TN93",
                              "JC69","F81","F84","upload"),
                  selected = "K80"),
      
      parameter_tabs,
      br(),
      
      ins_del_tabs,
      br(),

      actionButton("run", "run")),
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Deletion", 
                           helpText(" Deletion tab includes a plot for 
                                    deleltion size of distribution, and
                                    a detailed summary for deletion."),
                           plotOutput("d_plot"),
                           br(),
                           verbatimTextOutput("d_print")),
                  tabPanel("Insertion", 
                           helpText("Insertion tab is similar with 
                                    Deletion tab"),
                           plotOutput("i_plot"),
                           br(),
                           verbatimTextOutput("i_print")),
                  tabPanel("Parameters",
                           helpText("In parameter tab, there is a rate 
                                    of matrix plot and a plot for 
                                    equilibrium distribution. It also 
                                    has a summary for rate parameters."),
                           plotOutput("p_plot"),
                           br(),
                           verbatimTextOutput("p_print")),
                  tabPanel("Site Rates Plot", 
                           helpText("Site Rates Plot tab has a plot for 
                                    total site rate according to their positions"),
                           plotOutput("s_plot")),
                  tabPanel("Subtitions", 
                           helpText("A phylogenetic tree is included in 
                                    Subtition tab"),
                           
                           plotOutput("subst_plot"),
                           br(),
                           verbatimTextOutput("subst_print")),
                  tabPanel("Simulation", 
                           helpText("In the last tab Simulation, there is a large plot
                                    with colored sites for each sequence. Meanwhile, 
                                    the txt result of aligned sequences could be 
                                    downloaded from Simulation tab."),
                           plotOutput("sim_plot"),
                           br(),
                           downloadButton("Download", "Download alignment sequences"))
      )
    )
  )
)





# ---- Server logic ----
server <- function(input, output, session) {
  
  
  observeEvent(input$model_type,{
    updateTabsetPanel(session,"parameter_tabs",selected = input$model_type)
  })
  
  observeEvent(input$model_type,{
    if (input$model_type == 'upload'){
      hideTab('ins_del_tabs','root')
      hideTab('ins_del_tabs','ins_del')
    }
  })
  
  #model_type <- reactive(input$model_type)
  

  
  data <-  eventReactive(input$run, {
      
      if (input$model_type != 'upload'){
        if (input$model_type == 'HKY'){
          p <- HKY(rate.params=list("Alpha"=input$Alpha,"Beta"=input$Beta),
                   base.freqs=c(4,3,2,1)/10)
          
        }
        else if (input$model_type == 'K80'){
          p <- K80(rate.params=list("Alpha"=input$Alpha,"Beta"=input$Beta),
                   base.freqs=c(4,3,2,1)/10)
        }
        else if (input$model_type == 'GTR'){
          p <- GTR(rate.params=list(
            'a' = input$a, 'b' = input$b, 'c' = input$c,
            'd' = input$d, 'e' = input$e, 'f' = input$f
          ),
          base.freqs=c(2,2,1,1)/6)
        }
        else if (input$model_type == 'TN93'){
          p <- TN93(rate.params=list("Alpha1"=input$Alpha1,
                                     "Alpha2"=input$Alpha2,"Beta"=input$Beta),
                    base.freqs=c(4,3,2,1)/10)
        }
        else if (input$model_type == 'JC69'){
          p <- JC69()
        }
        else if (input$model_type == 'F81'){
          p<-F81(base.freqs=c(1,2,3,4)/10)
        }
        else if (input$model_type == 'F84'){
          p<-F84(rate.params=list("Kappa"=input$Kappa), base.freqs=c(1,2,3,4))
        }
        
        length = input$length
        deletion_rate = input$deletion_rate
        insertion_rate = input$insertion_rate
        trees =  input$trees
      }
      
      
      
      else if (input$model_type == 'upload'){
        file <- read.table(input$file1$datapath, 
                           header=T, fill = T, stringsAsFactors = FALSE,
                           comment.char = "#", sep="\t")
        
        if (file$model == 'HKY'){
          p <- HKY(rate.params=list("Alpha"=file$Alpha,"Beta"=file$Beta),
                   base.freqs=c(4,3,2,1)/10)
        }
        
        if (file$model == 'K80'){
          p <- K80(rate.params=list("Alpha"=file$Alpha,"Beta"=file$Beta),
                   base.freqs=c(4,3,2,1)/10)
        }
        
        if (file$model == 'TN93'){
          p <- TN93(rate.params=list("Alpha1"=file$Alpha1,
                                    "Alpha2"=file$Alpha2,
                                    "Beta"=file$Beta),
                   base.freqs=c(4,3,2,1)/10)
        }
        
        else if (file$model == 'GTR'){
          p <- GTR(rate.params=list(
            'a' = file$a, 'b' = file$b, 'c' = file$c,
            'd' = file$d, 'e' = file$e, 'f' = file$f
          ),
          base.freqs=c(2,2,1,1)/6)
        }
        else if (file$model == 'JC69'){
          p <- JC69()
        }
        else if (file$model == 'F81'){
          p<-F81(base.freqs=c(1,2,3,4)/10)
        }
        else if (input$model_type == 'F84'){
          p<-F84(rate.params=list("Kappa"=file$Kappa), base.freqs=c(1,2,3,4))
        }
        
        length = file$length
        deletion_rate = file$deletion_rate
        insertion_rate = file$insertion_rate
        trees =  file$trees
        
      }
      
      
      
      # create a sequence, attach process p
      s<-NucleotideSequence(length=length,processes=list(list(p)))
      # sample states
      sampleStates(s)
      # make the first five positions invariable
      setRateMultipliers(s,p,0,1:5)
      # get rate multipliers
      getRateMultipliers(s,p)
      
      ########## discrete deleltion ##########
      # construct deletion process object
      # proposing length in range 1:3
      d<-DiscreteDeletor(
        rate=deletion_rate,
        name="My_Deleletion",
        sizes=c(1:3),
        probs=c(3/6,2/6,1/6)
      )
      # attach d to s
      attachProcess(s,d)
      # create a region rejecting all deletions
      setDeletionTolerance(s,d,0,11:20)
      
      
      ########## discrete insertion ##########
      # construct insertion process object
      # proposing length in range 1:3
      i<-DiscreteInsertor(
        rate=insertion_rate,
        name="My_insertion",
        sizes=c(1:2),
        probs=c(1/2,1/2),
        template.seq=NucleotideSequence(length=1,processes=list(list(JC69())))
      ) 
      # states will be sampled from the JC69 equilibrium distribution
      # attach i to s
      attachProcess(s,i)
      # create a region rejecting all insertions
      setInsertionTolerance(s,i,0,11:20)
 
      # create a simulation object
      sim<-PhyloSim(root.seq=s,phylo=rcoal(trees))
      # run simulation
      Simulate(sim)
      
      # export the number of subtitions as a phylo object
      subst<-exportStatTree(sim,"substitution")
      
        return(
          list(d=d,i=i,p=p,s=s,subst=subst,sim=sim)
        )

      
  })
  
  
  
  output$d_print <- renderPrint({
    print(summary(data()$d))
  })
  
  output$d_plot <- renderPlot({
    plot(data()$d)
  })
  
  output$i_print <- renderPrint({
    print(summary(data()$i))
  })
  
  output$i_plot <- renderPlot({
    plot(data()$i)
  })
  
  output$p_print <- renderPrint({
    print(summary(data()$p))
  })
  
  output$p_plot <- renderPlot({
    plot(data()$p)
  })
  
  output$s_plot <- renderPlot({
    plot(data()$s)
  })
  
  output$subst_print <- renderPrint({
    print(summary(data()$subst))
  })
  
  output$subst_plot <- renderPlot({
    plot(data()$subst)
  })
  
  output$sim_plot <- renderPlot({
    plot(data()$sim)
  })
  
  output$Download <- downloadHandler(
    filename = 'alignment.txt',
    content = function(file) {
      write.table(data()$sim$alignment, file, row.names = T)
    }
  )
  
}
    


# ---- Run app ----
shinyApp(ui, server)
