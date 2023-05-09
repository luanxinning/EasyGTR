
# install.packages("shiny")
# package shiny is used to build interactive web applications (apps) straight from R
library(shiny)
library('phylosim')
library(ggplot2)

  
  

parameter_tabs <- tabsetPanel(
  id = "parameter_tabs",
  type = "hidden",
  tabPanel("GTR",
           sliderInput("a", "a", min = 1, max = 10, value = 1),
           sliderInput("b", "b", min = 1, max = 10, value = 2),
           sliderInput("c", "c", min = 1, max = 10, value = 3)
  ),
  tabPanel("K80", 
           label = "K80",
           sliderInput("Alpha",  "Alpha",
                       min = 1, max = 10, value = 6),
           sliderInput("Beta", "Beta",
                       min = 1, max = 10, value = 2)
  ),
  tabPanel("HYK", 
           label = "HYK",
           sliderInput("Alpha",  "Alpha",
                       min = 1, max = 10, value = 10),
           sliderInput("Beta", "Beta",
                       min = 1, max = 10, value = 2)
  ),
  tabPanel("TN93", 
           label = "TN93",
           sliderInput("Alpha1",  "Alpha1",
                       min = 1, max = 10, value = 4),
           sliderInput("Alpha2",  "Alpha2",
                       min = 1, max = 10, value = 3),
           sliderInput("Beta", "Beta",
                       min = 1, max = 10, value = 2)
  ),
  tabPanel("F84", 
           label = "F84",
           sliderInput("Kappa",  "Kappa",
                       min = 1, max = 10, value = 2)
  ),
  tabPanel("upload", 
           label = "Upload parameters file",
           fileInput("file1", "Choose txt File",
                     multiple = FALSE,
                     accept = c("text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv"))
    
  )
)



# ---- User interface ----
ui <- fluidPage(
  titlePanel('simulate genome by selected type of model.'),
  sidebarLayout(
    sidebarPanel(
      selectInput("model_type",
                  label = "Choose a model",
                  choices = c("GTR", "K80","HYK","TN93",
                              "JC69","F81","F84","upload"),
                  selected = "upload"),

      parameter_tabs,

      actionButton("run", "run")),
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("p_Plot", plotOutput("p_plot")),
                  tabPanel("summary", verbatimTextOutput("p_print")),
                  tabPanel("simulation alignment", verbatimTextOutput("sim_print")),
                  tabPanel("sim_Plot", plotOutput("sim_plot"))
      )
    )
  )
)





# ---- Server logic ----
server <- function(input, output, session) {
  
  
  observeEvent(input$model_type, {
    updateTabsetPanel(inputId = "parameter_tabs", selected = input$model_type)
  }) 
  
  model_type <- reactive(input$model_type)
  

  
  data <-  eventReactive(input$run, {
      input$model_type
      input$Alpha
      input$Alpha1
      input$Alpha2
      input$Beta
      input$a
      input$b
      input$c
      input$Kappa
      input$file1
      
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
          'd' = input$a, 'e' = input$b, 'f' = input$c
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
            'd' = file$a, 'e' = file$b, 'f' = file$c
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
      }
      
      
      
      # create a sequence, attach process p
      s<-NucleotideSequence(length=20,processes=list(list(p)))
      # sample states
      sampleStates(s)
      # make the first five positions invariable
      setRateMultipliers(s,p,0,1:5)
      # get rate multipliers
      getRateMultipliers(s,p)
      # create a simulation object
      sim<-PhyloSim(root.seq=s,phylo=rcoal(2))
      # run simulation
      Simulate(sim)
      
        return(
          list(p=p,sim=sim)
        )

      
  })
  
  
  
  output$p_print <- renderPrint({
    print(summary(data()$p))
  })
  
  output$p_plot <- renderPlot({
    plot(data()$p)
  })
  
  output$sim_print <- renderPrint({
    d <- data()$sim$alignment
    print(d)
    })
  
  output$sim_plot <- renderPlot({
   plot(data()$sim)
  })
  
}
    


# ---- Run app ----
shinyApp(ui, server)