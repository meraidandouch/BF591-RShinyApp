## Author: Taylor Falk
## tfalk@bu.edu
## BU BF591
## Assignment 7

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(shinyWidgets)
library(anicon)
library(shinydashboard)
library(fontawesome)
library(shinyjs)
library(tidyverse)


# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel(title=div(img(src="dna.png", height = 35, width = 35), "BF591 Assignment 7")),
  h5("To use this application, download the CSV deseq_res.csv from the data directory of this app's repository."),
  sidebarLayout(
    # Sidebar with a slider input
    sidebarPanel(
      tags$head(tags$style("
            .btn-file {
            background-color:#F575A2; 
            border-color: #F575A2;
            color: white;
            font-family: Arial, Helvetica, sans-serif; 
            font-size: 14px;
            font-weight: bold;
            box-shadow: 0 1px 3px rgba(0,0,0,0.12), 0 1px 2px rgba(0,0,0,0.24);
            }
            .progress-bar{
            background-color:#FFB27F;}")),
      fileInput("data", "Load differential expression results", accept = ".csv"),
      markdown ("A volcano plot can be generated with **log2 fold-change** on the x-axis and **p-adjusted** on the y-axis."),
      radioButtons("button1", "Choose the column for the x-axis",
                   c("baseMean", "log2FoldChange",'lfcSE', 'stat', 
                     'pvalue', 'padj')),
      radioButtons("button2", "Choose the column for the y-axis",
                   c("baseMean", "log2FoldChange",'lfcSE', 'stat', 
                     'pvalue', 'padj')),
      colourpicker::colourInput("col1", "Base point color", "#FFB27F"),
      colourpicker::colourInput("col2", "Highlight point color", "#9E0C7E"),
      tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background:#9107A9; 
                      border-top: 1px solid #9107A9; 
                      border-bottom: 1px solid #9107A9}" )),
      sliderInput("slider",
                  "Select the magnitude of the p adjusted coloring:",
                  min = -300,
                  max = 0,
                  value = -150),
      actionBttn(
        inputId = "Id114",
        label = "Plot", 
        style = "gradient",
        color = "danger",
        icon = icon("dna"))),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plot",
                 plotOutput("volcano")),
        tabPanel("Table",
                 tableOutput("table")),
       )
     )
   )
)
  

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    #' load_Data
    #'
    #' @details Okay this one is a little weird but bear with me here. This is 
    #' still a "function", but it will take no arguments. The `reactive({})` bit 
    #' says "if any of my inputs (as in, input$...) are changed, run me again". 
    #' This is useful when a user clicks a new button or loads a new file. In 
    #' our case, look for the uploaded file's datapath argument and load it with 
    #' read.csv. Return this data frame in the normal return() style.
    load_data <- reactive({
      req(input$data$datapath)
      data <- read.csv(input$data$datapath, header =TRUE)
      #data <- tibble::column_to_rownames(data)
      return(data)
    })
    
    #' Volcano plot
    #'
    #' @param dataf The loaded data frame.
    #' @param x_name The column name to plot on the x-axis
    #' @param y_name The column name to plot on the y-axis
    #' @param slider A negative integer value representing the magnitude of
    #' p-adjusted values to color. Most of our data will be between -1 and -300.
    #' @param color1 One of the colors for the points.
    #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
    #'
    #' @return A ggplot object of a volcano plot
    #' @details I bet you're tired of these plots by now. Me too, don't worry.
    #' This is _just_ a normal function. No reactivity, no bells, no whistles. 
    #' Write a normal volcano plot using geom_point, and integrate all the above 
    #' values into it as shown in the example app. The testing script will treat 
    #' this as a normal function.
    #' 
    #' !!sym() may be required to access column names in ggplot aes().
    #'
    #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
    volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
      vp <- dataf %>% as_tibble() %>% 
        mutate(significant = case_when( 
          padj < 1*10^(slider) ~ 'TRUE',
          padj > 1*10^(slider) ~ 'FALSE')) %>%
        ggplot() + 
        geom_point(mapping=aes(x=!!sym(x_name), y=-log10(!!sym(y_name)), color=significant)) + 
        theme_minimal() +
        scale_color_manual(values=c(color1, color2))  +
        theme(legend.position="bottom") + 
        ggtitle("Volcano plot of DESeq2 differential expression results")
      return(vp)
    }
    
    #' Draw and filter table
    #'
    #' @param dataf Data frame loaded by load_data()
    #' @param slider Negative number, typically from the slider input.
    #'
    #' @return Data frame filtered to p-adjusted values that are less than 
    #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits 
    #' displayed.
    #' @details Same as above, this function is a standard R function. Tests will 
    #' evaluate it normally. Not only does this function filter the data frame to 
    #' rows that are above the slider magnitude, it should also change the format 
    #' of the p-value columns to display more digits. This is so that it looks 
    #' better when displayed on the web page. I would suggest the function 
    #' `formatC()`
    #'
    #' @examples draw_table(deseq_df, -210)
    #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
    #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
    #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
    draw_table <- function(dataf, slider) {
      dataf$padj <- as.numeric(dataf$padj)
      dataf <- dplyr::filter(dataf, padj< (1*10^slider)) 
      dataf <- dataf %>% summarise(
        gene = X, 
        baseMean = formatC(baseMean),
        log2FC = formatC(log2FoldChange),
        lfcSE = formatC(lfcSE),
        stat = formatC(stat), 
        pvalue = formatC(pvalue),
        padj = formatC(padj)
      )
      return(dataf)
    }
    
    #' These outputs aren't really functions, so they don't get a full skeleton, 
    #' but use the renderPlot() and renderTabel() functions to return() a plot 
    #' or table object, and those will be displayed in your application.
    output$volcano <- renderPlot({
      volcano_plot(load_data(), input$button1, input$button2, input$slider, input$col1, input$col2)}) 
    
    # Same here, just return the table as you want to see it in the web page
    output$table <-renderTable({
      draw_table(load_data(), input$slider)
    }) 
}

# Run the application
shinyApp(ui = ui, server = server)
