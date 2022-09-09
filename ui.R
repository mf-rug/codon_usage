library(tidyverse)
library(shiny)
library(DT)
library(colourpicker)
library(geomtextpath)
library(colorspace)
library(shinyWidgets)
library(ggnewscale)
library(units)
library(shinycssloaders)

ui <- fluidPage(
  title = 'Codon Usage tables',
  titlePanel('Codon Usage tables'),
  sidebarLayout(
    position = 'right',
    sidebarPanel(
      width = 4, style = "overflow-y:scroll; max-height: 100%; position:relative;",
      tags$head(tags$script('var dimension = [0, 0];
                                $(document).on("shiny:connected", function(e) {
                                    dimension[0] = window.innerWidth;
                                    dimension[1] = window.innerHeight;
                                    Shiny.onInputChange("dimension", dimension);
                                });
                                $(window).resize(function(e) {
                                    dimension[0] = window.innerWidth;
                                    dimension[1] = window.innerHeight;
                                    Shiny.onInputChange("dimension", dimension);
                                });')),
      fluidRow(
        column(5, 
               radioGroupButtons('byaa', 'Color AA by', status = 'primary',
                                 choices = c('type', 'codon freq'), selected = 'type')
        ),
        column(4, align = 'center', HTML('<div style="line-height:1.5; margin-bottom:4px;"><strong>Color bases</strong></div>'),
               dropdownButton(
                 column(3, colourInput('g','Color A', value = '#9640FF',  closeOnClick = TRUE)),
                 column(3, colourInput('h','Color T', value = '#FF389F', closeOnClick = TRUE)),
                 column(3, colourInput('i','Color C', value = '#2B87FF', closeOnClick = TRUE)),
                 column(3, colourInput('j','Color G', value = '#5FFF54', closeOnClick = TRUE)),
                 width = '30vw', status = 'primary', label = 'Colors', icon = icon('palette'), tooltip = 'Customise genetic code sun colors', right = TRUE
               )
        ),
        column(3, 
               fluidRow(HTML('<div style="line-height:1.5; margin-bottom:4px;"><strong>Highlight AA</strong></div>')),
               fluidRow(textInput('highlight', NULL, width = '90%', placeholder = 'e.g. KR'))
        )
      ),
      withSpinner(plotOutput('plotsun', height = '60vh', hover = hoverOpts(id ="plot_hover"))),
      fluidRow(
        column(width = 12, align="center",
               htmlOutput("hover_text"),
               imageOutput("hover_info")
        )), br(),br(),br(),br(),br(),
    ),
    mainPanel(
      width = 8,
      hr(),
      div(style="display: inline-block;vertical-align:top", HTML('<p style="line-height: 2.5;"><strong>Organism&nbsp</strong></p>')),
      div(style="display: inline-block;vertical-align:top, padding-left:100px",
          selectizeInput('species', NULL, choices = NULL, selected = NULL, multiple = FALSE,
                         options = list(placeholder = 'Search by typing',
                                        maxOptions = 200))),
      div(style="display: inline-block;vertical-align:top", HTML("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp")),
      
      div(style="display: inline-block;vertical-align:top;line-height: 2.5;", 
          materialSwitch('all_species', 'Show all organisms (can be slow)', value = FALSE, status = 'primary', right = TRUE))
      ,
      div(style="display: inline-block;vertical-align:top", 
          tags$a(href = 'https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs',
                 target = "_blank",
                 tags$p(
                   style = "line-height: 2.5;",
                   tags$i(
                     class = "glyphicon glyphicon-info-sign",
                     style = "color:#0072B2;",
                     title = 'Data from https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs'
                   )
                 )
          )),
      HTML('<hr style="margin: 3px 0 10px" />'),
      div(style="display: inline-block;vertical-align:top",
          fluidRow(column(12,HTML('<strong>Plot&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp</strong>')))),
      div(style="display: inline-block;vertical-align:top",
          materialSwitch('by_aa', HTML('AA-wise&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp'), value = TRUE, status = 'primary', right = TRUE, width = '100%')),
      div(style="display: inline-block;vertical-align:top",
          materialSwitch('norm', 'normalised for each AA', value = TRUE, status = 'primary', right = TRUE)),
      withSpinner(plotOutput('plot')),
      HTML('<hr style="margin-bottom: 3px" />'),
      fluidRow(column(12,HTML('<strong>Table</strong>'))), br(),
      DTOutput('table')
    )
  )
)
