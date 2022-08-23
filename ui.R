library(tidyverse)
library(shiny)
library(DT)
library(colourpicker)
library(geomtextpath)
library(shinyWidgets)
library(ggnewscale)

ui <- fluidPage(
  title = 'Codon Usage tables',
  titlePanel('Codon Usage tables'),
  sidebarLayout(
    position = 'right',
    sidebarPanel(
      width = 4,
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
      fluidRow(column(12, HTML('<h3>Genetic code sun</h3>'))), 
      HTML('<hr style="margin: 0 0 10px" />'),
      fluidRow(
        column(3, colourInput('g','Color A', value = '#9640FF',  closeOnClick = TRUE)),
        column(3, colourInput('h','Color T', value = '#FF389F', closeOnClick = TRUE)),
        column(3, colourInput('i','Color C', value = '#2B87FF', closeOnClick = TRUE)),
        column(3, colourInput('j','Color G', value = '#5FFF54', closeOnClick = TRUE))
      ),
      radioGroupButtons('byaa', 'Color amino acid by', status = 'primary',
                        choices = c('AA property', 'Codon Frequency'), selected = 'AA property'),
      plotOutput('plotsun', height = '60vh')
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
                   ),
                 )
          )),
      HTML('<hr style="margin: 3px 0 10px" />'),
      div(style="display: inline-block;vertical-align:top",
          fluidRow(column(12,HTML('<strong>Plot&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp</strong>')))),
      div(style="display: inline-block;vertical-align:top",
          materialSwitch('by_aa', HTML('AA-wise&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp'), value = TRUE, status = 'primary', right = TRUE, width = '100%')),
      div(style="display: inline-block;vertical-align:top",
          materialSwitch('norm', 'normalised for each AA', value = TRUE, status = 'primary', right = TRUE)),
      plotOutput('plot', ),
      HTML('<hr style="margin-bottom: 3px" />'),
      fluidRow(column(12,HTML('<strong>Table</strong>'))), br(),
      DTOutput('table')
    )
  )
)
