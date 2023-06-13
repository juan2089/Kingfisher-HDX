
library(tidyverse, warn=FALSE)
library(readxl, warn=FALSE)
library(openxlsx, warn=FALSE)
library(writexl, warn=FALSE)
library(rlist, warn=FALSE)
library(seqinr, warn=FALSE)
library(gridExtra, warn=FALSE)
library(grid, warn=FALSE)
library(devtools,warn=FALSE)
library(shiny,warn=FALSE)
library(shinydashboard,warn=FALSE)
library(shinyjs,warn=FALSE)
library(NGLVieweR,warn=FALSE)
library(shinyWidgets,warn=FALSE)
library(shinydisconnect,warn=FALSE)
library(shinyBS,warn=FALSE)
library(shinydashboardPlus,warn=FALSE)

ui <- dashboardPage(
  
  dashboardHeader(title="Kingfisher HDX-MS analytics",
                  dropdownMenu( type = 'tasks',
                                icon=icon("twitter"),
                                headerText= tags$a(href="http://twitter.com/intent/tweet?text=I%20used%20Kingfisher%20HDX%20app %20from%20@GrossLab_WashU","Share us on Twitter",target="_blank"))
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introduction", tabName="intro", icon=icon("house")),
      menuItem("Input data", tabName="input_dash", icon=icon("table")),
      menuItem("Statistical analysis", tabName="stat_dash", icon=icon("arrow-right")),
      menuItem("Fun. data analysis",tabName="bayesian_dash",icon=icon("arrow-right"),badgeLabel = "NEW",badgeColor = "red"),
      menuItem("Global woodsplot", tabName="globwoods_dash", icon=icon("arrow-right")),
      menuItem("Woodsplot by timepoint", tabName="woodsbytimepoint_dash", icon=icon("arrow-right"),
               selectInput(inputId="timepoint",label = "Select labeling time:",choices ="", selected=NULL),
               menuSubItem("Show woods plot", tabName="woodsplotbytime_sub",
                           icon = icon("line-chart"))),
      menuItem("Digestion efficiency", tabName="digest_dash", icon=icon("arrow-right")),
      menuItem("Peptide map", tabName="pep_dash", icon=icon("arrow-right")),
      menuItem("Uptake plots", tabName="uptake_dash", icon=icon("arrow-right"),
               selectInput(inputId="selpeptide",label = "Select peptide:",choices ="", selected=NULL),
               menuSubItem("Show uptake plots", tabName="uptake_plotsub",
                           icon = icon("line-chart"))
      ),
      menuItem("3D Structure", tabName="3d_dash", icon=icon("cube")),
      menuItem("Export results", tabName="save_dash", icon=icon("save")),
      menuItem("FAQ", tabName="faq_dash",icon=icon("question")),
      menuItem("About", tabName="about_dash",icon=icon("info")),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
      div(style="display:inline-block;width:100%;text-align: center;",socialButton(
          href = "http://twitter.com/intent/tweet?text=I%20used%20Kingfisher%20HDX%20app %20from%20@GrossLab_WashU",
          icon = icon("twitter")
        ), socialButton(
          href = "https://github.com/juan2089/Kingfisher-HDX",
          icon = icon("github")
        ))
    
  )),
  
  ## Body content
  dashboardBody(
    tags$style(
      '
        @media (min-width: 768px){
          .sidebar-mini.sidebar-collapse .main-header .logo {
              width: 230px; 
          }
          .sidebar-mini.sidebar-collapse .main-header .navbar {
              margin-left: 230px;
          }
        }
        '
    ),
    tags$script(HTML("$('body').addClass('fixed');")),
    disconnectMessage(
      text = "Something went wrong. Please refresh the page and check your input data",
      refresh = "Refresh",
      background = "#FFFFFF",
      colour = "#000000",
      refreshColour = "#0461B3",
      overlayColour = "#000000",
      overlayOpacity = 0.75,
      width = 450,
      top = "center",
      size = 22,
      css = ""
    ),
    tags$head(tags$link(rel = "shortcut icon", href = "favicon.ico")),
    tabItems(
      #intro tab
      tabItem(tabName="intro",
              fluidRow(
                column(5,img (src="kingfisher.png",height=225,width=225,style=""),offset=5),
                column(1,img (src="universitylogo.png",height=100,style=""))),
              h2("  Kingfisher HDX-MS is an app developed by JuaNolan at Washington University in St. Louis",align="left"),
              fluidPage(
                tags$iframe(src='./introduction.html',
                            width = '100%', height = '320px',
                            frameborder = 0, scrolling = "no")
              ),
              h3(HTML("<b>1. Input data</b>")),
              tags$iframe(style="text-align:center; height:700px; width:100%", src="inputdatainstructions.pdf",align="middle"),
              h3(HTML("<b>2. Statistical analysis</b>")),
              tags$iframe(style="height:700px; width:100%", src="statisticalanalysisinstructions.pdf"),
              h3(HTML("<b>3. Global woodsplot</b>")),
              tags$iframe(style="height:700px; width:100%", src="globalwoodsplotinstructions.pdf"),
              h3(HTML("<b>4. Woodsplot by timepoint</b>")),
              tags$iframe(style="height:700px; width:100%", src="woodsplotbytimepointinstructions.pdf"),
              h3(HTML("<b>5. Digestion efficiency</b>")),
              tags$iframe(style="height:700px; width:100%", src="digestionefficiencyinstructions.pdf"),
              h3(HTML("<b>6. Peptide map</b>")),
              tags$iframe(style="height:700px; width:100%", src="peptidemapinstructions.pdf"),
              h3(HTML("<b>7. Uptake plots</b>")),
              tags$iframe(style="height:700px; width:100%", src="uptakeplotsinstructions.pdf"),
              h3(HTML("<b>8. 3D structure</b>")),
              tags$iframe(style="height:700px; width:100%", src="3dstructureinstructions.pdf"),
              h3(HTML("<b>9. Export results</b>")),
              tags$iframe(style="height:700px; width:100%", src="exportinstructions.pdf"),
              br(),
              tags$a(href="mailto:j.rinconpabon@wustl.edu", "Juan P Rincon"),
              br(),
              tags$a(href="mailto:nolanmclaughlin@wustl.edu", "Nolan McLaughlin")
      ),
      #1st tab
      tabItem(tabName="input_dash",
              h2("Data import and experiment details",align="center"),
              br(),
              # Input: Select a file ----
              fileInput("hdexaminerfile", "Select exported CSV file from HDExaminer",
                        multiple = FALSE,
                        accept = c(".csv")),
              fluidRow(style="margin-top: -30px; margin-left:20px;",
                       checkboxInput("nothdexaminer",
                                     label=tags$span("Check if you are NOT using an HDExaminer output file",
                                                     bsButton("testingbutton",label="",style="info",icon=icon("info"),size="extra-small")),value=FALSE)),
              bsPopover(id = "testingbutton",title = HTML("<b> Important </b>"),
                        content = paste0("If you are not using an HDExaminer output file, please format your data as in the test file available for download in the bottom right link."),
                        placement = "right",
                        trigger = "hover",
                        options = list(container = "body")
              ),
              
              
              # Input: Select a file ----
              fileInput("fastafile", "Select FASTA file",
                        multiple = FALSE,
                        accept = c(".FASTA")),
              
              fluidRow(
                       
                       column(4,
                              numericInput("states", 
                                           label=HTML("Number of protein states"), 
                                           value = 2,min=2)), 
                       
                       column(4,
                              numericInput("timepoints", 
                                           label=HTML("Number of labeling times"), 
                                           value = 4)), 
                       
                       
                       column(3,
                              numericInput("significancelevel", 
                                           label=HTML("Significance level"), 
                                           value = 0.01,min=0,max=0.5,step=0.01))
                       
              ),
              
              numericInput("replicates", 
                           "Number of replicates per labeling time", 
                           value = 4,min=2),
              radioButtons("cantdeutbutton", 
                           "Number of residues that can't withhold deuteration:", 
                           choices = c(1, 2),selected=2),
              textInput("labelingtimepoints", "Labeling times in seconds (separated by commas):", 
                        value = "30,60,800,2500"),
              uiOutput("analyze"),
              hr(style = "border-top: 1px solid #000000;"),
              h5('Download a test file', tags$a(href = "inputfile.csv","here"), align="right"),
              column(2,uiOutput("state1box")),
              column(2,uiOutput("state2box"))
      ),
      
      #2nd tab
      tabItem(tabName="stat_dash",
              fluidPage(shinyjs::useShinyjs(),
                fluidRow(
                  uiOutput("volcanotitle"),
                  column(7,plotOutput(outputId="volcanoplot",width="70%"),align="right"),br(),
                  shinyjs::hidden( div(id = "advanced",
                                       column(4,selectInput("colorstr", 
                                       label = "Select color for strong effects:",
                                       choices = c("forestgreen", 
                                                   "darkmagenta",
                                                   "deeppink",
                                                   "red"
                                       ),
                                       selected = "forestgreen"),
                         selectInput("colorint", 
                                     label = "Select color for intermediate effects:",
                                     choices = c("orange", 
                                                 "cyan3",
                                                 "aquamarine",
                                                 "palegreen2"
                                     ),
                                     selected = "orange"),offset=-4)))),
                fluidRow(column(7,
                       conditionalPanel(
                         condition="output.volcanoplot",
                         checkboxInput("clusteringsel","Check to cluster the data",value=FALSE)
                       ),align="right")
                ),
                br(),br(),
                column(1, uiOutput("volcanotype"),offset=3),
                column(2, uiOutput("volcanodownload"),offset=1),
                column(2, uiOutput("histdownload"),offset=0),
                uiOutput("volcanotext"),
                column(12,plotOutput(outputId="histogramclustplot",width="50%"),align="center")),
              br(),fluidPage(
                column(12,uiOutput("binsize"),align="center")
                
              )),
      
      tabItem(tabName="bayesian_dash",
              fluidPage(shinyjs::hidden(div(id="functional",
               h1("Functional data analysis",align="center"),br(),
               tags$iframe(src='./functionalanalysis.html',
                           width = '100%', height = '100px',
                           frameborder = 0, scrolling = "no")
              ,
               tabBox(title="", width="100%",id="bayesianbox", selected="Heat map",
                      tabPanel("Heat map", "additional description?",
                              fluidRow(column(10, plotOutput(outputId="heatmap"),align="center",offset=1)), 
                                fluidRow(column(3, selectInput("heattype", 
                                            label = NULL,
                                            choices = c(".pdf", 
                                                        ".svg",
                                                        ".eps"
                                            ),
                                            selected = ".pdf"),offset=3
                              ),column(3,
                                downloadButton("heatbutton", "Download Heat map")))
                              ),
                      tabPanel("Fitted kinetic plots", "additional description for kinetic plots?",
                               fluidRow(column(8, plotOutput(outputId="kineticplot"),align="center",offset=2)),br(),
                               column(8,align="center",selectInput(inputId="selpeptide2",label = "Select peptide:",choices =""),
                                       offset=2),br(),br(),
                               fluidRow(column(4, selectInput("kinetictype", 
                                                              label = NULL,
                                                              choices = c(".pdf", 
                                                                          ".svg",
                                                                          ".eps"
                                                              ),
                                                              selected = ".pdf"),offset=3
                               ),column(3,downloadButton("kineticbutton", "Download kinetic plot")))),
                      tabPanel("Manhattan plot","moreinfo?",
                               fluidRow(
                                 column(6, plotOutput(outputId="kineticplot2")),
                                 column(6, plotOutput(outputId="forestplot"))
                                 ),br(),
                               column(8,align="center",selectInput(inputId="selpeptide3",label = "Select peptide:",choices =""),offset=2),br(),br(),
                               fluidRow(column(4, selectInput("manhattantype", 
                                                              label = NULL,
                                                              choices = c(".pdf", 
                                                                          ".svg",
                                                                          ".eps"
                                                              ),
                                                              selected = ".pdf"),offset=3
                               ),column(3,downloadButton("manhattanbutton", "Download Manhattan plot"))))
                
               ))))) ,
      #3rd tab
      tabItem(tabName="globwoods_dash",
              fluidPage(
                fluidRow(
                  uiOutput("globwoodstitle"),
                  column(8,plotOutput(outputId="globwoodsplot",width="50%"),align="right"),
                  column(3,br(),shinyjs::hidden(div(id="advanced2",selectInput("colorflex", 
                                       label = "Select color for peptides with increased flexibility:",
                                       choices = c("firebrick2", 
                                                   "cyan",
                                                   "deeppink"
                                                   ),
                                       selected = "firebrick2"),
                         selectInput("colorprot", 
                                     label = "Select color for peptides with increased protection:",
                                     choices = c("blue", 
                                                 "darkcyan",
                                                 "chocolate"
                                                 ),
                                     selected = "blue"),offset=-2)))),
                br(),
                column(3,uiOutput("globwoodstype"),offset=3),
                column(2,uiOutput("globwoodsdownload"),offset=2),br(),
                uiOutput("woodsplottext")
              )
              
      ),
      #4th tab
      tabItem(tabName="woodsplotbytime_sub",
              fluidPage(
                fluidRow(
                  uiOutput("woodsbytimepointtitle"),
                  column(12,plotOutput(outputId="woodsbytimepointplot",width="60%"),align="center")),
                br(),
                column(3,uiOutput("woodsbytimepointtype"),offset=3),
                column(2,uiOutput("woodsbytimepointdownload"),offset=2),br(),
                uiOutput("woodsplotbytimetext")
              )
              
      ),
      #5th tab
      tabItem(tabName="digest_dash",
              fluidPage(
                fluidRow(
                  uiOutput("digestiontitle"),
                  column(12,plotOutput(outputId="digestionplot",width="60%"),align="center")),
                br(),
                column(3,uiOutput("digestiontype"),offset=3),
                column(2,uiOutput("digestiondownload"),offset=2),br(),
                uiOutput("digestiontext")                
              )
      ),
      #6th tab
      tabItem(tabName="pep_dash",
              fluidPage(
                fluidRow(
                  uiOutput("peptidetitle"),
                  column(12,plotOutput(outputId="peptideplot",width="60%"),align="center")),
                br(),
                column(3,uiOutput("peptidetype"),offset=3),
                column(2,uiOutput("peptidedownload"),offset=2),br(),
                uiOutput("peptidemaptext"),
                column(12,shinyjs::hidden(div(id="pep",selectInput("colorpep", 
                            label = "Select color for peptides:",
                            choices = c("blue", 
                                        "darkcyan",
                                        "chocolate"
                            ),
                            selected = "blue"),align="center")))
              )
      ),
      #7th tab
      tabItem(tabName="uptake_plotsub",
              fluidPage(
                fluidRow(
                  column(2,shinyjs::hidden(div(id="plotssettings",
                         radioButtons(inputId="coloring",label="Select the desired color palette:",
                                      choices=list("Solarized"='solarized',"Chalk"='chalk', "Dust"='dust')),
                         sliderInput(inputId="thickness",
                                     label="Line thickness:",
                                     min=0.5,
                                     max=5,
                                     value=1),
                         numericInput("fontsize",label="Font size",value=10,min=8,max=25)))),
                  column(10,plotOutput(outputId="uptakeplot",width="80%"),align="center")),
                br(),
                column(3,uiOutput("uptaketype"),offset=4),
                column(2,uiOutput("singledownload",offset=4)),
                column(2,uiOutput("uptakedownload")),
                br(),
                uiOutput("uptaketext")
              )
      ),
      #8th tab
      tabItem(tabName="3d_dash",
              fluidPage(
                fluidRow(
                  uiOutput("structuretitle"),
                  br(),
                  column(2,uiOutput("pdbselector")),br(),
                  column(2,
                         uiOutput("pdbtext"),
                         uiOutput("pdbinput")))),
              column(12,uiOutput("showactionbutton"),align="left"),
              br(),
              fluidRow(
                column(align="center",NGLVieweROutput("structure"),width=6,offset=3)),
              uiOutput("hrline"),
              uiOutput("representation"),
              uiOutput("spinbox"),
              uiOutput("structureoffset"),
              column(2,uiOutput("snapshot")),
              uiOutput("pymolbutton")
              
      ),
      #9th tab
      tabItem(tabName="save_dash",
              fluidPage(
                fluidRow(
                  uiOutput("exportheader"),br(),
                  uiOutput("exporttext"),
                  column(2,uiOutput("datadownload")),
                  column(3,uiOutput("allplotsdownload"))),
                br(),
                uiOutput("tableresults")
              )
      ),
      #faq tab
      tabItem(tabName="faq_dash",
              h1(" Frequently Asked Questions",align="center"),br(),
              fluidPage(
                tags$iframe(src='./faq.html',
                            width = '100%', height = '920px',
                            frameborder = 0, scrolling = "no")
                       )
              ),
      #about tab
      tabItem(tabName="about_dash",
              fluidPage(
                tags$iframe(src='./About.html',
                            width = '100%', height = '1000px',
                            frameborder = 0, scrolling = 'no')
              )
      )
    )
  )
)
