suppressPackageStartupMessages({
  library(shiny, warn = FALSE)
  library(shinydashboard, warn = FALSE)
  library(shinyjs, warn = FALSE)
  library(NGLVieweR, warn = FALSE)
  library(shinyWidgets, warn = FALSE)
  library(shinydisconnect, warn = FALSE)
  library(shinyBS, warn = FALSE)
  library(shinydashboardPlus, warn = FALSE)
  library(shinyalert, warn = FALSE)
  library(plotly, warn = FALSE)
  library(DT, warn = FALSE)
})

# ============================================================
# Small UI helper functions
# These improve readability and avoid repeating the same patterns
# ============================================================

section_description <- function(text) {
  tags$p(class = "section-description", text)
}

embedded_html <- function(src, height = "320px", scrolling = "no") {
  tags$iframe(
    src = src,
    width = "100%",
    height = height,
    frameborder = 0,
    scrolling = scrolling
  )
}

embedded_pdf <- function(src, height = "700px") {
  tags$iframe(
    style = paste0("height:", height, "; width:100%;"),
    src = src
  )
}

download_controls_row <- function(type_ui, download_ui, extra_ui = NULL,
                                  type_offset = 3, download_offset = 2) {
  fluidRow(
    column(3, type_ui, offset = type_offset),
    column(2, download_ui, offset = download_offset),
    if (!is.null(extra_ui)) column(3, extra_ui)
  )
}

# ============================================================
# MAIN UI
# ============================================================

ui <- dashboardPage(
  
  # ============================================================
  # HEADER
  # ============================================================
  dashboardHeader(
    title = "Kingfisher HDX-MS",
    titleWidth = 235,
    dropdownMenu(
      type = "tasks",
      icon = icon("twitter"),
      headerText = tags$a(
        href = "http://twitter.com/intent/tweet?text=I%20used%20Kingfisher%20HDX%20app%20from%20@GrossLab_WashU",
        "Share us on Twitter",
        target = "_blank"
      )
    )
  ),
  
  # ============================================================
  # SIDEBAR
  # ============================================================
  dashboardSidebar(
    width = 235,
    
    sidebarMenu(
      id = "tabs",
      
      menuItem("Introduction", tabName = "intro", icon = icon("house")),
      menuItem("Input data", tabName = "input_dash", icon = icon("table")),
      menuItem("Statistical analysis", tabName = "stat_dash", icon = icon("arrow-right")),
      menuItem(
        "Fun. data analysis",
        tabName = "bayesian_dash",
        icon = icon("arrow-right"),
        badgeLabel = "NEW",
        badgeColor = "red"
      ),
      menuItem("Global woodsplot", tabName = "globwoods_dash", icon = icon("arrow-right")),
      menuItem(
        "Woodsplot by timepoint",
        tabName = "woodsbytimepoint_dash",
        icon = icon("arrow-right"),
        selectInput(
          inputId = "timepoint",
          label = "Select labeling time:",
          choices = "",
          selected = NULL
        ),
        menuSubItem(
          "Show woods plot",
          tabName = "woodsplotbytime_sub",
          icon = icon("line-chart")
        )
      ),
      
      menuItem("Digestion efficiency", tabName = "digest_dash", icon = icon("arrow-right")),
      menuItem("Peptide map", tabName = "pep_dash", icon = icon("arrow-right")),
      menuItem(
        "Uptake plots",
        tabName = "uptake_dash",
        icon = icon("arrow-right"),
        selectInput(
          inputId = "selpeptide",
          label = "Select peptide:",
          choices = "",
          selected = NULL
        ),
        menuSubItem(
          "Show uptake plots",
          tabName = "uptake_plotsub",
          icon = icon("line-chart")
        )
      ),
      menuItem("3D Structure", tabName = "3d_dash", icon = icon("cube")),
      menuItem("Export results", tabName = "save_dash", icon = icon("save")),
      menuItem("FAQ", tabName = "faq_dash", icon = icon("question")),
      menuItem("About", tabName = "about_dash", icon = icon("info")),
      
      tags$li(
        class = "sidebar-footer-item",
        div(
          class = "sidebar-footer-stack",
          div(
            class = "sidebar-clear-wrap",
            actionButton(
              "refresh",
              label = tagList(
                icon("eraser"),
                span(class = "refresh-label", "Clear data")
              ),
              class = "btn btn-default"
            )
          ),
          div(
            class = "sidebar-social-wrap",
            socialButton(
              href = "http://twitter.com/intent/tweet?text=I%20used%20Kingfisher%20HDX%20app%20from%20@GrossLab_WashU",
              icon = icon("twitter")
            ),
            socialButton(
              href = "https://github.com/juan2089/Kingfisher-HDX",
              icon = icon("github")
            )
          )
        )
      )
    )
  ),
  
  # ============================================================
  # BODY
  # ============================================================
  dashboardBody(
    shinyjs::useShinyjs(),
    
    tags$head(
      tags$link(rel = "shortcut icon", href = "favicon.ico"),
      
      tags$script(HTML("
    $(document).on('shiny:connected', function() {
      fetch('https://ipapi.co/json/')
        .then(response => response.json())
        .then(data => {
          Shiny.setInputValue('visitor_country', {
            country_name: data.country_name || 'Unknown',
            country_code: data.country_code || 'UN'
          }, {priority: 'event'});
        })
        .catch(error => {
          Shiny.setInputValue('visitor_country', {
            country_name: 'Unknown',
            country_code: 'UN'
          }, {priority: 'event'});
        });
    });
  ")),
      tags$style(HTML("
  .main-header {
    position: fixed !important;
    top: 0;
    left: 0;
    right: 0;
    width: 100%;
    z-index: 1030;
  }

  .main-header .logo {
    position: fixed !important;
    top: 0;
    left: 0;
    z-index: 1031;
    width: 235px !important;
  }

  .main-header .navbar {
    position: fixed !important;
    top: 0;
    right: 0;
    left: 235px;
    z-index: 1030;
    margin-left: 0 !important;
  }

  @media (min-width: 768px){
    .sidebar-mini.sidebar-collapse .main-header .logo {
      width: 230px !important;
    }

    .sidebar-mini.sidebar-collapse .main-header .navbar {
      left: 230px !important;
      margin-left: 0 !important;
    }
  }

  .content-wrapper,
  .right-side,
  .main-footer {
    margin-top: 50px !important;
  }

  .content-wrapper,
  .right-side {
    background-color: #f7f7f7;
  }

  .section-description {
    color: #555;
    font-size: 14px;
    margin-top: -5px;
    margin-bottom: 15px;
  }

  .small-note {
    color: #666;
    font-size: 12px;
  }

  .shiny-input-checkboxgroup {
    margin-left: 0px !important;
    padding-left: 0px !important;
  }

  .shiny-input-checkboxgroup label.control-label {
    margin-left: 0px !important;
    padding-left: 0px !important;
  }

  .tight-space {
    margin-top: 8px;
    margin-bottom: 8px;
  }

  .instruction-card {
    background: #ffffff;
    border: 1px solid #dddddd;
    border-radius: 8px;
    padding: 12px 16px;
    margin-bottom: 15px;
  }

  .instruction-card h4 {
    margin-top: 0;
    margin-bottom: 10px;
  }

  .instruction-card .small-note {
    color: #666;
    font-size: 12px;
  }

  .uptake-action-inner {
    display: flex;
    justify-content: center;
    align-items: flex-end;
    gap: 8px;
    flex-wrap: wrap;
  }

  .uptake-type-wrap {
    width: 180px;
  }

  .uptake-action-bar .form-group,
  .uptake-action-bar .shiny-input-container {
    margin-bottom: 0 !important;
  }

  .uptake-btn-wrap {
    display: flex;
    align-items: flex-end;
  }

  .uptake-btn-wrap .btn,
  .uptake-btn-wrap .shiny-download-link {
    margin-top: 0 !important;
  }

  .main-sidebar,
  .left-side {
    position: fixed !important;
    top: 50px !important;
    left: 0;
    bottom: 0;
    width: 235px;
    height: auto !important;
    overflow: hidden !important;
    padding-top: 0 !important;
  }

  .main-sidebar .sidebar,
  .left-side .sidebar {
    height: 100% !important;
    display: flex;
    flex-direction: column;
    overflow: hidden !important;
    margin: 0 !important;
    padding: 0 !important;
  }

  .sidebar-menu,
  .main-sidebar .sidebar > .sidebar-menu,
  .left-side .sidebar > .sidebar-menu {
    flex: 1 1 auto;
    min-height: 0;
    overflow-y: auto !important;
    overflow-x: hidden !important;
    margin-top: 0 !important;
    margin-bottom: 0 !important;
    padding-top: 0 !important;
    padding-bottom: 8px !important;
    list-style: none;
  }

  .sidebar-footer-stack {
    flex: 0 0 auto;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 8px;
    padding: 10px 0 10px 0;
    margin: 0;
  }

  .sidebar-clear-wrap {
    width: 100%;
    display: flex;
    justify-content: center;
    align-items: center;
  }

  .sidebar-social-wrap {
    width: 100%;
    display: flex;
    justify-content: center;
    align-items: center;
    gap: 8px;
  }

  #refresh {
    display: inline-flex !important;
    align-items: center !important;
    justify-content: center !important;
    gap: 6px;
    margin: 0 auto !important;
    width: 85%;
    max-width: 220px;
    white-space: nowrap;
    transition: all 0.25s ease;
    background-color: #f4f4f4 !important;
    border: 1px solid #d2d6de !important;
    color: #444 !important;
  }

  #refresh:hover,
  #refresh:focus,
  #refresh:active {
    background-color: #e7e7e7 !important;
    border-color: #adadad !important;
    color: #333 !important;
  }

  #refresh i {
    font-size: 14px !important;
    line-height: 1 !important;
    display: inline-block !important;
  }

  .sidebar-mini.sidebar-collapse .main-sidebar,
  .sidebar-mini.sidebar-collapse .left-side {
    width: 50px !important;
  }

  .sidebar-mini.sidebar-collapse .sidebar-menu {
    padding-bottom: 8px !important;
  }

  .sidebar-mini.sidebar-collapse .sidebar-footer-stack {
    gap: 6px;
    padding: 8px 0 10px 0;
  }

  .sidebar-mini.sidebar-collapse .sidebar-social-wrap {
    flex-direction: column;
    gap: 6px;
  }

  .sidebar-mini.sidebar-collapse #refresh {
    width: 36px !important;
    min-width: 36px !important;
    height: 36px !important;
    padding: 0 !important;
    margin: 0 auto !important;
    border-radius: 6px;
    display: inline-flex !important;
    align-items: center !important;
    justify-content: center !important;
    gap: 0 !important;
  }

  .sidebar-mini.sidebar-collapse #refresh .refresh-label {
    display: none !important;
  }

  .sidebar-mini.sidebar-collapse #refresh i {
    font-size: 14px !important;
    margin: 0 !important;
    display: inline-block !important;
    line-height: 1 !important;
  }

  @media (min-width: 768px){
    .content-wrapper,
    .right-side,
    .main-footer {
      margin-left: 235px !important;
    }

    .sidebar-collapse .content-wrapper,
    .sidebar-collapse .right-side,
    .sidebar-collapse .main-footer {
      margin-left: 50px !important;
    }
  }

  .dashboard-footer {
    background: #ffffff !important;
    border-top: 1px solid #dcdcdc !important;
    padding: 8px 16px !important;
    font-size: 13px;
  }

  .visitor-item {
    display: inline-flex;
    align-items: center;
    gap: 6px;
    font-size: 13px;
  }

  .visitor-banner {
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    align-items: center;
    justify-content: center;
  }
 
  .stat-export-row {
    display: flex;
    justify-content: center;
    align-items: flex-end;
    gap: 16px;
    flex-wrap: wrap;
  }
  
  .stat-export-cell {
    display: flex;
    align-items: flex-end;
    justify-content: center;
  }
  
  .stat-export-type {
    width: 180px;
    min-width: 180px;
  }
  
  .stat-export-cell > div {
    width: 100%;
  }
  
  .stat-export-cell .form-group,
  .stat-export-cell .shiny-input-container {
    margin-bottom: 0 !important;
  }
  
  .stat-export-cell .btn,
  .stat-export-cell .shiny-download-link {
    margin-top: 0 !important;
    white-space: nowrap;
  }
  
  .pep-export-row {
    display: flex;
    justify-content: center;
    align-items: flex-end;
    gap: 16px;
    flex-wrap: wrap;
  }
  
  .pep-export-row .form-group,
  .pep-export-row .shiny-input-container {
    margin-bottom: 0 !important;
  }
  
  .pep-export-row {
    display: flex;
    justify-content: center;
    align-items: flex-end;
    gap: 16px;
    flex-wrap: wrap;
  }
  
  .pep-export-cell {
    display: flex;
    align-items: flex-end;
    justify-content: center;
  }
  
  .pep-export-type {
    width: 180px;
    min-width: 180px;
  }
  
  .pep-export-cell > div {
    width: 100%;
  }
  
  .pep-export-cell .form-group,
  .pep-export-cell .shiny-input-container {
    margin-bottom: 0 !important;
  }
  
  .pep-export-cell .btn,
  .pep-export-cell .shiny-download-link {
    margin-top: 0 !important;
    white-space: nowrap;
  }
")
  )
    ),
    
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
    uiOutput("app_error_banner"),
    tabItems(
      
      # ============================================================
      # INTRODUCTION TAB
      # ============================================================
      tabItem(
        tabName = "intro",
        
        # ------------------------------------------------------------
        # Header / branding banner
        # ------------------------------------------------------------
        div(
          class = "instruction-card",
          fluidRow(
            column(
              width = 2,
              align = "center",
              img(
                src = "kingfisher.png",
                height = 125,
                style = "max-width:100%; height:auto; margin-top: 10px;"
              )
            ),
            
            column(
              width = 8,
              
              h2(
                tagList(
                  "Kingfisher HDX-MS",
                  tags$span(
                    "v2.0",
                    title = "New in v2.0:
• Interactive Plotly zoom for volcano and Woods plots
• Peptides can now be identified interactively in the plots by hovering the mouse over each feature
• Data validation summary panel
• Improved 3D structure tab and sequence alignment preview
• Interactive peptide table and richer information
• Enhanced export and documentation layout
• Fixed bugs and enhanced GUI",
                    style = "display:inline-block; background:#3c8dbc; color:white; padding:3px 8px; border-radius:10px; font-size:12px; margin-left:8px; vertical-align:middle; cursor:help;"
                  )
                ),
                align = "left",
                style = "margin-top: 10px; margin-bottom: 8px;"
              ),
              
              tags$p(
                "Interactive HDX-MS analysis, visualization, and structural mapping.",
                style = "font-size:16px; color:#555; margin-top:-2px; margin-bottom:6px;"
              ),
              
              tags$p(
                "Developed by JuaNolan at Washington University in St. Louis",
                style = "font-size:13px; color:#777; margin-bottom:10px;"
              ),
              
              section_description(
                "This application supports HDX-MS data import, statistical comparison, functional fitting, visualization, structural mapping, and export of final results."
              ),
              
              tags$div(
                style = "margin-top:10px; margin-bottom:12px;",
                tags$span(
                  "HDX-MS",
                  style = "display:inline-block; background:#eef3f8; padding:4px 10px; border-radius:12px; margin-right:6px; margin-bottom:4px; font-size:12px;"
                ),
                tags$span(
                  "Statistics",
                  style = "display:inline-block; background:#eef3f8; padding:4px 10px; border-radius:12px; margin-right:6px; margin-bottom:4px; font-size:12px;"
                ),
                tags$span(
                  "Kinetics",
                  style = "display:inline-block; background:#eef3f8; padding:4px 10px; border-radius:12px; margin-right:6px; margin-bottom:4px; font-size:12px;"
                ),
                tags$span(
                  "Woods plots",
                  style = "display:inline-block; background:#eef3f8; padding:4px 10px; border-radius:12px; margin-right:6px; margin-bottom:4px; font-size:12px;"
                ),
                tags$span(
                  "3D mapping",
                  style = "display:inline-block; background:#eef3f8; padding:4px 10px; border-radius:12px; margin-right:6px; margin-bottom:4px; font-size:12px;"
                ),
                tags$span(
                  "Export",
                  style = "display:inline-block; background:#eef3f8; padding:4px 10px; border-radius:12px; margin-right:6px; margin-bottom:4px; font-size:12px;"
                )
              ),
              
              div(
                style = "margin-top:10px;",
                actionButton("go_input", "Get started", icon = icon("play")),
                tags$span(style = "margin-right:8px;"),
                tags$a(
                  href = "customcsvinputfileexample.csv",
                  download = NA,
                  class = "btn btn-default",
                  icon("download"),
                  " Example CSV"
                ),
                tags$span(style = "margin-right:8px;"),
                tags$a(
                  href = "NISTmAb.fasta",
                  download = NA,
                  class = "btn btn-default",
                  icon("download"),
                  " Example FASTA"
                ),
                tags$span(style = "margin-right:8px;"),
                actionButton("go_faq", "FAQ", icon = icon("question-circle")),
                tags$span(style = "margin-right:8px;"),
                actionButton("show_citation", "How to cite", icon = icon("quote-left"))
              )
            ),
            
            column(
              width = 2,
              style = "text-align:right; padding-right:10px;",
              img(
                src = "universitylogo.png",
                height = 72,
                style = "max-width:100%; height:auto; margin-top: 18px;"
              )
            )
          )
        ),
        
        br(),
        
        # ------------------------------------------------------------
        # Intro overview
        # ------------------------------------------------------------
        div(
          class = "instruction-card",
          h3("Overview"),
          embedded_html("./introduction.html", height = "320px", scrolling = "no")
        ),
        
        br(),
        
        # ------------------------------------------------------------
        # Documentation section
        # ------------------------------------------------------------
        div(
          class = "instruction-card",
          h3("Documentation by section"),
          section_description(
            "Select a section below to view the corresponding guide."
          ),
          
          selectInput(
            "intro_doc",
            "Choose a guide:",
            choices = c(
              "1. Input data" = "inputdatainstructions.pdf",
              "2. Statistical analysis" = "statisticalanalysisinstructions.pdf",
              "3. Global woodsplot" = "globalwoodsplotinstructions.pdf",
              "4. Woodsplot by timepoint" = "woodsplotbytimepointinstructions.pdf",
              "5. Digestion efficiency" = "digestionefficiencyinstructions.pdf",
              "6. Peptide map" = "peptidemapinstructions.pdf",
              "7. Uptake plots" = "uptakeplotsinstructions.pdf",
              "8. 3D structure" = "3dstructureinstructions.pdf",
              "9. Export results" = "exportinstructions.pdf"
            ),
            selected = "inputdatainstructions.pdf"
          ),
         tags$iframe(
            src = "",
            id = "intro_pdf_viewer",
            style = "height: 700px; width: 100%; border: none;"
          )
        ),
        
        br(),
        
        # ------------------------------------------------------------
        # Contact
        # ------------------------------------------------------------
        div(
          class = "instruction-card",
          h4("Contact"),
          tags$p("For questions, suggestions, or support:"),
          tags$a(href = "mailto:juanpablo.rinconpabon@manchester.ac.uk", "Juan P Rincon"),
          br(),
          tags$a(href = "mailto:nolanmclaughlin@wustl.edu", "Nolan McLaughlin")
        )
      ),
      
      # ============================================================
      # INPUT DATA TAB
      # ============================================================
      tabItem(
        tabName = "input_dash",
        
        h2("Data import and experiment details", align = "center"),
        tags$p(
          style = "color:#555; margin-bottom:20px;",
          "Upload your HDX-MS CSV file and FASTA file, then define the experimental design used for your dataset."
        ),
        
        # ----------------------------
        # File inputs and state boxes
        # ----------------------------
        div(
          style = "background:white; padding:15px; border-radius:8px; border:1px solid #ddd; margin-bottom:20px;",
          
          fluidRow(
            column(
              6,
              div(
                id = "hdexaminerfile_block",
                style = "margin-bottom: 0px; padding-bottom: 0px;",
                fileInput(
                  "hdexaminerfile",
                  "Select exported CSV file from HDExaminer",
                  multiple = FALSE,
                  accept = c(".csv")
                )
              )
            ),
            
            column(
              3,
              style = "padding-top: 25px;",
              uiOutput("state1box")
            ),
            
            column(
              3,
              style = "padding-top: 25px;",
              uiOutput("state2box")
            )
          ),
          
          div(
            id = "datacheck_block",
            style = "margin-top: -8px; margin-left: 0px; padding-left: 0px;",
            checkboxGroupInput(
              "datacheck",
              label = tags$span(
                "Check an option if you are NOT using HDExaminer output file",
                tags$span(style = "margin-left: 4px;"),
                bsButton(
                  "testingbutton",
                  label = "",
                  style = "info",
                  icon = icon("info"),
                  size = "extra-small"
                )
              ),
              c("WATERS/DynamX output file", "Custom .csv output file")
            )
          ),
          
          bsPopover(
            id = "testingbutton",
            title = HTML("<b>Important</b>"),
            content = paste0(
              "If you are not using an HDExaminer output file, please select the option that fits your data. ",
              "Check WATERS/DynamX if you are using a DynamX output file, or check Custom .csv for any other type of file. ",
              "Note: for custom .csv, you will need to manually format your data as in the example file available for download below."
            ),
            placement = "right",
            trigger = "hover",
            options = list(container = "body")
          ),
          
          fileInput(
            "fastafile",
            "Select FASTA file",
            multiple = FALSE,
            accept = c(".FASTA", ".fasta", ".fa")
          )
        ),
        
        # ----------------------------
        # Experiment settings
        # ----------------------------
        div(
          style = "background:white; padding:15px; border-radius:8px; border:1px solid #ddd; margin-bottom:20px;",
          
          fluidRow(
            column(
              4,
              numericInput(
                "states",
                label = HTML("Number of protein states"),
                value = 2,
                min = 2
              )
            ),
            
            column(
              4,
              numericInput(
                "timepoints",
                label = HTML("Number of labeling times"),
                value = 4
              )
            ),
            
            column(
              4,
              numericInput(
                "significancelevel",
                label = HTML("Significance level"),
                value = 0.01,
                min = 0,
                max = 0.5,
                step = 0.01
              )
            )
          ),
          
          fluidRow(
            column(
              6,
              numericInput(
                "replicates",
                "Number of replicates per labeling time",
                value = 4,
                min = 2
              )
            ),
            column(
              6,
              radioButtons(
                "cantdeutbutton",
                "Number of residues that can't withhold deuteration:",
                choices = c(1, 2),
                selected = 2,
                inline = TRUE
              )
            )
          ),
          
          textInput(
            "labelingtimepoints",
            "Labeling times in seconds (separated by commas):",
            value = "30,60,800,2500"
          )
        ),
        
        fluidRow(
          column(12, uiOutput("analyze"))
        ),
        
        br(),
        
        fluidRow(
          column(12, uiOutput("validation_summary"))
        ),
        
        hr(style = "border-top: 1px solid #000000;"),
        
        fluidRow(
          column(
            12,
            h5("Download a test .csv file", tags$a(href = "customcsvinputfileexample.csv", "here"), align = "right"),
            h5("Download a test FASTA file", tags$a(href = "NISTmAb.fasta", "here"), align = "right")
          )
        )
      ),
      
      # ============================================================
      # STATISTICAL ANALYSIS TAB
      # ============================================================
      tabItem(
        tabName = "stat_dash",
        
        fluidPage(
          
          # ------------------------------------------------------------
          # Header
          # ------------------------------------------------------------
          div(
            class = "instruction-card",
            h2("Statistical analysis", align = "center", style = "margin-top:0;"),
            tags$p(
              "This panel displays hybrid significance volcano plots and, when clustering is enabled, histograms of normalized differences.",
              style = "color:#555; text-align:center; margin-bottom:0;"
            )
          ),
          
          # ------------------------------------------------------------
          # Summary strip
          # ------------------------------------------------------------
          div(
            class = "instruction-card",
            tags$p(
              style = "color:#555; margin-bottom:10px; text-align:center;",
              "Quick overview of peptide-level statistical outcomes for the selected comparison."
            ),
            uiOutput("stat_summary_strip")
          ),
          
          # ------------------------------------------------------------
          # Volcano plot
          # ------------------------------------------------------------
          div(
            class = "instruction-card",
            h4("Volcano plot", style = "text-align:center; margin-top:0;"),
            tags$p(
              class = "small-note",
              style = "text-align:center;",
              "Use zoom and hover to inspect peptide-level differences interactively."
            ),
            fluidRow(
              column(
                width = 10,
                offset = 1,
                div(
                  style = "display:flex; justify-content:center;",
                  plotlyOutput(outputId = "volcanoplot", height = "600px", width = "100%")
                )
              )
            )
          ),
          
          # ------------------------------------------------------------
          # Clustering controls
          # ------------------------------------------------------------
          div(
            class = "instruction-card",
            h4("Optional clustering", style = "text-align:center; margin-top:0;"),
            
            fluidRow(
              column(
                width = 12,
                align = "center",
                conditionalPanel(
                  condition = "output.volcanoplot",
                  checkboxInput("clusteringsel", "Check to cluster the data", value = FALSE)
                )
              )
            ),
            
            shinyjs::hidden(
              div(
                id = "panel",
                div(
                  style = "display:flex; justify-content:center;",
                  div(
                    style = "width:100%; max-width:900px;",
                    fluidRow(
                      column(
                        6,
                        div(
                          style = "display:flex; justify-content:center;",
                          div(
                            style = "width:85%; max-width:320px;",
                            selectInput(
                              "colorstr",
                              label = "Select color for strong effects:",
                              choices = c("forestgreen", "darkmagenta", "deeppink", "red"),
                              selected = "forestgreen"
                            )
                          )
                        )
                      ),
                      column(
                        6,
                        div(
                          style = "display:flex; justify-content:center;",
                          div(
                            style = "width:85%; max-width:320px;",
                            selectInput(
                              "colorint",
                              label = "Select color for intermediate effects:",
                              choices = c("orange", "cyan3", "aquamarine", "palegreen2"),
                              selected = "orange"
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          ),
          
          # ------------------------------------------------------------
          # Downloads
          # ------------------------------------------------------------
          div(
            class = "instruction-card",
            h4("Export options", style = "text-align:center; margin-top:0;"),
            
            div(
              class = "stat-export-row",
              
              div(
                class = "stat-export-cell stat-export-type",
                uiOutput("volcanotype")
              ),
              
              div(
                class = "stat-export-cell",
                uiOutput("volcanodownload")
              ),
              
              conditionalPanel(
                condition = "input.clusteringsel == true",
                div(
                  class = "stat-export-cell",
                  uiOutput("histdownload")
                )
              )
            )
          ),
          
          div(
            class = "instruction-card",
            uiOutput("volcanotext")
          ),
          
          # ------------------------------------------------------------
          # Histogram (ONLY shown when clustering is checked)
          # ------------------------------------------------------------
          conditionalPanel(
            condition = "input.clusteringsel == true",
            div(
              class = "instruction-card",
              h4("Histogram of normalized differences", style = "text-align:center; margin-top:0;"),
              fluidRow(
                column(
                  width = 12,
                  align = "center",
                  plotOutput(outputId = "histogramclustplot", width = "60%")
                )
              ),
              br(),
              fluidRow(
                column(
                  width = 12,
                  align = "center",
                  uiOutput("binsize")
                )
              )
            )
          )
        )
      ),
      
      # ============================================================
      # FUNCTIONAL DATA ANALYSIS TAB
      # ============================================================
      tabItem(
        tabName = "bayesian_dash",
        fluidPage(
          shinyjs::hidden(
            div(
              id = "functional",
              
              h1("Functional data analysis", align = "center"),
              br(),
              
              section_description(
                "Fit kinetic models to peptide uptake curves, visualize fitted profiles, and compare states using forest and Manhattan plots."
              ),
              
              embedded_html("./functionalanalysis.html", height = "100px", scrolling = "no"),
              
              column(
                12,
                radioButtons(
                  inputId = "selectedeq",
                  label = h3("Select the desired equation to fit the data:"),
                  inline = TRUE,
                  choiceNames = list(
                    img(src = "equationwithd.png", width = 400, height = 45),
                    img(src = "equationwithoutd.png", width = 400, height = 45),
                    img(src = "equationwithdwithoutP.png", width = 400, height = 45),
                    img(src = "equationwithoutdorP.png", width = 400, height = 45)
                  ),
                  choiceValues = list(
                    "equationwithd",
                    "equationwithoutd",
                    "equationwithdwithoutP",
                    "equationwithoutdorP"
                  )
                ),
                align = "center"
              ),
              
              h4("Input initial parameters:"),
              section_description(
                "Provide starting values for non-linear fitting. Good starting values can improve convergence and reduce failed fits."
              ),
              
              fluidRow(
                column(
                  3,
                  numericInput("a_input", label = HTML("Parameter a:"), value = 0, min = 0, max = 100)
                ),
                column(
                  3,
                  numericInput("b_input", label = HTML("Parameter b:"), value = 0, min = 0, max = 100)
                ),
                column(
                  3,
                  numericInput("d_input", label = HTML("Parameter d:"), value = 0, min = 0, max = 100)
                ),
                column(
                  3,
                  numericInput("p_input", label = HTML("Parameter p:"), value = 1, min = 0, max = 1)
                )
              ),
              
              br(),
              actionButton("startfitting", "Calculate")
            )
          ),
          
          br(),
          
          shinyjs::hidden(
            div(
              id = "fdabox",
              tabBox(
                title = "",
                width = "100%",
                id = "bayesianbox",
                selected = "Heat map",
                
                # ----------------------------
                # Heat map tab
                # ----------------------------
                tabPanel(
                  "Heat map",
                  div(
                    class = "instruction-card",
                    "Displays the full qFeatures-derived exchange matrix across all samples and timepoints."
                  ),
                  fluidRow(column(10, plotOutput(outputId = "heatmap"), align = "center", offset = 1)),
                  fluidRow(
                    column(
                      3,
                      selectInput(
                        "heattype",
                        label = NULL,
                        choices = c(".pdf", ".svg", ".eps"),
                        selected = ".pdf"
                      ),
                      offset = 3
                    ),
                    column(3, downloadButton("heatbutton", "Download Heat map"))
                  )
                ),
                
                # ----------------------------
                # Fitted kinetic plots tab
                # ----------------------------
                tabPanel(
                  "Fitted kinetic plots",
                  div(
                    class = "instruction-card",
                    "Inspect fitted uptake curves and forest plots for individual peptides. Use the selector below to move between peptides."
                  ),
                  shinyjs::hidden(actionButton("reset", "Reset")),
                  fluidRow(
                    column(6, plotOutput(outputId = "kineticplot")),
                    column(6, plotOutput(outputId = "forestplot"))
                  ),
                  br(),
                  column(12, tableOutput("summary"), align = "center"),
                  column(
                    8,
                    align = "center",
                    selectInput(inputId = "selpeptide2", label = "Select peptide:", choices = ""),
                    offset = 2
                  ),
                  br(),
                  br(),
                  fluidRow(
                    column(
                      3,
                      selectInput(
                        "kinetictype",
                        label = NULL,
                        choices = c(".pdf", ".svg", ".eps"),
                        selected = ".pdf"
                      ),
                      offset = 2
                    ),
                    column(3, downloadButton("kineticbutton", "Download kinetic plot")),
                    column(3, actionButton("downloadallkinetic", "Download all kinetic plots", icon = icon("download")), align = "center")
                  )
                ),
                
                # ----------------------------
                # Manhattan plot tab
                # ----------------------------
                tabPanel(
                  "Manhattan plot",
                  div(
                    class = "instruction-card",
                    "Compares peptide-level statistical evidence across the selected labeling timepoint."
                  ),
                  br(),
                  column(
                    3,
                    selectInput(inputId = "seltimemanhattan", label = "Select labeling time:", choices = ""),
                    offset = 4,
                    align = "center"
                  ),
                  br(),
                  fluidRow(
                    column(12, plotOutput(outputId = "manhattanplot"), align = "center", offset = 0)
                  ),
                  br(),
                  br(),
                  fluidRow(
                    column(
                      4,
                      selectInput(
                        "manhattantype",
                        label = NULL,
                        choices = c(".pdf", ".svg", ".eps"),
                        selected = ".pdf"
                      ),
                      offset = 3
                    ),
                    column(3, downloadButton("manhattanbutton", "Download Manhattan plot"))
                  )
                )
              )
            )
          )
        )
      ),
      
      # ============================================================
      # GLOBAL WOODS PLOT TAB
      # ============================================================
      tabItem(
        tabName = "globwoods_dash",
        fluidPage(
          h2("Global Woods plot", align = "center"),
          section_description(
            "Shows residue coverage and significance across all peptides in a single global view."
          ),
          
          fluidRow(
            column(
              12,
              plotly::plotlyOutput("globwoodsplot", width = "85%", height = "600px"),
              align = "center"
            )
          ),
          
          fluidRow(
            br(),
            shinyjs::hidden(
              div(
                id = "panel2",
                fluidRow(
                  column(
                    6,
                    selectInput(
                      "colorflex",
                      label = "Select color for peptides with increased flexibility:",
                      choices = c("firebrick2", "cyan", "deeppink"),
                      selected = "firebrick2"
                    ),
                    align = "center"
                  ),
                  column(
                    6,
                    selectInput(
                      "colorprot",
                      label = "Select color for peptides with increased protection:",
                      choices = c("blue", "darkcyan", "chocolate"),
                      selected = "blue"
                    ),
                    align = "center"
                  )
                )
              )
            )
          ),
          
          br(),
          column(3, uiOutput("globwoodstype"), offset = 3),
          column(2, uiOutput("globwoodsdownload"), offset = 2),
          br(),
          uiOutput("woodsplottext")
        )
      ),
      
      # ============================================================
      # WOODS PLOT BY TIMEPOINT TAB
      # ============================================================
      tabItem(
        tabName = "woodsplotbytime_sub",
        fluidPage(
          h2("Woods plot by timepoint", align = "center"),
          section_description(
            "Displays peptide-level ΔHX values along the sequence for one labeling time at a time."
          ),
          
          fluidRow(
            column(12, plotlyOutput(outputId = "woodsbytimepointplot", height = "600px"), align = "center")
          ),
          
          br(),
          column(3, uiOutput("woodsbytimepointtype"), offset = 3),
          column(2, uiOutput("woodsbytimepointdownload"), offset = 2),
          br(),
          uiOutput("woodsplotbytimetext")
        )
      ),
      
      # ============================================================
      # DIGESTION EFFICIENCY TAB
      # ============================================================
      tabItem(
        tabName = "digest_dash",
        fluidPage(
          h2("Digestion efficiency", align = "center"),
          section_description(
            "Summarizes peptide coverage, average peptide length, and redundancy across the sequence."
          ),
          
          uiOutput("digestion_summary_panel"),
          
          br(),
          
          div(
            class = "instruction-card",
            h4("Peptide length distribution", style = "text-align:center; margin-top:0;"),
            tags$p(
              class = "small-note",
              style = "text-align:center;",
              "Distribution of peptide lengths observed in the current digestion dataset."
            ),
            fluidRow(
              column(
                12,
                align = "center",
                plotOutput(outputId = "digestionplot", width = "70%", height = "500px")
              )
            )
          ),
          
          br(),
          
          fluidRow(
            column(3, uiOutput("digestiontype"), offset = 3),
            column(2, uiOutput("digestiondownload"), offset = 2)
          ),
          
          br(),
          uiOutput("digestiontext")
        )
      ),
      
      # ============================================================
      # PEPTIDE MAP TAB
      # ============================================================
      tabItem(
        tabName = "pep_dash",
        
        fluidPage(
          
          div(
            class = "instruction-card",
            h2("Peptide map", align = "center", style = "margin-top:0;"),
            tags$p(
              "Visualizes peptide coverage across the protein sequence, with both a compact overview and the full stacked peptide layout.",
              style = "color:#555; text-align:center; margin-bottom:0;"
            )
          ),
          
          br(),
          
          # ------------------------------------------------------------
          # Summary peptide map
          # ------------------------------------------------------------
          shinyjs::hidden(
            div(
              id = "pep2",
              div(
                class = "instruction-card",
                h4("Summary peptide map", style = "text-align:center; margin-top:0;"),
                tags$p(
                  class = "small-note",
                  style = "text-align:center;",
                  "Compact overview of total sequence coverage and peptide stacking."
                ),
                
                fluidRow(
                  column(
                    width = 12,
                    align = "center",
                    plotOutput(
                      outputId = "summaryplot",
                      width = "80%",
                      height = "220px"
                    )
                  )
                ),
                
                br(),
                
                div(
                  style = "display:flex; justify-content:center;",
                  downloadButton("summarypepbutton", "Download summary plot")
                )
              )
            )
          ),
          
          br(),
          
          # ------------------------------------------------------------
          # Detailed peptide map
          # ------------------------------------------------------------
          div(
            class = "instruction-card",
            h4("Detailed peptide map", style = "text-align:center; margin-top:0;"),
            tags$p(
              class = "small-note",
              style = "text-align:center;",
              "Detailed stacked peptide map across the protein sequence."
            ),
            
            fluidRow(
              column(
                width = 12,
                align = "center",
                plotOutput(
                  outputId = "peptideplot",
                  width = "85%",
                  height = "620px"
                )
              )
            )
          ),
          
          br(),
          
          # ------------------------------------------------------------
          # Downloads button
          # ------------------------------------------------------------
          div(
            class = "instruction-card",
            h4("Display and export options", style = "text-align:center; margin-top:0;"),
            
            shinyjs::hidden(
              div(
                id = "pep",
                fluidRow(
                  column(
                    width = 12,
                    div(
                      style = "display:flex; justify-content:center;",
                      div(
                        style = "width:260px;",
                        selectInput(
                          "colorpep",
                          label = "Select color for peptides:",
                          choices = c("blue", "darkcyan", "chocolate"),
                          selected = "blue"
                        )
                      )
                    )
                  )
                )
              )
            ),
            
            div(
              class = "pep-export-row",
              
              div(
                class = "pep-export-cell pep-export-type",
                uiOutput("peptidetype")
              ),
              
              div(
                class = "pep-export-cell",
                uiOutput("peptidedownload")
              )
            )
          ),
          
          br(),
          
          div(
            class = "instruction-card",
            uiOutput("peptidemaptext")
          )
        )
      ),
      
      # ============================================================
      # UPTAKE PLOTS TAB
      # ============================================================
      tabItem(
        tabName = "uptake_plotsub",
        fluidPage(
          h2("Uptake plots", align = "center"),
          section_description(
            "Display uptake curves for the selected peptide, compare protein states, and customize the plot appearance."
          ),
          
          uiOutput("uptake_peptide_info"),
          
          br(),
          
          fluidRow(
            # ----------------------------------------------------------
            # Left settings panel
            # ----------------------------------------------------------
            column(
              width = 3,
              
              div(
                class = "instruction-card",
                h4("Plot settings"),
                shinyjs::hidden(
                  div(
                    id = "plotssettings",
                    
                    tags$label("Color palette"),
                    radioButtons(
                      inputId = "coloring",
                      label = NULL,
                      choices = list(
                        "Solarized" = "solarized",
                        "Chalk" = "chalk",
                        "Dust" = "dust"
                      )
                    ),
                    
                    sliderInput(
                      inputId = "thickness",
                      label = "Line thickness",
                      min = 0.5,
                      max = 5,
                      value = 1,
                      step = 0.1
                    ),
                    
                    numericInput(
                      "fontsize",
                      label = "Font size",
                      value = 10,
                      min = 8,
                      max = 25
                    )
                  )
                )
              ),
              
              div(
                class = "instruction-card",
                h4("Axis settings"),
                shinyjs::hidden(
                  div(
                    id = "limits",
                    fluidRow(
                      column(
                        6,
                        numericInput(
                          "yminvalue",
                          label = "Y min",
                          value = "",
                          min = 0,
                          max = 25
                        )
                      ),
                      column(
                        6,
                        numericInput(
                          "ymaxvalue",
                          label = "Y max",
                          value = "",
                          min = 0,
                          max = 25
                        )
                      )
                    )
                  )
                )
              ),
              
              div(
                class = "instruction-card",
                h4("Additional display options"),
                tags$p(
                  class = "small-note",
                  "Extra display controls will appear here when available."
                ),
                shinyjs::hidden(
                  div(
                    id = "allplotscheck",
                    checkboxInput(
                      inputId = "showallstates",
                      "Show all protein states",
                      value = FALSE
                    )
                  )
                )
              )
            ),
            
            # ----------------------------------------------------------
            # Right plot panel
            # ----------------------------------------------------------
            column(
              width = 9,
              
              div(
                class = "instruction-card",
                h4("Selected peptide uptake profile"),
                tags$p(
                  class = "small-note",
                  "Use the peptide selector in the sidebar to switch peptides. The plot updates automatically."
                ),
                
                div(
                  style = "display:flex; justify-content:center;",
                  plotOutput(
                    outputId = "uptakeplot",
                    width = "100%",
                    height = "520px"
                  )
                )
              )
            )
          ),
          
          br(),
          
          # ------------------------------------------------------------
          # Bottom action bar
          # ------------------------------------------------------------
          div(
            class = "instruction-card uptake-action-bar",
            style = "width: fit-content; margin: 0 auto; padding: 12px 18px;",
            
            div(
              class = "uptake-action-inner",
              
              div(
                class = "uptake-type-wrap",
                uiOutput("uptaketype")
              ),
              
              div(
                class = "uptake-btn-wrap",
                uiOutput("singledownload")
              ),
              
              div(
                class = "uptake-btn-wrap",
                uiOutput("uptakedownload")
              )
            )
          ),
          
          br(),
          
          uiOutput("uptaketext")
        )
      ),
      
      # ============================================================
      # 3D STRUCTURE TAB
      # ============================================================
      tabItem(
        tabName = "3d_dash",
        
        h2("3D Structure", align = "center"),
        tags$p(
          style = "color:#555; margin-bottom:20px;",
          "Visualize significant peptides on a protein structure, compare FASTA and PDB sequence alignment, and export screenshots or PyMOL scripts."
        ),
        
        fluidRow(
          # ---------------- LEFT PANEL: CONTROLS ----------------
          column(
            4,
            
            div(
              class = "instruction-card",
              h4("1. Load structure"),
              uiOutput("pdbselector"),
              uiOutput("pdbtext"),
              uiOutput("pdbinput"),
              br(),
              uiOutput("showactionbutton")
            ),
            
            div(
              class = "instruction-card",
              h4("2. Visualization settings"),
              uiOutput("representation"),
              fluidRow(
                column(6, uiOutput("spinbox")),
                column(
                  6,
                  shinyjs::hidden(
                    div(
                      id = "3dtoggles",
                      materialSwitch(
                        "clustercoloring",
                        label = "Cluster colors",
                        value = FALSE,
                        status = "success"
                      )
                    )
                  )
                )
              ),
              br(),
              uiOutput("structureoffset")
            ),
            
            div(
              class = "instruction-card",
              h4("3. Export"),
              tags$p(
                style = "font-size: 12px; color: #666;",
                "Save a screenshot of the current view or export a PyMOL coloring script."
              ),
              fluidRow(
                column(6, uiOutput("snapshot")),
                column(6, uiOutput("pymolbutton"))
              )
            ),
            
            div(
              class = "instruction-card",
              h4("Color legend"),
              tags$div(
                style = "font-size: 13px;",
                tags$p(HTML("<b>Default significance mode</b>")),
                tags$div(
                  style = "display:flex; align-items:center; gap:8px; margin-bottom:4px;",
                  tags$div(style = "width:14px; height:14px; background:gray; border:1px solid #333;"),
                  "No significant change"
                ),
                tags$div(
                  style = "display:flex; align-items:center; gap:8px; margin-bottom:4px;",
                  tags$div(style = "width:14px; height:14px; background:yellow; border:1px solid #333;"),
                  "Mixed behavior"
                ),
                tags$div(
                  style = "display:flex; align-items:center; gap:8px; margin-bottom:4px;",
                  tags$div(style = "width:14px; height:14px; background:#EE2C2C; border:1px solid #333;"),
                  "Increased flexibility"
                ),
                tags$div(
                  style = "display:flex; align-items:center; gap:8px;",
                  tags$div(style = "width:14px; height:14px; background:#0000FF; border:1px solid #333;"),
                  "Increased protection"
                ),
                br(),
                tags$p(
                  style = "font-size: 12px; color: #666;",
                  "If clustering mode is enabled, colors indicate effect strength rather than protection/deprotection."
                )
              )
            )
          ),
          
          # ---------------- RIGHT PANEL: VIEWER + SEQUENCES ----------------
          column(
            8,
            
            div(
              class = "instruction-card",
              uiOutput("structuretitle"),
              tags$p(
                style = "font-size: 12px; color: #666; margin-top: -5px;",
                "Tip: rotate the structure with the mouse, zoom with the wheel, and use the offset to align residue numbering."
              ),
              NGLVieweROutput("structure", height = "680px")
            ),
            
            uiOutput("hrline"),
            
            div(
              class = "instruction-card",
              h4("Sequence alignment preview"),
              tags$p(
                style = "font-size: 12px; color: #666;",
                "Compare the FASTA sequence with the shifted PDB sequence. The amino acid offset will move the PDB sequence accordingly."
              ),
              uiOutput("sequenceoffset"),
              uiOutput("pdbsequence")
            )
          )
        )
      ),
      
      # ============================================================
      # EXPORT TAB
      # ============================================================
      tabItem(
        tabName = "save_dash",
        fluidPage(
          h2("Export results", align = "center"),
          section_description(
            "Download the full analysis workbook and export the generated plots for reporting or publication."
          ),
          
          fluidRow(
            br(),
            column(2, uiOutput("datadownload")),
            column(3, uiOutput("allplotsdownload"))
          ),
          
          br(),
          uiOutput("tableresults")
        )
      ),
      
      # ============================================================
      # FAQ TAB
      # ============================================================
      tabItem(
        tabName = "faq_dash",
        h1("Frequently Asked Questions", align = "center"),
        section_description(
          "Common questions about Kingfisher"
        ),
        br(),
        fluidPage(
          embedded_html("./faq.html", height = "920px", scrolling = "no")
        )
      ),
      # ============================================================
      # ABOUT TAB
      # ============================================================
      tabItem(
        tabName = "about_dash",
        fluidPage(
          h1("About", align = "center"),
          section_description(
            "Project background, contributors, and supporting information."
          ),
          embedded_html("./About.html", height = "1450px", scrolling = "no")
        )
      )
    )
  ),
  footer = dashboardFooter(
    left = uiOutput("visitor_banner"),
    right = "Kingfisher HDX-MS analytics"
  )
)