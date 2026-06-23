
visitor_dir <- "logs"
if (!dir.exists(visitor_dir)) dir.create(visitor_dir, recursive = TRUE)

visitor_file <- file.path(visitor_dir, "visitor_counts.rds")

# ============================================================
# 1. Libraries and global options
# ============================================================

# ============================================================
# STEP 1 — Load libraries 
# ============================================================
suppressPackageStartupMessages({
  pkgs <- c(
    "readxl", "openxlsx", "writexl", "rlist", "seqinr", "gridExtra", "grid",
    "devtools", "shiny", "shinydashboard", "shinyjs", "NGLVieweR",
    "shinyWidgets", "shinydisconnect", "shinyBS", "hdxstats", "RColorBrewer",
    "pheatmap", "scales", "viridis", "patchwork", "Biostrings", "tidyverse",
    "qpdf", "ggpubr", "plyr"
  )
  invisible(lapply(pkgs, library, character.only = TRUE, warn.conflicts = FALSE))
})

options(shiny.maxRequestSize = 30 * 1024^2)

write_visitors <- function(df) {
  saveRDS(df, visitor_file)
}

# ============================================================
# 2. Helper functions
# ============================================================

# ------------------------------------------------------------
# 2.1 Build qFeatures object from long/wide HDX data
# ------------------------------------------------------------
.build_qf <- function(data, replicates, timepoints, states) {
  TF2 <- data %>%
    subset(select = -c(`pepnumber`)) %>%
    pivot_longer(
      cols = contains("_"),
      names_to = c("Time", "Replicate"),
      values_transform = list(Replicate = as.double, Time = as.double),
      names_pattern = "X?(.*)s_(.*)",
      values_to = "D"
    ) %>%
    dplyr::rename(
      "pep_sequence" = "Sequence",
      "pep_charge" = "Charge",
      "d" = "D",
      "hx_time" = "Time",
      "replicate_cnt" = "Replicate",
      "hx_sample" = "State"
    ) %>%
    mutate(
      hx_time = as.numeric(hx_time),
      replicate_cnt = as.numeric(replicate_cnt)
    ) %>%
    dplyr::group_by(pep_sequence, pep_charge, hx_time, replicate_cnt, hx_sample) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    filter(n == 1)
  
  HDX_wide <- pivot_wider(
    data.frame(TF2),
    values_from = d,
    names_from = c("hx_time", "replicate_cnt", "hx_sample"),
    id_cols = c("pep_sequence", "pep_charge")
  )
  
  new.colnames <- gsub("0_", "0rep", paste0("X", colnames(HDX_wide)[-c(1, 2)]))
  new.colnames <- gsub("_", "cond", new.colnames)
  new.colnames <- gsub(" ", "", new.colnames)
  new.colnames <- gsub(".", "", new.colnames, fixed = TRUE)
  
  HDXqf <- parseDeutData(
    object = DataFrame(HDX_wide),
    design = new.colnames,
    quantcol = 3:((replicates * timepoints * states) + 2)
  )
  
  return(HDXqf)
}

# ------------------------------------------------------------
# 2.2 Add mean and SD columns per labeling time
# ------------------------------------------------------------
.add_mean_sd <- function(HDXdatastates, timepoints, replicates, labelingtimepoints, start_state, end_state) {
  for (z in start_state:end_state) {
    HDXdatastates[[z]][["pepnumber"]] <- seq_len(nrow(HDXdatastates[[z]]))
    
    for (i in seq_len(timepoints)) {
      cols <- (5 + ((i - 1) * replicates) + 1):(5 + (i * replicates))
      vals <- as.matrix(HDXdatastates[[z]][, cols, drop = FALSE])
      
      mean1 <- rowMeans(vals, na.rm = TRUE)
      sd1 <- apply(vals, 1, sd, na.rm = TRUE)
      sd1[is.na(sd1)] <- 0
      
      HDXdatastates[[z]] <- cbind(HDXdatastates[[z]], mean1)
      colnames(HDXdatastates[[z]])[ncol(HDXdatastates[[z]])] <- paste("Average D at ", labelingtimepoints[i], "sec")
      
      HDXdatastates[[z]] <- cbind(HDXdatastates[[z]], sd1)
      colnames(HDXdatastates[[z]])[ncol(HDXdatastates[[z]])] <- paste("SD for ", labelingtimepoints[i], "sec")
    }
  }
  
  return(HDXdatastates)
}

# ------------------------------------------------------------
# 2.3 Calculate pooled standard deviation
# ------------------------------------------------------------
.calc_pooled_sd <- function(HDXdatastates, timepoints, replicates) {
  npooledSD <- 0
  dpooledSD <- 0
  
  for (z in 1:2) {
    for (i in seq_len(timepoints)) {
      cols <- (5 + ((i - 1) * replicates) + 1):(5 + (i * replicates))
      vals <- as.matrix(HDXdatastates[[z]][, cols, drop = FALSE])
      
      sds <- apply(vals, 1, sd, na.rm = TRUE)
      ns <- apply(vals, 1, function(x) sum(!is.na(x)))
      
      sds[is.na(sds)] <- 0
      valid <- ns > 1
      
      npooledSD <- npooledSD + sum((sds[valid]^2) * (ns[valid] - 1))
      dpooledSD <- dpooledSD + sum(ns[valid] - 1)
    }
  }
  
  return(sqrt(npooledSD / dpooledSD))
}

# ------------------------------------------------------------
# 2.4 Select kinetic fitting formula based on UI choice
# ------------------------------------------------------------
.fit_model <- function(HDXqf, all_peptides, starting_parameters, input) {
  if (input$selectedeq == "equationwithd") {
    hdxstats::analyse_kinetics(
      data = HDXqf,
      method = "fit",
      formula = value ~ a * (1 - exp(-b * (timepoint)^p)) + d,
      peptide_selection = all_peptides,
      start = starting_parameters,
      maxAttempts = 20
    )
  } else if (input$selectedeq == "equationwithoutd") {
    hdxstats::analyse_kinetics(
      data = HDXqf,
      method = "fit",
      formula = value ~ a * (1 - exp(-b * (timepoint)^p)),
      peptide_selection = all_peptides,
      start = starting_parameters,
      maxAttempts = 20
    )
  } else if (input$selectedeq == "equationwithdwithoutP") {
    hdxstats::analyse_kinetics(
      data = HDXqf,
      method = "fit",
      formula = value ~ a * (1 - exp(-b * (timepoint))) + d,
      peptide_selection = all_peptides,
      start = starting_parameters,
      maxAttempts = 20
    )
  } else {
    hdxstats::analyse_kinetics(
      data = HDXqf,
      method = "fit",
      formula = value ~ a * (1 - exp(-b * (timepoint))),
      peptide_selection = all_peptides,
      start = starting_parameters,
      maxAttempts = 20
    )
  }
}

# ------------------------------------------------------------
# 2.5 Compute delta D, p-values, and normalized values
# ------------------------------------------------------------
.compute_stats <- function(HDXdatastates, timepoints, replicates, cantdeut, labelingtimepoints, significancelevel) {
  listtestresults <- data.frame(
    HDXdatastates[[1]][2],
    HDXdatastates[[1]][3],
    HDXdatastates[[1]][4],
    HDXdatastates[[1]][5]
  )
  
  clusteringresults <- listtestresults
  
  for (m in seq_len(timepoints)) {
    cols <- (5 + ((m - 1) * replicates) + 1):(5 + (m * replicates))
    
    difference <- rep(NA_real_, nrow(HDXdatastates[[1]]))
    ttest <- rep(NA_real_, nrow(HDXdatastates[[1]]))
    normvalues <- rep(NA_real_, nrow(HDXdatastates[[1]]))
    
    for (n in seq_len(nrow(HDXdatastates[[1]]))) {
      dataset1 <- as.numeric(HDXdatastates[[1]][n, cols])
      dataset2 <- as.numeric(HDXdatastates[[2]][n, cols])
      
      if (sum(!is.na(dataset1)) < 2 || sum(!is.na(dataset2)) < 2) {
        ttest[n] <- NA
        difference[n] <- NA
        normvalues[n] <- NA
      } else {
        ttest[n] <- t.test(
          na.omit(dataset2),
          na.omit(dataset1),
          var.equal = FALSE,
          conf.level = significancelevel / 2,
          alternative = "two.sided"
        )$p.value
        
        difference[n] <- mean(dataset2, na.rm = TRUE) - mean(dataset1, na.rm = TRUE)
        
        exchangable <- str_length(HDXdatastates[[1]][n, 5]) -
          str_count(HDXdatastates[[1]][n, 5], "P") - cantdeut
        
        normvalues[n] <- difference[n] / exchangable
      }
    }
    
    listtestresults <- cbind(listtestresults, difference)
    colnames(listtestresults)[ncol(listtestresults)] <- paste0("Delta D at ", labelingtimepoints[m], "sec")
    
    listtestresults <- cbind(listtestresults, ttest)
    colnames(listtestresults)[ncol(listtestresults)] <- paste0("t test ", labelingtimepoints[m], "sec")
    
    clusteringresults <- cbind(clusteringresults, difference)
    colnames(clusteringresults)[ncol(clusteringresults)] <- paste0("Delta D at ", labelingtimepoints[m], "sec")
    
    clusteringresults <- cbind(clusteringresults, normvalues)
    colnames(clusteringresults)[ncol(clusteringresults)] <- paste0("Cluster", labelingtimepoints[m], "sec")
    
    clusteringresults <- cbind(clusteringresults, ttest)
    colnames(clusteringresults)[ncol(clusteringresults)] <- paste0("t test", labelingtimepoints[m], "sec")
  }
  
  return(list(listtestresults = listtestresults, clusteringresults = clusteringresults))
}

# ------------------------------------------------------------
# 2.6 Classify peptide significance across all timepoints
# ------------------------------------------------------------
.classify_significance <- function(listtestresults, timepoints, threshold, significancelevel) {
  significant <- rep(0, nrow(listtestresults))
  
  for (i in seq_len(nrow(listtestresults))) {
    var2 <- rep(0, timepoints)
    
    for (j in seq_len(timepoints)) {
      delta_col <- (2 * j) + 3
      p_col <- (2 * j) + 4
      
      if (is.na(listtestresults[i, delta_col]) || is.na(listtestresults[i, p_col])) {
        var2[j] <- 0
      } else if (listtestresults[i, delta_col] >= threshold && listtestresults[i, p_col] < significancelevel) {
        var2[j] <- 1
      } else if (listtestresults[i, delta_col] <= -threshold && listtestresults[i, p_col] < significancelevel) {
        var2[j] <- -1
      } else {
        var2[j] <- 0
      }
    }
    
    if (any(var2 > 0) && !any(var2 < 0)) {
      significant[i] <- 1
    } else if (any(var2 < 0) && !any(var2 > 0)) {
      significant[i] <- -1
    } else if (any(var2 > 0) && any(var2 < 0)) {
      significant[i] <- 5
    } else {
      significant[i] <- 0
    }
  }
  
  listtestresults <- cbind(listtestresults, significant)
  colnames(listtestresults)[ncol(listtestresults)] <- "SignificantResults"
  
  return(listtestresults)
}

# ------------------------------------------------------------
# 2.7 Build clustering input values
# ------------------------------------------------------------
.compute_allvalues <- function(listtestresults, HDXdatastates, timepoints, cantdeut) {
  exchangableamides <- rep(NA_real_, nrow(HDXdatastates[[1]]))
  allvalues <- c()
  
  for (i in seq_len(nrow(HDXdatastates[[1]]))) {
    peptideseq <- listtestresults[i, 4]
    exchangableamides[i] <- str_length(peptideseq) - str_count(peptideseq, "P") - cantdeut
    
    for (h in seq_len(timepoints)) {
      delta_col <- (h * 2) + 3
      if (!is.na(listtestresults[i, delta_col])) {
        allvalues <- c(allvalues, listtestresults[i, delta_col] / exchangableamides[i])
      }
    }
  }
  
  return(as.data.frame(allvalues))
}

# ------------------------------------------------------------
# 2.8 Utility helper for rounding numeric data frames
# ------------------------------------------------------------
round_df <- function(df, digits = 3) {
  df[] <- lapply(df, function(x) {
    if (is.numeric(x)) round(x, digits) else x
  })
  df
}

`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !is.na(a) && nzchar(a)) a else b
}

flag_img_url <- function(code) {
  code <- tolower(trimws(code))
  if (is.na(code) || nchar(code) != 2) {
    return("https://flagcdn.com/24x18/un.png")
  }
  paste0("https://flagcdn.com/24x18/", code, ".png")
}

to_hex_color <- function(x) {
  rgb_mat <- grDevices::col2rgb(x)
  grDevices::rgb(rgb_mat[1], rgb_mat[2], rgb_mat[3], maxColorValue = 255)
}

# ------------------------------------------------------------
# Progress / status helpers
# ------------------------------------------------------------
make_progress <- function(session, message = "Working...", max = 1, button_ids = NULL) {
  p <- shiny::Progress$new(session, min = 0, max = max)
  p$set(message = message, value = 0)
  
  if (!is.null(button_ids)) {
    for (id in button_ids) {
      try(shinyjs::disable(id), silent = TRUE)
    }
  }
  
  list(
    set = function(value = NULL, detail = NULL, message = NULL) {
      p$set(value = value, detail = detail, message = message)
    },
    inc = function(amount = 0.1, detail = NULL, message = NULL) {
      p$inc(amount, detail = detail, message = message)
    },
    close = function() {
      p$close()
      if (!is.null(button_ids)) {
        for (id in button_ids) {
          try(shinyjs::enable(id), silent = TRUE)
        }
      }
    }
  )
}

# ------------------------------------------------------------
# Error / status helpers
# ------------------------------------------------------------

stop_user <- function(...) {
  stop(paste0(...), call. = FALSE)
}

check_required_columns <- function(df, required, where = "input file") {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop_user(
      "The ", where, " is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

# ------------------------------------------------------------
# Lightweight error/status helpers
# ------------------------------------------------------------
format_user_error <- function(prefix, e) {
  msg <- conditionMessage(e)
  paste0(prefix, ": ", msg)
}

# ============================================================
# 3. Main server logic
# ============================================================
server <- function(input, output, session) {
  app_status <- reactiveVal("Ready")
  app_error  <- reactiveVal(NULL)
  last_error <- reactiveVal(NULL)
  
  
  output$app_error_banner <- renderUI({
    req(app_error())
    
    div(
      style = paste(
        "background:#f8d7da;",
        "color:#842029;",
        "border:1px solid #f5c2c7;",
        "border-radius:8px;",
        "padding:12px 14px;",
        "margin-bottom:15px;"
      ),
      tags$div(
        style = "display:flex; justify-content:space-between; align-items:flex-start; gap:12px;",
        tags$div(
          tags$strong("Please check your inputs"),
          tags$div(style = "margin-top:4px;", app_error())
        ),
        actionLink("dismiss_app_error", "Dismiss", style = "white-space:nowrap;")
      )
    )
  })
  
  observeEvent(input$dismiss_app_error, {
    app_error(NULL)
    last_error(NULL)
  })
  
  set_app_error <- function(msg) {
    # avoid spamming the same message repeatedly
    if (!identical(last_error(), msg)) {
      last_error(msg)
      app_error(msg)
    } else {
      app_error(msg)
    }
  }
  
  clear_app_error <- function() {
    app_error(NULL)
    last_error(NULL)
  }
  
  last_initialdata <- reactiveVal(NULL)
  last_variables <- reactiveVal(NULL)
  last_fitting <- reactiveVal(NULL)
  
  visitors_data <- reactiveFileReader(
    intervalMillis = 2000,
    session = session,
    filePath = visitor_file,
    readFunc = function(path) {
      if (file.exists(path)) {
        tryCatch(readRDS(path), error = function(e) {
          data.frame(
            country_code = character(),
            country_name = character(),
            visits = integer(),
            stringsAsFactors = FALSE
          )
        })
      } else {
        data.frame(
          country_code = character(),
          country_name = character(),
          visits = integer(),
          stringsAsFactors = FALSE
        )
      }
    }
  )
  # ------------------------------------------------------------
  # 3.1 Basic reactive wrappers
  # ------------------------------------------------------------
  cantdeut <- reactive({
    as.numeric(input$cantdeutbutton)
  })
  
  timepoints <- reactive({
    input$timepoints
  })
  
  replicates <- reactive({
    input$replicates
  })
  
  states <- reactive({
    input$states
  })
  
  significancelevel <- reactive({
    input$significancelevel
  })
  
  labelingtimepoints <- reactive({
    as.numeric(unlist(strsplit(input$labelingtimepoints, ",")))
  })
  
  visitor_logged <- reactiveVal(FALSE)
  
  observeEvent(input$visitor_country, {
    req(input$visitor_country)
    
    if (visitor_logged()) return()
    
    info <- input$visitor_country
    country_code <- toupper(info$country_code %||% "UN")
    country_name <- info$country_name %||% "Unknown"
    
    visitors <- visitors_data()
    
    if (country_code %in% visitors$country_code) {
      idx <- match(country_code, visitors$country_code)
      visitors$visits[idx] <- visitors$visits[idx] + 1
    } else {
      visitors <- rbind(
        visitors,
        data.frame(
          country_code = country_code,
          country_name = country_name,
          visits = 1,
          stringsAsFactors = FALSE
        )
      )
    }
    
    write_visitors(visitors)
    visitor_logged(TRUE)
  }, ignoreInit = TRUE)
  
  # ------------------------------------------------------------
  # 3.2 Input tab UI outputs and controls
  # ------------------------------------------------------------
  output$state1box <- renderUI({
    req(initialdata())
    selectInput(inputId = "state1", label = "Select reference state:", choices = "", selected = NULL)
  })
  
  output$state2box <- renderUI({
    req(initialdata())
    selectInput(inputId = "state2", label = "Select state 2:", choices = "", selected = NULL)
  })
  
  observeEvent(input$refresh, {
    session$reload()
  })
  
  output$analyze <- renderUI({
    req(
      input$hdexaminerfile,
      input$fastafile,
      input$timepoints,
      input$replicates,
      input$states,
      input$significancelevel,
      input$labelingtimepoints
    )
    
    actionButton("start", "Calculate")
  })
  
  # ------------------------------------------------------------
  # 3.3 Input data import and preprocessing
  # Reads uploaded HDX-MS files, standardizes column names,
  # filters invalid rows, and prepares the base analysis dataset
  # ------------------------------------------------------------
  initialdata <- eventReactive(input$start, {
    req(
      input$hdexaminerfile,
      input$fastafile,
      input$timepoints,
      input$replicates,
      input$states,
      input$significancelevel,
      input$labelingtimepoints
    )
    
    app_status("Loading input data...")
    
    tp <- timepoints()
    rp <- replicates()
    st <- states()
    lt <- labelingtimepoints()
    
    tryCatch({
      
      # --------------------------------------------------------
      # Input checks BEFORE reading files
      # --------------------------------------------------------
      if (length(lt) != tp) {
        stop_user(
          "You entered ", length(lt),
          " labeling times, but 'Number of labeling times' is set to ",
          tp, ". Please make these match."
        )
      }
      
      if (anyDuplicated(lt)) {
        stop_user("The labeling times contain duplicates. Please enter unique values.")
      }
      
      if (rp < 2) {
        stop_user("Number of replicates must be at least 2.")
      }
      
      if (st < 2) {
        stop_user("Number of protein states must be at least 2.")
      }
      
      # --------------------------------------------------------
      # 3.3.1 Read and normalize input format
      # --------------------------------------------------------
      if (length(input$datacheck) == 0) {
        rawdata <- read.csv(
          input$hdexaminerfile$datapath,
          header = FALSE,
          stringsAsFactors = FALSE
        )
        
        if (nrow(rawdata) < 3) {
          stop_user("The HDExaminer file appears to be empty or malformed.")
        }
        
        HDXdata <- rawdata[-1, ]
        colnames(HDXdata) <- HDXdata[1, ]
        HDXdata <- HDXdata[-1, ]
        colnames(HDXdata)[2] <- "pepnumber"
        HDXdata <- HDXdata[-c(6, 8)]
        
        for (x in seq_len(tp * rp)) {
          needed_col <- x + 6
          if (needed_col > ncol(HDXdata)) {
            stop_user(
              "The HDExaminer file does not contain enough uptake columns for ",
              tp, " timepoints and ", rp, " replicates."
            )
          }
          
          HDXdata <- HDXdata[-c(x + 6, x + 7, x + 9, x + 10, x + 11, x + 12, x + 13)]
          HDXdata[, x + 6] <- sapply(HDXdata[, x + 6], as.numeric)
        }
        
      } else if (input$datacheck == "WATERS/DynamX output file") {
        watersdata <- read.csv(
          input$hdexaminerfile$datapath,
          header = TRUE,
          stringsAsFactors = FALSE
        )
        
        check_required_columns(
          watersdata,
          c("Sequence", "Start", "End", "State", "Exposure", "Center", "z"),
          where = "WATERS/DynamX file"
        )
        
        watersdata <- watersdata[-c(5, 6, 7, 8, 11, 13, 14)]
        watersdata["Exposure"] <- watersdata["Exposure"] * 60
        watersdata$Exposure <- as.integer(watersdata$Exposure)
        
        watersdata <- watersdata %>%
          group_by(Sequence, State, z, Exposure) %>%
          dplyr::mutate(Exposure = paste0(Exposure, "s_", row_number()))
        
        watersdata <- watersdata %>%
          pivot_wider(names_from = Exposure, values_from = Center) %>%
          relocate(State)
        
        wdata <- subset(watersdata, select = c(1:6))
        wdata2 <- subset(watersdata, select = c(7:length(watersdata)))
        
        names <- c()
        for (i in seq_len(rp)) {
          names <- append(names, paste0(lt, "s_", i))
        }
        
        df <- data.frame(matrix(ncol = tp * rp, nrow = 0))
        names(df) <- names
        wdata2[setdiff(names(df), names(wdata2))] <- NA
        wdata2 <- wdata2[, order(as.numeric(sub("\\D*(\\d+).*", "\\1", colnames(wdata2))))]
        
        watersdata <- cbind(wdata, wdata2)
        
        names(watersdata)[2] <- "pepnumber"
        names(watersdata)[6] <- "Charge"
        watersdata <- data.frame(watersdata, check.names = FALSE)
        
        nonddata <- watersdata %>% select(starts_with("0s_"))
        nonddata$mean <- rowMeans(nonddata, na.rm = TRUE)
        
        watersdata <- cbind(watersdata, mean = nonddata$mean)
        initialdata <- watersdata %>% select(-starts_with("0s_"))
        initialdata <- initialdata %>% drop_na(mean)
        
        for (i in 7:(tp * rp + 6)) {
          initialdata[, i] <- (initialdata[, i] - initialdata$mean) * initialdata[, 6]
        }
        
        HDXdata <- initialdata %>% select(-mean)
        HDXdata$pepnumber <- ""
        
      } else {
        HDXdata <- read.csv(
          input$hdexaminerfile$datapath,
          header = TRUE,
          stringsAsFactors = FALSE
        )
        
        check_required_columns(
          HDXdata,
          c("State", "Start", "End", "Sequence"),
          where = "custom CSV file"
        )
        
        HDXdata[1:6] <- BiocGenerics::lapply(HDXdata[1:6], as.character)
        HDXdata$pepnumber <- ""
      }
      
      raw_rows <- nrow(HDXdata)
      
      # --------------------------------------------------------
      # 3.3.2 Remove rows with excessive missing values
      # --------------------------------------------------------
      HDXdata <- HDXdata[
        rowSums(is.na(HDXdata[, 7:((tp * rp) + 6)])) <= (tp * rp) - rp,
      ]
      
      # --------------------------------------------------------
      # 3.3.3 Keep only peptides observed in all states
      # --------------------------------------------------------
      HDXdata <- HDXdata %>%
        group_by(Sequence, Charge) %>%
        filter(n_distinct(State) == st)
      
      if (nrow(HDXdata) == 0) {
        stop_user(
          "After filtering, no peptides remain in all ", st,
          " states. Please check the file format and state count."
        )
      }
      
      pepcharge <- HDXdata[, 6]
      HDXdata <- HDXdata[, -6]
      HDXdata <- data.frame(HDXdata)
      
      # --------------------------------------------------------
      # 3.3.4 Rename uptake columns using labeling times/replicates
      # --------------------------------------------------------
      if ((5 + tp * rp) > ncol(HDXdata)) {
        stop_user(
          "The input file does not contain enough uptake columns for ",
          tp, " timepoints and ", rp, " replicates."
        )
      }
      
      colnames(HDXdata)[6:(5 + tp * rp)] <- unlist(lapply(seq_len(tp), function(i) {
        paste0(lt[i], "s_", seq_len(rp))
      }))
      
      shinyjs::show("functional")
      
      # --------------------------------------------------------
      # 3.3.5 Build simple validation diagnostics
      # --------------------------------------------------------
      uptake_cols <- grep("s_$|s_[0-9]+$", names(HDXdata), value = FALSE)
      if (length(uptake_cols) == 0) {
        uptake_cols <- which(grepl("s_", names(HDXdata), fixed = TRUE))
      }
      
      diagnostics <- list(
        rows_imported = raw_rows,
        rows_kept = nrow(HDXdata),
        rows_removed = raw_rows - nrow(HDXdata),
        detected_states = unique(as.character(HDXdata$State)),
        n_states_detected = length(unique(as.character(HDXdata$State))),
        n_unique_peptides = if (all(c("Sequence", "Start", "End") %in% colnames(HDXdata))) {
          dplyr::n_distinct(paste(HDXdata$Sequence, HDXdata$Start, HDXdata$End, sep = "_"))
        } else {
          NA_integer_
        },
        missing_uptake_values = if (length(uptake_cols) > 0) {
          sum(is.na(HDXdata[, uptake_cols, drop = FALSE]))
        } else {
          NA_integer_
        },
        total_uptake_cells = if (length(uptake_cols) > 0) {
          nrow(HDXdata) * length(uptake_cols)
        } else {
          NA_integer_
        },
        max_end = if ("End" %in% colnames(HDXdata)) {
          safe_max(HDXdata$End)
        } else {
          NA_real_
        }
      )
      
      result <- list(
        HDXdata = HDXdata,
        pepcharge = pepcharge,
        diagnostics = diagnostics
      )
      
      last_initialdata(result)
      clear_app_error()
      return(result)
      
    }, error = function(e) {
      set_app_error(
        format_user_error("Could not process input files", e)
       )
      
      # Fallback to last successful result if it exists
      last_good <- last_initialdata()
      if (!is.null(last_good)) {
        return(last_good)
      } else {
        return(NULL)
      }
    })
  })
  
  
  # ------------------------------------------------------------
  # 3.4 Data validation summary
  # Generates a quick QC summary after input preprocessing
  # ------------------------------------------------------------
  output$validation_summary <- renderUI({
    req(initialdata(), input$fastafile)
    
    diagnostics <- initialdata()$diagnostics
    
    fasta_length <- tryCatch({
      fasta_seq <- seqinr::read.fasta(
        input$fastafile$datapath,
        seqtype = "AA",
        whole.header = FALSE,
        seqonly = FALSE,
        as.string = TRUE
      )
      fasta_seq <- gsub(" ", "", as.character(fasta_seq))
      nchar(fasta_seq)
    }, error = function(e) {
      NA_integer_
    })
    
    missing_pct <- if (!is.na(diagnostics$total_uptake_cells) && diagnostics$total_uptake_cells > 0) {
      round(100 * diagnostics$missing_uptake_values / diagnostics$total_uptake_cells, 1)
    } else {
      NA_real_
    }
    
    fasta_ok <- if (!is.na(fasta_length) && !is.na(diagnostics$max_end)) {
      fasta_length >= diagnostics$max_end
    } else {
      NA
    }
    
    status_text <- if (isTRUE(fasta_ok)) {
      "FASTA length is compatible with the peptide residue range."
    } else if (identical(fasta_ok, FALSE)) {
      "FASTA appears shorter than the maximum peptide end position. Please check sequence numbering."
    } else {
      "FASTA compatibility could not be evaluated."
    }
    
    status_color <- if (isTRUE(fasta_ok)) {
      "#dff0d8"
    } else if (identical(fasta_ok, FALSE)) {
      "#f2dede"
    } else {
      "#fcf8e3"
    }
    
    state_tags <- lapply(diagnostics$detected_states, function(x) {
      tags$span(
        x,
        style = "display:inline-block; background:#3c8dbc; color:white; padding:3px 8px; border-radius:10px; margin-right:6px; margin-bottom:4px; font-size:12px;"
      )
    })
    
    box(
      title = "Data validation summary",
      width = 12,
      status = "primary",
      solidHeader = TRUE,
      collapsible = TRUE,
      
      fluidRow(
        column(
          3,
          tags$div(style = "font-weight:bold;", "Rows imported"),
          tags$div(style = "font-size:18px;", diagnostics$rows_imported)
        ),
        column(
          3,
          tags$div(style = "font-weight:bold;", "Rows kept"),
          tags$div(style = "font-size:18px;", diagnostics$rows_kept)
        ),
        column(
          3,
          tags$div(style = "font-weight:bold;", "Rows removed"),
          tags$div(style = "font-size:18px;", diagnostics$rows_removed)
        ),
        column(
          3,
          tags$div(style = "font-weight:bold;", "Unique peptides"),
          tags$div(style = "font-size:18px;", diagnostics$n_unique_peptides)
        )
      ),
      
      tags$hr(),
      
      fluidRow(
        column(
          4,
          tags$div(style = "font-weight:bold; margin-bottom:5px;", "Detected states"),
          tags$div(state_tags)
        ),
        column(
          4,
          tags$div(style = "font-weight:bold;", "Missing uptake values"),
          tags$div(
            style = "font-size:16px;",
            if (!is.na(missing_pct)) {
              paste0(diagnostics$missing_uptake_values, " (", missing_pct, "%)")
            } else {
              diagnostics$missing_uptake_values
            }
          )
        ),
        column(
          4,
          tags$div(style = "font-weight:bold;", "FASTA length / max peptide end"),
          tags$div(
            style = "font-size:16px;",
            paste0(fasta_length, " / ", diagnostics$max_end)
          )
        )
      ),
      
      tags$div(
        style = paste0(
          "margin-top:15px; padding:10px; border-radius:6px; background:",
          status_color,
          ";"
        ),
        tags$strong("Status: "),
        status_text
      )
    )
  })
  
  # ------------------------------------------------------------
  # 3.5 Input consistency checks and selector updates
  # ------------------------------------------------------------
  observe({
    if (length(input$datacheck) == 2) {
      showModal(modalDialog(
        title = "Please select only ONE option",
        footer = NULL,
        easyClose = TRUE
      ))
      updateCheckboxGroupInput(
        session,
        "datacheck",
        choices = c("WATERS/DynamX output file", "Custom .csv output file"),
        selected = NULL
      )
    }
  })
  
  observe({
    req(initialdata())
    HDXdata <- initialdata()$HDXdata
    states <- unique(HDXdata$State)
    updateSelectInput(session, "state1", label = "Select reference state:", choices = states, selected = states[1])
  })
  
  observe({
    req(initialdata())
    HDXdata <- initialdata()$HDXdata
    states <- unique(HDXdata$State)
    updateSelectInput(session, "state2", label = "Select state 2:", choices = states[!states == input$state1])
  })
  
  # ------------------------------------------------------------
  # 3.6 Main analysis reactive
  # ------------------------------------------------------------
  variables <- reactive({
    req(initialdata(), input$state1, input$state2, input$state1 != input$state2)
    
    app_status("Running statistical analysis...")
    
    tryCatch({
      
      tp <- timepoints()
      rp <- replicates()
      st <- states()
      lt <- labelingtimepoints()
      sl <- significancelevel()
      cd <- cantdeut()
      
      # ----------------------------
      # 3.6.1 Load FASTA sequence
      # ----------------------------
      sequence <- read.fasta(
        input$fastafile$datapath,
        seqtype = "AA",
        whole.header = FALSE,
        seqonly = FALSE,
        as.string = TRUE
      )
      sequence <- gsub(" ", "", sequence)
      
      if (nchar(sequence) == 0) {
        stop_user("The FASTA file could not be read correctly or is empty.")
      }
      
      withProgress(message = "Analyzing data", value = 0, {
        incProgress(1 / 10)
        
        # ----------------------------
        # 3.6.2 Split input table by selected states
        # ----------------------------
        HDXdata <- initialdata()$HDXdata
        HDXdata$State <- factor(HDXdata$State, levels = unique(HDXdata$State))
        HDXdatastates <- list()
        
        otherstates <- unique(HDXdata$State)
        otherstates <- otherstates[!otherstates == input$state1]
        otherstates <- otherstates[!otherstates == input$state2]
        
        HDXdatastates[[1]] <- filter(HDXdata, State == input$state1)
        HDXdatastates[[2]] <- filter(HDXdata, State == input$state2)
        
        if (nrow(HDXdatastates[[1]]) == 0 || nrow(HDXdatastates[[2]]) == 0) {
          stop_user("One of the selected states has no peptides after filtering.")
        }
        
        if (st > 2) {
          for (i in 3:st) {
            HDXdatastates[[i]] <- filter(HDXdata, State == otherstates[i - 2])
          }
        }
        
        # ----------------------------
        # 3.6.3 Build qFeatures objects
        # ----------------------------
        data <- add_column(initialdata()$HDXdata, initialdata()$pepcharge, .after = "Sequence")
        data <- data %>% dplyr::relocate(Charge, .after = Sequence)
        TF1 <- data
        
        HDXqf <- .build_qf(TF1, replicates = rp, timepoints = tp, states = st)
        
        incProgress(2 / 10)
        
        subHDXqf <- list()
        if (st > 2) {
          subTF2 <- TF1 %>% filter(State %in% c(input$state1, input$state2))
          subHDXqf <- .build_qf(subTF2, replicates = rp, timepoints = tp, states = 2)
        }
        
        incProgress(3 / 10)
        
        # ----------------------------
        # 3.6.4 Build heatmap input
        # ----------------------------
        heatmap <- pheatmap(
          t(assay(HDXqf)),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          color = brewer.pal(n = 9, name = "BuPu"),
          main = "Protein: Deuterium Incorporation",
          fontsize = 8,
          legend_breaks = c(0, 1, 2, 3, 4, 5, 6, max(assay(HDXqf))),
          legend_labels = c("0", "1", "2", "3", "4", "5", "6", "Incorporation")
        )
        
        # ----------------------------
        # 3.6.5 Add mean/SD columns per state
        # ----------------------------
        HDXdatastates <- .add_mean_sd(
          HDXdatastates = HDXdatastates,
          timepoints = tp,
          replicates = rp,
          labelingtimepoints = lt,
          start_state = 1,
          end_state = 2
        )
        
        if (isTRUE(input$showallstates) && st > 2) {
          HDXdatastates <- .add_mean_sd(
            HDXdatastates = HDXdatastates,
            timepoints = tp,
            replicates = rp,
            labelingtimepoints = lt,
            start_state = 3,
            end_state = st
          )
        }
        
        incProgress(4 / 10)
        
        # ----------------------------
        # 3.6.6 Compute pooled SD and significance threshold
        # ----------------------------
        pooledSD <- .calc_pooled_sd(HDXdatastates, timepoints = tp, replicates = rp)
        SEM <- sqrt(2 * (pooledSD^2) / rp)
        threshold <- qt(p = sl / 2, df = (2 * rp) - 2, lower.tail = FALSE) * SEM
        
        parameters <- data.frame(
          PooledSD = pooledSD,
          StrdErrorMean = SEM,
          alpha = sl,
          tValUsed = qt(p = sl / 2, df = (2 * rp) - 2, lower.tail = FALSE),
          StatisticalThreshold = threshold
        )
        
        parameters[5, 1] <- "Peptides presenting protection and deprotection at different labeling times are colored yellow (Check data in HDExaminer)"
        
        # ----------------------------
        # 3.6.7 Compute peptide-level statistics
        # ----------------------------
        stats_out <- .compute_stats(
          HDXdatastates = HDXdatastates,
          timepoints = tp,
          replicates = rp,
          cantdeut = cd,
          labelingtimepoints = lt,
          significancelevel = sl
        )
        
        listtestresults <- stats_out$listtestresults
        clusteringresults <- stats_out$clusteringresults
        
        incProgress(5 / 10)
        
        # ----------------------------
        # 3.6.8 Classify significance and build clustering inputs
        # ----------------------------
        listtestresults <- .classify_significance(
          listtestresults = listtestresults,
          timepoints = tp,
          threshold = threshold,
          significancelevel = sl
        )
        
        allvalues <- .compute_allvalues(
          listtestresults = listtestresults,
          HDXdatastates = HDXdatastates,
          timepoints = tp,
          cantdeut = cd
        )
        
        kmeansclust <- kmeans(abs(na.omit(allvalues)), 3, iter.max = 200)
        centers <- sort(kmeansclust$centers)
        
        limit1 <- -(((centers[3] - centers[2]) / 2) + centers[2])
        limit2 <- -(((centers[2] - centers[1]) / 2) + centers[1])
        limit3 <- -limit2
        limit4 <- -limit1
        
        incProgress(8 / 10)
        
        # ----------------------------
        # 3.6.9 Final clustering assignments and reporting objects
        # ----------------------------
        for (l in 1:nrow(HDXdatastates[[1]])) {
          for (m in 1:tp) {
            cluster_col <- (m * 3) + 3
            delta_col <- (m * 3) + 2
            p_col <- (m * 3) + 4
            
            if (is.na(clusteringresults[l, cluster_col])) {
              clusteringresults[l, cluster_col] <- 1
            } else if (clusteringresults[l, cluster_col] > limit4) {
              clusteringresults[l, cluster_col] <- 3
            } else if (clusteringresults[l, cluster_col] > limit3) {
              clusteringresults[l, cluster_col] <- 2
            } else if (clusteringresults[l, cluster_col] > limit2) {
              clusteringresults[l, cluster_col] <- 1
            } else if (clusteringresults[l, cluster_col] > limit1) {
              clusteringresults[l, cluster_col] <- 2
            } else {
              clusteringresults[l, cluster_col] <- 3
            }
            
            if (is.na(abs(clusteringresults[l, delta_col]))) {
              clusteringresults[l, cluster_col] <- 1
            } else if (abs(clusteringresults[l, delta_col]) < threshold) {
              clusteringresults[l, cluster_col] <- 1
            }
            
            if (is.na(clusteringresults[l, p_col])) {
              clusteringresults[l, cluster_col] <- 1
            } else if (-log10(clusteringresults[l, p_col]) < (-log10(sl))) {
              clusteringresults[l, cluster_col] <- 1
            }
          }
        }
        
        clustvec <- c()
        for (i in 1:nrow(clusteringresults)) {
          for (j in 1:tp) {
            delta_col <- (j * 3) + 2
            cluster_col <- (j * 3) + 3
            p_col <- (j * 3) + 4
            
            if (is.na(clusteringresults[i, delta_col])) clusteringresults[i, delta_col] <- 0
            if (is.na(clusteringresults[i, p_col])) clusteringresults[i, p_col] <- 0
            
            if (abs(clusteringresults[i, delta_col]) >= threshold && clusteringresults[i, p_col] <= sl) {
              clustvec <- append(clustvec, clusteringresults[i, cluster_col])
            } else {
              clustvec <- append(clustvec, 1)
            }
          }
          
          clusteringresults[i, (3 * tp) + 5] <- max(clustvec)
          clustvec <- c()
        }
        
        colnames(clusteringresults)[(3 * tp) + 5] <- "Highest significant cluster in peptide"
        
        a <- ncol(clusteringresults)
        clusteringresults[1, a + 1] <- "Cluster sizes"
        colnames(clusteringresults)[a + 1] <- ""
        clusteringresults[2, (a + 1):(a + 3)] <- kmeansclust$size
        colnames(clusteringresults)[a + 2] <- ""
        colnames(clusteringresults)[a + 3] <- ""
        clusteringresults[3, a + 1] <- "Cluster 1 (abs values)"
        clusteringresults[4, a + 1] <- "Cluster 2 (abs values)"
        clusteringresults[5, a + 1] <- "Cluster 3 (abs values)"
        clusteringresults[3, a + 2] <- paste("From", 0, "to", signif(limit3, 3))
        clusteringresults[4, a + 2] <- paste("From", signif(limit3, 3), "to", signif(limit4, 3))
        clusteringresults[5, a + 2] <- paste("From", signif(limit4, 3), "to max value")
        clusteringresults[6, a + 1] <- "k-means clustering iterations"
        clusteringresults[6, a + 2] <- kmeansclust$iter
        clusteringresults[7, a + 1] <- "Datapoints in cluster 1 are considered insignificant or negligible within your dataset"
        clusteringresults[8, a + 1] <- "Datapoints in cluster 2 are intermediate effects within your dataset"
        clusteringresults[9, a + 1] <- "Datapoints in cluster 3 are strong effects within your dataset"
        
        incProgress(9 / 10)
      })
      
      if (nchar(sequence) < max(as.numeric(HDXdatastates[[1]]$End))) {
        set_app_error(
          "Your FASTA file appears to have fewer residues than the maximum peptide end position. Please check the FASTA sequence or residue numbering.",
          status = "Input warning"
        )
      }
      
      result <- list(
        HDXdatastates = HDXdatastates,
        listtestresults = listtestresults,
        parameters = parameters,
        threshold = threshold,
        clusteringresults = clusteringresults,
        sequence = sequence,
        centers = centers,
        allvalues = allvalues,
        HDXqf = HDXqf,
        heatmap = heatmap,
        subHDXqf = subHDXqf
      )
      
      last_variables(result)
      clear_app_error()
      return(result)
      
    }, error = function(e) {
      set_app_error(
        format_user_error("Could not analyze the current inputs", e),
        status = "Analysis error"
      )
      
      last_good <- last_variables()
      if (!is.null(last_good)) {
        return(last_good)
      } else {
        return(NULL)
      }
    })
  })
  
  # ------------------------------------------------------------
  # 3.7 Heatmap outputs and download
  # ------------------------------------------------------------
  output$heatmap <- renderPlot({
    req(variables())
    return(variables()$heatmap)
  })
  
  output$heatbutton <- downloadHandler(
    filename = function() {
      paste0("Heat map", input$heattype)
    },
    content = function(file) {
      ggsave(file, variables()$heatmap, height = 7, width = 15)
    }
  )
  
  # ------------------------------------------------------------
  # 3.8 Functional fitting section
  # ------------------------------------------------------------
  realfit <- reactiveVal(value = 1)
  
  observe({
    if (input$selectedeq == "equationwithoutd") {
      shinyjs::disable("d_input")
      shinyjs::enable("p_input")
      updateNumericInput(session, "d_input", value = "NULL")
      updateNumericInput(session, "p_input", value = 1)
    }
    if (input$selectedeq == "equationwithd") {
      shinyjs::enable("d_input")
      shinyjs::enable("p_input")
      updateNumericInput(session, "d_input", value = 0)
      updateNumericInput(session, "p_input", value = 1)
    }
    if (input$selectedeq == "equationwithdwithoutP") {
      shinyjs::enable("d_input")
      shinyjs::disable("p_input")
      updateNumericInput(session, "d_input", value = 0)
      updateNumericInput(session, "p_input", value = "NULL")
    }
    if (input$selectedeq == "equationwithoutdorP") {
      shinyjs::disable("d_input")
      shinyjs::disable("p_input")
      updateNumericInput(session, "d_input", value = "NULL")
      updateNumericInput(session, "p_input", value = "NULL")
    }
  })
  
  counter <- reactiveVal(value = 0)
  
  observeEvent(input$reset, {
    counter(counter() + 1)
  })
  
  event_trigger <- reactive({
    list(input$state1, input$state2)
  })
  
  observeEvent(ignoreInit = TRUE, event_trigger(), {
    click("reset")
  })
  
  fitting <- eventReactive(input$startfitting, {
    req(variables())
    
    app_status("Running kinetic fitting...")
    
    tryCatch({
      
      st <- states()
      all_peptides <- rownames(variables()$HDXqf)[[1]]
      
      if (length(all_peptides) == 0) {
        stop_user("No peptides are available for fitting.")
      }
      
      if (st > 2) {
        suball_peptides <- rownames(variables()$subHDXqf)[[1]]
      }
      
      # ----------------------------
      # 3.8.1 Build starting parameter list
      # ----------------------------
      if (input$selectedeq == "equationwithd") {
        starting_parameters <- list(a = input$a_input, b = input$b_input, d = input$d_input, p = input$p_input)
      } else if (input$selectedeq == "equationwithdwithoutP") {
        starting_parameters <- list(a = input$a_input, b = input$b_input, d = input$d_input)
      } else if (input$selectedeq == "equationwithoutd") {
        starting_parameters <- list(a = input$a_input, b = input$b_input, p = input$p_input)
      } else {
        starting_parameters <- list(b = input$b_input, a = input$a_input)
      }
      
      subresults <- list()
      
      withProgress(message = "Fitting data", min = 0, max = 1, {
        incProgress(1 / 10)
        
        # ----------------------------
        # 3.8.2 Run hdxstats fitting
        # ----------------------------
        results <- .fit_model(
          HDXqf = variables()$HDXqf,
          all_peptides = all_peptides,
          starting_parameters = starting_parameters,
          input = input
        )
        
        if (st > 2) {
          subresults <- .fit_model(
            HDXqf = variables()$subHDXqf,
            all_peptides = suball_peptides,
            starting_parameters = starting_parameters,
            input = input
          )
        }
        
        incProgress(6 / 10)
        
        # ----------------------------
        # 3.8.3 Build fitted-plot outputs
        # ----------------------------
        graphics_kinetics <- hdxstats::visualise_hdx_data(results, type = "kinetics")
        graphics_forest <- if (st == 2) {
          hdxstats::visualise_hdx_data(results, type = "forest")
        } else {
          hdxstats::visualise_hdx_data(subresults, type = "forest")
        }
        
        summary <- BiocGenerics::lapply(results[["fitted_models"]]@statmodels, summary)
        
        incProgress(9 / 10)
      })
      
      shinyjs::show("fdabox")
      realfit(2)
      counter(0)
      
      result <- list(
        kinetics = graphics_kinetics,
        forest = graphics_forest,
        fittingparameters = results$functional_analysis@results,
        summary = summary,
        results = results,
        subresults = subresults
      )
      
      last_fitting(result)
      clear_app_error()
      return(result)
      
    }, error = function(e) {
      set_app_error(
        format_user_error("Could not fit the kinetic model", e),
        status = "Fitting error"
      )
      
      last_good <- last_fitting()
      if (!is.null(last_good)) {
        return(last_good)
      } else {
        return(NULL)
      }
    })
  })
  
  # ------------------------------------------------------------
  # 3.9 Statistical summary strip and fitted peptide outputs
  # ------------------------------------------------------------
  output$stat_summary_strip <- renderUI({
    req(variables(), input$state1, input$state2)
    
    listtestresults <- variables()$listtestresults
    threshold <- round(as.numeric(variables()$threshold), 3)
    
    sig <- as.character(listtestresults$SignificantResults)
    
    total_peptides <- nrow(listtestresults)
    n_protected <- sum(sig == "-1", na.rm = TRUE)
    n_deprotected <- sum(sig == "1", na.rm = TRUE)
    n_mixed <- sum(sig == "5", na.rm = TRUE)
    n_nochange <- sum(sig == "0", na.rm = TRUE)
    
    tagList(
      fluidRow(
        column(
          2,
          valueBox(
            value = total_peptides,
            subtitle = "Total peptides",
            icon = icon("table"),
            color = "light-blue",
            width = 12
          )
        ),
        column(
          2,
          valueBox(
            value = n_protected,
            subtitle = "Protected",
            icon = icon("shield"),
            color = "blue",
            width = 12
          )
        ),
        column(
          2,
          valueBox(
            value = n_deprotected,
            subtitle = "Deprotected",
            icon = icon("unlock"),
            color = "red",
            width = 12
          )
        ),
        column(
          2,
          valueBox(
            value = n_mixed,
            subtitle = "Mixed",
            icon = icon("shuffle"),
            color = "yellow",
            width = 12
          )
        ),
        column(
          2,
          valueBox(
            value = n_nochange,
            subtitle = "No change",
            icon = icon("minus"),
            color = "green",
            width = 12
          )
        ),
        column(
          2,
          valueBox(
            value = threshold,
            subtitle = "Threshold (Da)",
            icon = icon("sliders"),
            color = "purple",
            width = 12
          )
        )
      ),
      
      tags$div(
        style = "margin-top: -5px; margin-bottom: 12px; padding: 8px 12px; background: #f7f7f7; border-radius: 6px; font-size: 14px;",
        tags$strong("Comparison: "),
        paste0(input$state1, " vs ", input$state2)
      )
    )
  })
  
  # ------------------------------------------------------------
  # 3.9.1 Single fitted peptide plot
  # ------------------------------------------------------------
  singlefittedplot <- reactive({
    req(fitting())
    
    graphics_kinetics <- fitting()$kinetics
    graphics_forest <- fitting()$forest
    summary <- fitting()$summary
    wanted <- as.numeric(str_split(input$selpeptide2, ":")[[1]][1])
    
    null <- data.frame(coef(summary[[wanted]]$null), check.names = FALSE)
    null$variable <- substr(rownames(null), 1, 1)
    null <- null %>% select(variable, everything())
    null[, 6] <- "Null"
    colnames(null)[6] <- "State"
    
    for (i in 1:states()) {
      df <- data.frame(coef(summary[[wanted]][[i + 1]]), check.names = FALSE)
      df$variable <- substr(rownames(df), 1, 1)
      df <- df %>% select(variable, everything())
      df[, 6] <- paste0("State", i)
      colnames(df)[6] <- "State"
      null <- rbind(null, df)
    }
    
    if (is.null(graphics_forest[[wanted]])) {
      plot2 <- NULL
    } else {
      plot2 <- graphics_forest[[wanted]] +
        ggtitle(paste0("Forest plot ", input$state1, " vs. ", input$state2)) +
        theme(plot.title = element_text(hjust = 0.5))
    }
    
    return(list(plot1 = graphics_kinetics[[wanted]], plot2 = plot2, summary = null))
  })
  
  output$kineticplot <- renderPlot({
    print(singlefittedplot()$plot1)
  })
  
  output$forestplot <- renderPlot({
    test <- counter()
    text <- paste("\n Your reference/comparison states have changed. Please refit the data.\n")
    
    if (test != 0) {
      ggplot() +
        annotate("text", x = 4, y = 25, size = 7, label = text) +
        theme_void()
    } else {
      singlefittedplot()$plot2
    }
  })
  
  output$summary <- renderTable({
    singlefittedplot()$summary
  }, digits = 3)
  
  observe({
    req(fitting())
    labels <- c()
    for (i in 1:length(fitting()$kinetics)) {
      labels <- append(labels, paste0(i, ": ", fitting()$kinetics[[i]]$data$rowname[1]))
    }
    updateSelectInput(session, "selpeptide2", label = "Select peptide:", choices = labels)
  })
  
  output$kineticbutton <- downloadHandler(
    filename = function() {
      paste0("Fitted plot", input$kinetictype)
    },
    content = function(file) {
      pdf(file, width = 7, height = 7)
      print(singlefittedplot()$plot1)
      print(singlefittedplot()$plot2)
      dev.off()
    }
  )
  
  kineticdocument <- reactive({
    req(fitting())
    graphs <- list()
    a <- 1
    
    withProgress(message = "Creating document", min = 0.5, max = 1, {
      for (i in 1:length(fitting()$kinetics)) {
        graphs[[a]] <- fitting()$kinetics[[i]]
        graphs[[a + 1]] <- fitting()$forest[[i]]
        a <- a + 2
      }
      
      return(marrangeGrob(graphs, ncol = 2, nrow = 2, top = NULL, layout_matrix = matrix(1:4, 2, 2, TRUE)))
    })
  })
  
  observeEvent(input$downloadallkinetic, {
    showModal(modalDialog(
      title = "Please wait. Creating your document...",
      easyClose = FALSE,
      modalButton("Cancel"),
      downloadButton("allkineticbutton", "Download all kinetic plots"),
      footer = NULL
    ))
    
    shinyjs::disable(id = "allkineticbutton")
    a <- kineticdocument()
    shinyjs::enable(id = "allkineticbutton")
  })
  
  output$allkineticbutton <- downloadHandler(
    filename = "Allfittedplots.pdf",
    content = function(file) {
      withProgress(message = "Exporting", min = 0, max = 1, {
        ggsave(file, kineticdocument(), width = 14, height = 8)
      })
    }
  )
  
  # ------------------------------------------------------------
  # 3.10 Manhattan plot section
  # ------------------------------------------------------------
  observe({
    req(initialdata())
    updateSelectInput(
      session,
      "seltimemanhattan",
      label = "Select labeling time:",
      choices = as.numeric(unlist(strsplit(input$labelingtimepoints, ",")))
    )
  })
  
  manhattan <- reactive({
    req(fitting())
    
    if (states() == 2) {
      results <- fitting()$results
      HDXqf <- variables()$HDXqf
    } else {
      results <- fitting()$subresults
      HDXqf <- variables()$subHDXqf
    }
    
    data <- add_column(initialdata()$HDXdata, initialdata()$pepcharge, .after = "Sequence")
    labelingtimes <- as.numeric(unlist(strsplit(input$labelingtimepoints, ",")))
    time <- which(labelingtimes == input$seltimemanhattan)[1]
    
    diffdata <- rowMeans(
      assay(HDXqf)[, (((replicates() * (time - 1)) + 1):(replicates() * time))],
      na.rm = TRUE
    ) - rowMeans(
      assay(HDXqf)[, (((replicates() * timepoints()) + (replicates() * (time - 1)) + 1):((replicates() * timepoints()) + (replicates() * time)))],
      na.rm = TRUE
    )
    
    out <- hdxstats::processFunctional(object = HDXqf, params = results$fitted_models)
    
    data$names <- paste0(data$Sequence, "_", data$Charge)
    data$peptide <- paste0("[", data$Start, ":", data$End, "]")
    data <- data %>% distinct(names, .keep_all = TRUE)
    
    df <- data.frame(out@results)
    df$names <- rownames(df)
    
    diffdata <- data.frame(diffdata)
    diffdata$names <- rownames(diffdata)
    diffdata <- drop_na(diffdata)
    
    signs <- inner_join(diffdata, df)
    names <- inner_join(data, df)
    all <- inner_join(signs, names)
    all <- all %>% distinct(names, .keep_all = TRUE)
    
    manhplot <- ggplot(all, aes(x = 1:nrow(all), y = -log10(ebayes.fdr), color = factor(as.vector(sign(diffdata))))) +
      geom_point(size = 3.5) +
      scale_color_manual(
        breaks = c(-1, 1, 0),
        values = c(input$colorflex, input$colorprot, "gray"),
        labels = c("deprotected", "protected", "no difference")
      ) +
      scale_x_continuous(
        expand = c(0, 1.2),
        guide = guide_axis(n.dodge = 3),
        breaks = 1:length(all$fitcomplete),
        labels = all$peptide
      ) +
      theme_classic() +
      geom_hline(
        yintercept = -log10(input$significancelevel),
        linetype = "dashed",
        colour = "red",
        linewidth = 2
      ) +
      xlab("Peptide") +
      ylab("-log(ebayes FDR)") +
      theme(
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
      ) +
      labs(color = "Effect") +
      ggtitle(paste0("Manhattan plot ", input$state1, " vs. ", input$state2)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    return(manhplot)
  })
  
  output$manhattanplot <- renderPlot({
    test <- counter()
    text <- paste("\n Your reference/comparison states have changed. Please refit the data.\n")
    
    if (test != 0) {
      ggplot() +
        annotate("text", x = 4, y = 25, size = 10, label = text) +
        theme_void()
    } else {
      print(manhattan())
    }
  })
  
  output$manhattanbutton <- downloadHandler(
    filename = function() {
      paste0("ManhattanPlot", input$seltimemanhattan, "s", input$manhattantype)
    },
    content = function(file) {
      ggsave(
        file,
        (manhattan() + ggtitle("Manhattan plot") + theme(plot.title = element_text(hjust = 0.5, size = 20))),
        height = 7,
        width = 14
      )
    }
  )
  
  # ------------------------------------------------------------
  # 3.11 Volcano plot and clustering histogram
  # ------------------------------------------------------------
  observe({
    if (input$clusteringsel == TRUE) shinyjs::show("panel")
    if (input$clusteringsel == FALSE) shinyjs::hide("panel")
  })
  
  observeEvent(input$start, {
    shinyjs::show("panel2")
  })
  
  volcanoplot <- reactive({
    req(variables())
    
    threshold <- variables()$threshold
    listtestresults <- variables()$listtestresults
    clusteringresults <- variables()$clusteringresults
    
    maxvalue <- c()
    maxpvalue <- c()
    
    for (i in 1:timepoints()) {
      maxvalue <- append(maxvalue, listtestresults[, (2 * i) + 3])
      maxpvalue <- append(maxpvalue, listtestresults[, (2 * i) + 4])
    }
    
    volcanoplot <- ggplot(data = clusteringresults) +
      BiocGenerics::lapply(seq(timepoints()), function(f) {
        
        delta_vals <- as.numeric(clusteringresults[, (f * 3) + 2])
        p_vals <- as.numeric(clusteringresults[, (f * 3) + 4])
        cluster_vals <- as.numeric(clusteringresults[, (f * 3) + 3])
        
        effect_vals <- ifelse(
          is.na(delta_vals) | is.na(p_vals),
          "Missing",
          ifelse(
            delta_vals >= threshold & p_vals < significancelevel(),
            "Deprotected",
            ifelse(
              delta_vals <= -threshold & p_vals < significancelevel(),
              "Protected",
              "Not significant"
            )
          )
        )
        
        tooltip_vals <- paste0(
          "Peptide: ", listtestresults$pepnumber,
          "<br>Sequence: ", listtestresults$Sequence,
          "<br>Residues: ", listtestresults$Start, "-", listtestresults$End,
          "<br>Time: ", labelingtimepoints()[f], " s",
          "<br>ΔD: ", round(delta_vals, 3),
          "<br>p-value: ", signif(p_vals, 3),
          "<br>Class: ", effect_vals
        )
        
        if (input$clusteringsel == TRUE) {
          geom_point(
            aes(
              x = clusteringresults[, (f * 3) + 2],
              y = -log10(clusteringresults[, (f * 3) + 4]),
              color = factor(clusteringresults[, (f * 3) + 3]),
              text = tooltip_vals
            ),
            na.rm = TRUE
          )
        } else {
          geom_point(
            aes(
              x = clusteringresults[, (f * 3) + 2],
              y = -log10(clusteringresults[, (f * 3) + 4]),
              text = tooltip_vals
            ),
            na.rm = TRUE,
            color = "black"
          )
        }
      })
    
    ymaxlim <- if (-log10(min(abs(maxpvalue), na.rm = TRUE)) <= -log10(significancelevel())) {
      -log10(significancelevel()) + 0.1
    } else {
      -log10(min(abs(maxpvalue), na.rm = TRUE))
    }
    
    volcanoplot <- volcanoplot +
      geom_segment(x = threshold, xend = 100, y = -log10(significancelevel()), yend = -log10(significancelevel()), color = "red", linewidth = 1.1, linetype = "dashed") +
      geom_segment(x = -threshold, xend = -100, y = -log10(significancelevel()), yend = -log10(significancelevel()), color = "red", linewidth = 1.1, linetype = "dashed") +
      geom_segment(x = threshold, xend = threshold, y = -log10(significancelevel()), yend = 100, color = "red", linewidth = 1.1, linetype = "dashed") +
      geom_segment(x = -threshold, xend = -threshold, y = -log10(significancelevel()), yend = 100, color = "red", linewidth = 1.1, linetype = "dashed") +
      xlab("mass difference(Da)") +
      ylab("-log(p value)") +
      xlim((0 - max(abs(maxvalue), na.rm = TRUE)), (0 + max(abs(maxvalue), na.rm = TRUE))) +
      ylim(0, ymaxlim) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 12),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
      )
    
    if (input$clusteringsel == TRUE) {
      volcanoplot <- volcanoplot +
        scale_colour_manual(values = c("gray50", input$colorint, input$colorstr)) +
        theme(legend.position = "none")
    }
    
    return(volcanoplot)
  })
  
  output$volcanoplot <- renderPlotly({
    req(volcanoplot())
    
    ggplotly(volcanoplot(), tooltip = "text") %>%
      layout(dragmode = "zoom") %>%
      config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("select2d", "lasso2d", "autoScale2d")
      )
  })
  
  output$binsize <- renderUI({
    req(input$clusteringsel == TRUE)
    selectInput(inputId = "selbinsize", label = "Select number of bins:", c("little", "few", "many"), selected = "little")
  })
  
  histogram <- reactive({
    req(input$selbinsize)
    
    centers <- variables()$centers
    allvalues <- variables()$allvalues
    
    limit1 <- -(((centers[3] - centers[2]) / 2) + centers[2])
    limit2 <- -(((centers[2] - centers[1]) / 2) + centers[1])
    limit3 <- -limit2
    limit4 <- -limit1
    
    bins_n <- switch(input$selbinsize, little = 15L, few = 50L, many = 100L)
    bins <- seq(min(allvalues), max(allvalues), length.out = bins_n)
    
    hist <- ggplot(allvalues, aes(x = allvalues)) +
      geom_histogram(data = subset(allvalues, allvalues <= limit1), color = "black", fill = input$colorstr, breaks = bins) +
      geom_histogram(data = subset(allvalues, allvalues <= limit2 & allvalues > limit1), color = "black", fill = input$colorint, breaks = bins) +
      geom_histogram(data = subset(allvalues, allvalues <= limit3 & allvalues > limit2), color = "black", fill = "gray50", breaks = bins) +
      geom_histogram(data = subset(allvalues, allvalues <= limit4 & allvalues > limit3), color = "black", fill = input$colorint, breaks = bins) +
      geom_histogram(data = subset(allvalues, allvalues > limit4), color = "black", fill = input$colorstr, breaks = bins) +
      geom_vline(aes(xintercept = limit1), color = "red", linetype = "dashed") +
      geom_vline(aes(xintercept = limit2), color = "red", linetype = "dashed") +
      geom_vline(aes(xintercept = limit3), color = "red", linetype = "dashed") +
      geom_vline(aes(xintercept = limit4), color = "red", linetype = "dashed") +
      theme_bw() +
      xlab("normalized mass difference(Da)") +
      ylab("frequency")
    
    return(hist)
  })
  
  output$volcanotext <- renderUI({
    req(volcanoplot())
    tagList(
      tags$iframe(
        id = "volcano_iframe",
        src = "./volcano.html",
        width = "100%",
        height = "100px",
        frameborder = 0,
        scrolling = "no",
        align = "center",
        onload = "this.style.height = this.contentWindow.document.body.scrollHeight + 'px';"
      )
    )
  })
  
  output$histogramclustplot <- renderPlot({
    req(input$clusteringsel)
    if (input$clusteringsel == TRUE) print(histogram())
  })
  
  output$volcanotype <- renderUI({
    req(volcanoplot())
    selectInput("volcanotypes", label = NULL, choices = c(".pdf", ".svg", ".eps"), selected = ".pdf")
  })
  
  output$volcanodownload <- renderUI({
    req(volcanoplot())
    downloadButton("volcanobutton", "Download Volcano plot")
  })
  
  output$histdownload <- renderUI({
    req(histogram())
    req(input$clusteringsel == TRUE)
    downloadButton("histbutton", "Download Histogram")
  })
  
  output$histbutton <- downloadHandler(
    filename = function() {
      paste0("Histogram", input$volcanotypes)
    },
    content = function(file) {
      ggsave(
        file,
        histogram() + ggtitle("Histogram of normalized differences") +
          theme(plot.title = element_text(hjust = 0.5, size = 20))
      )
    }
  )
  
  output$volcanobutton <- downloadHandler(
    filename = function() {
      paste0("VolcanoPlot", input$volcanotypes)
    },
    content = function(file) {
      ggsave(
        file,
        volcanoplot() + ggtitle("Hybrid significance test") +
          theme(plot.title = element_text(hjust = 0.5, size = 20))
      )
    }
  )
 
  # ------------------------------------------------------------
  # 3.12 Global woods plot section
  # ------------------------------------------------------------
  woods_zoom <- reactiveValues(x = NULL, y = NULL)
  
  woodsplots <- reactive({  # calculation of global woods plot
    
    req(variables())
    
    listtestresults <- variables()$listtestresults
    sequence <- variables()$sequence
    
    listtestresults$SignificantResults <- as.numeric(listtestresults$SignificantResults)
    listtestresults$Start <- as.numeric(listtestresults$Start)
    listtestresults$End <- as.numeric(listtestresults$End)
    listtestresults$pepnumber <- as.numeric(listtestresults$pepnumber)
    
    colors <- c("0" = "gray43", "-1" = input$colorprot, "1" = input$colorflex, "5" = "yellow")
    
    p <- ggplot() +
      geom_rect(
        data = listtestresults,
        color = "black",
        linewidth = 0.1,
        mapping = aes(
          xmin = pepnumber - 0.40,
          xmax = pepnumber + 0.40,
          ymin = Start,
          ymax = End,
          fill = factor(SignificantResults)
        )
      ) +
      scale_fill_manual(values = colors, guide = "none") +
      scale_x_continuous(
        expand = c(0, 2),
        breaks = seq(0, 1000, by = 20)
      ) +
      theme_bw() +
      labs(x = "Peptide number", y = "Protein residue") +
      theme(
        legend.position = "",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      scale_y_continuous(
        expand = c(0, 2),
        breaks = seq(0, 1000, by = 20),
        limits = c(0, nchar(sequence) + 5)
      )
    
    return(p)
  })
  
  output$globwoodsplot <- renderPlotly({
    req(variables())
    
    df <- as.data.frame(variables()$listtestresults, stringsAsFactors = FALSE)
    sequence <- as.character(variables()$sequence)
    
    required_cols <- c("pepnumber", "Start", "End", "Sequence", "SignificantResults")
    validate(
      need(
        all(required_cols %in% colnames(df)),
        paste("Global Woods plot is missing required columns:",
              paste(setdiff(required_cols, colnames(df)), collapse = ", "))
      )
    )
    
    df$pepnumber <- suppressWarnings(as.numeric(df$pepnumber))
    df$Start <- suppressWarnings(as.numeric(df$Start))
    df$End <- suppressWarnings(as.numeric(df$End))
    df$SignificantResults <- as.character(df$SignificantResults)
    
    df <- df[!is.na(df$pepnumber) & !is.na(df$Start) & !is.na(df$End), , drop = FALSE]
    
    validate(
      need(nrow(df) > 0, "No valid peptide rectangles available for the Global Woods plot.")
    )
    
    color_map <- c(
      "0"  = to_hex_color("gray43"),
      "-1" = to_hex_color(input$colorprot),
      "1"  = to_hex_color(input$colorflex),
      "5"  = to_hex_color("yellow")
    )
    
    p <- plotly::plot_ly()
    
    for (i in seq_len(nrow(df))) {
      sig_i <- df$SignificantResults[i]
      if (!(sig_i %in% names(color_map))) sig_i <- "0"
      
      fill_col <- color_map[[sig_i]]
      
      hover_txt <- paste0(
        "Peptide: ", df$pepnumber[i],
        "<br>Sequence: ", df$Sequence[i],
        "<br>Residues: ", df$Start[i], "-", df$End[i]
      )
      
      p <- p %>%
        plotly::add_trace(
          x = c(
            df$pepnumber[i] - 0.40,
            df$pepnumber[i] + 0.40,
            df$pepnumber[i] + 0.40,
            df$pepnumber[i] - 0.40,
            df$pepnumber[i] - 0.40
          ),
          y = c(
            df$Start[i],
            df$Start[i],
            df$End[i],
            df$End[i],
            df$Start[i]
          ),
          type = "scatter",
          mode = "lines",
          fill = "toself",
          fillcolor = fill_col,
          line = list(color = "black", width = 0.5),
          text = hover_txt,
          hoverinfo = "text",
          hoveron = "fills",
          showlegend = FALSE
        )
    }
    
    p %>%
      plotly::layout(
        showlegend = FALSE,
        hovermode = "closest",
        dragmode = "zoom",
        margin = list(l = 60, r = 20, t = 20, b = 50),
        xaxis = list(
          title = "Peptide number",
          ticks = "outside",
          showgrid = FALSE,
          zeroline = FALSE,
          showline = TRUE,
          linecolor = "black",
          linewidth = 1,
          range = c(min(df$pepnumber, na.rm = TRUE) - 1, max(df$pepnumber, na.rm = TRUE) + 1)
        ),
        yaxis = list(
          title = "Protein residue",
          ticks = "outside",
          showgrid = FALSE,
          zeroline = FALSE,
          showline = TRUE,
          linecolor = "black",
          linewidth = 1,
          range = c(0, nchar(sequence) + 5)
        )
      ) %>%
      plotly::config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("select2d", "lasso2d", "autoScale2d")
      )
  })
  
  output$globwoodstype <- renderUI({
    req(woodsplots())
    selectInput("woodstypes", label = NULL, choices = c(".pdf", ".svg", ".eps"), selected = ".pdf")
  })
  
  output$woodsplottext <- renderUI({
    req(woodsplots())
    tags$iframe(
      id = "woods_iframe",
      src = "./woodstext.html",
      width = "100%",
      height = "100px",
      frameborder = 0,
      scrolling = "no",
      align = "center",
      onload = "this.style.height = (this.contentWindow.document.body.scrollHeight + 20) + 'px';"
    )
  })
  
  output$globwoodsdownload <- renderUI({
    req(woodsplots())
    downloadButton("woodsbutton", "Download")
  })
  
  output$woodsbutton <- downloadHandler(
    filename = function() {
      paste0("Global woods plots", input$woodstypes)
    },
    content = function(file) {
      ggsave(file, woodsplots() + ggtitle("Global woods plot") + theme(plot.title = element_text(hjust = 0.5, size = 20)))
    }
  )
  
  # ------------------------------------------------------------
  # 3.13 Woods-by-timepoint section
  # ------------------------------------------------------------
  observe({
    req(
      input$hdexaminerfile,
      input$fastafile,
      input$timepoints,
      input$replicates,
      input$states,
      input$significancelevel,
      input$labelingtimepoints
    )
    
    updateSelectInput(session, "timepoint", label = "Select labeling time:", choices = labelingtimepoints())
  })
  
  woodsbytimepoint <- reactive({ # wood's plots by time point
    req(variables())
    
    listtestresults <- variables()$listtestresults
    HDXdatastates <- variables()$HDXdatastates
    threshold <- variables()$threshold
    
    plotsbytimepoint <- list()
    
    listtestresults$Start <- as.numeric(HDXdatastates[[1]]$Start)
    listtestresults$End <- as.numeric(HDXdatastates[[1]]$End)
    
    maxvalue2 <- c()
    for (b in 1:timepoints()) {
      maxvalue2 <- append(maxvalue2, listtestresults[, (b * 2) + 3])
    }
    maxvalue2 <- as.numeric(maxvalue2)
    
    for (a in 1:timepoints()) {
      
      # Create fixed vectors for this specific timepoint
      delta_values <- as.numeric(listtestresults[, (2 * a) + 3])
      p_values <- as.numeric(listtestresults[, (2 * a) + 4])
      
      color1 <- c()
      for (i in 1:nrow(listtestresults)) {
        if (is.na(delta_values[i]) || is.na(p_values[i])) {
          color1 <- append(color1, 0)
        } else if ((delta_values[i] >= threshold) && (p_values[i] < significancelevel())) {
          color1 <- append(color1, 1)
        } else if ((delta_values[i] <= -threshold) && (p_values[i] < significancelevel())) {
          color1 <- append(color1, -1)
        } else {
          color1 <- append(color1, 0)
        }
      }
      
      color1 <- as.character(color1)
      
      plotdata <- data.frame(
        Start = listtestresults$Start,
        End = listtestresults$End,
        delta = delta_values,
        ymin = delta_values - 0.05,
        ymax = delta_values + 0.05,
        color1 = color1,
        Sequence = listtestresults$Sequence
      )
      
      plotdata$tooltip <- paste0(
        "Sequence: ", plotdata$Sequence,
        "<br>Residues: ", plotdata$Start, "-", plotdata$End,
        "<br>ΔHX: ", round(plotdata$delta, 3),
        "<br>Time: ", labelingtimepoints()[a], " s"
      )
      
      plotsbytimepoint[[a]] <- ggplot(plotdata) +
        geom_rect(
          aes(
            xmin = Start,
            xmax = End,
            ymin = ymin,
            ymax = ymax,
            fill = color1,
            text = tooltip
          ),
          color = "black",
          linewidth = 0.15,
          na.rm = TRUE
        ) +
        scale_x_continuous(expand = c(0, 2), breaks = seq(-1000, 1000, by = 35)) +
        scale_y_continuous(breaks = seq(-1000, 1000, by = 0.5)) +
        theme_bw() +
        ylim(
          -max(abs(maxvalue2), na.rm = TRUE) - 0.05,
          max(abs(maxvalue2), na.rm = TRUE) + 0.05
        ) +
        labs(
          x = "Residue",
          y = paste("delta HX (Da) ", labelingtimepoints()[a], "s", sep = "")
        ) +
        geom_hline(yintercept = threshold, color = "red", linewidth = 1.0, linetype = "dashed") +
        geom_hline(yintercept = -threshold, color = "red", linewidth = 1.0, linetype = "dashed") +
        theme(
          legend.position = "",
          panel.grid.major = element_blank()
        ) +
        scale_fill_manual(values = c("0" = "gray43", "-1" = input$colorprot, "1" = input$colorflex))
    }
    return(plotsbytimepoint)
  })
  
  output$woodsbytimepointplot <- renderPlotly({
    req(woodsbytimepoint(), input$timepoint)
    
    wanted <- as.numeric(input$timepoint)
    number <- match(wanted, labelingtimepoints())
    
    validate(
      need(!is.na(number), "Selected labeling time not found.")
    )
    
    plots <- woodsbytimepoint()
    p <- plotly::ggplotly(plots[[number]], tooltip = "text")
    
    # Hide legend if plotly recreates it
    for (i in seq_along(p$x$data)) {
      p$x$data[[i]]$showlegend <- FALSE
    }
    
    p %>%
      plotly::layout(
        showlegend = FALSE,
        dragmode = "zoom",
        margin = list(l = 60, r = 20, t = 20, b = 60),
        xaxis = list(
          ticks = "outside",
          showgrid = FALSE,
          zeroline = FALSE,
          showline = TRUE,
          linecolor = "black",
          linewidth = 1
        ),
        yaxis = list(
          ticks = "outside",
          showgrid = FALSE,
          zeroline = FALSE,
          showline = TRUE,
          linecolor = "black",
          linewidth = 1
        )
      ) %>%
      plotly::config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("select2d", "lasso2d", "autoScale2d")
      )
  })
  
  output$woodsbytimepointtype <- renderUI({
    req(woodsbytimepoint())
    selectInput("woodsbytimepointtypes", label = NULL, choices = c(".pdf", ".svg", ".eps"), selected = ".pdf")
  })
  
  output$woodsplotbytimetext <- renderUI({
    req(woodsbytimepoint())
    tags$iframe(
      id = "woodsbytime_iframe",
      src = "./woodsbytimetext.html",
      width = "100%",
      height = "100px",
      frameborder = 0,
      scrolling = "no",
      align = "center",
      onload = "
      var iframe = this;
      setTimeout(function() {
        iframe.style.height = (iframe.contentWindow.document.body.scrollHeight + 30) + 'px';
      }, 100);
    "
    )
  })
  
  output$woodsbytimepointdownload <- renderUI({
    req(woodsbytimepoint())
    downloadButton("woodsbytimepointbutton", "Download all")
  })
  
  output$woodsbytimepointbutton <- downloadHandler(
    filename = function() {
      paste0("Woodsplot by timepoint", input$woodsbytimepointtypes)
    },
    content = function(file) {
      ggsave(file, marrangeGrob(woodsbytimepoint(), nrow = timepoints(), ncol = 1, top = "Woods-plot by timepoint"), width = 8.5, height = 11)
    }
  )
  
  # ------------------------------------------------------------
  # 3.14 Digestion efficiency section
  # ------------------------------------------------------------
  digestion <- reactive({
    req(variables())
    
    HDXdatastates <- variables()$HDXdatastates
    sequence <- variables()$sequence
    
    coverage <- as.character(sequence)
    
    for (i in 1:2) {
      HDXdatastates[[i]]$Start <- as.numeric(HDXdatastates[[i]]$Start)
      HDXdatastates[[i]]$End <- as.numeric(HDXdatastates[[i]]$End)
    }
    
    for (i in 1:nrow(HDXdatastates[[1]])) {
      substr(
        coverage,
        HDXdatastates[[1]][i, 3],
        HDXdatastates[[1]][i, 4]
      ) <- strrep(
        "Z",
        HDXdatastates[[1]][i, 4] - HDXdatastates[[1]][i, 3] + 1
      )
    }
    
    digestion <- data.frame()
    digestion[1, 1] <- nrow(HDXdatastates[[1]])
    digestion[1, 2] <- signif((str_count(coverage, "Z") / str_length(coverage)) * 100, 4)
    digestion[1, 3] <- signif(mean(HDXdatastates[[1]]$End - HDXdatastates[[1]]$Start), 4)
    
    redundancy <- rep(0, str_length(sequence))
    for (i in 1:nrow(HDXdatastates[[1]])) {
      for (j in HDXdatastates[[1]][i, 3]:HDXdatastates[[1]][i, 4]) {
        redundancy[j] <- redundancy[j] + 1
      }
    }
    
    digestion[1, 4] <- signif(mean(redundancy), 4)
    names(digestion) <- c("Number of peptides", "% coverage", "Avg peptide length", "Redundancy")
    
    avglength <- ggplot(data = HDXdatastates[[1]], aes(x = End - Start)) +
      geom_histogram(binwidth = 1, color = "black", fill = "gray") +
      theme_bw() +
      xlab("Peptide length") +
      ylab("frequency") +
      geom_vline(
        aes(xintercept = mean(HDXdatastates[[1]]$End - HDXdatastates[[1]]$Start)),
        color = "red",
        linewidth = 1.5,
        linetype = "dashed"
      )
    
    efftable <- gridExtra::tableGrob(digestion, rows = NULL)
    
    download_plot <- gridExtra::arrangeGrob(
      avglength,
      efftable,
      ncol = 1,
      heights = c(3, 1),
      top = NULL
    )
    
    return(list(
      plot = avglength,              # used in the UI
      digestion = digestion,         # used in the summary panel
      download_plot = download_plot  # used in the download
    ))
  })
  
  
  output$digestionplot <- renderPlot({
    req(digestion())
    print(digestion()$plot)
  })
  
  output$digestion_summary_panel <- renderUI({
    req(digestion())
    
    dg <- digestion()$digestion
    
    metric_box <- function(title, value, bg = "#f7f7f7") {
      div(
        style = paste0(
          "background:", bg, ";",
          "border:1px solid #dddddd;",
          "border-radius:10px;",
          "padding:14px 16px;",
          "text-align:center;",
          "height:100%;"
        ),
        tags$div(
          style = "font-size:13px; color:#666; margin-bottom:6px; font-weight:600;",
          title
        ),
        tags$div(
          style = "font-size:24px; font-weight:700; color:#222;",
          value
        )
      )
    }
    
    div(
      class = "instruction-card",
      h4("Digestion summary", style = "text-align:center; margin-top:0;"),
      tags$p(
        class = "small-note",
        style = "text-align:center;",
        "Quick overview of the digestion quality metrics for the current dataset."
      ),
      
      fluidRow(
        column(
          3,
          metric_box("Number of peptides", dg[1, "Number of peptides"], "#eef3f8")
        ),
        column(
          3,
          metric_box("% coverage", paste0(dg[1, "% coverage"], "%"), "#eef8ee")
        ),
        column(
          3,
          metric_box("Avg peptide length", dg[1, "Avg peptide length"], "#fff8e8")
        ),
        column(
          3,
          metric_box("Redundancy", dg[1, "Redundancy"], "#f7eef8")
        )
      )
    )
  })
  outputOptions(output, "digestionplot", suspendWhenHidden = FALSE)
  
  outputOptions(output, "digestionplot", suspendWhenHidden = FALSE)
  
  output$digestiontype <- renderUI({
    req(digestion())
    selectInput("digestiontypes", label = NULL, choices = c(".pdf", ".svg", ".eps"), selected = ".pdf")
  })
  
  output$digestiondownload <- renderUI({
    req(digestion())
    downloadButton("digestionbutton", "Download")
  })
  
  output$digestionbutton <- downloadHandler(
    filename = function() {
      paste0("Digestion", input$digestiontypes)
    },
    content = function(file) {
      ggsave(
        file,
        plot = digestion()$download_plot,
        height = 8,
        width = 7
      )
    }
  )
  
  
  output$digestiontext <- renderUI({
    req(digestion())
    tags$iframe(src = "./digestiontext.html", width = "100%", height = "310px", frameborder = 0, scrolling = "auto", align = "center")
  })
  
  # ------------------------------------------------------------
  # 3.15 Peptide map section
  # ------------------------------------------------------------
  peptidemap <- reactive({
    req(variables())
    req(digestion())
    
    digestion <- digestion()$digestion
    shinyjs::show("pep2")
    
    HDXdatastates <- variables()$HDXdatastates
    sequence <- variables()$sequence
    
    HDXdatastates[[1]] <- HDXdatastates[[1]][, -1]
    pepmap <- HDXdatastates[[1]][, 1:3]
    pepmap[, 4] <- rep(1, nrow(pepmap))
    pepmap[, 1] <- as.numeric(pepmap[, 1])
    pepmap[, 2] <- as.numeric(pepmap[, 2])
    pepmap[, 3] <- as.numeric(pepmap[, 3])
    
    height <- list()
    height[[1]] <- rep(0, str_length(sequence))
    
    for (b in 1:nrow(pepmap)) {
      c <- 1
      if (all(height[[c]][pepmap[b, 2]:pepmap[b, 3]] == 0)) {
        pepmap[b, 4] <- c
        height[[c]][pepmap[b, 2]:pepmap[b, 3]] <- rep(1, pepmap[b, 3] - pepmap[b, 2] + 1)
      } else {
        while (is.element(1, height[[c]][pepmap[b, 2]:pepmap[b, 3]])) {
          c <- c + 1
          
          if ((c - 1) == length(height)) {
            height[[c]] <- rep(0, str_length(sequence))
          }
          
          if (all(height[[c]][pepmap[b, 2]:pepmap[b, 3]] == 0)) {
            pepmap[b, 4] <- c
            height[[c]][pepmap[b, 2]:pepmap[b, 3]] <- rep(1, pepmap[b, 3] - pepmap[b, 2] + 1)
            break
          }
        }
      }
    }
    
    coveragemap <- ggplot(data = pepmap) +
      geom_rect(
        aes(xmin = Start, xmax = End, ymin = pepmap[, 4] - 0.5, ymax = pepmap[, 4] + 0.5),
        fill = input$colorpep,
        color = "black",
        linewidth = 0.07
      ) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text = element_text(size = 8),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        axis.line.x.top = element_line(color = "white"),
        axis.text.x.top = element_text(face = "bold")
      )
    
    labels <- data.frame(
      x = c(
        seq(1, str_length(sequence), length.out = str_length(sequence)),
        c(1, (1:(str_length(sequence) / 25)) * 25)
      ),
      label = c(
        format(s2c(sequence), length.out = str_length(sequence)),
        paste0("\n", c(1, (1:(str_length(sequence) / 25)) * 25))
      )
    )
    
    plots <- list()
    for (d in 1:(ceiling(str_length(sequence) / 100))) {
      plots[[d]] <- coveragemap +
        scale_x_continuous(
          breaks = labels$x,
          labels = labels$label,
          minor_breaks = NULL,
          expand = c(0, 0),
          limits = c(0, NA)
        ) +
        coord_cartesian(xlim = c(((d - 1) * 100) + 1, (d * 100)))
    }
    
    coveragemap2 <- ggplot(data = pepmap) +
      geom_rect(aes(xmin = Start, xmax = End, ymin = 0, ymax = 1), fill = input$colorpep) +
      theme_void() +
      theme(panel.background = element_rect(fill = "gray")) +
      scale_x_continuous(expand = c(0, 0), limits = c(1, nchar(sequence))) +
      scale_y_continuous(expand = c(0, 0)) +
      ggtitle(paste0(digestion[1, 2], "% coverage with ", digestion[1, 1], " peptides")) +
      theme(plot.title = element_text(face = "bold"))
    
    coveragepep <- ggplot(data = pepmap) +
      geom_rect(aes(xmin = Start, xmax = End, ymin = -pepmap[, 4] - 0.44, ymax = -pepmap[, 4] + 0.44), fill = input$colorpep, color = "black", linewidth = 0.03) +
      scale_x_continuous(expand = c(0, 0), breaks = seq(0, nchar(sequence), by = 50), limits = c(1, nchar(sequence))) +
      scale_y_continuous(expand = c(0, 0)) +
      xlab("residue number") +
      theme(
        axis.text.x = element_text(face = "bold"),
        axis.line.x = element_line(color = "black", linewidth = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      )
    
    summary <- coveragemap2 / coveragepep + plot_layout(heights = c(1, 3))
    
    return(list(
      peptidemap = marrangeGrob(plots, nrow = ceiling(str_length(sequence) / 100), ncol = 1, top = NULL),
      summary = summary
    ))
  })
  
  output$summaryplot <- renderPlot({
    peptidemap()$summary
  })
  
  output$summarypepbutton <- downloadHandler(
    filename = function() {
      paste0("Summary Peptide map", input$peptidetypes)
    },
    content = function(file) {
      ggsave(file, peptidemap()$summary, height = 3.5, width = 14)
    }
  )
  
  output$peptideplot <- renderPlot({
    req(peptidemap())
    shinyjs::show("pep")
    print(peptidemap()$peptidemap)
  })
  
  output$peptidetype <- renderUI({
    req(peptidemap())
    selectInput("peptidetypes", label = NULL, choices = c(".pdf", ".svg", ".eps"), selected = ".pdf")
  })
  
  output$peptidedownload <- renderUI({
    req(peptidemap())
    downloadButton("peptidebutton", "Download")
  })
  
  output$peptidebutton <- downloadHandler(
    filename = function() {
      paste0("Peptide map", input$peptidetypes)
    },
    content = function(file) {
      ggsave(file, peptidemap()$summary)
    }
  )
  
  output$peptidemaptext <- renderUI({
    req(peptidemap())
    tags$iframe(
      id = "peptidemap_iframe",
      src = "./peptidemaptext.html",
      width = "100%",
      height = "100px",
      frameborder = 0,
      scrolling = "no",
      align = "center",
      onload = "
      var iframe = this;
      setTimeout(function() {
        iframe.style.height = (iframe.contentWindow.document.body.scrollHeight + 30) + 'px';
      }, 100);
    "
    )
  })
  
  # ------------------------------------------------------------
  # 3.16 Uptake plot section
  # ------------------------------------------------------------
  observe({
    req(variables())
    HDXdatastates <- variables()$HDXdatastates
    
    pepfield <- if ("pepnum" %in% colnames(HDXdatastates[[1]])) "pepnum" else "pepnumber"
    
    updateSelectInput(
      session,
      "selpeptide",
      label = "Select peptide:",
      choices = paste0(
        HDXdatastates[[1]][[pepfield]],
        ": (",
        HDXdatastates[[1]]$Start, "-",
        HDXdatastates[[1]]$End,
        ") ",
        HDXdatastates[[1]]$Sequence
      )
    )
  })
  
  output$uptaketype <- renderUI({
    req(variables())
    selectInput("uptaketypes", label = NULL, choices = c(".pdf", ".svg", ".eps"), selected = ".pdf")
  })
  
  # ------------------------------------------------------------
  # 3.16.1 Build single-peptide uptake plot
  # ------------------------------------------------------------
  singleplot <- reactive({
    req(variables())
    shinyjs::show("plotssettings")
    shinyjs::show("limits")
    req(input$coloring, input$thickness, input$fontsize)
    
    if (states() > 2) {
      shinyjs::show("allplotscheck")
    }
    
    wanted <- input$selpeptide
    selected <- as.numeric(str_split(wanted, ":")[[1]][1])
    
    HDXdatastates <- variables()$HDXdatastates
    listtestresults <- variables()$listtestresults
    tp <- timepoints()
    
    if (requireNamespace("ggthemr", quietly = TRUE)) {
      ggthemr::ggthemr(input$coloring)
    }
    
    deut <- c()
    sd <- c()
    state <- c()
    times <- c()
    
    end <- if (isTRUE(input$showallstates)) states() else 2
    
    for (h in 1:end) {
      rowdata <- HDXdatastates[[h]][HDXdatastates[[h]]$pepnumber == selected, , drop = FALSE]
      
      deut <- c(deut, as.numeric(rowdata %>% select(matches("Average"))))
      sd <- c(sd, as.numeric(rowdata %>% select(matches("SD"))))
      state <- c(state, rep(as.character(rowdata$State), each = tp))
      times <- c(times, labelingtimepoints())
    }
    
    uptake <- data.frame(Prostate = state, deuteration = deut, stddev = sd, time = times)
    
    sig_result <- listtestresults %>% filter(pepnumber == selected) %>% pull(SignificantResults)
    title_suffix <- ifelse(sig_result == 5, "**", ifelse(abs(sig_result) == 1, "*", ""))
    
    title <- paste(
      HDXdatastates[[1]][selected, 3], "-",
      HDXdatastates[[1]][selected, 4], ":",
      HDXdatastates[[1]][selected, 5], title_suffix,
      sep = ""
    )
    
    plot <- ggplot(data = uptake, aes(x = time, y = deuteration, color = Prostate)) +
      scale_x_continuous(trans = "log10") +
      scale_color_discrete(name = "Protein State") +
      xlab("Labeling time (sec)") +
      ylab("Average deuteration (Da)") +
      ggtitle(title) +
      theme_bw() +
      expand_limits(y = 0) +
      geom_line(linewidth = input$thickness) +
      geom_point(na.rm = TRUE, size = input$thickness * 1.5) +
      geom_errorbar(
        aes(ymin = deuteration - stddev, ymax = deuteration + stddev),
        width = 0.05,
        linewidth = input$thickness * 0.5,
        color = "black"
      ) +
      theme(
        legend.title = element_text(size = input$fontsize + 2),
        axis.text = element_text(size = input$fontsize + 2),
        axis.title = element_text(size = input$fontsize + 2),
        legend.text = element_text(size = input$fontsize),
        plot.title = element_text(hjust = 0.5, size = input$fontsize + 8)
      )
    
    if (is.na(input$yminvalue) | is.na(input$ymaxvalue)) {
      return(plot)
    } else {
      return(plot + ylim(input$yminvalue, input$ymaxvalue))
    }
  })
  
  output$uptakeplot <- renderPlot({
    singleplot()
  })
  output$uptake_peptide_info <- renderUI({
    req(variables(), input$selpeptide)
    
    wanted <- input$selpeptide
    selected <- suppressWarnings(as.numeric(stringr::str_split(wanted, ":")[[1]][1]))
    
    validate(
      need(!is.na(selected), "No peptide selected.")
    )
    
    HDXdatastates <- variables()$HDXdatastates
    listtestresults <- variables()$listtestresults
    
    rowdata <- HDXdatastates[[1]][HDXdatastates[[1]]$pepnumber == selected, , drop = FALSE]
    
    validate(
      need(nrow(rowdata) > 0, "Selected peptide not found.")
    )
    
    sig_result <- tryCatch({
      listtestresults %>%
        dplyr::filter(pepnumber == selected) %>%
        dplyr::pull(SignificantResults)
    }, error = function(e) {
      NA
    })
    
    sig_label <- "No significant change"
    sig_bg <- "#e8f5e9"
    
    if (length(sig_result) > 0 && !is.na(sig_result[1])) {
      if (as.character(sig_result[1]) == "1") {
        sig_label <- "Increased flexibility"
        sig_bg <- "#fdecea"
      } else if (as.character(sig_result[1]) == "-1") {
        sig_label <- "Increased protection"
        sig_bg <- "#e8f0fe"
      } else if (as.character(sig_result[1]) == "5") {
        sig_label <- "Mixed significance pattern"
        sig_bg <- "#fff8e1"
      }
    }
    
    tags$div(
      class = "instruction-card",
      fluidRow(
        column(
          3,
          tags$div(style = "font-weight:600; color:#666;", "Peptide number"),
          tags$div(style = "font-size:20px;", selected)
        ),
        column(
          3,
          tags$div(style = "font-weight:600; color:#666;", "Residue range"),
          tags$div(style = "font-size:20px;",
                   paste0(rowdata$Start[1], "–", rowdata$End[1]))
        ),
        column(
          4,
          tags$div(style = "font-weight:600; color:#666;", "Sequence"),
          tags$div(
            style = "font-size:18px; word-break: break-word;",
            as.character(rowdata$Sequence[1])
          )
        ),
        column(
          2,
          tags$div(style = "font-weight:600; color:#666;", "Status"),
          tags$div(
            style = paste0(
              "display:inline-block; padding:6px 10px; border-radius:10px; background:",
              sig_bg,
              "; font-weight:600;"
            ),
            sig_label
          )
        )
      )
    )
  })
  output$uptakedownload <- renderUI({
    req(variables())
    actionButton("uptakebutton", "Download all", icon = icon("download"))
  })
  
  observeEvent(input$selpeptide, {
    updateNumericInput(session, "yminvalue", value = "", min = 0, max = 25)
    updateNumericInput(session, "ymaxvalue", value = "", min = 0, max = 25)
  })
  
  output$singledownload <- renderUI({
    req(variables())
    downloadButton("singlebutton", "Download", icon = icon("download"))
  })
  
  x <- 3
  y <- 6
  
  observeEvent(input$uptakebutton, {
    showModal(modalDialog(
      title = "Please select number of rows and columns to export your uptakeplots",
      easyClose = TRUE,
      column(5, numericInput("columns", "Select number of columns:", value = 2, min = 1, max = 6)),
      column(4, numericInput("rows", "Select number of rows:", value = 4, min = 1, max = 6)),
      br(),
      footer = tagList(
        column(2, modalButton("Cancel"), offset = 6),
        downloadButton("downloadall", "Download uptake plots")
      )
    ))
    x <- input$columns
    y <- input$rows
  })
  
  output$singlebutton <- downloadHandler(
    filename = function() {
      paste0("Uptakeplot", input$uptaketypes)
    },
    content = function(file) {
      withProgress(message = "Creating plot", min = 0, max = 1, {
        ggsave(file, plot = singleplot(), width = 8.5)
      })
    }
  )
  
  # ------------------------------------------------------------
  # 3.16.2 Create all uptake plots for export
  # ------------------------------------------------------------
  allplots <- reactive({
    req(variables())
    
    app_status("Preparing all uptake plots...")
    
    pr <- make_progress(
      session = session,
      message = "Creating uptake plot collection...",
      max = 1
    )
    
    on.exit({
      pr$close()
      app_status("Ready")
    }, add = TRUE)
    
    HDXdatastates <- variables()$HDXdatastates
    
    if (requireNamespace("ggthemr", quietly = TRUE)) {
      ggthemr::ggthemr(input$coloring)
    }
    
    uptakeplots <- list()
    end <- if (input$showallstates) states() else 2
    
    for (i in 1:nrow(HDXdatastates[[1]])) {
      deut <- c()
      sd <- c()
      state <- c()
      times <- c()
      
      for (h in 1:end) {
        rowdata <- HDXdatastates[[h]][HDXdatastates[[h]]$pepnumber == i, , drop = FALSE]
        
        deut <- c(deut, as.numeric(rowdata %>% select(matches("Average"))))
        sd <- c(sd, as.numeric(rowdata %>% select(matches("SD"))))
        state <- c(state, rep(as.character(rowdata$State), each = timepoints()))
        times <- c(times, labelingtimepoints())
      }
      
      uptake <- data.frame(Prostate = state, deuteration = deut, stddev = sd, time = times)
      
      title <- paste(
        HDXdatastates[[1]][i, 3], "-",
        HDXdatastates[[1]][i, 4], ":",
        HDXdatastates[[1]][i, 5],
        sep = ""
      )
      
      uptakeplots[[i]] <- ggplot(data = uptake, aes(x = time, y = deuteration, color = Prostate)) +
        geom_point(na.rm = TRUE, size = 1.2) +
        geom_line(linewidth = 0.7) +
        geom_errorbar(
          aes(ymin = deuteration - stddev, ymax = deuteration + stddev),
          width = 0.05, linewidth = 0.1, color = "black"
        ) +
        scale_x_continuous(trans = "log10") +
        scale_color_discrete(name = "Protein State") +
        xlab("Labeling time (sec)") +
        ylab("Average deuteration (Da)") +
        ggtitle(title) +
        theme_bw() +
        expand_limits(y = 0) +
        theme(
          legend.title = element_text(size = 40 / (x + y)),
          axis.text = element_text(size = 64 / (x + y)),
          axis.title = element_text(size = 40 / (x + y / 2)),
          legend.text = element_text(size = 32 / (x + y)),
          plot.title = element_text(hjust = 0.5, size = 70 / (x + y)),
          legend.position = c(0.87, 0.26),
          legend.background = element_rect(fill = "white", color = "gray")
        )
      
      pr$set(
        value = i / nrow(HDXdatastates[[1]]),
        detail = paste("Preparing peptide", i, "of", nrow(HDXdatastates[[1]]))
      )
    }
    
    app_status("All uptake plots prepared")
    
    return(uptakeplots)
  })
  
  
  output$downloadall <- downloadHandler(
    filename = function() {
      paste0("All uptakeplots.pdf")
    },
    content = function(file) {
      withProgress(message = "Exporting", min = 0, max = 1, {
        on.exit(removeModal())
        ggsave(file, marrangeGrob(allplots(), nrow = input$rows, ncol = input$columns), width = 8.5, height = 11)
      })
    }
  )
  
  output$uptaketext <- renderUI({
    req(variables())
    tags$iframe(src = "./uptaketext.html", width = "100%", height = "310px", frameborder = 0, scrolling = "auto", align = "center")
  })
  
  # ------------------------------------------------------------
  # 3.17 3D structure section
  # ------------------------------------------------------------
  
  # ----------------------------
  # 3.17.1 PDB input and structure source selection
  # ----------------------------
  output$pdbselector <- renderUI({
    req(variables())
    radioButtons("pdbfiletype", "Please select how to import your PDB file:",
                 choices = c("Using PDB file", "Using PDB number"),
                 selected = "Using PDB file"
    )
  })
  
  output$pdbtext <- renderUI({
    req(input$pdbfiletype)
    req(variables())
    if (input$pdbfiletype == "Using PDB number") {
      textInput("pdbinputtext", "Type PDB code here:", value = "")
    } else {
      NULL
    }
  })
  
  output$pdbinput <- renderUI({
    req(input$pdbfiletype)
    req(variables())
    if (input$pdbfiletype == "Using PDB file") {
      fileInput("pdbfile", "Select PDB file:", multiple = FALSE, accept = c(".pdb"))
    } else {
      NULL
    }
  })
  
  output$showactionbutton <- renderUI({
    req(variables())
    actionButton("showplot", "Calculate")
  })
  
  pdbstructure <- eventReactive(input$showplot, {
    if (input$pdbfiletype == "Using PDB file") {
      req(input$pdbfile)
      pdbstructure <- input$pdbfile$datapath
    } else if (input$pdbfiletype == "Using PDB number") {
      req(input$pdbinputtext)
      pdbstructure <- input$pdbinputtext
    } else {
      pdbstructure <- NULL
    }
    
    return(pdbstructure)
  })
  
  # ----------------------------
  # 3.17.2 Significant residue mapping for structure coloring
  # ----------------------------
  significantresidues <- reactive({
    req(variables())
    
    sequence <- variables()$sequence
    listtestresults <- variables()$listtestresults
    clusteringresults <- variables()$clusteringresults
    
    listtestresults$Start <- as.numeric(listtestresults$Start)
    listtestresults$End <- as.numeric(listtestresults$End)
    
    if (input$clustercoloring == TRUE) {
      sigresidues <- as.numeric(rep(0, str_length(sequence)))
      
      for (i in 1:nrow(clusteringresults)) {
        for (j in clusteringresults[i, 2]:clusteringresults[i, 3]) {
          if (sigresidues[j] == 0) {
            sigresidues[j] <- clusteringresults[i, (timepoints() * 3) + 5]
          } else if (sigresidues[j] != 0 && clusteringresults[i, (timepoints() * 3) + 5] > sigresidues[j]) {
            sigresidues[j] <- clusteringresults[i, (timepoints() * 3) + 5]
          }
        }
      }
    } else {
      sigresidues <- as.numeric(rep(0, str_length(sequence)))
      
      for (i in 1:nrow(listtestresults)) {
        for (j in listtestresults[i, 2]:listtestresults[i, 3]) {
          if (sigresidues[j] == 0) {
            sigresidues[j] <- listtestresults[i, (timepoints() * 2) + 5]
          } else if (sigresidues[j] == 1 && (listtestresults[i, (timepoints() * 2) + 5] == -1 | listtestresults[i, (timepoints() * 2) + 5] == 5)) {
            sigresidues[j] <- 5
          } else if (sigresidues[j] == -1 && (listtestresults[i, (timepoints() * 2) + 5] == 1 | listtestresults[i, (timepoints() * 2) + 5] == 5)) {
            sigresidues[j] <- 5
          }
        }
      }
    }
    
    if (input$clustercoloring == FALSE) {
      return(list(
        all = paste((input$aaoffset + which(sigresidues == 0)), collapse = " or "),
        a = paste((input$aaoffset + which(sigresidues == 5)), collapse = " or "),
        b = paste((input$aaoffset + which(sigresidues == 1)), collapse = " or "),
        c = paste((input$aaoffset + which(sigresidues == -1)), collapse = " or ")
      ))
    } else {
      return(list(
        all = paste((input$aaoffset + which(sigresidues == 0)), collapse = " or "),
        a = paste((input$aaoffset + which(sigresidues == 1)), collapse = " or "),
        b = paste((input$aaoffset + which(sigresidues == 2)), collapse = " or "),
        c = paste((input$aaoffset + which(sigresidues == 3)), collapse = " or ")
      ))
    }
  })
  
  # ----------------------------
  # 3.17.3 NGL viewer rendering and controls
  # ----------------------------
  output$structure <- renderNGLVieweR({
    req(pdbstructure())
    
    if (input$colorflex == "firebrick2") {
      flex <- "#EE2C2C"
    } else if (input$colorflex == "cyan") {
      flex <- "#00FFFF"
    } else if (input$colorflex == "deeppink") {
      flex <- "#FF1493"
    } else {
      flex <- "#EE2C2C"
    }
    
    if (input$colorprot == "blue") {
      prot <- "#0000FF"
    } else if (input$colorprot == "darkcyan") {
      prot <- "#008B8B"
    } else if (input$colorprot == "chocolate") {
      prot <- "#D2691E"
    } else {
      prot <- "#0000FF"
    }
    
    a <- if (significantresidues()$a == "") "none" else significantresidues()$a
    b <- if (significantresidues()$b == "") "none" else significantresidues()$b
    c <- if (significantresidues()$c == "") "none" else significantresidues()$c
    
    if (input$clustercoloring == FALSE) {
      NGLVieweR(pdbstructure(), width = "50%") %>%
        stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
        addRepresentation(input$structuretype, param = list(sele = significantresidues()$all, color = "gray")) %>%
        addRepresentation(input$structuretype, param = list(sele = a, color = "yellow")) %>%
        addRepresentation(input$structuretype, param = list(sele = b, color = flex)) %>%
        addRepresentation(input$structuretype, param = list(sele = c, color = prot)) %>%
        setQuality("high") %>%
        setSpin(input$spinstructure)
    } else {
      NGLVieweR(pdbstructure(), width = "50%") %>%
        stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
        addRepresentation(input$structuretype, param = list(sele = significantresidues()$all, color = "gray")) %>%
        addRepresentation(input$structuretype, param = list(sele = a, color = "gray")) %>%
        addRepresentation(input$structuretype, param = list(sele = b, color = input$colorint)) %>%
        addRepresentation(input$structuretype, param = list(sele = c, color = input$colorstr)) %>%
        setQuality("high") %>%
        setSpin(input$spinstructure)
    }
  })
  
  output$hrline <- renderUI({
    req(pdbstructure())
    hr(style = "border-top: 1px solid #000000;")
  })
  
  output$representation <- renderUI({
    req(pdbstructure())
    shinyjs::show("3dtoggles")
    selectInput("structuretype", label = "Select Representation", choices = c("ball+stick", "surface", "cartoon"), selected = "cartoon")
  })
  
  output$spinbox <- renderUI({
    req(pdbstructure())
    materialSwitch("spinstructure", label = "Spin", value = TRUE, status = "primary")
  })
  
  observeEvent(input$structuretype, {
    NGLVieweR_proxy("structure", session = session) %>% removeSelection("allwaters")
  })
  
  output$structureoffset <- renderUI({
    req(pdbstructure())
    numericInput("aaoffset", label = "Amino acid offset", value = 0, min = -50, max = 50, step = 1)
  })
  
  observeEvent(input$spinstructure, {
    NGLVieweR_proxy("structure", session = session) %>% updateSpin(input$spinstructure)
  })
  
  output$snapshot <- renderUI({
    req(pdbstructure())
    actionButton("screenshot", "Save screenshot")
  })
  
  # ----------------------------
  # 3.17.4 FASTA/PDB sequence alignment preview
  # ----------------------------
  output$sequenceoffset <- renderUI({
    req(pdbstructure(), variables())
    
    sequence <- as.character(variables()$sequence)
    sequence <- gsub(" ", "", sequence)
    
    short_sequence <- substr(sequence, 1, 50)
    if (nchar(sequence) > 50) {
      short_sequence <- paste0(short_sequence, "...")
    }
    
    tags$div(
      style = "display: flex; align-items: flex-start; font-family: monospace; font-size: 13px;",
      tags$div(
        "FASTA:",
        style = "font-weight: bold; width: 60px; flex-shrink: 0;"
      ),
      tags$div(
        short_sequence,
        style = "white-space: pre-wrap; word-break: break-all;"
      )
    )
  })
  
  output$pdbsequence <- renderUI({
    req(pdbstructure(), variables(), input$structure_sequence, input$aaoffset)
    
    offset <- as.integer(input$aaoffset)
    if (is.na(offset)) offset <- 0
    
    pdb_seq <- paste(as.character(input$structure_sequence), collapse = "")
    pdb_seq <- gsub(" ", "", pdb_seq)
    
    # Shift the displayed PDB sequence according to the amino acid offset
    if (offset >= 0) {
      shifted_pdb_seq <- paste0(strrep(" ", offset), pdb_seq)
    } else {
      shift_left <- abs(offset)
      if (shift_left < nchar(pdb_seq)) {
        shifted_pdb_seq <- substr(pdb_seq, shift_left + 1, nchar(pdb_seq))
      } else {
        shifted_pdb_seq <- ""
      }
    }
    
    short_sequence <- substr(shifted_pdb_seq, 1, 50)
    if (nchar(shifted_pdb_seq) > 50) {
      short_sequence <- paste0(short_sequence, "...")
    }
    
    tags$div(
      style = "display: flex; align-items: flex-start; font-family: monospace; font-size: 13px;",
      tags$div(
        "PDB:",
        style = "font-weight: bold; width: 60px; flex-shrink: 0;"
      ),
      tags$div(
        short_sequence,
        style = "white-space: pre;"
      )
    )
  })
  
  # ----------------------------
  # 3.17.5 Screenshot and PyMOL export helpers
  # ----------------------------
  observeEvent(input$screenshot, {
    NGLVieweR_proxy("structure", session = session) %>%
      setQuality("high") %>%
      snapShot("Structurescreenshot.PNG", param = list(
        antialias = TRUE,
        trim = TRUE,
        transparent = TRUE,
        scale = 1
      ))
  })
  
  output$pymolbutton <- renderUI({
    req(pdbstructure())
    downloadButton("pymolscript", "Download pymol script")
  })
  
  observe({
    if (input$clustercoloring == TRUE) {
      showModal(modalDialog(
        title = "IMPORTANT!",
        "When coloring using clustering results, colors indicate if the peptides have an intermediate or strong effect. Colors DO NOT indicate protection or deprotection",
        footer = NULL,
        easyClose = TRUE
      ))
    }
  })
  
  pymoloutput <- reactive({
    req(pdbstructure())
    req(significantresidues())
    
    if (input$clustercoloring == TRUE) {
      colorstr <- col2rgb(input$colorstr)
      colorstr <- paste0("[", colorstr[1], ",", colorstr[2], ",", colorstr[3], "]")
      colorint <- col2rgb(input$colorint)
      colorint <- paste0("[", colorint[1], ",", colorint[2], ",", colorint[3], "]")
    }
    
    significantresidues <- significantresidues()
    export <- data.frame()
    
    mixed <- gsub(" or ", "+", significantresidues$a)
    deprotected <- gsub(" or ", "+", significantresidues$b)
    protected <- gsub(" or ", "+", significantresidues$c)
    
    export[1, 1] <- "hide everything"
    export[2, 1] <- "show cartoon"
    export[3, 1] <- "color gray,all"
    
    if (input$clustercoloring == TRUE) {
      export[4, 1] <- paste0("set_color strong, ", colorstr)
      export[5, 1] <- paste0("select strongeffects, (i. ", protected, ")")
      export[6, 1] <- "color strong, strongeffects"
      export[7, 1] <- paste0("set_color inter, ", colorint)
      export[8, 1] <- paste0("select intermediateeffects, (i. ", deprotected, ")")
      export[9, 1] <- "color inter, intermediateeffects"
      export[10, 1] <- paste0("select negligibleeffects, (i. ", mixed, ")")
      export[11, 1] <- "color gray, negligibleeffects"
      export[12, 1] <- "deselect"
    } else {
      export[4, 1] <- paste0("select mixedpep, (i. ", mixed, ")")
      export[5, 1] <- "color yellow, mixedpep"
      
      if (input$colorflex == "firebrick2") {
        flexcode <- "[238,44,44]"
      } else if (input$colorflex == "cyan") {
        flexcode <- "[0,255,255]"
      } else if (input$colorflex == "deeppink") {
        flexcode <- "[255,20,147]"
      } else {
        flexcode <- "[238,44,44]"
      }
      
      export[6, 1] <- paste0("set_color deprotected, ", flexcode)
      export[7, 1] <- paste0("select deprotectedpep, (i. ", deprotected, ")")
      export[8, 1] <- "color deprotected, deprotectedpep"
      
      if (input$colorprot == "blue") {
        protcode <- "[0,0,255]"
      } else if (input$colorprot == "darkcyan") {
        protcode <- "[0,139,139]"
      } else if (input$colorprot == "chocolate") {
        protcode <- "[210,105,30]"
      } else {
        protcode <- "[0,0,255]"
      }
      
      export[9, 1] <- paste0("set_color protected, ", protcode)
      export[10, 1] <- paste0("select protectedpep, (i. ", protected, ")")
      export[11, 1] <- "color protected, protectedpep"
      export[12, 1] <- "deselect"
    }
    
    return(export)
  })
  
  output$pymolscript <- downloadHandler(
    filename = function() {
      paste0("Pymol script.pml")
    },
    content = function(file) {
      write.table(pymoloutput(), file, na = "", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  )
  
  # ------------------------------------------------------------
  # 3.18 Export section
  # ------------------------------------------------------------
  output$datadownload <- renderUI({
    req(variables())
    downloadButton("exportbutton", "Export results")
  })
  
  exportdata <- reactive({
    req(variables())
    req(realfit())
    
    app_status("Exporting results workbook...")
    
    pr <- make_progress(
      session = session,
      message = "Exporting data...",
      max = 1
    )
    
    on.exit({
      pr$close()
      app_status("Ready")
    }, add = TRUE)
    
    pepcharge <- initialdata()$pepcharge
    HDXdatastates <- variables()$HDXdatastates
    clusteringresults <- variables()$clusteringresults
    listtestresults <- variables()$listtestresults
    parameters <- variables()$parameters
    threshold <- variables()$threshold
    
    parameters[3, 1] <- paste0("Positive deltaD values (State2-State ref) means increase in flexibility by state 2 compared to the reference (Peptides colored ", input$colorflex, ")")
    parameters[4, 1] <- paste0("Negative deltaD values (State2-State ref) means increase in protection by state 2 (Peptides colored ", input$colorprot, ")")
    
    pr$set(value = 0.20, detail = "Preparing state tables")
    
    for (i in 1:2) {
      HDXdatastates[[i]] <- cbind(
        HDXdatastates[[i]][, 1:5],
        pepcharge,
        HDXdatastates[[i]][, 6:((timepoints() * replicates()) + (timepoints() * 2) + 5)]
      )
    }
    
    wb <<- openxlsx::createWorkbook()
    
    pr$set(value = 0.40, detail = "Writing worksheets")
    
    openxlsx::addWorksheet(wb, paste("State 1", HDXdatastates[[1]][1, 1], "(Ref)"))
    openxlsx::writeData(wb, paste("State 1", HDXdatastates[[1]][1, 1], "(Ref)"), HDXdatastates[[1]])
    
    openxlsx::addWorksheet(wb, paste("State 2", as.character(HDXdatastates[[2]][1, 1])))
    openxlsx::writeData(wb, paste("State 2", as.character(HDXdatastates[[2]][1, 1])), HDXdatastates[[2]])
    
    openxlsx::addWorksheet(wb, "Statistical parameters")
    openxlsx::writeData(wb, "Statistical parameters", parameters)
    
    openxlsx::addWorksheet(wb, "Statistics 2 vs ref")
    listtestresults$SignificantResults[listtestresults$SignificantResults == 5] <- "check"
    openxlsx::writeData(wb, "Statistics 2 vs ref", listtestresults)
    
    openxlsx::addWorksheet(wb, paste("Clustering", 2, "vs ref"))
    openxlsx::writeData(wb, paste("Clustering", 2, "vs ref"), clusteringresults)
    
    style1 <- openxlsx::createStyle(fontColour = input$colorflex)
    style2 <- openxlsx::createStyle(fontColour = input$colorprot)
    style3 <- openxlsx::createStyle(fontColour = "yellow3")
    style4 <- openxlsx::createStyle(fontColour = input$colorprot, textDecoration = "bold")
    style5 <- openxlsx::createStyle(fontColour = input$colorflex, textDecoration = "bold")
    
    pr$set(value = 0.60, detail = "Writing fitting outputs")
    
    if (realfit() == 2) {
      peptide <- c()
      for (i in 1:length(fitting()$kinetics)) {
        peptide <- append(peptide, fitting()$kinetics[[i]]$data$rowname[1])
      }
      
      openxlsx::addWorksheet(wb, "Fitting parameters")
      openxlsx::writeData(wb, "Fitting parameters", cbind(fitting()$fittingparameters, peptide))
      
      table <- data.frame()
      summary <- fitting()$summary
      
      for (k in 1:length(fitting()$kinetics)) {
        null <- data.frame(coef(summary[[k]]$null), check.names = FALSE)
        null$variable <- substr(rownames(null), 1, 1)
        null <- null %>% select(variable, everything())
        null[, 6] <- "Null"
        colnames(null)[6] <- "State"
        
        for (l in 1:states()) {
          df <- data.frame(coef(summary[[k]][[l + 1]]), check.names = FALSE)
          df$variable <- substr(rownames(df), 1, 1)
          df <- df %>% select(variable, everything())
          df[, 6] <- paste0("State", l)
          colnames(df)[6] <- "State"
          null <- rbind.fill(null, df)
        }
        
        null <- null %>% unite("variable", c(variable, State), remove = TRUE)
        null <- pivot_wider(null, names_from = variable, values_from = c(Estimate, "Std. Error", "t value", "Pr(>|t|)"))
        table <- rbind.fill(table, null)
      }
      
      openxlsx::addWorksheet(wb, "Calculated parameters")
      openxlsx::writeData(wb, "Calculated parameters", cbind(peptide, table))
    }
    
    pr$set(value = 0.80, detail = "Applying workbook styling")
    
    openxlsx::addStyle(wb, sheet = "Statistical parameters", style = style1, rows = 4, cols = 1, gridExpand = TRUE)
    openxlsx::addStyle(wb, sheet = "Statistical parameters", style = style2, rows = 5, cols = 1, gridExpand = TRUE)
    openxlsx::addStyle(wb, sheet = "Statistical parameters", style = style3, rows = 6, cols = 1, gridExpand = TRUE)
    
    for (j in 1:nrow(listtestresults)) {
      if (listtestresults$SignificantResults[j] == "-1") {
        openxlsx::addStyle(wb, sheet = paste("Statistics", 2, "vs ref"), style = style2, rows = j + 1, cols = 1:ncol(listtestresults), gridExpand = TRUE)
      } else if (listtestresults$SignificantResults[j] == "1") {
        openxlsx::addStyle(wb, sheet = paste("Statistics", 2, "vs ref"), style = style1, rows = j + 1, cols = 1:ncol(listtestresults), gridExpand = TRUE)
      } else if (listtestresults$SignificantResults[j] == "check") {
        openxlsx::addStyle(wb, sheet = paste("Statistics", 2, "vs ref"), style = style3, rows = j + 1, cols = 1:ncol(listtestresults), gridExpand = TRUE)
      }
    }
    
    for (i in 1:nrow(listtestresults)) {
      for (j in 1:timepoints()) {
        if (!is.na(listtestresults[i, (2 * j) + 3]) && !is.na(listtestresults[i, (2 * j) + 4])) {
          if (listtestresults[i, (2 * j) + 3] >= threshold && listtestresults[i, (2 * j) + 4] < significancelevel()) {
            openxlsx::addStyle(wb, sheet = paste("Statistics", 2, "vs ref"), style = style5, rows = i + 1, cols = (((2 * j) + 3):((2 * j) + 4)))
          } else if (listtestresults[i, (2 * j) + 3] <= -threshold && listtestresults[i, (2 * j) + 4] < significancelevel()) {
            openxlsx::addStyle(wb, sheet = paste("Statistics", 2, "vs ref"), style = style4, rows = i + 1, cols = (((2 * j) + 3):((2 * j) + 4)))
          }
        }
      }
    }
    
    pr$set(value = 1, detail = "Done")
    app_status("Workbook export prepared")
    
    return(wb)
  })
  
  output$tableresults <- renderUI({
    req(variables())
    DT::DTOutput("outputresults")
  })
  
  output$outputresults <- DT::renderDT({
    req(variables())
    
    listtestresults <- variables()$listtestresults
    listtestresults$SignificantResults[listtestresults$SignificantResults == 5] <- "check"
    
    tbl <- round_df(listtestresults, digits = 3)
    
    DT::datatable(
      tbl,
      rownames = FALSE,
      filter = "top",
      selection = list(mode = "single", target = "row"),
      extensions = c("Buttons", "Scroller"),
      options = list(
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel"),
        pageLength = 10,
        lengthMenu = c(10, 25, 50, 100),
        scrollX = TRUE,
        deferRender = TRUE,
        scrollY = 400,
        scroller = TRUE,
        autoWidth = TRUE
      )
    )
  }, server = FALSE)
  
  output$exportbutton <- downloadHandler(
    filename = function() {
      paste0("HDX Analysis results.xlsx")
    },
    content = function(file) {
      openxlsx::saveWorkbook(exportdata(), file)
    }
  )
  
  # ------------------------------------------------------------
  # 3.18.1 Export all plots
  # ------------------------------------------------------------
  output$allplotsdownload <- renderUI({
    req(variables())
    actionButton("exportallaction", "Export all plots", icon = icon("download"))
  })
  
  document2 <- reactive({
    volcano <- volcanoplot()
    histogram_plot <- NULL
    if (input$clusteringsel == TRUE) histogram_plot <- histogram()
    woodsplots_obj <- woodsplots()
    woods <- woodsbytimepoint()
    digestion_obj <- digestion()
    peptide <- peptidemap()
    
    return(list(
      volcano = volcano,
      histogram = histogram_plot,
      woodsplots = woodsplots_obj,
      digestion = digestion_obj,
      peptide = peptide,
      woods = woods
    ))
  })
  
  observeEvent(input$exportallaction, {
    showModal(modalDialog(
      title = "Please wait. Creating your document...",
      easyClose = FALSE,
      modalButton("Cancel"),
      downloadButton("exportallplotsbutton", "Export all plots"),
      footer = NULL
    ))
    
    shinyjs::disable(id = "exportallplotsbutton")
    if (realfit() == 2) a <- kineticdocument()
    b <- allplots()
    c <- document2()
    shinyjs::enable(id = "exportallplotsbutton")
  })
  
  output$exportallplotsbutton <- downloadHandler(
    filename = "All plots.pdf",
    content = function(file) {
      app_status("Exporting all plots...")
      
      pr <- make_progress(
        session = session,
        message = "Exporting all plots...",
        max = 1,
        button_ids = c("exportallplotsbutton")
      )
      
      on.exit({
        pr$close()
        app_status("Ready")
      }, add = TRUE)
      
      pdf(paste0(file, ".7x7"), width = 7, height = 7)
      print(document2()$volcano)
      pr$set(value = 1 / 7, detail = "Exported volcano plot")
      
      if (input$clusteringsel == TRUE && !is.null(document2()$histogram)) {
        print(document2()$histogram)
      }
      
      print(document2()$woodsplots)
      print(marrangeGrob(document2()$woods, nrow = timepoints(), ncol = 1, top = "Woods-plot by timepoint"))
      pr$set(value = 2 / 7, detail = "Exported Woods plots")
      
      print(document2()$peptide$summary)
      pr$set(value = 3 / 7, detail = "Exported peptide map")
      
      print(digestion()$plot)
      dev.off()
      
      pdf(paste0(file, ".8x11"), width = 8, height = 11)
      pr$set(value = 4 / 7, detail = "Exporting uptake plots")
      print(marrangeGrob(allplots(), nrow = 4, ncol = 2, top = NULL))
      dev.off()
      
      if (realfit() == 2) {
        pdf(paste0(file, ".14x8"), width = 14, height = 8)
        print(kineticdocument())
        dev.off()
      }
      
      pr$set(value = 5 / 7, detail = "Combining PDF outputs")
      
      qpdf::pdf_subset(
        paste0(file, ".8x11"),
        pages = 2:(qpdf::pdf_length(paste0(file, ".8x11"))),
        output = paste0(file, ".8")
      )
      
      if (realfit() == 1) {
        qpdf::pdf_combine(paste0(file, c(".7x7", ".8")), output = file)
        file.remove(paste0(file, c(".7x7", ".8")))
      } else {
        qpdf::pdf_combine(paste0(file, c(".7x7", ".8", ".14x8")), output = file)
        file.remove(paste0(file, c(".7x7", ".8", ".14x8")))
      }
      
      pr$set(value = 1, detail = "Done")
      app_status("All plots exported")
    }
  )
  
  output$visitor_banner <- renderUI({
    visitors <- visitors_data()
    
    if (nrow(visitors) == 0) {
      return(tags$div(class = "visitor-banner"))
    }
    
    visitors <- visitors[order(-visitors$visits), ]
    
    tags$div(
      class = "visitor-banner",
      lapply(seq_len(nrow(visitors)), function(i) {
        tags$span(
          class = "visitor-item",
          tags$img(
            src = flag_img_url(visitors$country_code[i]),
            height = "16px",
            style = "border:1px solid #ddd; border-radius:2px;"
          ),
          tags$span(
            style = "font-size: 12px; font-weight: 600;",
            visitors$visits[i]
          )
        )
      })
    )
  })
  observeEvent(input$go_input, {
    updateTabItems(session, "tabs", selected = "input_dash")
  })
  
  observeEvent(input$go_faq, {
    updateTabItems(session, "tabs", selected = "faq_dash")
  })
  
  observe({
    req(input$intro_doc)
    
    shinyjs::runjs(sprintf(
      "document.getElementById('intro_pdf_viewer').src = '%s';",
      input$intro_doc
    ))
  })
  
  output$download_citation_endnote <- downloadHandler(
    filename = function() {
      "Kingfisher_HDX_MS_Citation.ris"
    },
    content = function(file) {
      ris_lines <- c(
        "TY  - JOUR",
        "T1  - Kingfisher: An open-sourced web-based platform for the analysis of hydrogen exchange mass spectrometry data",
        "AU  - McLaughlin, Nolan K.",
        "AU  - Rincon Pabon, Juan P.",
        "AU  - Gies, Samantha",
        "AU  - Dastvan, Reza",
        "AU  - Gross, Michael L.",
        "Y1  - 2025/04/01",
        "PY  - 2025",
        "DA  - 2025/04/01",
        "DO  - 10.1002/pro.70096",
        "T2  - Protein Science",
        "JF  - Protein Science",
        "JO  - Protein Science",
        "JA  - Protein Science",
        "SP  - e70096",
        "VL  - 34",
        "IS  - 4",
        "KW  - HDX statistics",
        "KW  - hydrogen–deuterium exchange mass spectrometry",
        "KW  - open-source software",
        "KW  - protein",
        "KW  - R programming",
        "KW  - structural proteomics",
        "PB  - John Wiley & Sons, Ltd",
        "SN  - 0961-8368",
        "UR  - https://doi.org/10.1002/pro.70096",
        "Y2  - 2026/06/22",
        paste0(
          "N2  - Abstract Hydrogen-deuterium exchange mass spectrometry (HDX-MS) is now a critical tool in molecular biology and structural proteomics. ",
          "It is routinely used to probe protein and conformational dynamics through a well-established experiment where amide hydrogens exchange with deuterium atoms in a buffer containing D2O. ",
          "Although there have been numerous advances in the field, data analysis still poses challenges mainly due to the need for manual curation of the data and the lack of standardized statistics and accessible software. ",
          "In response, we developed Kingfisher, an open-source, user-friendly, web-based solution that facilitates downstream analysis using well-established statistics and provides advanced high-resolution representations of the HDX results. ",
          "Kingfisher is able to read data directly as exported from common software packages and usually takes less than a minute to run the analysis, without the need to download the raw code or install any software. ",
          "We foresee Kingfisher as a valuable tool for both newcomers and experts in the field of Hydrogen Exchange Mass Spectrometry. ",
          "Kingfisher is available to all users as an interactive web application at https://kingfisher.wustl.edu/."
        ),
        "ER  -"
      )
      writeLines(ris_lines, file, useBytes = TRUE)
    }
  )
  
  observeEvent(input$show_citation, {
    showModal(
      modalDialog(
        title = "How to cite Kingfisher HDX-MS",
        easyClose = TRUE,
        size = "l",
        footer = tagList(
          downloadButton("download_citation_endnote", "Download for EndNote"),
          modalButton("Close")
        ),
        
        tags$p(
          "If you use Kingfisher HDX-MS in your work, please cite the manuscript below:"
        ),
        
        tags$div(
          class = "instruction-card",
          tags$p(tags$b("Plain-text citation")),
          tags$pre(
            "McLaughlin NK, Rincon Pabon JP, Gies S, Dastvan R, Gross ML. Kingfisher: An open-sourced web-based platform for the analysis of hydrogen exchange mass spectrometry data. Protein Science. 2025;34(4):e70096. https://doi.org/10.1002/pro.70096"
          )
        ),
        
        tags$div(
          class = "instruction-card",
          tags$p(tags$b("DOI")),
          tags$a(
            href = "https://doi.org/10.1002/pro.70096",
            target = "_blank",
            "https://doi.org/10.1002/pro.70096"
          )
        ),
        
        tags$div(
          class = "instruction-card",
          tags$p(tags$b("BibTeX")),
          tags$pre(
            "@article{https://doi.org/10.1002/pro.70096,
author = {McLaughlin, Nolan K. and Rincon Pabon, Juan P. and Gies, Samantha and Dastvan, Reza and Gross, Michael L.},
title = {Kingfisher: An open-sourced web-based platform for the analysis of hydrogen exchange mass spectrometry data},
journal = {Protein Science},
volume = {34},
number = {4},
pages = {e70096},
keywords = {HDX statistics, hydrogen–deuterium exchange mass spectrometry, open-source software, protein, R programming, structural proteomics},
doi = {https://doi.org/10.1002/pro.70096},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/pro.70096},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/pro.70096},
abstract = {Abstract Hydrogen–deuterium exchange mass spectrometry (HDX-MS) is now a critical tool in molecular biology and structural proteomics. It is routinely used to probe protein and conformational dynamics through a well-established experiment where amide hydrogens exchange with deuterium atoms in a buffer containing D2O. Although there have been numerous advances in the field, data analysis still poses challenges mainly due to the need for manual curation of the data and the lack of standardized statistics and accessible software. In response, we developed Kingfisher, an open-source, user-friendly, web-based solution that facilitates downstream analysis using well-established statistics and provides advanced high-resolution representations of the HDX results. Kingfisher is able to read data directly as exported from common software packages and usually takes less than a minute to run the analysis, without the need to download the raw code or install any software. We foresee Kingfisher as a valuable tool for both newcomers and experts in the field of Hydrogen Exchange Mass Spectrometry. Kingfisher is available to all users as an interactive web application at https://kingfisher.wustl.edu/.},
year = {2025}
}

"
          )
        )
      )
    )
  })
  
  
}