library(shiny)
library(shinyWidgets)
library(dplyr)
library(tidyverse)
library(gsubfn)
library(zip)

source('nanostring-normalizer.R')

ui <- fluidPage(
  
  titlePanel('Nanostring Normalization Tool'), 
  
  sidebarLayout(
    
    sidebarPanel(
      plotOutput('Ratatouille', width = '100%', height = '100%'),
      helpText(tags$p(
        'This tool is intended to be a replacement for nSolver in analyzing PlexSets.',
        'It is intended to be used directly with the raw RCC files.',
        'Please contact its owner, Aaron Yu, if there are any issues or questions.'
      )),
      hr(),
      textInput('exptname', 'Experiment/Export File Name', placeholder = 'Enter Experiment Name...'),
      textAreaInput('hkgs', 'Housekeeping Genes', 
                    placeholder = paste('Enter each HKG on a new line.', 
                    'HKGs are case-insensitive.', 
                    'All probes containing any of these names will be included.'), 
                    height = '150px', 
                    resize = 'vertical'),
      div(style = 'display:block',
          actionButton('submitbutton', 'Analyze', class = 'btn btn-primary', width = '78%'),
          actionButton('advancedSettings', '', icon = icon("cog", lib = "glyphicon"), width = '20%')
      )
    ), 
    
    mainPanel(
      tabsetPanel(id = 'tabset',
        tabPanel('File Uploads', 
                 helpText('Please upload the raw RCC files by clicking the buttons below.'),
                 br(),
                 fileInput(
                   'uploadCalibrator', 
                   'Calibrator RCC File:', 
                   accept = c('.RCC'),
                   multiple = FALSE, 
                   buttonLabel = 'Upload File'
                 ), 
                 tableOutput('fileTable_calibrator'), 
                 # uiOutput('platefiles'), #not using separate plates anymore
                 fileInput(
                   'uploadData', 
                   'Data RCC Files:', 
                   accept = c('.RCC'),
                   multiple = TRUE, 
                   buttonLabel = 'Upload Files'
                 ), 
                 tableOutput('fileTable_data')
        ),
        tabPanel('Platemaps',
                 helpText(
                   tags$p(
                     'Please upload a ', tags$strong('.csv'), 'file containing the platemap using 
                     the templates and convention shown. Do not include the information for the 
                     calibrator file in the platemap file.'
                   )
                 ), 
                 div(
                   style = 'display:block', 
                   downloadButton('downloadTemplate', 'Download Template'),
                   downloadButton('downloadExample', 'Download Example')
                 ),
                 br(),
                 fileInput(
                   'uploadPlatemap',
                   'Platemap CSV File:',
                   accept = c('.csv'),
                   multiple = FALSE,
                   buttonLabel = 'Upload File'
                 ), 
                 hr(), 
                 helpText(tags$strong('Example Platemap Information')),
                 tableOutput('platemap_example')
        ), 
        tabPanel('Output',
                 helpText(tags$strong('Final Output Table Below')),
                 # div(
                 #   style = 'display:block', 
                 #   downloadButton('downloadzip', 'Download Data'),
                 #   downloadButton('downloadData', 'Download Data'),
                 #   downloadButton('downloadInfo', 'Download Info')
                 # ),
                 downloadButton('downloadzip', 'Download Data'),
                 hr(),
                 tableOutput('dataTable'))
      )
    )
    
  )
  
)


server <- function(input, output, session) {
  
  # Initialization
  downloadable <- FALSE
  output.data.raw <- data.frame()
  output.data.normalized <- data.frame()
  output.data.normtype <- ''
  output.info <- list()
  use_highest_of_neg_ctrls <- TRUE
  threshold_value <- 20
  separate_tissues <- TRUE
  tech_norm_geomean <- 8000
  content_norm_geomean <- 2000
  use_technormstatic <- TRUE
  use_contentnormstatic <- TRUE
  
  # Amazing Image
  output$Ratatouille <- renderImage({
    list(
      src = 'quote-edited.jpg',
      contentType = 'image/jpeg',
      height = session$clientData$output_Ratatouille_width,
      width = session$clientData$output_Ratatouille_width
    )
  }, deleteFile = FALSE)
  
  # Initial defining functions
  getRCCfiles <- function() {
    validate(
      need(
        input$uploadCalibrator,
        message = 'Please upload a valid calibrator file!'
      )
    )
    files.calibrator <- data.frame(rcc.files = input$uploadCalibrator$name, 
                        tempname = '0.RCC')
    files.dirs.calibrator <- dirname(input$uploadCalibrator$datapath)
    validate(
      need(
        input$uploadCalibrator,
        message = 'Please upload a valid calibrator file!'
      )
    )
    files <- data.frame(
      rcc.files = input$uploadData$name, 
      tempname = paste(
        1:length(input$uploadData$name) - 1,
        rep('RCC', times = length(input$uploadData$name)), 
        sep = '.'
      )
    )
    files.dirs <- dirname(input$uploadData$datapath)
    files <- list(files.calibrator, files)
    files.dirs <- list(files.dirs.calibrator, files.dirs)
    return(list(files, files.dirs))
  }
  
  convertRCCs <- function(rccfiledir, rccfilenames) {
    out <- NACHO::load_rcc(data_directory = rccfiledir, ssheet_csv = rccfilenames, id_colname = 'tempname')
    return(get_rawCounts_from_NACHO(out$nacho))
  }
  
  getwellnamefrominfo <- function(row, column, cdf.filename, column.name) {
    col <- str_pad(gsub('\\s', '', column), 2, pad = '0')
    paste0(
      paste(
        'Set',
        row, 
        paste0('(', row, col, ')'),
        paste(cdf.filename, column.name, col, sep = '_')
      ),
      '.RCC'
    )
  }
  
  platemaps.convert <- function() {
    validate(
      need(
        input$uploadPlatemap,
        message = 'Please upload a valid platemap file!'
      )
    )
    out <- read.csv(input$uploadPlatemap$datapath, check.names = FALSE) %>% as.data.frame(.)
    rownames(out) <- apply(out, 1, function(info) {
      getwellnamefrominfo(info['Row'], info['Column'], info['CDF_Filename'], info['Column_Name'])
    })
    return(out)
  }
  
  parseHKGs <- function(hkgtext) {
    hkgtext.processed <- str_split(hkgtext, pattern = '(\\n|\\s)') %>% unlist(.)
    if(hkgtext.processed == '' 
       || is.null(hkgtext.processed) 
       || is.na(hkgtext.processed)) {
      return(NULL)
    } else {
      return(hkgtext.processed[hkgtext.processed != '' 
                               && !is.null(hkgtext.processed) 
                               && !is.na(hkgtext.processed)])
    }
  }
  
  constructcalibratorplatemap <- function(calibrator.names) {
    samples <- length(calibrator.names)
    Plate <- rep('0', times = samples)
    Tissue <- rep('Calibrator', times = samples)
    Animal <- rep('0', times = samples)
    Group <- rep('Calibrator', times = samples)
    Time <- rep('0', times = samples)
    Well <- sapply(calibrator.names, function(i) {
      start <- str_locate(i, '\\(')[1] + 1
      end <- str_locate(i, '\\)')[1] - 1
      substr(i, start, end)
    })
    Row <- sapply(Well, function(i) substr(i, 1, 1))
    Column <- sapply(Well, function(i) as.integer(substr(i, 2, 3)))
    Column_Name <- rep('0', times = samples)
    CDF_Filename <- rep('0', times = samples)
    RNA_loaded <- rep('0', times = samples)
    out <- data.frame(Plate, CDF_Filename, Row, Column, Column_Name, Group, Animal, Tissue, Time, RNA_loaded)
    rownames(out) <- calibrator.names
    return(out)
  }
  
  settingsModal <- function() {
    modalDialog(
      title = 'Advanced Settings (for Nanostring experts ONLY)',
      footer = tagList(
        actionButton('settingsOK', 'Save', class = 'btn btn-primary'),
        modalButton('Cancel')
      ),
      helpText('Default settings are all checked.'),
      checkboxInput('separatetissues', 'Normalize by tissue type?', value = separate_tissues),
      div(
        style = 'display:flex',
        checkboxInput('negCtrls', 'Threshold based on negative controls?', value = use_highest_of_neg_ctrls),
        uiOutput('numericthreshold')
      ),
      div(
        style = 'display:flex',
        checkboxInput('technormmeanstatic', 'Set all positive control geomeans to a constant value? (Recommended)', 
                      value = use_technormstatic),
        uiOutput('technormvalue')
      ),
      div(
        style = 'display:flex',
        checkboxInput('contentnormmeanstatic', 'Set all HKG geomeans to a constant value? (Recommended)', 
                      value = use_contentnormstatic),
        uiOutput('contentnormvalue')
      )
    )
  }
  
  # Outputs
  output$fileTable_calibrator <- renderTable({
    req(input$uploadCalibrator)
    input$uploadCalibrator$name
  })
  
  output$fileTable_data <- renderTable({
    req(input$uploadData)
    input$uploadData$name
  })
  
  output$platemap_example <- renderTable({
    Plate <- c('Plate # for organization', '1', '2', '1')
    Tissue <- c('Tissue/Organ - use consistent naming', 'Gastroc', 'Gastroc', 'VL')
    Animal <- c('#', '1', '3', '2')
    Group <- c('#', '1001', '3', '15')
    Time <- c('days/hrs', '24', '24', '48')
    Row <- c('Letter', 'A', 'D', 'G')
    Column <- c('#', '1', '11', '7')
    Column_Name <- c('Column name from .cdf file', 'Col1', 'Gas4', '7')
    CDF_Filename <- c('Name of .cdf file used to generate data (excluding \'.cdf\')', '20211106_2021-10-26-CNM', 
                      '20220406_2022-04-05-FOP-titration', '20220708_2022-07-07-704-DM2-Pset2')
    RNA_loaded <- c('ng (optional)', '105.8', '82.67', '201.2')
    data.frame(Plate, CDF_Filename, Row, Column, Column_Name, Group, Animal, Tissue, Time, RNA_loaded)
  })
  
  output$numericthreshold <- renderUI({
    if(!input$negCtrls) {
      numericInput('threshold', 'Threshold Value', min = 1, value = threshold_value)
    }
  })
  
  output$technormvalue <- renderUI({
    if(input$technormmeanstatic) {
      numericInput('posCtrlsgeomean', 'Value', min = 1, value = tech_norm_geomean)
    } else {
      helpText('Using average of geomeans instead (not recommended)...')
    }
  })
  
  output$contentnormvalue <- renderUI({
    if(input$contentnormmeanstatic) {
      numericInput('hkgsgeomean', 'Value', min = 1, value = content_norm_geomean)
    } else {
      helpText('Using average of geomeans instead (not recommended)...')
    }
  })
  
  # Action reactive events
  observeEvent(input$advancedSettings, {
    showModal(settingsModal())
  })
  
  observeEvent(input$settingsOK, {
    use_highest_of_neg_ctrls <<- input$negCtrls
    threshold_value <<- ifelse(use_highest_of_neg_ctrls, threshold_value, input$threshold)
    use_technormstatic <<- input$technormmeanstatic
    if(use_technormstatic) {
      tech_norm_geomean <<- input$posCtrlsgeomean
    }
    use_contentnormstatic <<- input$contentnormmeanstatic
    if(use_contentnormstatic) {
      content_norm_geomean <<- input$hkgsgeomean
    }
    removeModal()
  })
  
  observeEvent(input$submitbutton, {
    withProgress(message = 'Calculating normalization:', value = 0, {
      n <- 7
      incProgress(1/n, detail = 'converting platemaps...')
      platemap.all <- platemaps.convert()
      incProgress(1/n, detail = 'converting RCC files...')
      list[rcc.files, rcc.files.dirs] <- getRCCfiles()
      counts.nacho <- lapply(1:length(rcc.files.dirs), function(i) {
        convertRCCs(rccfiledir = rcc.files.dirs[[i]],
                    rccfilenames = rcc.files[[i]])
      })
      genelist <- names(counts.nacho[[1]])
      start_gene <- genelist[5] ### based on NACHo parser function that adds additional information columns in the beginning
      validate(
        need(
          all(sort(Reduce(intersect, lapply(counts.nacho, names))) == sort(names(counts.nacho[[1]]))), 
          message = tags$p('The RCC files must created using the same RCF! 
                         Please ensure the correct RCC files were selected.')
        )
      )
      incProgress(1/n, detail = 'combining raw counts')
      ### Remove calibrator files from main files
      counts.nacho[[2]] <- counts.nacho[[2]][!(rownames(counts.nacho[[2]]) %in% rownames(counts.nacho[[1]])), ]
      ###
      counts.raw.only.all <- Reduce(rbind, counts.nacho)
      platemap.calibrator <- constructcalibratorplatemap(rownames(counts.nacho[[1]]))
      platemap.all <- rbind(platemap.calibrator, platemap.all)
      samples <- intersect(rownames(counts.raw.only.all), rownames(platemap.all))
      validate(
        need(
          !is.null(samples),
          message = tags$p('The platemap does not seem to contain the same samples as the RCC files!
                         Please ensure that the platemap was entered correctly.')
        )
      )
      counts.raw.only.subset <- counts.raw.only.all[rownames(counts.raw.only.all) %in% samples, ]
      platemap.subset <- platemap.all[rownames(platemap.all) %in% samples, ]
      counts.raw <- cbind(platemap.subset, counts.raw.only.subset)
      output.data.raw <<- counts.raw
      incProgress(1/n, detail = 'technial normalization')
      tech_norm_geomean <<- ifelse(use_technormstatic, tech_norm_geomean, 0)
      counts.technorm <- technicalNormalize(counts.raw, 
                                            calibrator.sample.names = rownames(platemap.calibrator),
                                            read_startCol_gene = start_gene,
                                            highest_of_neg_ctrls = use_highest_of_neg_ctrls,
                                            threshold_value = threshold_value,
                                            normalize_value = tech_norm_geomean)
      incProgress(1/n, detail = 'content normalization')
      content_norm_geomean <<- ifelse(use_contentnormstatic, content_norm_geomean, 0)
      HKGs.input <- parseHKGs(input$hkgs)
      HKGs.grep <- unique(
        grep(
          paste(HKGs.input, collapse = '|'), 
          genelist, 
          ignore.case = TRUE, 
          value = TRUE
        )
      )
      output.data.normalized <<- counts.technorm
      output.data.normtype <<- 'TechNormCounts'
      information <- list(
        c('HKGs', 'NULL'),
        c('Threshold', ifelse(use_highest_of_neg_ctrls, 'Highest of negative controls', threshold_value)),
        c('Normalize by tissue type', separate_tissues),
        c('Positive controls average geomean value', tech_norm_geomean),
        c('HKGs average geomean value', 'NULL')
      )
      if(!is_empty(HKGs.grep)) {
        counts.norm <- contentNormalize(counts.technorm,
                                        tissue_colName = 'Tissue',
                                        HK_genes = HKGs.grep,
                                        TOIs = 'All',
                                        read_startCol_gene = start_gene,
                                        highest_neg_ctrl = use_highest_of_neg_ctrls,
                                        threshold_value = threshold_value,
                                        separate_tissues = separate_tissues,
                                        normalize_value = content_norm_geomean)
        information <<- list(
          c('HKGs', HKGs.grep),
          c('Threshold', ifelse(use_highest_of_neg_ctrls, 'Highest of negative controls', threshold_value)),
          c('Normalize by tissue type', separate_tissues),
          c('Positive controls average geomean value', tech_norm_geomean),
          c('HKGs average geomean value', content_norm_geomean)
        )
        output.data.normalized <<- counts.norm
        output.data.normtype <<- 'NormalizedCounts'
      }
      output.info <<- paste(paste(information,
                                  sep = '\t'),
                            collapse = '\n')
      incProgress(1/n, detail = 'outputing data')
      output$dataTable <- renderTable(output.data.normalized, rownames = TRUE)
      downloadable <<- TRUE
      incProgress(1/n, detail = 'done!')
      updateTabsetPanel(session, 'tabset', selected = 'Output')
    })
  })
  
  output$downloadTemplate <- downloadHandler(
    filename = 'platemaps_template.csv',
    content = function(file) {
      file.copy('platemaps_template.csv', file)
    },
    contentType = 'text/csv'
  )
  
  output$downloadExample <- downloadHandler(
    filename = 'platemaps_example.csv',
    content = function(file) {
      file.copy('platemaps_example.csv', file)
    },
    contentType = 'text/csv'
  )
  
  # output$downloadData <- downloadHandler(
  #   filename = function() {
  #     req(input$exptname)
  #     paste0(input$exptname, '_', output.data.normtype, '.csv')
  #   },
  #   content = function(file) {
  #     write.csv(output.data, file)
  #   },
  #   contentType = 'text/csv'
  # )
  # 
  # output$downloadInfo <- downloadHandler(
  #   filename = function() {
  #     req(input$exptname)
  #     paste0(input$exptname, '_', 'Info', '.txt')
  #   },
  #   content = function(file) {
  #     writeLines(output.info, file)
  #   },
  #   contentType = 'text/plain'
  # )
  
  output$downloadzip <- downloadHandler(
    filename = function() {
      paste0(input$exptname, '_', 'NanostringData', '.zip')
    },
    content = function(filename) {
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      files <- c(
        paste0(input$exptname, '_', 'RawCounts', '.csv'),
        paste0(input$exptname, '_', output.data.normtype, '.csv'),
        paste0(input$exptname, '_', 'info', '.txt')
      )
      files.data <- list(
        output.data.raw,
        output.data.normalized,
        output.info
      )
      # lapply(1:length(files), function(i) {
      #   write(files.data[[i]], files[i])
      # })
      write.csv(files.data[[1]], files[1])
      write.csv(files.data[[2]], files[2])
      writeLines(files.data[[3]], files[3])
      zip(zipfile = filename, files = files)
    },
    contentType = 'application/zip'
  )
  
}

shinyApp(ui, server)