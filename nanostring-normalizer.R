suppressWarnings({
  library(NACHO)
  library(tidyr)
  library(tidyverse)
  library(dplyr)
  library(stringr)
  library(psych)
})

options(scipen = 999)

##function to identify first gene column index
getStartColumn <- function(data, gene) {
  read_startCol <- unique(grep(gene, colnames(as.data.frame(data)), ignore.case = TRUE, value = FALSE))
  return(read_startCol)
}

ctrls.neg <- c("NEG_1","NEG_2","NEG_3","NEG_4","NEG_5","NEG_6","NEG_7","NEG_8")
ctrls.pos <- c("POS_1","POS_2","POS_3","POS_4","POS_5","POS_6","POS_7","POS_8")
ctrls.all <- c(ctrls.neg, ctrls.pos)

#thresholding function
threshold <- function(df, cutoff = 20, kill = FALSE) {
  df[df <= cutoff | df == Inf | df == -Inf | is.na(df)] <- ifelse(kill, 0, cutoff)
  return(df)
}

build_calibrator_matrix <- function(pos_norm_calibrator_column, default_row = 'A') {
  # grab the default row that will be used for calibration
  # calibrator should be in the format with rownames(calibrator) = sample_lanes and colnames(calibrator) = gene_names
  # calibrator should not contain any control genes
  ##pattern <- paste0('\\(', default_row)
  ##basis <- pos_norm_calibrator_column[grepl(pattern, rownames(pos_norm_calibrator_column)), ]
  basis <- pos_norm_calibrator_column[rownames(pos_norm_calibrator_column) == default_row, ]
  scale <- apply(pos_norm_calibrator_column, 1, function(row) basis / row) %>% Reduce(rbind, .)
  rownames(scale) <- rownames(pos_norm_calibrator_column)
  
  return(scale)
}

technicalNormalize <- function(rawData, use_all_pos_ctrls = FALSE, calibrator.sample.names, highest_of_neg_ctrls = FALSE, read_startCol_gene = '', threshold_value = 20, normalize_value = 8000) {
  read_startCol <- getStartColumn(rawData, gene = read_startCol_gene)
  rawData.samplenames <- rawData[, 1:(read_startCol-1)]
  rawData[, 1:(read_startCol-1)] <- NULL
  
  # positive ctrl normalization
  ctrls <- rawData[, colnames(rawData) %in% ctrls.all] 
  ctrls.mean <- ctrls %>% lapply(., mean) %>% as.data.frame(.) %>% t(.) %>% as.data.frame(.)
  colnames(ctrls.mean) <- c('reads')
  num_ctrls <- ifelse(use_all_pos_ctrls, 8, 3)
  top_ctrls <- ctrls.mean %>% arrange(desc(reads)) %>% slice(1:num_ctrls) %>% rownames(.)
  geomeans <- apply(ctrls[, colnames(ctrls) %in% top_ctrls], 1, geometric.mean)
  avg_geomean <- ifelse(normalize_value == 0, 
                        mean(geomeans),
                        normalize_value)
  well.scale <- sapply(geomeans, function(x) avg_geomean / x)
  pos_norm_data <- rawData[, !colnames(rawData) %in% ctrls.all]
  samples <- nrow(pos_norm_data)
  genes <- ncol(pos_norm_data)
  scale <- matrix(rep(well.scale, times = genes), nrow = samples, ncol = genes)
  pos_norm_data <- scale*as.matrix(pos_norm_data) %>% as.data.frame()
  
  well.scale <- as.data.frame(well.scale)
  colnames(well.scale) <- "Positive Control Normalization Factor"
  # rownames(well.scale) <- rownames(ctrls)
  
  
  # calibration matrix
  cutoff <- ifelse(highest_of_neg_ctrls, max(ctrls[, colnames(ctrls) %in% ctrls.neg]), threshold_value)
  pos_norm_calibrator_column <- pos_norm_data[rownames(pos_norm_data) %in% calibrator.sample.names, ]
  pos_norm_calibrator_column <- threshold(pos_norm_calibrator_column, cutoff = cutoff)
  ##
  rownames(pos_norm_calibrator_column) <- rawData.samplenames$Row[rownames(pos_norm_data) %in% rownames(pos_norm_calibrator_column)]
  ##
  calibration_matrix <- build_calibrator_matrix(pos_norm_calibrator_column)
  ##calibratotion.rows_in_matrix <- rownames(calibration_matrix) %>% gsub("^.*Set\\s*|\\s*\\(.*$", "", .) %>% gsub(" ", "", ., fixed = TRUE)
  calibratotion.rows_in_matrix <- rownames(calibration_matrix)
  tech_norm_data <- lapply(rownames(pos_norm_data), function(sample.name) {
    ##row <- gsub("^.*Set \\s*|\\s* \\(.*$", "", sample.name)
    row <- rawData.samplenames[sample.name, 'Row']
    sample.row <- pos_norm_data[sample.name, ]
    calibration.row <- calibration_matrix[calibratotion.rows_in_matrix == row, ] %>% as.matrix()
    sample.row <- as.matrix(sample.row)
    return(as.data.frame(calibration.row * sample.row))
  }) %>% Reduce(rbind, .) %>% as.data.frame(.)
  tech_norm_data <- as.data.frame(lapply(tech_norm_data, function(x) round(x, digits = 2)))
  tech_norm_data <- threshold(tech_norm_data, cutoff = cutoff)
  tech_norm_data <- cbind(rawData.samplenames, tech_norm_data, ctrls, well.scale) %>% arrange(., Row)
  
  return(tech_norm_data)
}

##function for content normalization
contentNormalize <- function(TechNormData, tissue_colName, HK_genes, TOIs = 'All', read_startCol_gene = '', separate_tissues = TRUE, highest_neg_ctrl = FALSE, threshold_value = 20, normalize_value = 2000) {
  ##subset tissue types and get tissues of interest (TOI)
  TechNormData <- as.data.frame(TechNormData)
  tissue_colName.grep <- unique(grep(paste(tissue_colName, collapse = '|'), colnames(TechNormData), ignore.case = TRUE, value = TRUE))
  if (any(TOIs == "All")) {
    TOIs <- unique(TechNormData[[tissue_colName.grep]]) %>% .[nzchar(.)] %>% .[. != 'Calibrator']
  } else {
    TOIs <- unique(grep(paste(TOIs, collapse = '|'), TechNormData[[tissue_colName.grep]], ignore.case = TRUE, value = TRUE))
  }
  TOIs_data <- TechNormData[TechNormData[[tissue_colName.grep]] %in% TOIs,]
  read_startCol <- getStartColumn(TOIs_data, gene = read_startCol_gene)
  
  ##organize data into subsets for each TOI (optional) and get HKG geomean
  if (separate_tissues) {
    calibrator <- TechNormData[TechNormData[[tissue_colName.grep]] == 'Calibrator',]
    calibrator[, "Content Normalization Factor"] <- c(rep(1, times = nrow(calibrator)))
    TOIs_data_norm <- lapply(TOIs, function(TOI) contentNormalize(TOIs_data, tissue_colName = tissue_colName.grep, 
                                                                  HK_genes = HK_genes, TOIs = TOI, 
                                                                  read_startCol_gene = read_startCol_gene, 
                                                                  separate_tissues = FALSE)) %>% 
      compact(.) %>% Reduce(rbind, .) %>% rbind(., calibrator) %>% arrange(., Row)
    return(TOIs_data_norm)
  } else {
    TOIs_ctrls <- TOIs_data[, colnames(TOIs_data) %in% ctrls.all]
    TOIs_data <- TOIs_data[, !colnames(TOIs_data) %in% ctrls.all]
    
    HKGsearch <- paste(HK_genes, collapse = '|')
    HKGs <- unique(grep(HKGsearch, colnames(TechNormData), ignore.case = TRUE, value = TRUE))
    cutoff <- ifelse(highest_neg_ctrl, max(TOIs_ctrls[, colnames(TOIs_ctrls) %in% ctrls.neg]), threshold_value)
    TOIs_data_HKG <- threshold(TOIs_data, cutoff = cutoff, kill = TRUE)
    TOIs_data_HKG[, 1:(read_startCol-1)] <- NULL     #remove unnecessary columns
    needsNormalizing <- colnames(TOIs_data_HKG) %in% HKGs
    if (!any(needsNormalizing)) {
      out <- cbind(TOIs_data, TOIs_ctrls) %>% rbind(TechNormData[TechNormData[[tissue_colName.grep]] == 'Calibrator',])
      return(out)
    } else {
      positive_ctrl_norm <- TOIs_data[, "Positive Control Normalization Factor"] %>% as.data.frame(.)
      colnames(positive_ctrl_norm) <- "Positive Control Normalization Factor"
      rownames(positive_ctrl_norm) <- rownames(TOIs_data)
      TOIs_data[, "Positive Control Normalization Factor"] <- NULL
      
      TOIs_data_HKG <- TOIs_data_HKG[, needsNormalizing]
      TOIs_data_HKG_geomeans <- apply(TOIs_data_HKG, 1, geometric.mean)
      TOIs_data_HKG_geomeanAVG <- ifelse(normalize_value == 0, 
                                         mean(TOIs_data_HKG_geomeans),
                                         normalize_value)
      TOIs_data_HKG_scale.well <- sapply(TOIs_data_HKG_geomeans, function(x) TOIs_data_HKG_geomeanAVG / x)
      
      ##scale the data by the HKG scaling factor (content normalization)
      TOIs_data_norm <-TOIs_data
      TOIs_data_samplenames <- TOIs_data_norm[, 1:(read_startCol-1)]
      TOIs_data_norm[, 1:(read_startCol-1)] <- NULL
      TOIs_data_norm <- threshold(TOIs_data_norm, cutoff = cutoff, kill = TRUE)
      samples <- nrow(TOIs_data_norm)
      genes <- ncol(TOIs_data_norm)
      TOIs_data_HKG_scale <- matrix(rep(TOIs_data_HKG_scale.well, times = genes), nrow = samples, ncol = genes)
      
      TOIs_data_HKG_scale.well <- as.data.frame(TOIs_data_HKG_scale.well)
      colnames(TOIs_data_HKG_scale.well) <- "Content Normalization Factor"
      rownames(TOIs_data_HKG_scale.well) <- rownames(TOIs_data_norm)
      
      TOIs_data_norm <- TOIs_data_HKG_scale*as.matrix(TOIs_data_norm) %>% as.data.frame(.)
      TOIs_data_norm <- as.data.frame(lapply(TOIs_data_norm, function(x) round(x, digits = 2)))
      TOIs_data_norm <- threshold(TOIs_data_norm, cutoff = cutoff, kill = TRUE)
      TOIs_data_norm <- cbind(TOIs_data_samplenames, TOIs_data_norm, TOIs_ctrls, positive_ctrl_norm, TOIs_data_HKG_scale.well) %>% arrange(., Row)
      
      return(TOIs_data_norm)
    }
  }
}

get_rawCounts_from_NACHO <- function(nachoCounts) {
  nacho_columns <- c(
    'rcc.files',
    'Sample_Attributes.sample_GeneRLF',
    'Lane_Attributes.lane_ID', 
    'plexset_id', 
    'Name', 
    'Count',
    'BD',
    'CartridgeID',
    'FoV'
  )
  counts.raw <- nachoCounts %>% as.data.frame() %>%  .[, colnames(.) %in% nacho_columns]
  colnames(counts.raw) <- c(
    'Filename',
    'CartridgeID',
    'Column',
    'Row',
    'Gene',
    'Raw Counts',
    'Binding Density',
    'RLF Filename',
    'FOV %'
  )
  counts.raw$Row <- utf8ToInt('A') - 1 + as.integer(gsub('S', '', counts.raw$Row))
  counts.raw$Row <- lapply(counts.raw$Row, intToUtf8)
  counts.raw$Column <- str_pad(counts.raw$Column, 2, pad = "0")
  counts.raw$Readname <- paste0(rep('Set ', nrow(counts.raw)), counts.raw$Row, rep(' (', nrow(counts.raw)), 
                                counts.raw$Row, counts.raw$Column, rep(') ', nrow(counts.raw)), counts.raw$Filename)
  samples <- unique(counts.raw$Readname)
  counts.raw <- lapply(samples, function (sample.name) {
    well.samples <- counts.raw[counts.raw$Readname == sample.name, ]
    well.samples.counts <- well.samples[, colnames(counts.raw) %in% c('Gene', 'Raw Counts')]
    rownames(well.samples.counts) <- well.samples.counts[, 'Gene']
    well.samples.counts[, 'Gene'] <- NULL
    colnames(well.samples.counts) <- sample.name
    well.samples.counts <- t(well.samples.counts) %>% as.data.frame(.) %>% select(sort(names(.)))
    well.info <- data.frame(
      row.names = sample.name,
      CartridgeID = well.samples[1, 'CartridgeID'],
      'RLF Filename' = well.samples[1, 'RLF Filename'],
      'FOV %' = well.samples[1, 'FOV %'],
      'Binding Density' = well.samples[1, 'Binding Density'],
      check.names = FALSE
    )
    return(cbind(well.info, well.samples.counts))
  }) %>% Reduce(rbind, .) %>% select(which(!colnames(.) %in% ctrls.all), which(colnames(.) %in% ctrls.all))
  return(counts.raw)
}

