# MAIN SCRIPTS
library(dplyr)
##### Set working directory ####
dir <- "/Users/vinhthuantran/Desktop/ME/QFAB/Project/"
setwd(dir)

# Select data directory -------------

#data.folder <- paste(dir, "Data", sep = "/")
data.folder <- paste0(dir, "Data@QFAB/HCSRawData")
file.names = list.files(data.folder, pattern = "\\.csv$")


for (ind in 1 : length(file.names)) {
  
  file.wo.ext = paste0(substr(file.names[ind], 10, nchar(file.names[ind]) - 4))
  data.df <- paste(file.wo.ext, "_data", sep = '')
  assign(data.df, read.csv(file=file.path(data.folder, file.names[ind]),skip=9))
  
  # This If condition is exclusively used for Intensity_data.
  if (any(grepl("Cell.Area", colnames(get(data.df))))) {
    area_col_index  <- grep("Cell.Area", colnames(get(data.df))) # Search for Cell.Area col
    nuclei_col_index <- grep("^(?=.*Nuclei.Selected)(?!.*Cell.Area)", 
                             colnames(get(data.df)), 
                             ignore.case = T, perl =T) # Search for columns of interest except Cell.Area col
    
    
    intensity.data <- get(data.df)
    # Implement division between Nuclei..Columns and Cell Area Column using sweep()
    scaled_data <- sweep(intensity.data[ ,nuclei_col_index], 
                         MARGIN = 1, FUN = "/", STATS = intensity.data[ ,area_col_index])
    
    intensity.data[,nuclei_col_index] <- scaled_data
    intensity.data <- select(intensity.data, -area_col_index) # Drop the Cell Area Column
    
    # Replace orginal data with scaled one
    assign(data.df, intensity.data)
    
  }
  # Load data --------------------
  mydata <- select(get(data.df), 
                   Row, Column, Compound, Cell.Type, starts_with("Nuclei.Selected"))
  
  # Mean, Median and SD ----------------------------
  # mydata_stats (Texture_stats or Intensity_stats)
  # stores calculated mean, median and standard deviation for each well.
  mydata_stats <- paste(file.wo.ext, "stats", sep = '_')
  assign(mydata_stats, mydata %>% 
           group_by(Row,Column)  %>% 
           summarise_at(vars(starts_with("Nuclei.Selected")), 
                        funs(mean,median,sd)) %>% 
           ungroup())
    
  # Calculate number of wells (or distinct (row,col))
  # mydata_wells (Texture_wells or Intensity_wells)
  mydata_wells <- paste(file.wo.ext, "wells", sep = '_')
  assign(mydata_wells, mydata %>% 
           select(-starts_with("Nuclei.Selected")) %>% 
           group_by(Row,Column, Compound, Cell.Type) %>% 
           summarise(Total.of.Well = n()) %>%
           ungroup())
  
  # Join datasets 
  # mydata_cells (Texture_cells or Intensity_cells)
  # stores all mydata_stats and mydata_wells
  
 mydata_cells <- paste(file.wo.ext, "cells", sep = '_')
 assign(mydata_cells, 
        left_join(get(mydata_wells), get(mydata_stats), by = c("Row", "Column")))
  
  # Calculate Coefficient of variance (CV) --------------
  # CV = (SD/Mean) * 100
  
  meanCols <- grep('_mean', colnames(get(mydata_cells)))
  meanColNames <- names((get(mydata_cells))[, meanCols])
    
  sdCols <- grep('_sd', colnames(get(mydata_cells)))
  sdColNames <- names((get(mydata_cells))[, sdCols])
  
  cvColNames <- paste0(substr(sdColNames, 
                              1, nchar(colnames((get(mydata_cells))[meanCols])) - 5), '_cv')
  
  for (ind in 1 : length(cvColNames)) {
    # SD_numerator <- sdColNames[ind]
    # Mean_denominator <- meanColNames[ind]
    CV <- paste0(sdColNames[ind], '/', meanColNames[ind], '* 100')
    assign(mydata_cells, get(mydata_cells) %>% 
      mutate_(.dots = setNames(CV, cvColNames[ind])))
  }
  
 
  # SAVE Mean, Median and Standard Deviation for each well---------------------------
  # cell.stats.file <- paste(file.wo.ext, "stats.csv", sep = '_')
  # write.csv(get(cell.stats), file = cell.stats.file)
  
  # A quick look at number of cell types in dataset
  cellTypes <- mydata %>% 
    group_by(Cell.Type) %>%
    summarise(Frequency = n())
  
  type_names <- cellTypes$Cell.Type  # list of 14 type names:
  # 1403  1705  1803  1809  1813  1815  1816  2141  2261  2361  2509  2704 21551 22551
  
  by_DMSO <- get(mydata_cells) %>% 
    filter(Compound == "DMSO") # Sort the data by DMSO compound
  
  # test <- by_DMSO %>% dplyr::group_by(.,Cell.Type, Compound) %>%
  # mutate_at(vars(ends_with('_median')), funs(Mean = mean)) %>%
  #   rename_at(vars(ends_with('_median_Mean')),
  #             funs(paste0('Mean ', substr(., 1, nchar(.) - 5))))
  
  
  
  log2_summary <- {}
  
  for (eachType in 1:length(type_names)) {
    
    # DMSO.eachType will be sth like "Texture_DSMO1403", 
    # "Texture_DMSO1705", "Texture_DMSO1803", etc.
    
    DMSO.eachType <- paste("DMSO", type_names[eachType], sep = "")
    DMSO.eachType <- paste(file.wo.ext, DMSO.eachType, sep = "_")
    
    assign(DMSO.eachType, by_DMSO %>%
             select(Cell.Type, Compound, ends_with("_median")) %>%
             filter(Cell.Type == type_names[eachType]))
    
    # DMSO MEAN for each type ----------------------------
    # DMSO.eachType.mean looks like Texture_DMSO1403Mean
    
    DMSO.eachType.mean <- paste(DMSO.eachType, "Mean", sep = "")
    
    assign(DMSO.eachType.mean, get(DMSO.eachType) %>%
             dplyr::group_by(.,Cell.Type, Compound) %>%
             summarise_all(funs(Mean = mean)) %>%
             rename_at(vars(ends_with('_Mean')),
                       funs(paste0('Mean ', substr(., 1, nchar(.) - 5)))) %>%
             ungroup())
    
    
    
    ####### Log transformation for DMSO Mean of each type #####
    assign(DMSO.eachType.mean, get(DMSO.eachType.mean) %>% 
             mutate_at(vars(-Cell.Type, -Compound), funs(Log2 = log2)) %>%
             rename_at(vars(ends_with('_Log2')),
                       funs(paste0('Log2 ', substr(., 1, nchar(.) - 5)))))
    
    
    # # Left join Mean and log transformation of DMSO to mydata_cells which is Texture_cells or Intensity_cells
    # assign(mydata_cells, left_join(get(mydata_cells), Texture_DMSO1403Mean, by = c("Compound","Cell.Type")))
    
    # SAVE DMSO Mean ------------------------------
    # DMSO_mean_file_wExt = paste0(DMSO.eachType.mean, ".csv")
    # write.csv(get(DMSO.eachType.mean), DMSO_mean_file_wExt)
    
    
    # CALCULATE OVERALL MEAN -------------------------------
    #Texture_mydata_1403
    mydata.eachType <- paste("mydata_",type_names[eachType], sep='')
    mydata.eachType <- paste(file.wo.ext, mydata.eachType, sep = "_")
    
    assign(mydata.eachType, get(mydata_cells) %>% 
             filter(Cell.Type == type_names[eachType])) # Sort the data by type.
    
    #Texture_mydata_1403_mean
    mydata.eachType.mean <- paste(mydata.eachType, "_mean", sep = '')
    
    assign(mydata.eachType.mean, get(mydata.eachType) %>%
             dplyr::select(Cell.Type, ends_with('_median')) %>%
             group_by(Cell.Type) %>%
             summarise_all(funs(Mean = mean)) %>%
             rename_at(vars(ends_with('_Mean')),
                       funs(paste0('Mean ', substr(., 1, nchar(.) - 5)))))
    
    #Log Transfromation
    assign(mydata.eachType.mean, get(mydata.eachType.mean) %>% 
             mutate_at(vars(-Cell.Type), funs(Log2 = log2)) %>%
             rename_at(vars(ends_with('_Log2')),
                       funs(paste0('Log2 ', substr(., 1, nchar(.) - 5)))))
    
    
    # Join data
    assign(mydata.eachType, left_join(get(mydata.eachType), get(mydata.eachType.mean), by = "Cell.Type"))
    
    #SAVE overall mean file ----------------
    # overall_mean_file_wExt = paste0(mydata.eachType.mean, ".csv")
    # write.csv(get(mydata.eachType.mean), overall_mean_file_wExt)
    
    # Log2(avg(DMSO)/avg(overall)) -----------------------------
    M <- left_join(get(mydata.eachType), get(DMSO.eachType.mean), by = c("Compound", "Cell.Type"))
    
    M <- M %>%
      dplyr::select(Row, Column, Compound, Cell.Type, starts_with("Log2"))
    
    S <- M[, grepl("*\\.x$", names(M))] - M[, grepl("*\\.y$", names(M))]
    
    mydata_centred_combined <- cbind(get(mydata.eachType)[,c(1:4)],S)
    log2_summary <- rbind(log2_summary, mydata_centred_combined)
    
  }
  
  #####SAVE Log2 Summary#####
  mydata.log2 <- paste0(file.wo.ext, "_log2_summary")
  assign(mydata.log2, log2_summary)
  
  log2_file <- paste0(file.wo.ext, "_log2.csv")
  # assign(log2_file, log2_summary)
  write.csv(log2_summary, log2_file)
  
}