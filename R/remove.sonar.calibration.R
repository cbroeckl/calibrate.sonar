#' generate.sonar.calibration
#'
#' Create calibration parameters from a low collision energy SONAR scan - traditional or 'hybrid' Sonar
#' @details This function uses Apex3D64.exe, distributed by Waters Corp, to perform 4d peak detection in SONAR data, reads that output data, filters, if desired, and calibrates the linear relationships between sonar drift bin and m/z value.  An R-object is returned, and a .csv file exported.  These data contain details on the instrumentation and acquisition method used to collect the calibration data, to ensure that the calibration is only applied to compatible data files.  Masslynx Version, Instrument Serial number, and sonar settings must all match.
#'
#' @param raw.file full path (or path within working directory) to the .RAW file(s) to which sonar calibrations have previously been applied  An entire directory can be input instead, and if a directory is chosen all appropriate .raw files in that directory will be uncalibrated.
#' @param recursive logical.  if 'raw.file' is a path to a directory of raw files, should the function also look for raw files in all subdirectories? see ?list.files
#' @return returns nothing to the R environment
#' @return selected .RAW files will have been modified to remove the _sonar.INF and optionally some file copying and renaming if using hybrid sonar.
#' @concept Waters SONAR
#' @concept mass spectrometry
#' @concept data-independent MS/MS
#' @author Corey Broeckling
#' @export

remove.sonar.calibration <- function(
  raw.file = "R:/RSTOR-PMF/Projects/Broeckling_Corey/SONAR/20220124-hemp-cal/20210805_MSIMM_LT_1080.PRO/Data/20220124_sonarcaltest_011.raw",
  recursive = FALSE
) {

  if(!file.exists(raw.file)) {
    stop("raw file: ", raw.file, "does not exist")
  }
  
  if(tolower(substring(raw.file, first = nchar(raw.file)-2, last = nchar(raw.file))) != 'raw') {
    raw.file <- list.files(raw.file, pattern = ".raw", ignore.case = TRUE, recursive = recursive)
  }
  
  for(i in 1:length(raw.file)) {
    h <- readLines(paste0(raw.file[i], "/_HEADER.TXT"))

    
    ## get SONAR parameters:
    method <- readLines(paste0(raw.file[i], "/_extern.inf"))
    
    f1 <- method[grep("Function 1", method):(grep("Function 2", method)-1)]
    f1.sonar <- unlist(strsplit(f1[grep("UseSONARMode", f1)], "\t", fixed = TRUE))
    f1.sonar <- as.logical(f1.sonar[length(f1.sonar)])
    f2 <- method[grep("Function 2", method):(grep("Function 3", method)-1)]
    f2.sonar <- unlist(strsplit(f2[grep("UseSONARMode", f2)], "\t", fixed = TRUE))
    f2.sonar <- as.logical(f2.sonar[length(f2.sonar)])
    
    if(!(f1.sonar | f2.sonar)) {
      next(raw.file[i], "did not use SONAR", '\n')
    }
    
    hybrid <-  !(f1.sonar && f2.sonar)
    
    file.remove(paste0(raw.file[i], "/_sonar.INF"))
                
    ## if hybrid sonar, we need to calibrate of function 2, which requires some hacking
    if(hybrid) {
      # you could move to _func003.cdt and _func003.ind to a separate folder
      suppressWarnings(dir.create(paste0(dirname(raw.file[i]), "/SonarCalTmp")))
      file.rename(from = paste0(raw.file[i], c("/_func003_.cdt", "/_func003_.ind")),
                to = paste0(raw.file[i], c("/_func003.cdt", "/_func003.ind")))
      
      # Create a copy of _TYPES.INF and change the content so it mimics an ordinary SONAR/HDMSE file as illustrated below
      #  #Function 1 : 1
      #  #Function 2 : 1
      #  #Function 3 : 1
      types <- suppressWarnings(readLines(paste0(raw.file[i], "/_TYPES.INF")))
      types[1] <- "#Function 1 : 2"
      sink(paste0(raw.file[i], "/_TYPES.INF"))
      cat(paste(types, collapse = '\n'))
      sink()
      
      
      # create copies of _func002.cdt and _func002.ind and rename them to _func001.cdt and _func001.ind.
      file.remove(paste0(raw.file, c("/_func001.cdt", "/_func001.ind")))
    }
    
  }
}






