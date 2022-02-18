#' generate.sonar.calibration
#'
#' Create calibration parameters from a low collision energy SONAR scan - traditional or 'hyrbid' Sonar
#' @details This function uses Apex3D64.exe, distributed by Waters Corp, to perform 4d peak detection in SONAR data, reads that output data, filters, if desired, and calibrates the linear relationships between sonar drift bin and m/z value.  An R-object is returned, and a .csv file exported.  These data contain details on the instrumentation and acquisition method used to collect the calibration data, to ensure that the calibration is only applied to compatible data files.  Masslynx Version, Instrument Serial number, and sonar settings must all match.
#'
#' @param raw.file full path (or path within working directory) to the .RAW file from which the calibration is to be made
#' @param apex3d.dir full path to the directory containing Apex3d64.exe executable file.
#' @param low.energy.function integer value. The mass lynx function containing the low collision energy sonar data. For traditional sonar, set to 1.  For hybrid sonar, set to 2.
#' @param rt.range vector of length two.  i.e. c(2,8).  Retention time range, in minutes, for calibration data to be generated on.
#' @param cal.out ile name or path to export CSV file to.  If you have an existing file and want to append the data, enter full file path/name to that file, and new data will be appended to the same csv file.
#' @param filter logical.  If TRUE (default), calculate approximate sonar bin to mass relationship based on method settings.  Calculate ratio of actual vs predicted bin and express as a proportion.  1 = perfect prediction. Ratio values outside 1 +/ filter.range value are removed before calibrating.
#' @param filter.range numeric. default = 0.1.  Range (as a proportion) of drift bin assignment for filtering.  
#' @param allow.high.energy.calibration logical If FALSE (default), trying generating a calibration using sonar at high collision energy will generate an error.  If TRUE, only a warning is issued.  Calibrating on high collision energy data comes with risk, and is more likely to be inaccurate.  
#' @return returns an R object which is a data from containing calibration information.
#' @return Additionally exports a .csv file with name defined using the 'cal.out' option.
#' @concept Waters SONAR
#' @concept mass spectrometry
#' @concept data-independent MS/MS
#' @author Corey Broeckling
#' @export

generate.sonar.calibration <- function(
  raw.file = "C:/Users/cbroeckl/Documents/temp/hybridSonarCal.raw",
  apex3d.dir = "C:/Users/cbroeckl/Documents/SONAR/lib",
  low.energy.function = 1,
  rt.range = c(2,12),
  cal.out = "cal.out.csv",
  filter = TRUE,
  filter.range = 0.1,
  allow.high.energy.calibration = FALSE
) {

  if(!file.exists(raw.file)) {
    stop("raw file: ", raw.file, "does not exist")
  }

  if(!file.exists(paste0(apex3d.dir, "/", "Apex3d64.exe"))) {
    stop("'Apex3D64.exe' is missing from directory: ", apex3d.dir)
  }

  # ## SONAR INF file to be exported.
  # ## file named _sonar.INF
  # sonar.inf <-
  #   c("[VERSION#1]",
  #     "MassLynx v4.2 SCN1018",
  #     "", "[DATE]",
  #     "Acquired Date=dd-mmmm-yyyy",   ## must update for each output file
  #     "<---- comment",
  #     "Acquired Time=hh:mm:ss",       ## must update for each output file
  #     "[INSTRUMENT]",
  #     "Quad Resolution=15.0",
  #     "Collision Cell Pressure=3.0",
  #     "Collision Cell Gradient=5.0",
  #     "Ion Energy=0.2",
  #     "Post Filter=2",
  #     "Resolution Mode=Resolution",
  #     "Pusher Period=76.25",         ## must update for each output file
  #     "[SONAR#1]",
  #     "Scan Time=scan.time",         ## must update for each output file
  #     "Start M/Z=q.start",           ## must update for each output file
  #     "End M/Z=q.start",             ## must update for each output file
  #     "Quadrupole Window=q.width",  ## must update for each output file
  #     "Calibration=q.start,y.low1,y.mid.1,y.high.1",  ## must update for each output file
  #     "Calibration=q.mid,y.low2,y.mid.2,y.high.2",    ## must update for each output file
  #     "Calibration=q.stop,y.low3,y.mid.3,y.high.3"    ## must update for each output file
  #   )

  ## get inst and run info from _HEADER.TXT
  h <- suppressWarnings(readLines(paste0(raw.file, "/_HEADER.TXT")))
  ms.model <- h[grep("Instrument: ", h)]
  ms.serial.number <- strsplit(ms.model, "#", fixed = TRUE)[[1]][2]
  ms.model <- strsplit(strsplit(ms.model, "Instrument:")[[1]][2], "#", fixed = TRUE)[[1]][1]

  ## get SONAR parameters:
  method <- suppressWarnings(readLines(paste0(raw.file, "/_extern.inf")))
  mass.lynx.version <- as.numeric(strsplit(strsplit(method[2], " v")[[1]][2], " ")[[1]][1])
  mass.lynx.scn <- as.numeric(strsplit(method[2], " SCN")[[1]][2])
  # "Function Parameters - Function 2 - TOF MS FUNCTION"
  # "SONARQuadrupoleStartMass\t\t\t60.0"
  # "SONARQuadrupoleStopMass\t\t\t\t440.0"
  # "Scan Time (sec)\t\t\t\t\t0.300"
  # "SONARQuadrupolePeakWidth12"
  method <- method[grep(paste("Function", low.energy.function), method):(grep(paste("Function", (low.energy.function+1)), method)-1)]
  method <- gsub('\t', "", method)
  
  ##check that this is a sonar method
  is.sonar <- method[grep("UseSONARMode", method)]
  is.sonar <- substr(is.sonar, nchar(is.sonar)-4+1, nchar(is.sonar))
  if(!is.sonar == "TRUE") {
    cat(is.sonar, '\n')
    stop(paste(" - function", low.energy.function, "is not a sonar function", '\n'))
  }
  
  ## check collision energy
  CE <- method[which(startsWith(method, "Collision Energy"))]
  CE <- as.numeric(unlist(strsplit(CE, ")"))[2])
  if(CE >= 10) {
    cat(CE, '\n')
    if(allow.high.energy.calibration) {
      stop(paste("Calibration aborted - function", low.energy.function, "used a collion energy of", CE, '\n'))
    }
    warning(paste(" - function", low.energy.function, "used a collion energy of", CE, '\n'))
  }
  
  scan.time <- method[grep("Scan Time", method)]
  scan.time <- as.numeric(unlist(strsplit(scan.time, ")"))[2])
  q.start <- method[grep("SONARQuadrupoleStartMass", method)]
  q.start <- as.numeric(gsub("SONARQuadrupoleStartMass", "", q.start))
  q.stop <- method[grep("SONARQuadrupoleStopMass", method)]
  q.stop <- as.numeric(gsub("SONARQuadrupoleStopMass", "", q.stop))
  q.width <- method[grep("SONARQuadrupolePeakWidth", method)]
  q.width <- as.numeric(gsub("SONARQuadrupolePeakWidth", "", q.width))
  interscan.delay <- method[grep("Interscan Time", method)]
  interscan.delay <- as.numeric(gsub("Interscan Time (sec)", "", interscan.delay, fixed = TRUE))
  pusher.period <- method[grep("ADC Pusher Frequency", method)]
  pusher.period <- as.numeric(strsplit(pusher.period, "s)")[[1]][2])

  ## if hybrid sonar, we need to calibrate of function 2, which requires some hacking
  if(low.energy.function == 2) {
    # you could move to _func003.cdt and _func003.ind to a separate folder
    # suppressWarnings(dir.create(paste0(dirname(raw.file), "/SonarCalTmp")))
    file.rename(from = paste0(raw.file, c("/_func003.cdt", "/_func003.ind")),
                to = paste0(raw.file, c("/_func003_.cdt", "/_func003_.ind")))

    # Create a copy of _TYPES.INF and change the content so it mimics an ordinary SONAR/HDMSE file
    types <- suppressWarnings(readLines(paste0(raw.file, "/_TYPES.INF")))
    types[1] <- "#Function 1 : 1"
    sink(paste0(raw.file, "/_TYPES.INF"))
    cat(paste(types, collapse = '\n'))
    sink()

    # create copies of _func002.cdt and _func002.ind and rename them to _func001.cdt and _func001.ind.
    file.copy(from = paste0(raw.file, c("/_func002.cdt", "/_func002.ind")),
              to = paste0(raw.file, c("/_func001.cdt", "/_func001.ind")))
  }


  ## run Apex3d via system call

  system(
    paste(
      normalizePath(paste0(apex3d.dir, "/Apex3d64.exe")),
      "-pRawDirName", normalizePath(raw.file),
      "-outputDirName", normalizePath(dirname(raw.file)),
      "-outputUserDirName", normalizePath(raw.file),
      "-lockMassZ1 556.2771",
      "-lockmassToleranceAMU 0.25",
      "-writeBinary 0 -bCSVOutput 1",
      "-leThresholdCounts 2000",
      "-startingRTMin", rt.range[1], "-endingRTMin", rt.range[2],
      "-function", low.energy.function
    ), show.output.on.console = FALSE
  )
  
  d <- read.csv(
    paste0(raw.file, "/", gsub(".raw", "_Apex3DIons.csv", basename(raw.file))),
    header = TRUE,
    check.names = FALSE
  )
  
  if(nrow(d) < 200) {
    cat("rerunning apex3d with increase sensitivity to detect more peaks for calibration", '\n')
    system(
      paste(
        normalizePath(paste0(apex3d.dir, "/Apex3d64.exe")),
        "-pRawDirName", normalizePath(raw.file),
        "-outputDirName", normalizePath(dirname(raw.file)),
        "-outputUserDirName", normalizePath(raw.file),
        "-lockMassZ1 556.2771",
        "-lockmassToleranceAMU 0.25",
        "-writeBinary 0 -bCSVOutput 1",
        "-leThresholdCounts 200",
        "-startingRTMin", rt.range[1], "-endingRTMin", rt.range[2],
        "-function", low.energy.function-1
      ), show.output.on.console = FALSE
    )
    
    d <- read.csv(
      paste0(raw.file, "/", gsub(".raw", "_Apex3DIons.csv", basename(raw.file))),
      header = TRUE,
      check.names = FALSE
    )
  }
  
  ## if hybrid sonar, unhack hacking

  if(low.energy.function == 2) {
    # you could move to _func003.cdt and _func003.ind to a separate folder
    # suppressWarnings(dir.create(paste0(dirname(raw.file), "/SonarCalTmp")))
    file.rename(from = paste0(raw.file, c("/_func003_.cdt", "/_func003_.ind")),
                to = paste0(raw.file, c("/_func003.cdt", "/_func003.ind")))
    
    # Create a copy of _TYPES.INF and change the content so it mimics an ordinary SONAR/HDMSE file
    types <- suppressWarnings(readLines(paste0(raw.file, "/_TYPES.INF")))
    types[1] <- "#Function 1 : 2"
    sink(paste0(raw.file, "/_TYPES.INF"))
    cat(paste(types, collapse = '\n'))
    sink()
    
    # create copies of _func002.cdt and _func002.ind and rename them to _func001.cdt and _func001.ind.
    file.remove(paste0(raw.file, c("/_func001.cdt", "/_func001.ind")))
  }
  
  d <- read.csv(
    paste0(raw.file, "/", gsub(".raw", "_Apex3DIons.csv", basename(raw.file))),
    header = TRUE,
    check.names = FALSE
  )
  
  
  par(mfrow = c(2,1))
  
  
  if(filter) {

    ## will estimate a and b values to provide an approximate filter and remove noise. 
    q.range <- q.stop - q.start
    bins.per.mass <- q.width / (q.range/200)
    est.a <- 200/ q.range ## - bins.per.mass
    est.b <- median(d$mobility - (est.a*d$m_z))
    
    pred.mobility <- est.a*d$m_z + est.b
    cols <- rep(1, nrow(d))
    rm.pts <- which(d$mobility/pred.mobility < (1 - filter.range) | d$mobility/pred.mobility > (1 + filter.range) )
    par(mfrow = c(2,2))
    if(length(rm.pts) > 0) {
      cols[rm.pts] <- 2 
    }
    plot(d$m_z, (d$mobility/pred.mobility), col = cols, xlab = "m/z", ylab = "actual vs predicted sonar bin", pch = 19, cex = 0.5, main = "pre-calibration filtering")
    
    if(length(rm.pts) > 0) {
      d <- d[-rm.pts,]
    }
    plot(d$m_z, (d$mobility/pred.mobility[-rm.pts]), xlab = "m/z", ylab = "actual vs predicted sonar bin", pch = 19, cex = 0.5, main = "outliers removed")
    
  }

  
  fit.mobility <- lm(d$mobility~d$m_z)
  fit.start <- lm(d$atDrift1Start ~ d$m_z)
  fit.stop <- lm(d$atDrift1Stop ~ d$m_z)

  plot(d$m_z, d$mobility, pch = 19, col = 1, cex = 0.5,
         xlab = "mz", ylab = "sonar.bin", main = "sonar calibration")
  legend(x = "topleft", legend = c("low.mass.limit", "center.mass", "high.mass.limit"), text.col = c(4,1,2), bty = "n", cex = 1)
  points(abline(fit.mobility), col = 1)
  intercept.mid <- round(as.numeric(fit.mobility$coefficients[1]), digits =5)
  slope.mid <- round(as.numeric(fit.mobility$coefficients[2]), digits = 5)


  points(d$m_z, d$atDrift1Start, pch = 19, col = 2, cex = 0.5)
  points(abline(fit.start), col = 2)
  intercept.start <- round(as.numeric(fit.start$coefficients[1]), digits =5)
  slope.start <- round(as.numeric(fit.start$coefficients[2]), digits = 5)


  points(d$m_z, d$atDrift1Stop, pch = 19, col = 4, cex = 0.5)
  points(abline(fit.stop), col = 4)
  intercept.stop <- round(as.numeric(fit.stop$coefficients[1]), digits =5)
  slope.stop <- round(as.numeric(fit.stop$coefficients[2]), digits = 5)

  mtext(paste0(q.start, ":", q.stop, ", width = ", q.width), side = 3, outer = TRUE)
  
  ## check for consistent slope
  slope.cv <- sd(c(slope.mid, slope.start, slope.stop))/mean(c(slope.mid, slope.start, slope.stop))
  if(slope.cv > 0.05) {
    warning("the slope values have more variation that expected, please consider refining peak detection", '\n')
  }

  
  ## check for strong r^2
  if(min(
    summary(fit.mobility)$r.squared,
    summary(fit.start)$r.squared,
    summary(fit.stop)$r.squared
  ) < 0.995
  ) {
    warning("the minimum r-squared values is less than 0.995, please consider refining peak detection", '\n')
  }

  plot(fit.mobility$model$'d$mobility', fit.mobility$residuals, main = "residuals (center mass)", pch = 19, cex = 0.5, 
       xlab = "sonar bin", ylab = "residual")
  
  out <- c(
    "date.time" = paste(Sys.time()),
    "ms.model" = ms.model,
    "ms.serial.number" = ms.serial.number,
    "mass.lynx.version" = mass.lynx.version,
    "mass.lynx.scn" = mass.lynx.scn,
    "interscan.delay" = interscan.delay,
    "pusher.period" = pusher.period,
    "scan.time" = scan.time,
    "q.start" = round(q.start, digits = 1),
    "q.stop" = round(q.stop, digits = 1),
    "q.width" = round(q.width, digits = 1),
    "intercept.start" = round(intercept.start, digits = 4),
    'intercept.mid' = round(intercept.mid, digits = 4),
    "intercept.stop" = round(intercept.stop, digits = 4),
    "slope.start" = round(slope.start, digits = 4),
    "slope.mid" = round(slope.mid, digits = 4),
    "slope.stop" = round(slope.stop, digits = 4)
  )

  out <- as.data.frame(t(as.matrix(out)))


  if(file.exists(paste0(getwd(), "/", cal.out))) {
    existing <- read.csv(cal.out)
    if(!all(names(existing) == names(out))) {
      warning("column names in cal.out file provided are not as expected, new cal.out writtend to:", '\n',
              paste0(getwd(), "/new.calibration.output.csv"))
      write.csv(out, file = paste0(getwd(), "/new.calibration.output.csv"))
    } else {
      out <- rbind(out, existing)
      write.csv(out, file = cal.out, row.names = FALSE)
      cat("output calibration file supplemented with new calibrations and written: ", '\n',
          cal.out)
    }
  } else {
    write.csv(out, file = cal.out, row.names = FALSE)
    cat("new output calibration file written: ", '\n',
        cal.out)
  }

  class(out) <- c(class(out), "sonar.calibration")
  return(out)
}






