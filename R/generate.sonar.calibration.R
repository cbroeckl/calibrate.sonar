#' generate.sonar.calibration
#'
#' Create calibration parameters from a low collision energy SONAR scan - traditional or 'hyrbid' Sonar
#' @details This function uses Apex3D64.exe, distributed by Waters Corp, to perform 4d peak detection in SONAR data, reads that output data, filters, if desired, and calibrates the linear relationships between sonar drift bin and m/z value.  An R-object is returned, and a .csv file exported.  These data contain details on the instrumentation and acquisition method used to collect the calibration data, to ensure that the calibration is only applied to compatible data files.  Masslynx Version, Instrument Serial number, and sonar settings must all match.
#'
#' @param raw.file full path (or path within working directory) to the .RAW file from which the calibration is to be made
#' @param apex3d.dir full path to the directory containing Apex3d64.exe executable file.
#' @param low.energy.function integer value. The mass lynx function containing the low collision energy sonar data. For traditional sonar, set to 1.  For hybrid sonar, set to 2.
#' @param rt.range vector of length two.  i.e. c(2,8).  Retention time range, in minutes, for calibration data to be generated on.
#' @param cal.out numeric.  the ratio of isotope signal (all isotopes) divided by total spectrum signal * 100 much be greated than minPercentSignal to evaluate charge state. Value should be between 0 and 100.
#' @param filter file name or path to export CSV file to.  If you have an existing file and want to append the data, enter full file path/name to that file, and new data will be appended to the same csv file.
#' @return returns an R object which is a data from containing calibration information.
#' @return Additionally exports a .csv file with name defined using the 'cal.out' option.
#' @concept Waters SONAR
#' @concept mass spectrometry
#' @concept data-indepdenent MS/MS
#' @author Corey Broeckling
#' @export

generate.sonar.calibration <- function(
  raw.file = "R:/RSTOR-PMF/Projects/Broeckling_Corey/SONAR/20220124-hemp-cal/20210805_MSIMM_LT_1080.PRO/Data/20220124_sonarcaltest_016.raw",
  apex3d.dir = "C:/Users/cbroeckl/Documents/SONAR/lib",
  low.energy.function = 2,
  rt.range = c(2,12),
  cal.out = "cal.out.csv",
  filter = TRUE
) {

  if(!file.exists(raw.file)) {
    stop("raw file: ", raw.file, "does not exist")
  }

  if(!file.exists(paste0(apex3d.dir, "/", "Apex3d64.exe"))) {
    stop("'Apex3D64.exe' is missing from directory: ", apex3d.dir)
  }

  ## SONAR INF file to be exported.
  ## file named _sonar.INF
  sonar.inf <-
    c("[VERSION#1]",
      "MassLynx v4.2 SCN1018",
      "", "[DATE]",
      "Acquired Date=dd-mmmm-yyyy",   ## must update for each output file
      "<---- comment",
      "Acquired Time=hh:mm:ss",       ## must update for each output file
      "[INSTRUMENT]",
      "Quad Resolution=15.0",
      "Collision Cell Pressure=3.0",
      "Collision Cell Gradient=5.0",
      "Ion Energy=0.2",
      "Post Filter=2",
      "Resolution Mode=Resolution",
      "Pusher Period=76.25",         ## must update for each output file
      "[SONAR#1]",
      "Scan Time=scan.time",         ## must update for each output file
      "Start M/Z=q.start",           ## must update for each output file
      "End M/Z=q.start",             ## must update for each output file
      "Quadrupole Window=q.width",  ## must update for each output file
      "Calibration=q.start,y.low1,y.mid.1,y.high.1",  ## must update for each output file
      "Calibration=q.mid,y.low2,y.mid.2,y.high.2",    ## must update for each output file
      "Calibration=q.stop,y.low3,y.mid.3,y.high.3"    ## must update for each output file
    )

  ## get inst and run info from _HEADER.TXT
  h <- readLines(paste0(raw.file, "/_HEADER.TXT"))
  ms.model <- h[grep("Instrument: ", h)]
  ms.serial.number <- strsplit(ms.model, "#", fixed = TRUE)[[1]][2]
  ms.model <- strsplit(strsplit(ms.model, "Instrument:")[[1]][2], "#", fixed = TRUE)[[1]][1]

  ## get SONAR parameters:
  method <- readLines(paste0(raw.file, "/_extern.inf"))
  mass.lynx.version <- as.numeric(strsplit(strsplit(method[2], " v")[[1]][2], " ")[[1]][1])
  mass.lynx.scn <- as.numeric(strsplit(method[2], " SCN")[[1]][2])
  # "Function Parameters - Function 2 - TOF MS FUNCTION"
  # "SONARQuadrupoleStartMass\t\t\t60.0"
  # "SONARQuadrupoleStopMass\t\t\t\t440.0"
  # "Scan Time (sec)\t\t\t\t\t0.300"
  # "SONARQuadrupolePeakWidth12"
  method <- method[grep("Function 2", method):(grep("Function 3", method)-1)]
  method <- gsub('\t', "", method)
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

  ## if hybrid sonar, we need to calibrate of function 2, which requires some hacking
  if(low.energy.function == 2) {
    # you could move to _func003.cdt and _func003.ind to a separate folder
    dir.create(paste0(dirname(raw.file), "/SonarCalTmp"))
    file.copy(from = paste0(raw.file, c("/_func003.cdt", "/_func003.ind")),
              to = paste0(dirname(raw.file), "/SonarCalTmp", c("/_func003.cdt", "/_func003.ind")))
    file.remove(paste0(raw.file, c("/_func003.cdt", "/_func003.ind")))

    # Create a copy of _TYPES.INF and change the content so it mimics an ordinary SONAR/HDMSE file as illustrated below
    #  #Function 1 : 1
    #  #Function 2 : 1
    #  #Function 3 : 1
    types <- readLines(paste0(raw.file, "/_TYPES.INF"))
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
      "-function", low.energy.function-1
    )
  )

  if(low.energy.function == 2) {
    # you could move to _func003.cdt and _func003.ind to a separate folder
    dir.create(paste0(dirname(raw.file), "/SonarCalTmp"))
    file.copy(to = paste0(raw.file, c("/_func003.cdt", "/_func003.ind")),
              from = paste0(dirname(raw.file), "/SonarCalTmp", c("/_func003.cdt", "/_func003.ind")))

    # Create a copy of _TYPES.INF and change the content so it mimics an ordinary SONAR/HDMSE file as illustrated below
    #  #Function 1 : 1
    #  #Function 2 : 1
    #  #Function 3 : 1
    types <- readLines(paste0(raw.file, "/_TYPES.INF"))
    types[1] <- "#Function 1 : 2"
    sink(paste0(raw.file, "/_TYPES.INF"))
    cat(paste(types, collapse = '\n'))
    sink()

    # create copies of _func002.cdt and _func002.ind and rename them to _func001.cdt and _func001.ind.
    file.remove(paste0(raw.file, c("/_func001.cdt", "/_func001.ind")))
    file.remove(paste0(dirname(raw.file), "/SonarCalTmp"))
  }


  d <- read.csv(
    paste0(raw.file, "/", gsub(".raw", "_Apex3DIons.csv", basename(raw.file))),
    header = TRUE,
    check.names = FALSE
  )

  fit.mobility <- lm(d$mobility~d$m_z)
  fit.start <- lm(d$atDrift1Start ~ d$m_z)
  fit.stop <- lm(d$atDrift1Stop ~ d$m_z)
  if(filter & min(summary(fit.mobility)$r.squared,
         summary(fit.start)$r.squared,
         summary(fit.stop)$r.squared) < 0.99
  ) {
    cat(" -- filtering", '\n')
    fitted <- predict(fit.mobility, list(d$mobility))
    plot(log2(d$mobility/fitted))
    abline(h = 0.1, col = 2)
    abline(h = -0.1, col = 2)
    rem.pts <- which(abs(log2(d$mobility/fitted))> 0.1)

    if((length(rem.pts)/nrow(d)) > 0.1) {
      stop('filtering more than 10% of points may indicate the calibration is poor', '\n')
    }

    d <- d[-rem.pts,]
    fit.mobility <- lm(d$mobility~d$m_z)
    fit.start <- lm(d$atDrift1Start ~ d$m_z)
    fit.stop <- lm(d$atDrift1Stop ~ d$m_z)

  }

  plot(d$m_z, d$mobility, pch = 19, col = 1, cex = 0.5,
         xlab = "mz", ylab = "sonar.bin", main = "sonar calibration")
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

  out <- c(
    "ms.model" = ms.model,
    "ms.serial.number" = ms.serial.number,
    "mass.lynx.version" = mass.lynx.version,
    "mass.lynx.scn" = mass.lynx.scn,
    "interscan.delay" = interscan.delay,
    "scan.time" = scan.time,
    "q.start" = q.start,
    "q.stop" = q.stop,
    "intercept.start" = intercept.start,
    'intercept.mid' = intercept.mid,
    "intercept.stop" = intercept.stop,
    "slope.start" = slope.start,
    "slope.mid" = slope.mid,
    "slope.stop" = slope.stop
  )

  out <- as.data.frame(t(as.matrix(out)))


  if(file.exists(cal.out)) {
    existing <- read.csv(cal.out)
    if(!all(names(existing) == names(out))) {
      warning("column names in cal.out file provided are not as expected, new cal.out writtend to:", '\n',
              paste0(getwd(), "/new.calibration.output.csv"))
      write.csv(out, file = paste0(getwd(), "/new.calibration.output.csv"))
    } else {
      existing[nrow(existing)+1,] <- out[1,]
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






