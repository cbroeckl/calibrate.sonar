% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate.sonar.calibration.R
\name{generate.sonar.calibration}
\alias{generate.sonar.calibration}
\title{generate.sonar.calibration}
\usage{
generate.sonar.calibration(
  raw.file = "C:/Users/cbroeckl/Documents/temp/hybridSonarCal.raw",
  apex3d.dir = "C:/Users/cbroeckl/Documents/SONAR/lib",
  low.energy.function = 1,
  rt.range = c(2, 12),
  cal.out = "cal.out.csv",
  filter = TRUE,
  filter.range = 0.1,
  allow.high.energy.calibration = FALSE
)
}
\arguments{
\item{raw.file}{full path (or path within working directory) to the .RAW file from which the calibration is to be made}

\item{apex3d.dir}{full path to the directory containing Apex3d64.exe executable file.}

\item{low.energy.function}{integer value. The mass lynx function containing the low collision energy sonar data. For traditional sonar, set to 1.  For hybrid sonar, set to 2.}

\item{rt.range}{vector of length two.  i.e. c(2,8).  Retention time range, in minutes, for calibration data to be generated on.}

\item{cal.out}{ile name or path to export CSV file to.  If you have an existing file and want to append the data, enter full file path/name to that file, and new data will be appended to the same csv file.}

\item{filter}{logical.  If TRUE (default), calculate approximate sonar bin to mass relationship based on method settings.  Calculate ratio of actual vs predicted bin and express as a proportion.  1 = perfect prediction. Ratio values outside 1 +/ filter.range value are removed before calibrating.}

\item{filter.range}{numeric. default = 0.1.  Range (as a proportion) of drift bin assignment for filtering.}

\item{allow.high.energy.calibration}{logical If FALSE (default), trying generating a calibration using sonar at high collision energy will generate an error.  If TRUE, only a warning is issued.  Calibrating on high collision energy data comes with risk, and is more likely to be inaccurate.}
}
\value{
returns an R object which is a data from containing calibration information.

Additionally exports a .csv file with name defined using the 'cal.out' option.
}
\description{
Create calibration parameters from a low collision energy SONAR scan - traditional or 'hyrbid' Sonar
}
\details{
This function uses Apex3D64.exe, distributed by Waters Corp, to perform 4d peak detection in SONAR data, reads that output data, filters, if desired, and calibrates the linear relationships between sonar drift bin and m/z value.  An R-object is returned, and a .csv file exported.  These data contain details on the instrumentation and acquisition method used to collect the calibration data, to ensure that the calibration is only applied to compatible data files.  Masslynx Version, Instrument Serial number, and sonar settings must all match.
}
\author{
Corey Broeckling
}
\concept{Waters SONAR}
\concept{data-independent MS/MS}
\concept{mass spectrometry}
