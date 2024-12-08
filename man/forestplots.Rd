% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Forestplots_tikka41124.R
\name{forestplots}
\alias{forestplots}
\title{Forestplottings}
\usage{
forestplots(NAFLD, Outcome, Group, name, ordera, oute, first, e, xlim)
}
\arguments{
\item{NAFLD}{The dataset with steroids and cases}

\item{Outcome}{The case (e.g. steatosis). Typically a column name in the dataset}

\item{Group}{All, Female, Male (typically a gender column)}

\item{name}{An header add (of the case) for the boxplot's header}

\item{oute}{The label title. Here, almost the same as case.}
}
\value{
The boxplots. Yes.
}
\description{
Converts data to forestplots.
Written by Pauli Tikka. University of Turku. Updated 1.11.24.
}
\examples{
ie=tv_half_log22;Outcome='Steatosis Grade';Out='Steatosis'; oute='Steatosis Grade';Group='All';
boxplots(ie,Group,Outcome,Out,oute)
}
