% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing.R
\name{calculate_sample_snp}
\alias{calculate_sample_snp}
\title{calculate_sample_snp}
\usage{
calculate_sample_snp(x, sample = "sample", method = "CPM")
}
\arguments{
\item{x}{dataframe of JH4 pipeline output}

\item{sample}{column name indicating sample}

\item{method}{normalization method; either by CPM (method = "CPM") or percentages (method = "\%"), defaults to "CPM"}
}
\description{
calculate SNP count matrix per sample
}
