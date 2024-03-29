#' nnspat: A package for NN Methods and Their Use in Testing Spatial Patterns
#'
#' \code{nnspat} is a package for computation of spatial pattern tests based on NN relations and 
#' generation of various spatial patterns.
#'
#' The \code{nnspat} package contains the functions for segregation/association tests based on nearest neighbor contingency
#' tables (NNCTs), and tests for species correspondence, NN symmetry and reflexivity based on the corresponding
#' contingency tables and functions for generating patterns of segregation, association, uniformity and various
#' non-random labeling (RL) patterns for disease clustering for data in two (or more) dimensions.
#' See (\insertCite{dixon:1994,ceyhan:eest-2010,ceyhan:jkss-posthoc-2017;textual}{nnspat}).
#'
#' @section The \code{nnspat} functions:
#' The \code{nnspat} functions can be grouped as Auxiliary Functions, NNCT Functions, SCCT Functions, RCT Functions
#' NN-Symmetry Functions and the Pattern (Generation) Functions.
#'
#' @section Auxiliary Functions:
#' Contains the auxiliary functions used in NN methods, such as indices of NNs, number of shared NNs, Q, R and T 
#' values, and so on. 
#' In all these functions the data sets are either matrices or data frames.
#'
#' @section NNCT Functions:
#' Contains the functions for testing segregation/association using the NNCT. The types of the tests are cell-
#' specific tests, class-specific tests and overall tests of segregation.
#' See (\insertCite{ceyhan:stat-neer-class2009,ceyhan:exact-NNCT,ceyhan:eest-2010,ceyhan:SJScorrected2010,ceyhan:jnps-NNCT-2010,ceyhan:jkss-posthoc-2017;textual}{nnspat}).
#'
#' @section SCCT Functions:
#' Contains the functions used for testing species correspondence using the NNCT. The types are NN self and self-sum tests
#' and the overall test of species correspondence.
#' See (\insertCite{ceyhan:NNCorrespond2018;textual}{nnspat}).
#'
#' @section RCT Functions:
#' Contains the functions for testing reflexivity using the reflexivity contingency table (RCT). The types are NN 
#' self reflexivity and NN mixed-non reflexivity.
#' See (\insertCite{ceyhan:NNreflexivity2017,ceyhan:NNreflex1D2018;textual}{nnspat}).
#'
#' @section Symmetry Functions:
#' Contains the functions for testing NN symmetry using the NNCT and \eqn{Q}-symmetry contingency table. The types are NN 
#' symmetry and symmetry in shared NN structure.
#' See (\insertCite{ceyhan:SWJ-spat-sym2014;textual}{nnspat}).
#' 
#' @section Pattern Functions:
#' Contains the functions for generating and visualization of spatial patterns of segregation, association, uniformity
#' clustering and non-RL.
#' See (\insertCite{ceyhan:SiM-seg-ind2014,ceyhan:serra-2014;textual}{nnspat}).
#'
#' @references
#' \insertAllCited{}
#'
#' @docType package
#' @name nnspat
NULL
