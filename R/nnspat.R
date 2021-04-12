#' nnspat: A package for NN Methods and Their Use in Testing Spatial Patterns
#'
#' nnspat is a package for computation of spatial pattern tests based on NN relations and 
#' generation of various spatial patterns.
#'
#' The nnspat package contains the functions for segregation/association tests based on nearest neighbor contingency
#' tables (NNCTs), and tests for species correspondence, NN symmetry and reflexivity based on the corresponding
#' contingency tables and functions for generating patterns of segregation, association, uniformness and various
#' non-random labeling (RL) patterns for disease clustering for data in two (or more) dimensions.
#' (see ???(\insertCite{ceyhan:Phd-thesis;textual}{nnspat}).
#'
#' @section The nnspat functions:
#' The nnspat functions can be grouped as Auxiliary Functions, NNCT Functions, SCCT Functions, RCT Functions
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
#'
#' @section SCCT Functions:
#' Contains the functions used for testing species correspondence using the NNCT. The types are NN self and self-sum tests
#' and the overall test of species correspondence.
#'
#' @section RCT Functions:
#' Contains the functions for testing reflexivity using the reflexivity contingency table (RCT). The types are NN 
#' self reflexivity and NN mixed-non reflexivity.
#'
#' @section Symmetry Functions:
#' Contains the functions for testing NN symmetry using the NNCT and Qsymmetry contingency table. The types are NN 
#' symmetry and symmetry in shared NN structure.
#' 
#' @section Pattern Functions:
#' Contains the functions for generating and visualization of spatial patterns of segregation, association, uniformness
#' clustering and non-RL.
#'
#' @references
#' \insertAllCited{}
#'
#' @docType package
#' @name nnspat
NULL
