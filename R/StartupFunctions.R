## Startup functions ------------------------------------

#' .onAttach start message
#'
#' @param libname defunct
#' @param pkgname defunct
#'
#' @return invisible()
.onAttach <- function(libname, pkgname) {
  start_message <- c( "      nnspat: functions for NN methods and\n\n",
                      "     their application for spatial patterns\n\n",
                      "    by Dr. Elvan Ceyhan <elvanceyhan@gmail.com>\n\n"
  )
  packageStartupMessage(start_message)
  invisible()
}
#'

################################################################
#'
#' .onLoad getOption package settings
#'
#' @param libname defunct
#' @param pkgname defunct
#'
#' @return invisible()
#'
#' @examples
#' getOption("nnspat.name")
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.nnspat <- list(
    #nnspat.path = "~/R-dev",
    nnspat.install.args  = "",
    nnspat.name          = "Elvan Ceyhan",
    nnspat.desc.author   = "Elvan Ceyhan <elvanceyhan@gmail.com> [aut, cre]",
    nnspat.desc.license  = "GPL-2",
    nnspat.desc.suggests = NULL,
    nnspat.desc          = list()
  )
  toset <- !(names(op.nnspat) %in% names(op))
  if (any(toset)) options(op.nnspat[toset])

  invisible()
}
