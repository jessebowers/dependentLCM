##############
############## UTILITIES
##############

CAST_FUN <- list("integer"=as.integer, "double"=as.numeric)

#' Calculate mode (most common value)
#' @keywords internal
getMode <- function(x) {
  xtype <- typeof(x)
  xcounts <- table(x)
  xmode <- names(xcounts)[which.max(xcounts)]
  xmode <- CAST_FUN[[xtype]](xmode) 
  return(xmode)
}


#' Same as log(sum(exp(x))), but adjusted to improve precision.
#' Assumes x are logged values. Calculates sum(e^x) and then converts back to log scale
#' @param x numeric vector in log scale
#' @keywords internal
expSumLog <- function(x) {
  # Handle precision: log(e^a+e^b+e^c) = log(a * (1+e^(b-a)+e^(c-1))) = log(a) + log(1+e^(b-a)+e^(c-1)))
  xmax <- max(x)
  return(xmax + log(sum(exp(x-xmax))))
}

#' @name set_dimnames
#' @title set_dimnames
#' @description Names each axis of your array by naming each row/column axisName#
#' @param xarray The array you wish to name the axis of
#' @param axis_names The name of each axis of this array
#' @keywords internal
set_dimnames <- function(xarray, axis_names) {
  dimnames(xarray) <- lapply(
    seq_along(axis_names)
    , function(iaxis) {
      iaxis_name <- axis_names[iaxis]
      if (is.null(iaxis_name)) {
        return(null)
      }
      paste0(
        iaxis_name
        , seq_len(dim(xarray)[iaxis])
      )
    }
  )
  return(xarray)
}


#' @name sample.df
#' @title sample.df
#' @description Samples rows from a dataframe or matrix
#' @param x The dataframe or matrix we are sampling from
#' @param size Number of rows to sample
#' @param ... Other options per sample.int()
#' @export
sample.df <- function(x, size, ...) {
  return(x[sample.int(n=nrow(x), size=size, ...),])
}

#' @name sample.integers
#' @title sample.integers
#' @description Samples values from an provided integer vector.
#' Safer version of sample() which has bad behavior if length(x)==1.
#' @param x integer vector we are sampling from
#' @param size Number of rows to sample
#' @param ... Other options per sample.int()
#' @keywords internal
sample.integers <- function(x, size, ...) {
  return(x[sample.int(n=length(x), size=size, ...)])
}

`%>%` <- magrittr::`%>%`
