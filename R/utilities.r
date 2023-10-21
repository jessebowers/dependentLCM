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

#' @name set_dimnames
#' @title set_dimnames
#' @description Names each axis of your array by naming each row/column axisName#
#' @param xarray The array you wish to name the axis of
#' @param axis_names The name of each axis of this array
#' @keywords internal
set_dimnames <- function(xarray, axis_names) {
  
  
  axis_names[dim(xarray)==0] <- NA
  
  dimnames(xarray) <- lapply(
    seq_along(axis_names)
    , function(iaxis) {
      iaxis_name <- axis_names[iaxis]
      if (is.na(iaxis_name)) {
        return(NULL)
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

#' @name match_lookup
#' @title match_lookup
#' @description As base::match(), but with an additional lookup argument.
#' Where 'x' values match 'table' keys, return the correspond values from 'lookup'.
#' @keywords export
match_lookup <- function(lookup, ...) {
  indexes <- match(...)
  return(lookup[indexes])
}

`%>%` <- dplyr::`%>%`
