### subset.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 12 2025 (10:35) 
## Version: 
## Last-Updated: jun 18 2025 (16:52) 
##           By: Brice Ozenne
##     Update #: 31
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * subset (documentation)
##' @title Interim Analysis Data
##' @description Extract the data corresponding to an interim analysis.
##' 
##' @param object output of \code{simTrial}
##' @param interim [numeric, >0, <1] timepoint up to which the data is restriected.

## * subset (code)
subset.simTrial <- function(object, interim){

    tol <- 1e-12

    ## ** normalize user input
    admin.censoring <- attr(object, "admin.censoring")

    if(is.null(interim)){
        if(!is.null(attr(object,"interim"))){
            interim <- attr(object,"interim")
        }else{
            interim <- Inf
        }
    }else if(length(interim)!=1){
        stop("Argument \'interim\' should have length 1. \n")
    }
    if(!is.numeric(interim) || interim <= 0){
        stop("Argument \'interim\' should be a strictly positive numeric value. \n")
    }

    if(!is.null(attr(object,"interim"))){
        if(attr(object,"interim")<interim){
            stop("Cannot apply the subset function when the data has already been subset at an earlier interim. \n")
        }
    }

    ## ** subset
    out <- object[object$timeInclusion < interim,,drop=FALSE]
    out$pipeline <- ((out$timeInclusion + out$timeSurv) > interim)
    out$statusSurv <- out$statusSurv * ((out$timeInclusion + out$timeSurv) <= interim)
    out$timeSurv <- pmin(out$timeInclusion + out$timeSurv, interim) - out$timeInclusion
    out$statusTox <- out$statusTox * ((out$timeInclusion + out$timeTox) <= interim)
    out$timeTox <- pmin(out$timeInclusion + out$timeTox, interim) - out$timeInclusion
    attr(out,"interim") <- interim

    ## ** add type
    out$state <- factor("alive", levels = c("alive","alive & tox","dead","dead & tox","tox & dropout","dropout"))
    out$state[out$statusSurv==1 & out$statusTox==1] <- "dead & tox"
    out$state[out$statusSurv==1 & out$statusTox==0] <- "dead"
    out$state[out$statusSurv==0 & out$statusTox==1 & out$timeSurv + tol >= admin.censoring] <- "alive & tox"
    out$state[out$statusSurv==0 & out$statusTox==1 & out$timeSurv + tol < admin.censoring] <- "tox & dropout"
    out$state[out$statusSurv==0 & out$statusTox==0 & out$timeSurv + tol < admin.censoring] <- "dropout"

    ## ** export
    return(out)

}


##----------------------------------------------------------------------
### subset.R ends here
