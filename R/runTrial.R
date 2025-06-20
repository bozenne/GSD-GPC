### runTrial.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 17 2025 (17:43) 
## Version: 
## Last-Updated: jun 20 2025 (10:05) 
##           By: Brice Ozenne
##     Update #: 150
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * runTrial (documentation)
##' @title Run a Group Sequential Trial
##' @description Based on a dataset and a statistical model, fit the model at interim/decision/final, evaluate the boundaries, and report the results of the trail.
##'
##' @param object output of \code{simTrial}.
##' @param data [data.frame] dataset containing the survival, toxicity, and group membership for each participant.
##' @param interim [numeric, >0, <1] timepoint at which to perform the interim analysis.
##' @param typeOfDesign method used to compute the boundaries:
##' \code{"fixed"}: fixed design, i.e., no interim analysis and only a final analysis
##' \code{"none"}: no adjustment of the boundaries, i.e., use the fixed design boundary at all stages (interim/decision/final)
##' \code{"bonferroni"}: use the same boundaries at all stages defined as the fixed design boundary bonferonni adjusted for 2 tests.
##' \code{"delayed"}: use GSD methodology for delayed response to define the boundaries: method 1, non-binding futility, specific error-spending function.
##' \code{"asP", "asOF", "asKD", "asHSD"}: use GSD methodology with deletion and no futility boundaries.
##' @param maxInfo [numeric, >0] Planned maximum information.
##' @param id.var [character] name of the variable indexing the participants.

## * runTrial (example)
##' @examples
##' library(rpact)
##' library(BuyseTest)
##' library(survival)
##' 
##' library(pbapply)
##' library(data.table)
##' 
##' ## simulate data
##' df.sim <- simTrial(n.E = 100, n.C = 200, seed = 1)
##'
##' ## analysis via Cox
##' e.Cox <- coxph(Surv(timeSurv, statusSurv) ~ group, data = df.sim)
##' runTrial(e.Cox, data = df.sim, interim = 0.5, typeOfDesign = "asOF", maxInfo = 30)
##'
##' n.sim <- 1e3
##' ls.simCox <- pblapply(1:n.sim, function(iSeed){ ## iSeed <- 1
##'   iData <- simTrial(n.E = 100, n.C = 200, seed = iSeed)
##'   iGSD <- runTrial(e.Cox, data = iData, interim = 0.5, typeOfDesign = "asOF", maxInfo = 30)
##'   cbind(seed = iSeed, design = "asOF", iGSD)
##' }, cl = 10)
##' dt.simCox <- as.data.table(do.call(rbind,ls.simCox))
##'
##' table(dt.simCox$decision)
##' dt.simCox$se[dt.simCox$]
##' dt.simCox[,.(futility=sum(decision=="futility",na.rm=TRUE),efficacy=sum(decision=="efficacy",na.rm=TRUE))]
##' 
##' ## analysis via GPC
##' e.BT <- BuyseTest(group ~ tte(timeSurv, statusSurv), data = df.sim, scoring.rule = "Gehan", trace = FALSE)
##' runTrial(e.BT, data = df.sim, interim = 0.5, typeOfDesign = "asOF", maxInfo = 550)


## * runTrial (code)
runTrial <- function(object, data, interim, typeOfDesign, maxInfo,
                     id.var = "id"){

    alpha <- 0.025

    ## ** normalize user input
    valid.gsd <- c("asP", "asOF", "asKD", "asHSD")
    typeOfDesign <- match.arg(typeOfDesign, choices = c("fixed","none","bonferroni","delayed",valid.gsd), several.ok = TRUE)
    if("delayed" %in% typeOfDesign){
        if(requireNamespace("DelayedGSD") == FALSE){
            stop("Argument method set to \"delayed\" requires to have installed the package DelayedGSD. \n",
                 "Consider running remotes::install_github(\"bozenne/DelayedGSD\"). \n")
        }
    }

    ## ** re-create data
    data.interim <- subset(data, interim =  interim)
    data.decision <- data[data[[id.var]] %in% data.interim[[id.var]],]
    data.final <- data
    
    ## ** fit and test
    if(inherits(object,"S4BuyseTest")){
        ls.args <- object@call[setdiff(names(object@call),"data")]
        ls.model <- list(interim = do.call(BuyseTest, args = c(ls.args, list(data = data.interim))),
                         decision = do.call(BuyseTest, args = c(ls.args, list(data = data.decision))),
                         final = do.call(BuyseTest, args = c(ls.args, list(data = data.final))))

        ## confint(ls.model[[1]]) ## 2*(1-pnorm(abs(atanh(-0.1479111)/(0.0649747/(1-(-0.1479111)^2)))))
        ## confint(ls.model[[1]], transform = FALSE) ## 2*(1-pnorm(abs(-0.1479111/0.0649747)))
        GSD.confint <- do.call(rbind, lapply(ls.model, function(iO){
            iO.CI <- confint(iO)
            ## add atanh transformation
            iO.CI[,"se"] <- iO.CI[,"se"]/(1-iO.CI[,"estimate"]^2)
            iO.CI[,"estimate"] <- atanh(iO.CI[,"estimate"])
            ## export
            return(iO.CI[NROW(iO.CI),c("estimate","se","lower.ci","upper.ci","p.value"),drop=FALSE])
        }))
        ## 2*(1-pnorm(abs(GSD.confint$estimate/GSD.confint$se)))
    }else{
        ls.model <- list(interim = update(object, data = data.interim),
                         decision = update(object, data = data.decision),
                         final = update(object, data = data.final))

        ## take minus the coefficient so high values are good
        GSD.confint <- as.data.frame(do.call(rbind,lapply(ls.model, function(iM){ ## iM <- ls.model[[1]]
            iOut <- c(-coef(iM)[1], sqrt(vcov(iM)[1,1]))
            c(iOut, iOut[1] + qnorm(alpha) * iOut[2], iOut[1] + qnorm(1-alpha) * iOut[2], 2*(1-pnorm(abs(iOut[1]/iOut[2]))))            
        })))
        
    }
    names(GSD.confint) <- c("estimate","se","lower","upper","p.value")
         
    GSD.confint <- cbind(stage = c("interim","decision","final"),
                         n.obs = c(NROW(data.interim),NROW(data.decision),NROW(data.final)),
                         n.event = NA,
                         n.pipeline = c(sum(data.interim$pipeline),sum(data.decision$pipeline),sum(data.final$pipeline)),
                         info.pc = (1/GSD.confint$se^2)/maxInfo,
                         GSD.confint)
    GSD.confint$n.event <- list(c(sum(data.interim$statusSurv), sum(data.interim$statusTox)),
                                c(sum(data.decision$statusSurv), sum(data.decision$statusTox)),
                                c(sum(data.final$statusSurv), sum(data.final$statusTox)))

    ## ** update boundaries and take decision
    out <- do.call(rbind,lapply(typeOfDesign, function(iType){ ## iType <- "none"
        iOut <- cbind(design = iType, GSD.confint, bound = .updateBound(GSD.confint, alpha = alpha, typeOfDesign = iType, maxInfo = maxInfo))
        return(.updateDecision(iOut, typeOfDesign = iType))
    }))

    ## ** export
    rownames(out) <- NULL
    return(out)
}

## * helper
## ** updateBound
.updateBound <- function(object, alpha, typeOfDesign, maxInfo){

    ## output
    bound <- setNames(rep(as.numeric(NA),3), c("interim","decision","final"))
    
    ## update
    if(typeOfDesign == "fixed"){
        bound["interim"] <- Inf
        bound["decision"] <- NA
        bound["final"] <- qnorm(1-alpha)
    }else if(typeOfDesign == "none"){
        bound[] <- qnorm(1-alpha)
    }else if(typeOfDesign == "bonferroni"){
        bound[] <- qnorm(1-alpha/2)
    }else if(object["interim","info.pc"]>1){ ## over-run: go to decision
        bound["interim"] <- -Inf
        bound["decision"] <- qnorm(1-alpha)
        bound["final"] <- NA
    }else if(typeOfDesign == "delayed"){

        if(object["interim","info.pc"] < object["decision","info.pc"]){

            ## spend all alpha at final (evaluated regardless to the information at final, to over or under-running)
            gsd <- DelayedGSD::CalcBoundaries(kMax = 2,  
                                              alpha = alpha, 
                                              beta = 0.2,  
                                              InfoR.i = object["interim","info.pc"],  
                                              rho_alpha = 2,  
                                              rho_beta = 2,  
                                              method = 1, 
                                              cNotBelowFixedc = FALSE, 
                                              InfoR.d = c(object["decision","info.pc"],1),   
                                              bindingFutility = FALSE,
                                              alternative = "greater")        
            
        }else{ ## decreasing information between interim and decision

            ## put decision after interim in term of information while keeping the information ratio between the two
            ## should ensure balanced reversal probability
            ## spend all alpha at final (evaluated regardless to the information at final, to over or under-running)
            gsd <- DelayedGSD::CalcBoundaries(kMax = 2,  
                                              alpha = alpha, 
                                              beta = 0.2,  
                                              InfoR.i = object["interim","info.pc"],  
                                              rho_alpha = 2,  
                                              rho_beta = 2,  
                                              method = 1, 
                                              cNotBelowFixedc = FALSE, 
                                              InfoR.d = c(max(object["interim","info.pc"]*(object["interim","info.pc"]/object["decision","info.pc"]),1),1),   
                                              bindingFutility = FALSE,
                                              alternative = "greater")

            ## Sigma <- diag(1, 3)
            ## Sigma[1,2] <- Sigma[2,1] <- sqrt(object["decision","info.pc"]/object["interim","info.pc"])
            ## Sigma[2,3] <- Sigma[3,2] <- sqrt(object["interim","info.pc"]/1)
            ## Sigma[1,3] <- Sigma[3,1] <- sqrt(object["decision","info.pc"]/1)
            ## reversal1 <- pmvnorm(lower = c(gsd$planned$uk, -Inf), upper = c(Inf, gsd$planned$ck[1]), mean = c(0,0), sigma = Sigma[1:2,1:2]) 
            ## reversal2 <- pmvnorm(lower = c(-Inf, gsd$planned$ck[1]), upper = c(gsd$planned$lk, Inf), mean = c(0,0), sigma = Sigma[1:2,1:2])
            ## reversal1 - reversal2
        }

        bound["interim"] <- gsd$planned$uk
        bound[c("decision","final")] <- gsd$planned$ck
    
        ## bound["decision"] - calc_ck(uk=bound["interim"],
        ##                             lk=gsd$planned$lk,
        ##                             Info.i = object["interim","info.pc"],
        ##                             Info.d = object["decision","info.pc"],
        ##                             Info.max = maxInfo,
        ##                             ImaxAnticipated = object["interim","info.pc"]>1,
        ##                             rho_alpha=NA,
        ##                             alpha=alpha,
        ##                             bindingFutility = FALSE)
    
    }else{

        ## spend all alpha at final (evaluated regardless to the information at final, to over or under-running)
        gsd <- getDesignGroupSequential(kMax = 2,
                                        alpha = alpha,
                                        sided = 1,
                                        informationRates = c(object["interim","info.pc"],1),
                                        typeOfDesign = typeOfDesign)
        bound[c("interim","final")] <- gsd$criticalValues

        if(is.infinite(bound["interim"]) && bound["interim"]>0){ ## information close to 0 at interim: continue at interim and no decision
            bound["decision"] <- as.numeric(NA)
        }else{
            ## Whitehead (Cont. Clin.Trials, 1992) proposed the "deletion method".
            ## The analysis k at which termination occurs is deleted
            ## and one behaves as if analysis k had occurred with the information level \tilde{I}k arising from the final set of responses.
            if(object["decision","info.pc"] >= 1){
                bound["decision"] <- qnorm(1-alpha)
            }else{
                bound["decision"] <- getDesignGroupSequential(kMax = 2,
                                                              alpha = alpha,
                                                              sided = 1,
                                                              informationRates = c(object["decision","info.pc"],1),
                                                              typeOfDesign = typeOfDesign)$criticalValues[1]
            }
        }
    }

    ## typeBetaSpending <- switch(typeOfDesign,
    ##                            "P" = "bsP", "asP" = "bsP",
    ##                            "OF" = "bsOF", "asOF" = "bsOF",
    ##                            "asKD" = "bsKD",
    ##                            "asHSD" = "bsHSD",
    ##                            NA)
    ## if(is.na(typeOfDesign)){
    ##     stop("Argument \'typeOfDesign\' does not match available types for non-binding futility bounds. \n",
    ##          "Possible values: \"P\", \"asP\", \"OF\", \"asOF\", \"asKD\", \"asHSD\". \n")
    ## }
    ## if(object["interim","info.pc"]>object["final","info.pc"] || object["final","info.pc"] >= 1){
    ##     gsd <- getDesignGroupSequential(kMax = 2,
    ##                                     alpha = alpha,
    ##                                     sided = 1,
    ##                                     beta = 0.2, bindingFutility = FALSE, typeBetaSpending = typeBetaSpending,
    ##                                     informationRates = c(object["interim","info.pc"],1),
    ##                                     typeOfDesign = typeOfDesign)
    ## }else{
    ##     gsd <- getDesignGroupSequential(kMax = 2,
    ##                                     alpha = alpha,
    ##                                     sided = 1,
    ##                                     beta = 0.2, bindingFutility = FALSE, typeBetaSpending = typeBetaSpending,
    ##                                     informationRates = object[c("interim","final"),"info.pc"],
    ##                                     typeOfDesign = typeOfDesign)
    ## }
    ## bound[c("interim","final")] <- gsd$criticalValues
    ## bound["decision"] <- calc_ck(uk=bound["interim"],
    ##                              lk=gsd$futilityBounds,
    ##                              Info.i = object["interim","info.pc"],
    ##                              Info.d = object["decision","info.pc"],
    ##                              Info.max = maxInfo,
    ##                              ImaxAnticipated = object["interim","info.pc"]>1,
    ##                              rho_alpha=NA,
    ##                              alpha=alpha,
    ##                              bindingFutility = FALSE)
                

    ## export
    return(bound)

}

## ** .updateDecision
.updateDecision <- function(object, typeOfDesign){

    object$statistic <- object$estimate/object$se
    object$decision <- NA
    object$reversal <- NA
    if(typeOfDesign == "fixed"){
        object["interim","decision"] <- "continue"
        object["final","decision"] <- c("futility","efficacy")[1+(object["final","statistic"] >= object["final","bound"])]
    }else{
        if(typeOfDesign %in% c("none","bonferroni")){
            object["interim","decision"] <- c("continue","stop")[1+(object["interim","statistic"] >= object["interim","bound"])]
        }else{
            object["interim","decision"] <- ifelse(object["interim","info.pc"]>1,
                                                   "infoMax",
                                                   c("continue","stop")[1+(object["interim","statistic"] >= object["interim","bound"])])
        }
        object["decision","decision"] <- ifelse(object["interim","decision"] == "continue",
                                                NA,
                                                c("futility","efficacy")[1+(object["decision","statistic"] >= object["decision","bound"])])
        object["final","decision"] <- ifelse(object["interim","decision"] == "continue",
                                             c("futility","efficacy")[1+(object["final","statistic"] >= object["final","bound"])],
                                             NA)
        object["decision","reversal"] <- ifelse(object["interim","decision"] == "continue",
                                                NA,
                                                (object["interim","decision"]=="stop") & (object["decision","decision"]=="futility"))
    }

    return(object)

}
##----------------------------------------------------------------------
### runTrial.R ends here
