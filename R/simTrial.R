### simTrial.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun  3 2025 (09:38) 
## Version: 
## Last-Updated: jun 18 2025 (16:52) 
##           By: Brice Ozenne
##     Update #: 352
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * simTrial (documentation)
##' @title Simulate Bivariate Time-to-event for 2-stage GSD
##' @description Simulate data corresponding to a two-arm group sequential design with a single interim analysis and two time to event outcomes: one is survival and the other is toxicity.
##' The time to event outcomes may be right-censored due to patient drop-out.
##' 
##' @param n.E [integer] sample size in the experimental group
##' @param n.C [integer] sample size in the control group
##' @param shape.recruitment.C,shape.recruitment.E [numeric vector of length 2] shape of the distribution used to model the inclusion of the patients in the control/experimental group.
##' When \code{NA},  deterministic and equally spread recruitment time between 0 and 1.
##' Otherwise beta distributed recruitment time - \code{c(1,1)} corresponds to stochastic inclusion based on a uniform distribution.
##' @param scale.C,scale.E [numeric vector/list of length 2] first parameter of the survival and toxicity distribution in the control/experimental group.
##' @param shape.C,shape.E [numeric vector/list of length 2] second parameter of the survival and toxicity distribution in the control/experimental group.
##' @param dist.C,dist.E [character] distribution used to generate the survival and toxicity time to events in the control/experimental group.
##' Can be a Weibull distribution (\code{"weibull"}), a uniform distribution (\code{"uniform"}), or based on piecewise constant hazards (\code{"piecewiseExp"}).
##' @param rho.C,rho.E [numeric] correlation coefficient used to parametrize the Gaussian copula inducing correlation between the survival and toxicity time to events in the control/experimental group.
##' @param scale.censoring.E,scale.censoring.E [numeric] first parameter of the censoring distribution in the control/experimental group
##' @param shape.censoring.C,shape.censoring.E [numeric] second parameter of the censoring distribution in the control/experimental group
##' @param dist.censoring.C,dist.censoring.E [character] distribution used to generate the censoring time in the control/experimental group.
##' @param  [character] distribution used to generate the censoring time in the experimental group.
##' Can be a Weibull distribution (\code{"weibull"}), a uniform distribution (\code{"uniform"}), or based on piecewise constant hazards (\code{"piecewiseExp"}).
##' @param admin.censoring [numeric, >0] maximum possible follow-up time for each participant.
##' @param seed [integer, >0] set the initial state of the random number generator (for reproducibility).
##' @param latent [logical] should the censoring times and uncensored event time be output?
##'
##' @details Time to event can either be Weibull distributed (scale,shape), uniform distributed (min,max), or with piecewise constant hazard (hazard, cutpoints).

## * simTrial (example)
##' @examples
##' 
##' ## simulate data
##' df.sim <- simTrial(n.E = 100, n.C = 200, seed = 1, latent = TRUE)
##'
##' ## display data
##' df.sim
##' summary(df.sim)
##' plot(df.sim, type = "lexis")
##' plot(df.sim, facet = ~group, labeller = "label_both")
##' plot(df.sim, facet = ~state, labeller = "label_both")
##' plot(df.sim, facet = group~sntate, labeller = "label_both")
##' 
##' plot(df.sim, type = "recruitment")
##' plot(df.sim, type = "survival")
##' plot(df.sim, type = "toxicity")
##'
##' ## interim analysis
##' subset(df.sim, interim = 0.75)
##' plot(df.sim, interim = 0.75)
##' plot(df.sim, interim = 0.75, facet =  ~state)
##' plot(df.sim, interim = 0.75, facet =  pipeline~state)
##' plot(df.sim, interim = 0.75, type = "recruitment")
##'
##' ## same data without latent times
##' dfObs.sim <- simTrial(n.E = 50, n.C = 100, seed = 1, latent = FALSE)
##' 
##' plot(dfObs.sim, type = "survival")
##' plot(dfObs.sim, type = "toxicity")
##'
##'
##' ## add correlation between endpoints
##' dfCor.sim <- simTrial(n.E = 100, n.C = 200, seed = 1, latent = TRUE, rho.C = -0.5, rho.E = 0.5)
##'
##' ggCor <- ggplot(dfCor.sim, aes(x = eta.timeTox, y = eta.timeSurv))
##' ggCor <- ggCor + geom_point() + facet_wrap(~group)
##' ggCor
##' 
##' plot(dfCor.sim, facet = group~state, labeller = label_both)

## * simTrial (code)
simTrial <- function(n.C, shape.recruitment.C = c(NA,NA),
                     scale.C = c(0.75,0.5), shape.C = c(3,0.75), dist.C = c("weibull", "weibull"), rho.C = 0,
                     scale.censoring.C = 1, shape.censoring.C = 1, dist.censoring.C = "weibull",
                     n.E = NULL, shape.recruitment.E = NULL,
                     scale.E = NULL, shape.E = NULL, dist.E = NULL, rho.E = NULL,
                     scale.censoring.E = NULL, shape.censoring.E = NULL, dist.censoring.E = NULL,
                     admin.censoring = 0.75, seed = NULL, latent = FALSE){

    tol <- 1e-12

    ## ** normalize user input
    name.inclusion <- "timeInclusion"
    name.time <- c("timeSurv","timeTox")
    name.etatime <- paste0("eta.",name.time)
    name.censtime <- "timeCens"
    name.status <- c("statusSurv","statusTox")
    
    ## set experimental arm parameter to control arm parameter when missing
    if(is.null(n.E)){n.E <- n.C}
    if(is.null(shape.recruitment.E)){shape.recruitment.E <- shape.recruitment.C}
    if(is.null(scale.E)){scale.E <- scale.C}
    if(is.null(shape.E)){shape.E <- shape.C}
    if(is.null(dist.E)){dist.E <- dist.C}
    if(is.null(rho.E)){rho.E <- rho.C}
    if(is.null(scale.censoring.E)){scale.censoring.E <- scale.censoring.C}
    if(is.null(shape.censoring.E)){shape.censoring.E <- shape.censoring.C}
    if(is.null(dist.censoring.E)){dist.censoring.E <- dist.censoring.C}

    ## distribution of the recruitment times
    shape.recruitment.C <- check_recruitmentDist(value = shape.recruitment.C, name = "shape.recruitment.C")
    shape.recruitment.E <- check_recruitmentDist(value = shape.recruitment.E, name = "shape.recruitment.E")

    ## outcome distribution
    scale.C <- check_outcomeScale(scale.C, name = "scale.C")
    shape.C <- check_outcomeShape(shape.C, name = "shape.C")
    dist.C <- check_outcomeDist(dist.C, name = "dist.C")
    rho.C <- check_outcomeCor(rho.C, name = "rho.C")

    scale.E <- check_outcomeScale(scale.E, name = "scale.E")
    shape.E <- check_outcomeShape(shape.E, name = "shape.E")
    dist.E <- check_outcomeDist(dist.E, name = "dist.E")
    rho.E <- check_outcomeCor(rho.E, name = "rho.E")

    ## censoring distribution
    scale.censoring.C <- check_censoringScale(scale.censoring.C, name = "scale.censoring.C")
    shape.censoring.C <- check_censoringShape(shape.censoring.C, name = "shape.censoring.C")
    dist.censoring.C <- check_censoringDist(dist.censoring.C, name = "dist.censoring.C")

    scale.censoring.E <- check_censoringScale(scale.censoring.E, name = "scale.censoring.E")
    shape.censoring.E <- check_censoringShape(shape.censoring.E, name = "shape.censoring.E")
    dist.censoring.E <- check_censoringDist(dist.censoring.E, name = "dist.censoring.E")

    ## ** data generating mechanism
    if(!is.null(seed)){
        set.seed(seed)
    }
    data.C <- data.frame(group = rep("control",n.C))
    data.E <- data.frame(group = rep("experimental",n.E))

    ## inclusion time
    if(is.null(shape.recruitment.C)){
        data.C[[name.inclusion]] <- seq(from = 0, to = 1, length.out = n.C)
    }else{
        data.C[[name.inclusion]] <- sort(rbeta(n.C, shape1 = shape.recruitment.C[1], shape2 = shape.recruitment.C[2]))
    }
    if(is.null(shape.recruitment.E)){
        data.E[[name.inclusion]] <- seq(from = 0, to = 1, length.out = n.E)
    }else{
        data.E[[name.inclusion]] <- sort(rbeta(n.E, shape1 = shape.recruitment.E[1], shape2 = shape.recruitment.E[2]))
    }

    ## censoring distribution
    if(all(is.infinite(scale.censoring.C))){
        data.C[[name.censtime]] <- Inf
    }else{
        data.C[[name.censtime]] <- switch(dist.censoring.C,
                                          "uniform" = runif(n.C, min = scale.censoring.C, max = shape.censoring.C),
                                          "weibull" = rweibull(n.C, scale = scale.censoring.C, shape = shape.censoring.C),
                                          "piecewiseExp" = rexppiecewise(n.C, rate = 1/scale.censoring.C, breaks = shape.censoring.C))
    }
    if(all(is.infinite(scale.censoring.E))){
        data.E[[name.censtime]] <- Inf
    }else{
        data.E[[name.censtime]] <- switch(dist.censoring.E,
                                          "uniform" = runif(n.E, min = scale.censoring.E, max = shape.censoring.E),
                                          "weibull" = rweibull(n.E, scale = scale.censoring.E, shape = shape.censoring.E),
                                          "piecewiseExp" = rexppiecewise(n.E, rate = 1/scale.censoring.E, breaks = shape.censoring.E))
    }

    ## outcome distribution
    if(abs(rho.C)<0){
        for(iterO in 1:2){
            data.C[[name.etatime[iterO]]] <- switch(dist.C[iterO],
                                                    "uniform" = runif(n.C, min = scale.C[[iterO]], max = shape.C[[iterO]]),
                                                    "weibull" = rweibull(n.C, scale = scale.C[[iterO]], shape = shape.C[[iterO]]),
                                                    "piecewiseExp" = rexppiecewise(rate = 1/scale.C[[iterO]], breaks = shape.C[[iterO]]))
        }
    }else{ ## Gaussian copula
        requireNamespace("mvtnorm")        
        Sigma.C <- diag(1-rho.C,2,2) + rho.C
        Mnorm.C <- mvtnorm::rmvnorm(n.C, mean = rep(0,2), sigma = Sigma.C)
        Munif.C <- matrix(pnorm(Mnorm.C), nrow = n.C, ncol = 2)
        ## qnorm(Munif.C[1,]) - Mnorm.C[1,]
        for(iterO in 1:2){
            data.C[[name.etatime[iterO]]] <- switch(dist.C[iterO],
                                                    "uniform" = qunif(Munif.C[,iterO], min = scale.C[[iterO]], max = shape.C[[iterO]]),
                                                    "weibull" = qweibull(Munif.C[,iterO], scale = scale.C[[iterO]], shape = shape.C[[iterO]]),
                                                    "piecewiseExp" = qexppiecewise(Munif.C[,iterO], rate = 1/scale.C[[iterO]], breaks = shape.C[[iterO]]))
        }
    }
   
    if(abs(rho.E)<0){
        for(iterO in 1:2){
            data.E[[name.etatime[iterO]]] <- switch(dist.E[iterO],
                                                    "uniform" = runif(n.E, min = scale.E[[iterO]], max = shape.E[[iterO]]),
                                                    "weibull" = rweibull(n.E, scale = scale.E[[iterO]], shape = shape.E[[iterO]]),
                                                    "piecewiseExp" = rexppiecewise(rate = 1/scale.E[[iterO]], breaks = shape.E[[iterO]]))
        }
    }else{
        requireNamespace("mvtnorm")        
        Sigma.E <- diag(1-rho.E,2,2) + rho.E
        Mnorm.E <- mvtnorm::rmvnorm(n.E, mean = rep(0,2), sigma = Sigma.E)
        Munif.E <- matrix(pnorm(Mnorm.E), nrow = n.E, ncol = 2)
        ## qnorm(Munif.E[1,]) - Mnorm.E[1,]
        for(iterO in 1:2){
            data.E[[name.etatime[iterO]]] <- switch(dist.E[iterO],
                                                    "uniform" = qunif(Munif.E[,iterO], min = scale.E[[iterO]], max = shape.E[[iterO]]),
                                                    "weibull" = qweibull(Munif.E[,iterO], scale = scale.E[[iterO]], shape = shape.E[[iterO]]),
                                                    "piecewiseExp" = qexppiecewise(Munif.E[,iterO], rate = 1/scale.E[[iterO]], breaks = shape.E[[iterO]]))
        }
    }

    ## observed (right-censored) outcome distribution
    data.C[[name.time[1]]] <- pmin(data.C[[name.etatime[1]]], data.C[[name.censtime]], admin.censoring)
    data.C[[name.status[1]]] <- as.numeric((data.C[[name.etatime[1]]] <= data.C[[name.censtime]]) & (data.C[[name.etatime[1]]] <= admin.censoring))

    data.C[[name.time[2]]] <- pmin(data.C[[name.etatime[1]]], data.C[[name.etatime[2]]], data.C[[name.censtime]], admin.censoring)
    term1 <- (data.C[[name.etatime[2]]] <= data.C[[name.etatime[1]]]) & (data.C[[name.etatime[2]]] <= data.C[[name.censtime]]) & (data.C[[name.etatime[2]]] <= admin.censoring)
    term2 <- (data.C[[name.etatime[1]]] < data.C[[name.etatime[2]]]) & (data.C[[name.etatime[1]]] <= data.C[[name.censtime]]) & (data.C[[name.etatime[1]]] <= admin.censoring)
    data.C[[name.status[2]]] <- term1 + 2*term2

    data.E[[name.time[1]]] <- pmin(data.E[[name.etatime[1]]], data.E[[name.censtime]], admin.censoring)
    data.E[[name.status[1]]] <- as.numeric((data.E[[name.etatime[1]]] <= data.E[[name.censtime]]) & (data.E[[name.etatime[1]]] <= admin.censoring))

    data.E[[name.time[2]]] <- pmin(data.E[[name.etatime[1]]], data.E[[name.etatime[2]]], data.E[[name.censtime]], admin.censoring)
    term1 <- (data.E[[name.etatime[2]]] <= data.E[[name.etatime[1]]]) & (data.E[[name.etatime[2]]] <= data.E[[name.censtime]]) & (data.E[[name.etatime[2]]] <= admin.censoring)
    term2 <- (data.E[[name.etatime[1]]] < data.E[[name.etatime[2]]]) & (data.E[[name.etatime[1]]] <= data.E[[name.censtime]]) & (data.E[[name.etatime[1]]] <= admin.censoring)
    data.E[[name.status[2]]] <- term1 + 2*term2

    ## ** collect data
    out <- rbind(data.C, data.E)
    out.order <- cbind(id = 1:NROW(out),out[order(out$timeInclusion),])
    if(latent==FALSE){
        out.order[name.censtime] <- NULL
        out.order[name.etatime[1]] <- NULL
        out.order[name.etatime[2]] <- NULL
    }
    
    ## ** export
    class(out.order) <- append("simTrial",class(out.order))
    attr(out.order,"seed") <- seed
    attr(out.order,"admin.censoring") <- admin.censoring
    return(out.order)
}


## * helper
## ** check_recruitment
check_recruitmentDist <- function(value, name){

    if(length(value)==0 || any(is.na(value))){
        value <- NULL
    }else if(length(value)==1){
        if(is.numeric(value)){
            value <- rep(value,2)
        }else{
            stop("Argument \'",name,"\' should either be a numeric vector of length 2 (beta distributed inclusion time), \n",
                 "or NULL (deterministic equally spread inclusion times). \n")
        }
    }else if(length(value)==2){
        if(!is.numeric(value)){
            stop("Argument \'",name,"\' should either be a numeric vector of length 2 (beta distributed inclusion time), \n",
                 "or NULL (deterministic equally spread inclusion times). \n")
        }
    }else if(length(value)>2){
        stop("Argument \'",name,"\' should have length at most 2. \n")
    }

    return(value)
}

## ** check_outcome
check_outcomeScale <- function(value, name){
    if(length(value) != 2){
        stop("Argument \'",name,"\' should have length 2. \n")
    }
    return(value)
}
check_outcomeShape <- function(value, name){
    if(length(value) != 2){
        stop("Argument \'",name,"\' should have length 2. \n")
    }
    return(value)
}
check_outcomeDist <- function(value, name){

    if(length(value) != 2){
        stop("Argument \'",name,"\' should have length 2. \n")
    }
    value[1] <- match.arg(value[1], choices = c("weibull", "uniform","piecewiseExp"))
    value[2] <- match.arg(value[2], choices = c("weibull", "uniform","piecewiseExp"))

    return(value)
}
check_outcomeCor <- function(value, name){
    if(length(value) != 1){
        stop("Argument \'",name,"\' should have length 1. \n")
    }
    if(abs(value) >= 1){
        stop("Argument \'",name,"\' should correspond to a correlation: numeric, strictly between -1 and 1. \n")
    }
    return(value)
}

## ** check_censoring
check_censoringScale <- function(value, name){
    if(length(value) != 1){
        stop("Argument \'",name,"\' should have length 1. \n")
    }
    return(value)
}
check_censoringShape <- function(value, name){
    if(length(value) != 1){
        stop("Argument \'",name,"\' should have length 1. \n")
    }
    return(value)
}
check_censoringDist <- function(value, name){

    if(length(value) != 1){
        stop("Argument \'",name,"\' should have length 1. \n")
    }
    value <- match.arg(value, choices = c("weibull", "uniform","piecewiseExp"))

    return(value)
}

## ** rexppiecewise
##' @examples
##' ## EXAMPLE 1
##' rate1 <- c(0.1,0.2)
##' 
##' set.seed(1) 
##' res1 <- rexppiecewise(5, rate = rate1, breaks = 10)
##' res1
##'
##' set.seed(1)
##' E <- rexp(5)
##' (E <= 10*rate1[1])*E/rate1[1] + (E > 10*rate1[1])*(E+10*rate1[1])/rate1[2]
##'
##' ## EXAMPLE 2 
##' set.seed(2)
##' res2 <- rexppiecewise(1e4, rate = c(1,0.01,0.5,2), breaks = 1:3, method = "single")
##' library(survival)
##' e.KM <- survfit(Surv(time,status)~1, data = data.frame(time=res2,status=1))
##' df.KM <- data.frame(cumhaz = e.KM$cumhaz,
##'                     time = e.KM$time,
##'                     time0 = pmin(e.KM$time,1),
##'                     time1 = pmin(pmax(e.KM$time-1,0),1),
##'                     time2 = pmin(pmax(e.KM$time-2,0),1),
##'                     time3 = pmax(e.KM$time-3,0)
##' )
##' plot(x = df.KM$time, y = df.KM$cumhaz, xlim = c(0,4))
##' lm(cumhaz ~ 0 + time0 + time1 + time2 + time3, data = df.KM)
##' 
##' ## EXAMPLE 3
##' rate3 <- c(0.1,0.1,0.1)
##' 
##' set.seed(3) 
##' res3 <- rexppiecewise(5, rate = rate3, breaks = c(1,10))
##' res3
##'
##' set.seed(3)
##' E <- rexp(5, rate = unique(rate3))
##' E
rexppiecewise <- function(n, rate, breaks, method = "single"){

    method <- match.arg(method, c("multiple","single"))
    n.rate <- length(rate)
    diff.breaks <- diff(c(0,breaks,Inf))
        
    if(method == "multiple"){
        ## idea: simulate an exponential with each of the rate 
        ##     : if the simulated value exceed the interval width move to the next time interval
        ##     : stop at the first simulated value below the interval width and add the previous interval widths to get the output

        M.simExp <- matrix(Inf, ncol = n.rate, nrow = n)
        M.simExp[,rate>0] <- do.call(cbind,lapply(rate[rate>0], rexp, n = n))

        ls.out <- apply(M.simExp, MARGIN = 1, FUN = function(iRow){ ## iRow <- M.simExp[3,]
            iIndex <- which(iRow<diff.breaks)[1]
            c(0,breaks)[iIndex] + iRow[iIndex]
        }, simplify = FALSE)

        out <- unlist(ls.out)
    }else if(method == "single"){
        ## idea: An exponential variable has survival function S(t)=P[E>t]=exp(-t)
        ##     : so P[A^-1(E)>t] = P[E>A(t)] = exp(-A(t))
        ##     : so applying the reciprocal of the cumulative hazard function to simulated exponential function
        ##     : will lead to a random variable with the appropriate survival function
        ## A(t) = \sum_{k, t_k < t} \alpha_k \Delta_k + \alpha_K (t-t_K)
        ## therefore t = (A(t) - \sum_{k, t_k < t} \alpha_k \Delta_k) / alpha_K + t_K

        E <- rexp(n)

        A.min <- cumsum(c(0,rate*diff.breaks)[1:n.rate])
        UA.min <- unique(A.min)
        index.A <- cut(E, breaks = c(UA.min,Inf))

        out <- (E - UA.min[as.numeric(index.A)])/rate[rate>0][index.A] + c(0,breaks)[rate>0][index.A]
    }

    return(out)
}

## ** pexppiecewise
##' @examples
##' ## EXAMPLE 1
##' rate1 <- c(0.1,0.2)
##' q1 <- seq(0,15,by=0.1)
##' 
##' res1 <- pexppiecewise(q1, rate = rate1, breaks = 10)
##' res1
##' all(diff(res1)>0)
##' res1 - (pexp(q1, rate = rate1[1])*(q1<=10) + (1 - pexp(10, rate = rate1[1], lower.tail = FALSE)*pexp(q1 - 10, rate = rate1[2], lower.tail = FALSE))*(q1>10))
##' 
##' ## EXAMPLE 3
##' rate3 <- c(0.1,0.1,0.1)
##' q3 <- c(0,0.1,1.1,10.1)
##' 
##' res3 <- pexppiecewise(q = q3, rate = rate3, breaks = c(1,10))
##' res3
##' res3 - pexp(q = q3, rate = unique(rate3))
##' 
pexppiecewise <- function(q, rate, breaks){

    tol <- 1e-12    
    n.rate <- length(rate)
    diff.breaks <- diff(c(0,breaks,Inf))

    ## evaluate the cumulative hazard at the end of each interval (except the last interval, Inf)
    A.min <- cumsum(c(0,rate*diff.breaks)[1:n.rate])
    ## find the interval corresponding to the timepoint
    index.A <- cut(q, breaks = c(-tol,breaks,Inf))
    indexNum.A <- as.numeric(index.A)
    ## evaluate the cumulative hazard for each timepoint
    A.q <- A.min[indexNum.A] + rate[indexNum.A] * (q - c(0,breaks)[indexNum.A])
    out <- 1 - exp(-A.q)
    return(out)
}

## ** qexppiecewise
##' @examples
##' ## EXAMPLE 1
##' rate1 <- c(0.1,0.2)
##' p1 <- seq(0,1,by=0.1)
##' 
##' res1 <- qexppiecewise(p1, rate = rate1, breaks = 10)
##' res1
##' pexppiecewise(res1, rate = rate1, breaks = 10) - p1
##' 
##' ## EXAMPLE 3
##' rate3 <- c(0.1,0.1,0.1)
##' p3 <- seq(0,1,0.1)
##' 
##' res3 <- qexppiecewise(p = p3, rate = rate3, breaks = c(1,10))
##' res3
##' pexppiecewise(res3, rate = rate3, breaks = c(1,10)) - p3
qexppiecewise <- function(p, rate, breaks){

    tol <- 1e-12    
    n.rate <- length(rate)
    diff.breaks <- diff(c(0,breaks,Inf))

    ## evaluate the cumulative hazard at the end of each interval (except the last interval, Inf)
    A.min <- cumsum(c(0,rate*diff.breaks)[1:n.rate])
    ## find the interval corresponding to the timepoint
    index.A <- cut(p, breaks = c(-tol,pexppiecewise(breaks, rate = rate, breaks = breaks),1+tol))
    indexNum.A <- as.numeric(index.A)
    ## revert the cumulative hazard for each timepoint
    out <- c(0,breaks)[indexNum.A] + (- log(1 - p) - A.min[indexNum.A])/rate[indexNum.A]
    return(out)
}


##----------------------------------------------------------------------
### simTrial.R ends here
