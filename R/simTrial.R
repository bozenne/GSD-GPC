### simTrial.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun  3 2025 (09:38) 
## Version: 
## Last-Updated: jun 12 2025 (12:11) 
##           By: Brice Ozenne
##     Update #: 171
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * simTrial (documentation)
##' @param n.E [integer] sample size in the experimental group
##' @param n.C [integer] sample size in the control group
##' @param dist.recruitment [character] distribution used to model the inclusion of the patients.
##' Can be \code{"uniform"} for deterministic equally spread inclusion between 0 and 1
##' or \code{"runif"} for stochastic inclusion based on a uniform distribution.
##' 
##' @examples
##'
##' ## simulate data
##' df.sim <- simTrial(n.E = 50, n.C = 100, seed = 1, latent = TRUE)
##'
##' ## display data
##' df.sim
##' summary(df.sim)
##' plot(df.sim, type = "lexis")
##' plot(df.sim, facet = ~group, labeller = label_both)
##' plot(df.sim, facet = ~state, labeller = label_both)
##' plot(df.sim, facet = group~state, labeller = label_both)
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

## * simTrial (code)
simTrial <- function(n.E, n.C = NULL, dist.recruitment = "uniform",
                     scale.E = c(0.5,0.5), shape.E = c(0.5,2), dist.E = c("weibull", "weibull"),
                     scale.C = NULL, shape.C = NULL, dist.C = NULL,
                     scale.censoring.E = 1, shape.censoring.E = 1, dist.censoring.E = "weibull",
                     scale.censoring.C = NULL, shape.censoring.C = NULL, dist.censoring.C = NULL,
                     rho = rep(0,2), admin.censoring = 0.5, seed = NULL, latent = FALSE
                     ){
    requireNamespace("lava")

    ## ** normalize user input
    name.time <- c("timeSurv","timeTox")
    name.etatime <- paste0("eta.",name.time)
    name.censtime <- paste0("cens.",name.time)
    name.status <- c("statusSurv","statusTox")
    
    ## sample size per arm
    if(is.null(n.C)){
        n.C <- n.E
    }

    ## distribution of the recruitment times
    dist.recruitment <- match.arg(dist.recruitment, choices = c("uniform","runif"))

    ## outcome distribution
    if(length(scale.E) != 2){
        stop("Argument \'scale.E\' should have length 2. \n")
    }
    if(length(shape.E) != 2){
        stop("Argument \'scale.E\' should have length 2. \n")
    }
    if(length(dist.E) != 2){
        stop("Argument \'scale.E\' should have length 2. \n")
    }else{
        dist.E[1] <- match.arg(dist.E[1], choices = c("weibull", "uniform","piecewiseExp"))
        dist.E[2] <- match.arg(dist.E[2], choices = c("weibull", "uniform","piecewiseExp"))
    }
    if(is.null(scale.C)){
        scale.C <- scale.E
    }else if(length(scale.C) != 2){
        stop("Argument \'scale.C\' should have length 2. \n")
    }
    if(is.null(shape.C)){
        shape.C <- shape.E
    }else if(length(shape.C) != 2){
        stop("Argument \'scale.C\' should have length 2. \n")
    }
    
    if(is.null(dist.C)){
        dist.C <- dist.E
    }else if(length(dist.C) != 2){
        stop("Argument \'dist.C\' should have length 2. \n")
    }else{
        dist.C[1] <- match.arg(dist.C[1], choices = c("weibull", "uniform","piecewiseExp"))
        dist.C[2] <- match.arg(dist.C[2], choices = c("weibull", "uniform","piecewiseExp"))
    }

    if(length(scale.censoring.E) != 1){
        stop("Argument \'scale.censoring.E\' should have length 1. \n")
    }
    if(length(shape.censoring.E) != 1){
        stop("Argument \'scale.censoring.E\' should have length 1. \n")
    }
    if(length(dist.censoring.E) != 1){
        stop("Argument \'scale.censoring.E\' should have length 1. \n")
    }else{
        dist.censoring.E <- match.arg(dist.censoring.E, choices = c("weibull", "uniform","piecewiseExp"))
    }
    if(is.null(scale.censoring.C)){
        scale.censoring.C <- scale.censoring.E
    }else if(length(scale.censoring.C) != 1){
        stop("Argument \'scale.censoring.C\' should have length 1. \n")
    }
    if(is.null(shape.censoring.C)){
        shape.censoring.C <- shape.censoring.E
    }else if(length(shape.censoring.C) != 1){
        stop("Argument \'scale.censoring.C\' should have length 1. \n")
    }
    
    if(is.null(dist.censoring.C)){
        dist.censoring.C <- dist.censoring.E
    }else if(length(dist.censoring.C) != 1){
        stop("Argument \'dist.censoring.C\' should have length 1. \n")
    }else{
        dist.censoring.C <- match.arg(dist.censoring.C, choices = c("weibull", "uniform","piecewiseExp"))
    }

    ## ** data generating mechanism
    model.C <- lava::lvm()
    ## lava::distribution(model.C,~id) <- lava::Sequence.lvm(1,n.C)
    lava::categorical(model.C,labels="control") <- "group"
    model.E <- lava::lvm()
    ## lava::distribution(model.E,~id) <- lava::Sequence.lvm(n.C+1,n.C+n.E)
    lava::categorical(model.E,labels="experimental") <- "group"
    
    ## inclusion time
    if(dist.recruitment == "uniform"){
        lava::distribution(model.C, "timeInclusion") <- lava::Sequence.lvm(0,1)
        lava::distribution(model.E, "timeInclusion") <- lava::Sequence.lvm(0,1)
    }else if(dist.recruitment == "runif"){
        lava::distribution(model.C, "timeInclusion") <- lava::uniform.lvm(0,1)
        lava::distribution(model.E, "timeInclusion") <- lava::uniform.lvm(0,1)
    }
    
    ## censoring distribution
    lava::distribution(model.C, "timeCens") <- switch(dist.censoring.C,
                                                      "uniform" = lava::uniform.lvm(a = scale.censoring.C, b = shape.censoring.C),
                                                      "weibull" = lava::weibull.lvm(scale = scale.censoring.C, shape = 1/shape.censoring.C),
                                                      "piecewiseExp" = lava::coxExponential.lvm(scale = scale.censoring.C, timecut = shape.censoring.C))
    lava::distribution(model.E, "timeCens") <- switch(dist.censoring.E,
                                                      "uniform" = lava::uniform.lvm(a = scale.censoring.E, b = shape.censoring.E),
                                                      "weibull" = lava::weibull.lvm(scale = scale.censoring.E, shape = 1/shape.censoring.E),
                                                      "piecewiseExp" = lava::coxExponential.lvm(scale = scale.censoring.E, timecut = shape.censoring.E))

    ## outcome distribution
    for(iterO in 1:2){
        lava::distribution(model.C, name.etatime[iterO]) <- switch(dist.C[iterO],
                                                                   "uniform" = lava::uniform.lvm(a = scale.C[iterO], b = shape.C[iterO]),
                                                                   "weibull" = lava::weibull.lvm(scale = scale.C[iterO], shape = 1/shape.C[iterO]),
                                                                   "piecewiseExp" = lava::coxExponential.lvm(scale = scale.C[iterO], timecut = shape.C[iterO]))
        lava::distribution(model.E, name.etatime[iterO]) <- switch(dist.E[iterO],
                                                                   "uniform" = lava::uniform.lvm(a = scale.E[iterO], b = shape.E[iterO]),
                                                                   "weibull" = lava::weibull.lvm(scale = scale.E[iterO], shape = 1/shape.E[iterO]),
                                                                   "piecewiseExp" = lava::coxExponential.lvm(scale = scale.E[iterO], timecut = shape.E[iterO]))
    }
    model.C <- lava::latent(model.C, reformulate(termlabels = c(name.etatime,"timeCens",name.censtime)))
    model.E <- lava::latent(model.E, reformulate(termlabels = c(name.etatime,"timeCens",name.censtime)))

    ## observed (right-censored) outcome distribution
    for(iterO in 1:2){
        if(iterO==1){
            iFF <- reformulate("timeCens", name.censtime[iterO])
        }else{
            iFF <- reformulate(c("timeCens",name.etatime[1:(iterO-1)]), name.censtime[iterO])
        }
        lava::transform(model.C, iFF) <- function(x){apply(cbind(admin.censoring,x),1,min)}
        lava::transform(model.E, iFF) <- function(x){apply(cbind(admin.censoring,x),1,min)}

        txtSurv <- paste0(name.time[iterO], "~min(",name.etatime[iterO],"=1,",name.censtime[iterO],"=0)")
        model.C <- lava::eventTime(model.C, stats::as.formula(txtSurv), name.status[iterO])
        model.E <- lava::eventTime(model.E, stats::as.formula(txtSurv), name.status[iterO])
    }

    ## ** simulate data
    if(!is.null(seed)){
        set.seed(seed)
    }
    out <- rbind(lava::sim(model.E, n = n.E, latent = latent),
                 lava::sim(model.C, n = n.C, latent = latent))
    out.order <- cbind(id = 1:NROW(out),out[order(out$timeInclusion),])
    
    ## ** export
    class(out.order) <- append("simTrial",class(out.order))
    attr(out.order,"seed") <- seed
    attr(out.order,"admin.censoring") <- admin.censoring
    return(out.order)
}




##----------------------------------------------------------------------
### simTrial.R ends here
