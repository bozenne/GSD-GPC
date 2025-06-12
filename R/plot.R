### plot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 12 2025 (10:34) 
## Version: 
## Last-Updated: jun 12 2025 (12:07) 
##           By: Brice Ozenne
##     Update #: 51
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * plot (code)
plot.simTrial <- function(object, interim = NULL, id = NULL,
                          type = "lexis", size = 2, linewidth = 1, alpha = 1, bins = NULL,
                          facet = NULL, facet_type = ggplot2::facet_grid, labeller = "label_value"){

    tol <- 1e-12
    requireNamespace("ggplot2")
    requireNamespace("prodlim")

    ## ** normalize user input
    admin.censoring <- attr(object, "admin.censoring")
    type <- match.arg(type, choices = c("lexis","recruitment","survival","toxicity"))

    object <- subset(object, interim = interim) ## add state variable
    if(!is.null(id)){
        object <- object[object$id %in% id,,drop=FALSE]
    }

    ## ** graphical display
    gg <- ggplot2::ggplot()
    if(type == "recruitment"){
        if(is.null(interim)){
            seqTime <- sort(unique(c(0, object$timeInclusion, object$timeInclusion+object$timeSurv, 1+admin.censoring)))
        }else{
            seqTime <- sort(unique(c(0, object$timeInclusion, setdiff(object$timeInclusion+object$timeSurv, interim), interim-10*tol)))
        }
        table.n <- do.call(rbind,lapply(seqTime, function(iTime){
            iObject <- object[object$timeInclusion < iTime,,drop=FALSE]
            c(n = NROW(iObject), pipeline = sum((iObject$timeInclusion + iObject$timeSurv) > iTime))
        }))
        dfW.n <- data.frame(time = seqTime, n = table.n[,"n"], completed = table.n[,"n"] - table.n[,"pipeline"], pipeline = table.n[,"pipeline"])

        gg <- gg + ggplot2::geom_step(data = dfW.n, ggplot2::aes(x = time, y = n, color = "recruited", group = "recruited"), linewidth = linewidth)
        gg <- gg + ggplot2::geom_step(data = dfW.n, ggplot2::aes(x = time, y = completed, color = "completed", group = "completed"), linewidth = linewidth)
        gg <- gg + ggplot2::geom_step(data = dfW.n, ggplot2::aes(x = time, y = pipeline, color = "pipeline", group = "pipeline"), linewidth = linewidth)
        gg <- gg + ggplot2::labs(x = "Time since start of the trial", y = "Number of participants",
                                 color = "", shape = "")

    }else if(type == "lexis"){
        gg <- gg + ggplot2::geom_segment(data = object,
                                         mapping = ggplot2::aes(x = timeInclusion, y = 0, xend = timeInclusion + timeSurv, yend = timeSurv), linewidth = linewidth)
        gg <- gg + ggplot2::geom_point(data = object[object$statusTox==1,,drop=FALSE],
                                       mapping = ggplot2::aes(x = timeInclusion + timeTox, y = timeTox, shape = "toxicity", color = "toxicity"), size = size)
        gg <- gg + ggplot2::geom_point(data = object[object$statusSurv==1,,drop=FALSE],
                                       mapping = ggplot2::aes(x = timeInclusion + timeSurv, y = timeSurv, shape = "death", color = "death"), size = size)
        gg <- gg + ggplot2::geom_point(data = object[(object$statusSurv==0) & (object$timeSurv + tol < admin.censoring),,drop=FALSE],
                                       mapping = ggplot2::aes(x = timeInclusion + timeSurv, y = timeSurv, shape = "drop-out", color = "drop-out"), size = size)
        gg <- gg + ggplot2::geom_point(data = object[(object$statusSurv==0) & (object$timeSurv + tol >= admin.censoring),,drop=FALSE],
                                       mapping = ggplot2::aes(x = timeInclusion + timeSurv, y = timeSurv, shape = "alive", color = "alive"), size = size)
        gg <- gg + ggplot2::labs(x = "Time since start of the trial", y = "Time since inclusion of the participant",
                                 color = "Type of event", shape = "Type of event")
        gg <- gg + ggplot2::scale_colour_manual(breaks = c("alive","toxicity","death","drop-out"), values = c("green","orange","red","gray"))
        gg <- gg + ggplot2::scale_shape_manual(breaks = c("alive","toxicity","death","drop-out"), values = c(19,18,15,8))
        gg <- gg + ggplot2::coord_cartesian(xlim = c(0,1+admin.censoring), ylim = c(0,admin.censoring))
        if(!is.null(facet)){
            gg <- gg + do.call(facet_type, args = list(facet, labeller = labeller))
        }        
    }else if(type == "toxicity"){
        if("eta.timeTox" %in% names(object) == FALSE){
            object$time <- pmin(object$timeSurv, object$timeTox)
            object$event <- c("alive" = 0,"alive & tox" = 1, "dead" = 2, "dead & tox" = 1, "tox & dropout" = 1, "dropout" = 0)[as.character(object$state)]
            e.prodlim <- prodlim::prodlim(prodlim::Hist(time,event) ~ group, data = object)
            gg <- prodlim::ggprodlim(e.prodlim, cause = 1)
        }else{
            gg <- gg + ggplot2::geom_histogram(data = object, mapping = ggplot2::aes(x = eta.timeTox, group = group, fill = group),
                                               alpha = alpha, bins = bins, position = "dodge")
            gg <- gg + ggplot2::labs(x = "Time to toxicity (uncensored)")
        }
    }else if(type == "survival"){
        if("eta.timeSurv" %in% names(object) == FALSE){
            e.prodlim <- prodlim::prodlim(prodlim::Hist(timeSurv,statusSurv) ~ group, data = object)
            gg <- prodlim::ggprodlim(e.prodlim)
        }else{
            gg <- gg + ggplot2::geom_histogram(data = object, mapping = ggplot2::aes(x = eta.timeSurv, group = group, fill = group),
                                               alpha = alpha, bins = bins, position = "dodge")
            gg <- gg + ggplot2::labs(x = "Time to death (uncensored)")
        }
    }
    
    ## ** export
    return(gg)
    
}


##----------------------------------------------------------------------
### plot.R ends here
