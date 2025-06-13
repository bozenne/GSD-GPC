### summary.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 12 2025 (10:34) 
## Version: 
## Last-Updated: jun 13 2025 (15:15) 
##           By: Brice Ozenne
##     Update #: 22
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary (code)
summary.simTrial <- function(object, interim = NULL){

    admin.censoring <- attr(object, "admin.censoring")

    ## ** apply interim and find states
    data.interim <- subset(object, interim = interim)
    interim.state <- table(data.interim$group, data.interim$state)

    ## ** count final state
    M.event <- cbind(Alive = rowSums(interim.state[,c("alive","alive & tox")]),
                     Dead = rowSums(interim.state[,c("dead","dead & tox")]),
                     Tox = rowSums(interim.state[,c("alive & tox","dead & tox","tox & dropout")]),
                     Dropout = rowSums(interim.state[,c("dropout","tox & dropout")])
                     )

    M.event[] <- paste0(M.event[]," (",formatC(100*M.event[]/matrix(table(object$group), nrow = 2, ncol = NCOL(M.event), byrow = FALSE), digits = 2, format = "f"),"%)")

    ## ** count final state
    M.state <- interim.state
    M.state[] <- paste0(M.state[]," (",formatC(100*M.state[]/matrix(table(object$group), nrow = 2, ncol = NCOL(M.state), byrow = FALSE), digits = 2, format = "f"),"%)")

    ## ** interim
    if(!is.null(interim)){
        n.interim <- c("interim" = NROW(data.interim), "interimCC" = sum(data.interim$pipeline == FALSE), "pipeline" = sum(data.interim$pipeline))
        npc.interim <- paste0(n.interim," (",formatC(100*n.interim/NROW(object), digits = 2, format = "f"),"%)")
    }

    ## ** display
    if(is.null(interim)){        
        cat("  - At final: t=",1,"+",admin.censoring,"=",1 + admin.censoring,"\n",sep="")
        cat("      ",NROW(object)," individuals with opportunity for complete follow-up \n\n",sep="")
    }else{
        cat("  - At interim: t=",interim,"+",admin.censoring,"=",interim + admin.censoring,"\n",sep="")
        cat("      ",npc.interim[1]," individuals: ",npc.interim[2]," with opportunity for complete follow-up and ",npc.interim[3]," pipeline patients\n\n",sep="")
    }

    cat("  - Marginal frequency of the events:\n")
    rownames(M.event) <- paste0("  ", rownames(M.event))
    print(M.event, quote = FALSE)
    
    cat("\n")
    cat("  - Possible trajectories:\n")
    rownames(M.state) <- paste0("  ", rownames(M.state))
    print(M.state, quote = FALSE)

    

    ## ** export
    return(invisible(NULL))
}


##----------------------------------------------------------------------
### summary.R ends here
