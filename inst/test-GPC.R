### test-GPC.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 12 2025 (12:13) 
## Version: 
## Last-Updated: jun 13 2025 (17:21) 
##           By: Brice Ozenne
##     Update #: 12
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * source files
library(BuyseTest)
library(pbapply)
library(parallel)
library(rpact)

if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    path <- "~/Github/GSD-GPC/R/"
}else if(system("whoami",intern=TRUE)=="hpl802"){  
    path <- "./R"
}

sapply(list.files(path,full.names=TRUE), source)

## * planning
n.E <- 200
n.C <- 200

Mplan <- do.call(rbind,pblapply(1:100, function(iSeed){
    iDf <- simTrial(n.E = n.E, n.C = n.C, seed = iSeed, scale.censoring.C = Inf)
    iDf$statusTox.bin <- (iDf$statusTox == 1) 
    iBT <- BuyseTest(group ~ tte(timeSurv, statusSurv,threshold = 0.1) + bin(statusTox.bin, operator = "<0"),
                     data = iDf, trace = FALSE)
    return(cbind(seed = iSeed, endpoint = rownames(confint(iBT)), confint(iBT)))
}))
rownames(Mplan) <- NULL
Infoplan <- 1/tapply(Mplan$se^2, Mplan$endpoint, mean)

designOF <- getDesignGroupSequential(kMax = 2,
                                     alpha = 0.025,
                                     sided = 1,
                                     informationRates = c(0.6,1),
                                     typeOfDesign = "asOF")
designP <- getDesignGroupSequential(kMax = 2,
                                    alpha = 0.025,
                                    sided = 1,
                                    informationRates = c(0.6,1),
                                    typeOfDesign = "asP")

critical.threshold <- rbind(fixed = c(Inf,Inf,qnorm(1-0.05/2)),
                            gsd.none = c(qnorm(1-0.05/2),qnorm(1-0.05/2),qnorm(1-0.05/2)),
                            gsd.of = c(designP$criticalValues[1],designP$criticalValues),
                            gsd.p = c(designOF$criticalValues[1],designOF$criticalValues),
                            gsd.bonferroni = c(qnorm(1-0.05/4),qnorm(1-0.05/4),qnorm(1-0.05/4)))
colnames(critical.threshold) <- c("interim","decision","final")

## * warper
warper <- function(iSeed){ ## iSeed <- 1

    ## ** simulate data
    df.final <- simTrial(n.E = n.E, n.C = n.C, seed = iSeed, scale.censoring.C = Inf)
    df.final$statusTox.bin <- (df.final$statusTox == 1) 
    ## plot(df.final)
    df.interim <- subset(df.final, interim =  0.75)
    df.interim$statusTox.bin <- (df.interim$statusTox == 1) 
    ## summary(df.interim)
    ## plot(df.interim)
    df.decision <- df.final[df.final$id %in% df.interim$id,]

    ## ** run GPC
    eBT.interim <- BuyseTest(group ~ tte(timeSurv, statusSurv,threshold = 0.1) + bin(statusTox.bin, operator = "<0"),
                             data = df.interim, trace = FALSE)
    eBTinfo.interim <- (1/confint(eBT.interim)[2,"se"]^2)/Infoplan[2]
    eBTstat.interim <- confint(eBT.interim)[2,"estimate"]/confint(eBT.interim)[2,"se"]

    eBT.decison <- BuyseTest(group ~ tte(timeSurv, statusSurv,threshold = 0.1) + bin(statusTox.bin, operator = "<0"),
                           data = df.decision, trace = FALSE)
    eBTinfo.decision <- (1/confint(eBT.decision)[2,"se"]^2)/Infoplan[2]
    eBTstat.decision <- confint(eBT.decision)[2,"estimate"]/confint(eBT.decision)[2,"se"]

    eBT.final <- BuyseTest(group ~ tte(timeSurv, statusSurv,threshold = 0.1) + bin(statusTox.bin, operator = "<0"),
                           data = df.final, trace = FALSE)
    eBTinfo.final <- (1/confint(eBT.final)[2,"se"]^2)/Infoplan[2]
    eBTstat.final <- confint(eBT.final)[2,"estimate"]/confint(eBT.final)[2,"se"]

    ## ** retrieve boundaries
    bound <- critical.threshold
    bound["gsd.of",1] <- getDesignGroupSequential(kMax = 2,
                                                 alpha = 0.025,
                                                 sided = 1,
                                                 informationRates = c(eBTinfo.interim,1),
                                                 typeOfDesign = "asOF")$criticalValues
    bound["gsd.p",1] <- getDesignGroupSequential(kMax = 2,
                                                alpha = 0.025,
                                                sided = 1,
                                                informationRates = c(eBTinfo.interim,1),
                                                typeOfDesign = "asP")$criticalValues

    
    ## ** retrieve boundaries
    decision.interim <- eBTstat.interim >= bound[,"interim"]
    decision.final <- ifelse(decision.interim, eBTstat.decision >= bound[,"decision"], NA)
    decision.final <- ifelse(decision.interim, NA, eBTstat.final >= bound[,"final"])

    out <- data.frame(seed = iSeed,
                      method = rownames(bound),
                      n.interim = NROW(df.interim),
                      n.pipeline = sum(df.interim$pipeline),
                      estimate.interim = confint(eBT.interim)[2,"estimate"],
                      se.interim = confint(eBT.interim)[2,"se"],
                      pvalue.interim = confint(eBT.interim)[2,"p.value"],
                      bound.interim = bound[,"interim"],
                      decision.interim = decision.interim,
                      estimate.final = confint(eBT.final)[2,"estimate"],
                      se.final = confint(eBT.final)[2,"se"],
                      pvalue.final = confint(eBT.final)[2,"p.value"],
                      bound.final = bound[,"final"],
                      decision.final = decision.final,
                      decision = ifelse(decision.interim,decision.interim,decision.final))

    return(out)
}

## * simulation study
cl <- makeCluster(50)
parallel::clusterExport(cl, varlist = c("path","n.E","n.C","Infoplan","critical.threshold"))
xxx <- parallel::clusterCall(cl = cl, function(x){
    library(BuyseTest)
    library(rpact)
    sapply(list.files(path,full.names=TRUE), source)
})

ls.sim <- pblapply(1:1000, warper, cl = cl)

library(data.table)
dt.sim <- as.data.table(do.call(rbind,ls.sim))
dt.sim[,.(.N,
          decision.interim = mean(decision.interim),
          decision.final = sum(decision.final, na.rm=TRUE)/.N,
          decision = mean(decision)),
       by = "method"]


##----------------------------------------------------------------------
### test-GPC.R ends here
