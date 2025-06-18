### test-GPC.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 12 2025 (12:13) 
## Version: 
## Last-Updated: jun 18 2025 (17:21) 
##           By: Brice Ozenne
##     Update #: 94
##----------------------------------------------------------------------
## 
### Commentary: 
## cd /projects/biostat01/people/hpl802/GSD-GPC/R
## source("../inst/test-GPC.R")
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

options(width = 150)


## * source files
library(BuyseTest)
library(survival)
library(pbapply)
library(parallel)
library(rpact)
library(data.table)

if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    path <- "~/Github/GSD-GPC/R/"
}else if(system("whoami",intern=TRUE)=="hpl802"){  
    path <- "."
}

sapply(list.files(path,full.names=TRUE), source)

## * planning
n.E <- 200
n.C <- 200

## * warper
warper <- function(iSeed, maxInfo = c(GPC = 300, HR = 63), typeOfDesign = c("fixed","none","bonferroni","asOF","asP","delayed"),
                   formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox), rho  = 0, ...){ ## iSeed <- 84

    

    ## ** simulate data
    df.final <- simTrial(n.E = n.E, n.C = n.C, seed = iSeed, scale.censoring.C = Inf, 
                         rho.C = rho, rho.E = rho, ...)

    ## ** run model
    e.GPC <- BuyseTest(formula, data = df.final, scoring.rule = "Gehan", trace = FALSE)
    e.Cox <- coxph(Surv(timeSurv, statusSurv) ~ group, data = df.final)

    ## ** extact information and test statistic
    if(is.null(maxInfo)){
        names(maxInfo) <- c("GPC","HR")
    }
    eGSD.GPC <- cbind(estimand = "NTB", runTrial(e.GPC, data = df.final, interim = 0.5, typeOfDesign = typeOfDesign, maxInfo = maxInfo["GPC"]))
    eGSD.Cox <- cbind(estimand = "HR", runTrial(e.Cox, data = df.final, interim = 0.5, typeOfDesign = typeOfDesign, maxInfo = maxInfo["HR"]))

    ## ** export
    out <- cbind(seed = iSeed, rho = rho, rbind(eGSD.GPC, eGSD.Cox))
    return(out)
}

warper(1)

## * simulation study
cl <- makeCluster(100)
parallel::clusterExport(cl, varlist = c("path","n.E","n.C","warper"))
tempo <- parallel::clusterCall(cl = cl, function(x){
    library(BuyseTest)
    library(rpact)
    library(survival)
    sapply(list.files(path,full.names=TRUE), source)
})

## ** planned information
ls.planInfo <- pblapply(1:1000, function(iSeed){
    warper(iSeed, typeOfDesign = "fixed")
}, cl = cl)
dtL.planInfo <- as.data.table(do.call(rbind,ls.planInfo))
dtS.planInfo <- dtL.planInfo[stage == "final", .(.N, info = mean(1/se^2), info.min = min(1/se^2), info.max = max(1/se^2)), by = "estimand"]
dtS.planInfo
##    estimand     N      info  info.min  info.max
##      <char> <int>     <num>     <num>     <num>
## 1:      NTB  1000 300.38569 295.64606 309.47063
## 2:       HR  1000  62.97779  53.68651  70.16954

## ** type 1 error 
ls.simH0 <- pblapply(1:20000, function(iSeed){
    rbind(warper(iSeed, rho = -0.5, maxInfo = c(GPC = 300, HR = 63)),
          warper(iSeed, rho = -0.25, maxInfo = c(GPC = 300, HR = 63)),
          warper(iSeed, rho = 0, maxInfo = c(GPC = 300, HR = 63)),
          warper(iSeed, rho = 0.25, maxInfo = c(GPC = 300, HR = 63)),
          warper(iSeed, rho = 0.5, maxInfo = c(GPC = 300, HR = 63)))
}, cl = cl)
dt.simH0 <- as.data.table(do.call(rbind,ls.simH0))



## *** information
dtL.simH0.info <- dt.simH0[design=="fixed",.(.N, info = mean(1/se^2), info.min = min(1/se^2), info.max = max(1/se^2), info.pc = mean(info.pc)), by = c("estimand","stage","rho")]
dtL.simH0.info[, info2 := paste0(round(info,2)," (",round(100*info.pc,2),"%) [",round(info.min,2),";",round(info.max,2),"]")]
dtW.simH0.info <- reshape(dtL.simH0.info[,.(N, estimand,stage,rho,info2)], direction = "wide", timevar = "stage", idvar = c("estimand","rho"), v.names = "info2")
dtW.simH0.info[order(dtW.simH0.info$estimand,),]
dtW.simH0.info <- reshape(dtL.simH0.info[,.(N, estimand,stage,rho,info2)], direction = "wide", timevar = "stage", idvar = c("estimand","rho"), v.names = "info2")
dtW.simH0.info[order(dtW.simH0.info$estimand,),]
##        N estimand   rho                    info2.interim                  info2.decision                      info2.final
##     <int>   <char> <num>                           <char>                          <char>                           <char>
##  1: 20000       HR -0.50      14.29 (22.69%) [8.67;20.67]    47.26 (75.02%) [38.96;55.15]    63.08 (100.12%) [53.87;71.67]
##  2: 20000       HR -0.25      14.29 (22.69%) [8.91;20.91]    47.27 (75.03%) [39.47;56.34]     63.08 (100.13%) [53.24;72.2]
##  3: 20000       HR  0.00       14.3 (22.69%) [8.64;20.68]    47.27 (75.03%) [38.97;55.95]    63.08 (100.13%) [52.99;73.22]
##  4: 20000       HR  0.25      14.31 (22.71%) [8.42;20.95]    47.27 (75.04%) [38.66;55.73]    63.08 (100.13%) [53.99;72.73]
##  5: 20000       HR  0.50      14.31 (22.71%) [8.68;20.49]    47.26 (75.02%) [38.98;55.16]    63.07 (100.11%) [52.72;71.46]
##  6: 20000      NTB -0.50 413.24 (137.75%) [349.87;512.18] 225.05 (75.02%) [215.88;243.36]  300.07 (100.02%) [292.34;317.9]
##  7: 20000      NTB -0.25 401.53 (133.84%) [332.84;495.87] 225.12 (75.04%) [216.46;244.36] 300.15 (100.05%) [292.53;319.88]
##  8: 20000      NTB  0.00  391.06 (130.35%) [323.25;479.9] 225.27 (75.09%) [218.23;242.71] 300.34 (100.11%) [292.53;317.33]
##  9: 20000      NTB  0.25  381.31 (127.1%) [324.88;474.98] 225.56 (75.19%) [218.43;243.46] 300.72 (100.24%) [292.97;319.72]
## 10: 20000      NTB  0.50  371.39 (123.8%) [318.04;465.51]  226.1 (75.37%) [218.72;249.76] 301.42 (100.47%) [293.84;326.52]

## *** bound
dtL.simH0.bound <- dt.simH0[rho==0,.(rho = rho[1], .N, bound = mean(bound,na.rm=TRUE), bound.min = min(bound,na.rm=TRUE), bound.max = max(bound,na.rm=TRUE)), by = c("estimand","stage","design")]
dtL.simH0.bound[, bound2 := paste0(round(bound,2)," [",round(bound.min,2),";",round(bound.max,2),"]")]
reshape(dtL.simH0.bound[,.(rho,.N,estimand,stage,design,bound2)], direction = "wide", timevar = "stage", idvar = c("estimand","design"), v.names = "bound2")
##       rho     N estimand     design   bound2.interim  bound2.decision     bound2.final
##     <num> <int>   <char>     <char>           <char>           <char>           <char>
##  1:     0    36      NTB      fixed    Inf [Inf;Inf]      NaN [NA;NA] 1.96 [1.96;1.96]
##  2:     0    36      NTB       none 1.96 [1.96;1.96] 1.96 [1.96;1.96] 1.96 [1.96;1.96]
##  3:     0    36      NTB bonferroni 2.24 [2.24;2.24] 2.24 [2.24;2.24] 2.24 [2.24;2.24]
##  4:     0    36      NTB       asOF -Inf [-Inf;-Inf] 1.96 [1.96;1.96]      NaN [NA;NA]
##  5:     0    36      NTB        asP -Inf [-Inf;-Inf] 1.96 [1.96;1.96]      NaN [NA;NA]
##  6:     0    36      NTB    delayed -Inf [-Inf;-Inf] 1.96 [1.96;1.96]      NaN [NA;NA]
##  7:     0    36       HR      fixed    Inf [Inf;Inf]      NaN [NA;NA] 1.96 [1.96;1.96]
##  8:     0    36       HR       none 1.96 [1.96;1.96] 1.96 [1.96;1.96] 1.96 [1.96;1.96]
##  9:     0    36       HR bonferroni 2.24 [2.24;2.24] 2.24 [2.24;2.24] 2.24 [2.24;2.24]
## 10:     0    36       HR       asOF 4.58 [3.74;5.94] 2.01 [1.98;2.05] 1.98 [1.96;2.18]
## 11:     0    36       HR        asP  2.4 [2.28;2.56] 2.26 [2.24;2.26] 2.09 [2.05;2.15]
## 12:     0    36       HR    delayed 3.02 [2.78;3.31] 1.55 [1.45;1.64] 1.98 [1.97;1.99]

## *** Decision
dt.simH0.decision <- dt.simH0[rho==0,.(nSim = (.N/3),
                                       earlyStop = sum(decision %in% c("stop","infoMax"),na.rm=TRUE)/(.N/3),
                                       earlyH0 = sum(!is.na(decision) & decision == "futility" & stage == "decision", na.rm=TRUE)/(.N/3),
                                       earlyH1 = sum(!is.na(decision) & decision == "efficacy" & stage == "decision", na.rm=TRUE)/(.N/3),
                                       H0 = sum(decision == "futility",na.rm=TRUE)/(.N/3),
                                       H1 = sum(decision == "efficacy",na.rm=TRUE)/(.N/3),
                                       maxInfo = sum(decision == "infoMax",na.rm=TRUE)/(.N/3),
                                       reversal = sum(reversal, na.rm=TRUE)/(.N/3)),
                              by = c("estimand","design")]
dt.simH0.decision
##     estimand     design  nSim earlyStop earlyH0 earlyH1      H0      H1 maxInfo reversal
##       <char>     <char> <num>     <num>   <num>   <num>   <num>   <num>   <num>    <num>
##  1:      NTB      fixed 20000   1.00000 0.00000 0.00000 0.00000 0.00000       1  0.00000
##  2:      NTB       none 20000   1.00000 0.97255 0.02745 0.97255 0.02745       1  0.00000
##  3:      NTB bonferroni 20000   1.00000 0.98555 0.01445 0.98555 0.01445       1  0.00000
##  4:      NTB       asOF 20000   1.00000 0.97255 0.02745 0.97255 0.02745       1  0.00000
##  5:      NTB        asP 20000   1.00000 0.97255 0.02745 0.97255 0.02745       1  0.00000
##  6:      NTB    delayed 20000   1.00000 0.97255 0.02745 0.97255 0.02745       1  0.00000
##  7:       HR      fixed 20000   0.00000 0.00000 0.00000 0.97475 0.02525       0  0.00000
##  8:       HR       none 20000   0.02320 0.01870 0.00450 0.97390 0.02610       0  0.01870
##  9:       HR bonferroni 20000   0.01085 0.00925 0.00160 0.98675 0.01325       0  0.00925
## 10:       HR       asOF 20000   0.00000 0.00000 0.00000 0.97605 0.02395       0  0.00000
## 11:       HR        asP 20000   0.00640 0.00510 0.00130 0.98185 0.01815       0  0.00510
## 12:       HR    delayed 20000   0.00045 0.00020 0.00025 0.97550 0.02450       0  0.00020
quantile(sapply(1:1e3, function(iSeed){set.seed(iSeed);mean(rbinom(1e4, size = 1, prob = 0.025))}))
##     0%    25%    50%    75%   100% 
## 0.0199 0.0239 0.0250 0.0262 0.0297 

warper(1, rho=0)


## ** power [TO BE DONE]
## ls.simH1 <- pblapply(1:1000, function(iSeed){
##     rbind(warper(iSeed, rho = -0.5, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox),
##                  scale.C = c(1,0.8), scale.E = c(0.9,5), shape.C = c(1,1), shape.E = c(1,1)),
##           warper(iSeed, rho = -0.25, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox),
##                  scale.C = c(1,0.8), scale.E = c(0.9,5), shape.C = c(1,1), shape.E = c(1,1)),
##           warper(iSeed, rho = 0, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox),
##                  scale.C = c(1,0.8), scale.E = c(0.9,5), shape.C = c(1,1), shape.E = c(1,1)),
##           warper(iSeed, rho = 0.25, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox),
##                  scale.C = c(1,0.8), scale.E = c(0.9,5), shape.C = c(1,1), shape.E = c(1,1)),
##           warper(iSeed, rho = 0.5, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox),
##                  scale.C = c(1,0.8), scale.E = c(0.9,5), shape.C = c(1,1), shape.E = c(1,1)))
## }, cl = cl)

## dt.simH1 <- as.data.table(do.call(rbind,ls.simH1))
## dt.simH1[,.(.N,
##             info = paste(round(100*mean(info.interim),2),round(100*mean(info.decision),2),round(100*mean(info.final),2),sep=";"),
##             stop.interim = 100*mean(!is.na(decision.interim)),
##             reversal.interim = 100*sum(reversal.interim, na.rm=TRUE)/.N,
##             decision.interim = 100*sum(decision.interim, na.rm=TRUE)/.N,
##             decision.final = 100*sum(decision.final, na.rm=TRUE)/.N,
##             decision = mean(decision)),
##          by = c("rho","method")]


 

## ## example of reversal
## df.final <- simTrial(n.E = 10*n.E, n.C = 10*n.C, seed = 1, scale.censoring.C = Inf, scale.E = c(0.6,1))
## df.final$statusTox.bin <- (df.final$statusTox == 1) 

## df.interim <- subset(df.final, interim = 0.75)

## plot(df.final, facet = ~group)
## plot(df.interim, facet = ~group)

## eBT.interim <- BuyseTest(group ~ tte(timeSurv, statusSurv, threshold = 0.1) + bin(statusTox.bin, operator = "<0"),
##                          data = df.interim, scoring.rule = "Gehan", trace = FALSE)
## confint(eBT.interim)

## eBT.final <- BuyseTest(group ~ tte(timeSurv, statusSurv, threshold = 0.1) + bin(statusTox.bin, operator = "<0"),
##                        data = df.final, scoring.rule = "Gehan", trace = FALSE)
## confint(eBT.final)

## eBT.final <- BuyseTest(group ~ cont(timeSurv) + cont(timeTox, operator = "<0"),
##                        data = df.final, trace = FALSE)
## confint(eBT.final)

## confint(BuyseTest(group ~ cont(timeSurv, threshold = 0.1) + bin(statusTox.bin, operator = "<0"), data = df.final, trace = FALSE))

## ** export
if(FALSE){
    saveRDS(dtL.simH0.info, file = file.path("../inst/results","H0-info.rds"))
    saveRDS(dtL.simH0.bound, file = file.path("../inst/results","H0-bound.rds"))
    saveRDS(dt.simH0.decision, file = file.path("../inst/results","H0-decision.rds"))
}

##----------------------------------------------------------------------
### test-GPC.R ends here



