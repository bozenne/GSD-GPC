### test-GPC.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 12 2025 (12:13) 
## Version: 
## Last-Updated: jun 19 2025 (16:46) 
##           By: Brice Ozenne
##     Update #: 190
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
library(ggplot2)
library(data.table)

if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    path <- "~/Github/GSD-GPC/R/"
}else if(system("whoami",intern=TRUE)=="hpl802"){  
    path <- "."
}

## file.remove("Rplots.pdf")
sapply(list.files(path,full.names=TRUE,pattern="*.R"), source)

## * planning
n.E <- 200
n.C <- 200
seqRho <- c(-0.5,0.25,0,0.25,0.5)
n.sim <- 2e4
n.cpus <- 100

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
    eGSD.GPC <- cbind(estimand = "NTB", runTrial(e.GPC, data = df.final, interim = 0.75, typeOfDesign = typeOfDesign, maxInfo = maxInfo["GPC"]))
    eGSD.Cox <- cbind(estimand = "HR", runTrial(e.Cox, data = df.final, interim = 0.75, typeOfDesign = typeOfDesign, maxInfo = maxInfo["HR"]))

    ## ** export
    out <- cbind(seed = iSeed, rho = rho, rbind(eGSD.GPC, eGSD.Cox))
    return(out)
}

warper(1)

## * initialize cluster 
cl <- makeCluster(n.cpus)
parallel::clusterExport(cl, varlist = c("path","n.E","n.C","warper","seqRho","n.sim"))
tempo <- parallel::clusterCall(cl = cl, function(x){
    library(BuyseTest)
    library(rpact)
    library(survival)
    sapply(list.files(path,full.names=TRUE), source)
})

## * simulation study (no delay) 
cat("No Delay \n")

## ** under the null hypothesis
## *** example 
dfH0.nodelay <- simTrial(n.E = n.E, n.C = n.C, seed = 1, scale.censoring.C = Inf, scale.C = c(0.01,0.01), admin.censoring = 0.01)
## plot(dfH0.nodelay, facet =~group)
## plot(dfH0.nodelay, type = "recruitment")

coxH0.nodelay <- coxph(Surv(timeSurv, statusSurv) ~ group, data = dfH0.nodelay)
runTrial(coxH0.nodelay, data = dfH0.nodelay, interim = 0.75, typeOfDesign = "delayed", maxInfo = 63)
GPCH0.nodelay <- BuyseTest(group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox), data = dfH0.nodelay, scoring.rule = "Gehan", trace = FALSE)
runTrial(GPCH0.nodelay, data = dfH0.nodelay, interim = 0.75, typeOfDesign = "delayed", maxInfo = 301)

## warper(iSeed, rho = 0, scale.C = c(0.01,0.01), admin.censoring = 0.01, maxInfo = c(GPC = 301, HR = 63))

runTrial(coxH0.nodelay, data = dfH0.nodelay, interim = 0.75, typeOfDesign = "asP", maxInfo = 63)

## *** planned information
ls.planInfoH0.nodelay <- pblapply(1:min(n.sim,1000), function(iSeed){
    warper(iSeed, typeOfDesign = "fixed", scale.C = c(0.01,0.01), admin.censoring = 0.01)
}, cl = cl)
dtL.planInfoH0.nodelay <- as.data.table(do.call(rbind,ls.planInfoH0.nodelay))
planInfoH0.nodelay <- dtL.planInfoH0.nodelay[stage == "final", .(.N, info = mean(1/se^2), info.min = min(1/se^2), info.max = max(1/se^2)), by = "estimand"]
planInfoH0.nodelay
##    estimand     N      info  info.min  info.max
##      <char> <int>     <num>     <num>     <num>
## 1:      NTB  1000 300.92834 295.91965 311.18833
## 2:       HR  1000  62.97779  53.68651  70.16954

## *** simulation
cat(" - type 1 error \n")
ls.nodelayH0 <- pblapply(1:n.sim, function(iSeed){ ## iSeed <- 1
    do.call(rbind,lapply(seqRho, function(iRho){
        warper(iSeed, rho = iRho, scale.C = c(0.01,0.01), admin.censoring = 0.01, maxInfo = c(GPC = 301, HR = 63))
    }))
}, cl = cl)
dt.nodelayH0 <- as.data.table(do.call(rbind,ls.nodelayH0))

## *** information
dtL.nodelayH0.info <- dt.nodelayH0[design=="fixed",.(.N, info = mean(1/se^2), info.min = min(1/se^2), info.max = max(1/se^2), info.pc = mean(info.pc)), by = c("estimand","stage","rho")]
dtL.nodelayH0.info[, info2 := paste0(round(info,2)," (",round(100*info.pc,2),"%) [",round(info.min,2),";",round(info.max,2),"]")]
dtW.nodelayH0.info <- reshape(dtL.nodelayH0.info[,.(N, estimand,stage,rho,info2)], direction = "wide", timevar = "stage", idvar = c("estimand","rho"), v.names = "info2")
dtW.nodelayH0.info[order(dtW.nodelayH0.info$estimand,),]
dtW.nodelayH0.info <- reshape(dtL.nodelayH0.info[,.(N, estimand,stage,rho,info2)], direction = "wide", timevar = "stage", idvar = c("estimand","rho"), v.names = "info2")
dtW.nodelayH0.info[order(dtW.nodelayH0.info$estimand,),]
##         N estimand   rho                   info2.interim                  info2.decision                      info2.final
##     <int>   <char> <num>                          <char>                          <char>                           <char>
## 1: 20000       HR -0.50    46.74 (74.19%) [38.72;54.37]    47.26 (75.02%) [38.96;55.15]   63.08 (100.12%) [53.87;71.67]
## 2: 40000       HR  0.25    46.75 (74.21%) [38.16;54.97]    47.27 (75.04%) [38.66;55.73]   63.08 (100.13%) [53.99;72.73]
## 3: 20000       HR  0.00      46.74 (74.2%) [38.2;54.96]    47.27 (75.03%) [38.97;55.95]   63.08 (100.13%) [52.99;73.22]
## 4: 20000       HR  0.50    46.74 (74.19%) [38.73;54.39]    47.26 (75.02%) [38.98;55.16]   63.07 (100.11%) [52.72;71.46]
## 5: 20000      NTB -0.50  228.81 (76.02%) [214.4;233.45]  223.09 (74.11%) [210.5;225.06]  298.1 (99.04%) [287.23;300.02]
## 6: 40000      NTB  0.25 228.69 (75.98%) [216.46;234.21]  224.26 (74.51%) [213.19;228.6]  299.64 (99.55%) [286.55;303.6]
## 7: 20000      NTB  0.00   228.47 (75.9%) [215.3;233.19]   223.64 (74.3%) [210.96;227.2] 298.82 (99.27%) [285.83;302.21]
## 8: 20000      NTB  0.50 229.35 (76.19%) [216.91;235.49] 225.32 (74.86%) [213.13;230.81] 301.03 (100.01%) [287.58;306.9]

## *** bound
dtL.nodelayH0.bound <- dt.nodelayH0[rho==0,.(rho = rho[1], .N, bound = mean(bound,na.rm=TRUE), bound.min = min(bound,na.rm=TRUE), bound.max = max(bound,na.rm=TRUE)), by = c("estimand","stage","design")]
dtL.nodelayH0.bound[, bound2 := paste0(round(bound,2)," [",round(bound.min,2),";",round(bound.max,2),"]")]
reshape(dtL.nodelayH0.bound[,.(rho,N,estimand,stage,design,bound2)], direction = "wide", timevar = "stage", idvar = c("estimand","design"), v.names = "bound2")
##       rho     N estimand     design   bound2.interim  bound2.decision     bound2.final
##     <num> <int>   <char>     <char>           <char>           <char>           <char>
##  1:     0  1000      NTB      fixed    Inf [Inf;Inf]      NaN [NA;NA] 1.96 [1.96;1.96]
##  2:     0  1000      NTB       none 1.96 [1.96;1.96] 1.96 [1.96;1.96] 1.96 [1.96;1.96]
##  3:     0  1000      NTB bonferroni 2.24 [2.24;2.24] 2.24 [2.24;2.24] 2.24 [2.24;2.24]
##  4:     0  1000      NTB       asOF 2.94 [2.91;3.04] 2.99 [2.96;3.08] 1.97 [1.97;1.97]
##  5:     0  1000      NTB        asP 2.15 [2.15;2.17] 2.16 [2.16;2.18]  2.2 [2.19;2.21]
##  6:     0  1000      NTB    delayed 2.49 [2.47;2.53] 1.77 [1.76;1.78] 2.02 [2.01;2.02]
##  7:     0  1000       HR      fixed    Inf [Inf;Inf]      NaN [NA;NA] 1.96 [1.96;1.96]
##  8:     0  1000       HR       none 1.96 [1.96;1.96] 1.96 [1.96;1.96] 1.96 [1.96;1.96]
##  9:     0  1000       HR bonferroni 2.24 [2.24;2.24] 2.24 [2.24;2.24] 2.24 [2.24;2.24]
## 10:     0  1000       HR       asOF    3 [2.74;3.35] 2.97 [2.71;3.32] 1.97 [1.96;1.98]
## 11:     0  1000       HR        asP 2.16 [2.12;2.22] 2.16 [2.11;2.22]  2.2 [2.17;2.22]
## 12:     0  1000       HR    delayed  2.51 [2.4;2.65] 1.43 [0.62;1.79]    2.02 [2;2.03]

## *** Decision
dt.nodelayH0.decision <- dt.nodelayH0[,.(nSim = (.N/3),
                                       earlyStop = sum(decision %in% c("stop","infoMax"),na.rm=TRUE)/(.N/3),
                                       earlyH0 = sum(!is.na(decision) & decision == "futility" & stage == "decision", na.rm=TRUE)/(.N/3),
                                       earlyH1 = sum(!is.na(decision) & decision == "efficacy" & stage == "decision", na.rm=TRUE)/(.N/3),
                                       H0 = sum(decision == "futility",na.rm=TRUE)/(.N/3),
                                       H1 = sum(decision == "efficacy",na.rm=TRUE)/(.N/3),
                                       maxInfo = sum(decision == "infoMax",na.rm=TRUE)/(.N/3),
                                       reversal = sum(reversal, na.rm=TRUE)/(.N/3)),
                              by = c("estimand","rho","design")]
dt.nodelayH0.decision[rho==0]
##     estimand   rho     design  nSim earlyStop earlyH0 earlyH1      H0      H1 maxInfo reversal
##       <char> <num>     <char> <num>     <num>   <num>   <num>   <num>   <num>   <num>    <num>
##  1:      NTB     0      fixed 20000   0.00000 0.00000 0.00000 0.97345 0.02655       0  0.00000 Net Treatment Benefit (GPC)
##  2:      NTB     0       none 20000   0.02620 0.00320 0.02300 0.96495 0.03505       0  0.00320 Net Treatment Benefit (GPC)
##  3:      NTB     0 bonferroni 20000   0.01320 0.00160 0.01160 0.98105 0.01895       0  0.00160 Net Treatment Benefit (GPC)
##  4:      NTB     0       asOF 20000   0.00995 0.00000 0.00995 0.97330 0.02670       0  0.00000 Net Treatment Benefit (GPC)
##  5:      NTB     0        asP 20000   0.02150 0.00870 0.01280 0.98240 0.01760       0  0.00870 Net Treatment Benefit (GPC)
##  6:      NTB     0    delayed 20000   0.01475 0.00005 0.01470 0.97390 0.02610       0  0.00005 Net Treatment Benefit (GPC)
##  7:       HR     0      fixed 20000   0.00000 0.00000 0.00000 0.97370 0.02630       0  0.00000          Hazard Ratio (Cox)
##  8:       HR     0       none 20000   0.02645 0.00280 0.02365 0.96410 0.03590       0  0.00280          Hazard Ratio (Cox)
##  9:       HR     0 bonferroni 20000   0.01360 0.00130 0.01230 0.98160 0.01840       0  0.00130          Hazard Ratio (Cox)
## 10:       HR     0       asOF 20000   0.00980 0.00000 0.00980 0.97330 0.02670       0  0.00000          Hazard Ratio (Cox)
## 11:       HR     0        asP 20000   0.02220 0.00875 0.01345 0.98250 0.01750       0  0.00875          Hazard Ratio (Cox)
## 12:       HR     0    delayed 20000   0.01435 0.00000 0.01435 0.97295 0.02705       0  0.00000          Hazard Ratio (Cox)
quantile(sapply(1:1e3, function(iSeed){set.seed(iSeed);mean(rbinom(1e4, size = 1, prob = 0.025))}))
##     0%    25%    50%    75%   100% 
## 0.0199 0.0239 0.0250 0.0262 0.0297 

dt.nodelayH0.decision$estimand2 <- sapply(dt.nodelayH0.decision$estimand, switch,
                                          "HR" = "Hazard Ratio (Cox)",
                                          "NTB" = "Net Treatment Benefit (GPC)")
gg.nodelayH0.type1error <- ggplot(dt.nodelayH0.decision, aes(x = rho, y = H1, group = design, shape = design, color = design))
gg.nodelayH0.type1error <- gg.nodelayH0.type1error + geom_point(size = 3) + geom_line(linewidth = 1)
gg.nodelayH0.type1error <- gg.nodelayH0.type1error + facet_wrap(~estimand2)
gg.nodelayH0.type1error <- gg.nodelayH0.type1error + geom_hline(yintercept = 0.025, color = "black", linetype = 2)
gg.nodelayH0.type1error <- gg.nodelayH0.type1error + labs(y = "Type 1 error")
gg.nodelayH0.type1error <- gg.nodelayH0.type1error + theme(text = element_text(size=15), 
                                                           axis.line = element_line(linewidth = 1.25),
                                                           axis.ticks = element_line(linewidth = 2),
                                                           axis.ticks.length=unit(.25, "cm"),
                                                           legend.key.size = unit(2,"line"))
## gg.nodelayH0.type1error


## ** power 
## *** example
dfH1.nodelay <- simTrial(n.E = n.E, n.C = n.C, seed = 1, scale.censoring.C = Inf, scale.C = c(0.01,0.01), scale.E = c(0.011,0.02), admin.censoring = 0.01)
## plot(dfH1.nodelay, facet = ~group)
## plot(dfH1.nodelay, type = "recruitment")

coxH1.nodelay <- coxph(Surv(timeSurv, statusSurv) ~ group, data = dfH1.nodelay)
runTrial(coxH1.nodelay, data = dfH1.nodelay, interim = 0.75, typeOfDesign = "delayed", maxInfo = 49)
GPCH1.nodelay <- BuyseTest(group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox), data = dfH1.nodelay, scoring.rule = "Gehan", trace = FALSE)
runTrial(GPCH1.nodelay, data = dfH1.nodelay, interim = 0.75, typeOfDesign = "delayed", maxInfo = 280)

## warper(iSeed, rho = 0, scale.C = c(0.01,0.01), admin.censoring = 0.01, maxInfo = c(GPC = 301, HR = 63))

## *** planned information
ls.planInfoH1.nodelay <- pblapply(1:min(n.sim,1000), function(iSeed){
    warper(iSeed, rho = 0, scale.C = c(0.01,0.01), scale.E = c(0.011,0.02), admin.censoring = 0.01, typeOfDesign = "fixed")
}, cl = cl)
dtL.planInfoH1.nodelay <- as.data.table(do.call(rbind,ls.planInfoH1.nodelay))
planInfoH1.nodelay <- dtL.planInfoH1.nodelay[stage == "final", .(.N, info = mean(1/se^2), info.min = min(1/se^2), info.max = max(1/se^2)), by = "estimand"]
planInfoH1.nodelay
##    estimand     N      info  info.min  info.max
##      <char> <int>     <num>     <num>     <num>
## 1:      NTB  1080 294.43307 276.77529 302.09842
## 2:       HR  1080  57.28441  49.36116  65.71514

## *** simulation
cat(" - power \n")
ls.nodelayH1 <- pblapply(1:(n.sim/4), function(iSeed){ ## iSeed <- 1
    do.call(rbind,lapply(seqRho, function(iRho){
        warper(iSeed, rho = iRho, scale.C = c(0.01,0.01), scale.E = c(0.011,0.02), admin.censoring = 0.01, maxInfo = c(GPC = 295, HR = 58))
    }))
}, cl = cl)
dt.nodelayH1 <- as.data.table(do.call(rbind,ls.nodelayH1))

## *** information
dtL.nodelayH1.info <- dt.nodelayH1[design=="fixed",.(.N, info = mean(1/se^2), info.min = min(1/se^2), info.max = max(1/se^2), info.pc = mean(info.pc)), by = c("estimand","stage","rho")]
dtL.nodelayH1.info[, info2 := paste0(round(info,2)," (",round(100*info.pc,2),"%) [",round(info.min,2),";",round(info.max,2),"]")]
dtW.nodelayH1.info <- reshape(dtL.nodelayH1.info[,.(N, estimand,stage,rho,info2)], direction = "wide", timevar = "stage", idvar = c("estimand","rho"), v.names = "info2")
dtW.nodelayH1.info[order(dtW.nodelayH1.info$estimand,),]
dtW.nodelayH1.info <- reshape(dtL.nodelayH1.info[,.(N, estimand,stage,rho,info2)], direction = "wide", timevar = "stage", idvar = c("estimand","rho"), v.names = "info2")
dtW.nodelayH1.info[order(dtW.nodelayH1.info$estimand,),]
##         N estimand   rho                   info2.interim                  info2.decision                      info2.final
##     <int>   <char> <num>                          <char>                          <char>                           <char>
## 1:  5000       HR -0.50    42.44 (73.17%) [34.51;50.61]       42.92 (74%) [34.51;50.84]        57.3 (98.8%) [47.74;66.5]
## 2: 10000       HR  0.25     42.47 (73.22%) [33.54;50.4]    42.95 (74.06%) [33.54;50.89]     57.33 (98.85%) [48.23;67.48]
## 3:  5000       HR  0.00    42.47 (73.22%) [34.08;50.37]    42.95 (74.06%) [34.07;51.34]     57.32 (98.83%) [48.57;67.54]
## 4:  5000       HR  0.50      42.46 (73.2%) [32.38;51.3]    42.94 (74.03%) [32.58;52.02]     57.31 (98.81%) [47.64;66.27]
## 5:  5000      NTB -0.50  223.82 (75.87%) [201.4;232.06] 218.04 (73.91%) [196.26;225.02]  291.41 (98.78%) [269.11;299.63]
## 6: 10000      NTB  0.25 226.43 (76.76%) [204.41;236.64] 221.65 (75.14%) [199.74;231.01]  296.12 (100.38%) [271.81;306.9]
## 7:  5000      NTB  0.00 225.28 (76.37%) [204.45;234.43] 220.19 (74.64%) [200.32;228.24]   294.2 (99.73%) [272.25;304.02]
## 8:  5000      NTB  0.50 228.21 (77.36%) [208.83;239.65] 223.75 (75.85%) [203.16;234.34] 298.87 (101.31%) [273.29;312.15]

## *** bound
dtL.nodelayH1.bound <- dt.nodelayH1[rho==0,.(rho = rho[1], .N, bound = mean(bound,na.rm=TRUE), bound.min = min(bound,na.rm=TRUE), bound.max = max(bound,na.rm=TRUE)), by = c("estimand","stage","design")]
dtL.nodelayH1.bound[, bound2 := paste0(round(bound,2)," [",round(bound.min,2),";",round(bound.max,2),"]")]
reshape(dtL.nodelayH1.bound[,.(rho, N,estimand,stage,design,bound2)], direction = "wide", timevar = "stage", idvar = c("estimand","design"), v.names = "bound2")
##       rho     N estimand     design   bound2.interim  bound2.decision     bound2.final
##     <num> <int>   <char>     <char>           <char>           <char>           <char>
##  1:     0  5000      NTB      fixed    Inf [Inf;Inf]      NaN [NA;NA] 1.96 [1.96;1.96]
##  2:     0  5000      NTB       none 1.96 [1.96;1.96] 1.96 [1.96;1.96] 1.96 [1.96;1.96]
##  3:     0  5000      NTB bonferroni 2.24 [2.24;2.24] 2.24 [2.24;2.24] 2.24 [2.24;2.24]
##  4:     0  5000      NTB       asOF 2.31 [2.26;2.45]    2.01 [2;2.02]    2.02 [2;2.02]
##  5:     0  5000      NTB        asP 2.03 [2.02;2.06] 2.26 [2.25;2.26] 2.26 [2.25;2.26]
##  6:     0  5000      NTB    delayed 2.18 [2.15;2.26] 1.88 [1.85;1.89] 2.07 [2.06;2.08]
##  7:     0  5000       HR      fixed    Inf [Inf;Inf]      NaN [NA;NA] 1.96 [1.96;1.96]
##  8:     0  5000       HR       none 1.96 [1.96;1.96] 1.96 [1.96;1.96] 1.96 [1.96;1.96]
##  9:     0  5000       HR bonferroni 2.24 [2.24;2.24] 2.24 [2.24;2.24] 2.24 [2.24;2.24]
## 10:     0  5000       HR       asOF  2.38 [2.14;2.7] 2.01 [1.98;2.05] 2.01 [1.98;2.04]
## 11:     0  5000       HR        asP    2.05 [2;2.11] 2.26 [2.23;2.26] 2.26 [2.23;2.26]
## 12:     0  5000       HR    delayed 2.22 [2.08;2.38] 1.75 [1.56;1.93] 2.07 [1.82;2.09]

## *** Decision
dt.nodelayH1[seed == 1 & design == "fixed"]
dt.nodelayH1.decision <- dt.nodelayH1[,.(nSim = (.N/3),
                                       earlyStop = sum(decision %in% c("stop","infoMax"),na.rm=TRUE)/(.N/3),
                                       earlyH0 = sum(!is.na(decision) & decision == "futility" & stage == "decision", na.rm=TRUE)/(.N/3),
                                       earlyH1 = sum(!is.na(decision) & decision == "efficacy" & stage == "decision", na.rm=TRUE)/(.N/3),
                                       H0 = sum(decision == "futility",na.rm=TRUE)/(.N/3),
                                       H1 = sum(decision == "efficacy",na.rm=TRUE)/(.N/3),
                                       maxInfo = sum(decision == "infoMax",na.rm=TRUE)/(.N/3),
                                       reversal = sum(reversal, na.rm=TRUE)/(.N/3)),
                              by = c("estimand","rho","design")]
dt.nodelayH1.decision[rho==0]
##     estimand   rho     design  nSim earlyStop earlyH0 earlyH1     H0     H1 maxInfo reversal
##       <char> <num>     <char> <num>     <num>   <num>   <num>  <num>  <num>   <num>    <num>
##  1:      NTB     0      fixed  5000    0.0000  0.0000  0.0000 0.2352 0.7648       0   0.0000
##  2:      NTB     0       none  5000    0.6344  0.0186  0.6158 0.2294 0.7706       0   0.0186
##  3:      NTB     0 bonferroni  5000    0.5302  0.0186  0.5116 0.3258 0.6742       0   0.0186
##  4:      NTB     0       asOF  5000    0.5042  0.0004  0.5038 0.2456 0.7544       0   0.0004
##  5:      NTB     0        asP  5000    0.6050  0.0870  0.5180 0.3754 0.6246       0   0.0870
##  6:      NTB     0    delayed  5000    0.5526  0.0002  0.5524 0.2560 0.7440       0   0.0002
##  7:       HR     0      fixed  5000    0.0000  0.0000  0.0000 0.4058 0.5942       0   0.0000
##  8:       HR     0       none  5000    0.4698  0.0150  0.4548 0.3842 0.6158       0   0.0150
##  9:       HR     0 bonferroni  5000    0.3550  0.0148  0.3402 0.4976 0.5024       0   0.0148
## 10:       HR     0       asOF  5000    0.3052  0.0000  0.3052 0.4144 0.5856       0   0.0000
## 11:       HR     0        asP  5000    0.4350  0.0824  0.3526 0.5394 0.4606       0   0.0824
## 12:       HR     0    delayed  5000    0.3658  0.0000  0.3658 0.4242 0.5758       0   0.0000


dt.nodelayH1.decision$estimand2 <- sapply(dt.nodelayH1.decision$estimand, switch,
                                          "HR" = "Hazard Ratio (Cox)",
                                          "NTB" = "Net Treatment Benefit (GPC)")
gg.nodelayH1.power <- ggplot(dt.nodelayH1.decision, aes(x = rho, y = H1, group = design, color = design))
gg.nodelayH1.power <- gg.nodelayH1.power + geom_point(size = 2) + geom_line(linewidth = 1)
gg.nodelayH1.power <- gg.nodelayH1.power + facet_wrap(~estimand2)
gg.nodelayH1.power <- gg.nodelayH1.power + labs(y = "Power")
gg.nodelayH1.power <- gg.nodelayH1.power + theme(text = element_text(size=15), 
                                                 axis.line = element_line(linewidth = 1.25),
                                                 axis.ticks = element_line(linewidth = 2),
                                                 axis.ticks.length=unit(.25, "cm"),
                                                 legend.key.size = unit(2,"line"))
## gg.nodelayH1.power
cat("\n")

## * simulation study (delay) 
cat("With Delay \n")

## ** under the null hypothesis
## *** example 
dfH0.delay <- simTrial(n.E = n.E, n.C = n.C, seed = 1, scale.censoring.C = Inf)
plot(dfH0.delay, facet =~group)
plot(dfH0.delay, type = "recruitment")

coxH0.delay <- coxph(Surv(timeSurv, statusSurv) ~ group, data = dfH0.delay)
runTrial(coxH0.delay, data = dfH0.delay, interim = 0.75, typeOfDesign = "delayed", maxInfo = 63)
GPCH0.delay <- BuyseTest(group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox), data = dfH0.delay, scoring.rule = "Gehan", trace = FALSE)
runTrial(GPCH0.delay, data = dfH0.delay, interim = 0.75, typeOfDesign = "delayed", maxInfo = 300)

## warper(iSeed, rho = 0, maxInfo = c(GPC = 300, HR = 63))

## *** planned information
ls.planInfoH0.delay <- pblapply(1:min(n.sim,1000), function(iSeed){
    warper(iSeed, typeOfDesign = "fixed")
}, cl = cl)
dtL.planInfoH0.delay <- as.data.table(do.call(rbind,ls.planInfoH0.delay))
planInfoH0.delay <- dtL.planInfoH0.delay[stage == "final", .(.N, info = mean(1/se^2), info.min = min(1/se^2), info.max = max(1/se^2)), by = "estimand"]
planInfoH0.delay
##    estimand     N      info  info.min  info.max
##      <char> <int>     <num>     <num>     <num>
## 1:      NTB  1000 300.92834 295.91965 311.18833
## 2:       HR  1000  62.97779  53.68651  70.16954

## *** simulation
cat(" - type 1 error \n")
ls.delayH0 <- pblapply(1:n.sim, function(iSeed){ ## iSeed <- 1
    do.call(rbind,lapply(seqRho, function(iRho){
        warper(iSeed, rho = iRho, maxInfo = c(GPC = 300, HR = 63))
    }))
}, cl = cl)
dt.delayH0 <- as.data.table(do.call(rbind,ls.delayH0))

## *** information
dtL.delayH0.info <- dt.delayH0[design=="fixed",.(.N, info = mean(1/se^2), info.min = min(1/se^2), info.max = max(1/se^2), info.pc = mean(info.pc)), by = c("estimand","stage","rho")]
dtL.delayH0.info[, info2 := paste0(round(info,2)," (",round(100*info.pc,2),"%) [",round(info.min,2),";",round(info.max,2),"]")]
dtW.delayH0.info <- reshape(dtL.delayH0.info[,.(N, estimand,stage,rho,info2)], direction = "wide", timevar = "stage", idvar = c("estimand","rho"), v.names = "info2")
dtW.delayH0.info[order(dtW.delayH0.info$estimand,),]
dtW.delayH0.info <- reshape(dtL.delayH0.info[,.(N, estimand,stage,rho,info2)], direction = "wide", timevar = "stage", idvar = c("estimand","rho"), v.names = "info2")
dtW.delayH0.info[order(dtW.delayH0.info$estimand,),]
##         N estimand   rho                   info2.interim                  info2.decision                      info2.final
##     <int>   <char> <num>                          <char>                          <char>                           <char>
## 1: 20000       HR -0.50      14.29 (22.69%) [8.67;20.67]    47.26 (75.02%) [38.96;55.15]   63.08 (100.12%) [53.87;71.67]
## 2: 40000       HR  0.25      14.31 (22.71%) [8.42;20.95]    47.27 (75.04%) [38.66;55.73]   63.08 (100.13%) [53.99;72.73]
## 3: 20000       HR  0.00       14.3 (22.69%) [8.64;20.68]    47.27 (75.03%) [38.97;55.95]   63.08 (100.13%) [52.99;73.22]
## 4: 20000       HR  0.50      14.31 (22.71%) [8.68;20.49]    47.26 (75.02%) [38.98;55.16]   63.07 (100.11%) [52.72;71.46]
## 5: 20000      NTB -0.50 411.26 (137.09%) [342.14;509.48]  223.02 (74.34%) [210.39;224.8] 298.03 (99.34%) [287.21;299.82]
## 6: 40000      NTB  0.25   379.3 (126.43%) [321.5;471.11]  223.52 (74.51%) [212.31;226.8]  298.67 (99.56%) [285.2;301.59]
## 7: 20000      NTB  0.00 389.06 (129.69%) [322.89;469.36] 223.23 (74.41%) [210.86;225.37] 298.29 (99.43%) [285.06;300.57]
## 8: 20000      NTB  0.50 369.36 (123.12%) [314.59;460.89] 224.06 (74.69%) [212.21;227.67] 299.36 (99.79%) [286.11;303.23]

## *** bound
dtL.delayH0.bound <- dt.delayH0[rho==0,.(rho = rho[1], .N, bound = mean(bound,na.rm=TRUE), bound.min = min(bound,na.rm=TRUE), bound.max = max(bound,na.rm=TRUE)), by = c("estimand","stage","design")]
dtL.delayH0.bound[, bound2 := paste0(round(bound,2)," [",round(bound.min,2),";",round(bound.max,2),"]")]
reshape(dtL.delayH0.bound[,.(rho,.N,estimand,stage,design,bound2)], direction = "wide", timevar = "stage", idvar = c("estimand","design"), v.names = "bound2")
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
## 10:     0    36       HR       asOF 4.58 [3.74;5.94] 2.01 [1.98;2.05] 1.96 [1.96;1.96]
## 11:     0    36       HR        asP  2.4 [2.28;2.56] 2.26 [2.24;2.26] 2.09 [2.04;2.14]
## 12:     0    36       HR    delayed 3.02 [2.78;3.31] 1.55 [1.45;1.64] 1.98 [1.97;1.99]

## *** Decision
dt.delayH0.decision <- dt.delayH0[,.(nSim = (.N/3),
                                       earlyStop = sum(decision %in% c("stop","infoMax"),na.rm=TRUE)/(.N/3),
                                       earlyH0 = sum(!is.na(decision) & decision == "futility" & stage == "decision", na.rm=TRUE)/(.N/3),
                                       earlyH1 = sum(!is.na(decision) & decision == "efficacy" & stage == "decision", na.rm=TRUE)/(.N/3),
                                       H0 = sum(decision == "futility",na.rm=TRUE)/(.N/3),
                                       H1 = sum(decision == "efficacy",na.rm=TRUE)/(.N/3),
                                       maxInfo = sum(decision == "infoMax",na.rm=TRUE)/(.N/3),
                                       reversal = sum(reversal, na.rm=TRUE)/(.N/3)),
                              by = c("estimand","rho","design")]
dt.delayH0.decision[rho==0]
##     estimand   rho     design  nSim earlyStop earlyH0 earlyH1      H0      H1 maxInfo reversal                   estimand2
##       <char> <num>     <char> <num>     <num>   <num>   <num>   <num>   <num>   <num>    <num>                      <char>
##  1:      NTB     0      fixed 20000   1.00000 0.00000 0.00000 0.00000 0.00000       1  0.00000 Net Treatment Benefit (GPC)
##  2:      NTB     0       none 20000   1.00000 0.97375 0.02625 0.97375 0.02625       1  0.00000 Net Treatment Benefit (GPC)
##  3:      NTB     0 bonferroni 20000   1.00000 0.98705 0.01295 0.98705 0.01295       1  0.00000 Net Treatment Benefit (GPC)
##  4:      NTB     0       asOF 20000   1.00000 0.97375 0.02625 0.97375 0.02625       1  0.00000 Net Treatment Benefit (GPC)
##  5:      NTB     0        asP 20000   1.00000 0.97375 0.02625 0.97375 0.02625       1  0.00000 Net Treatment Benefit (GPC)
##  6:      NTB     0    delayed 20000   1.00000 0.97375 0.02625 0.97375 0.02625       1  0.00000 Net Treatment Benefit (GPC)
##  7:       HR     0      fixed 20000   0.00000 0.00000 0.00000 0.97370 0.02630       0  0.00000          Hazard Ratio (Cox)
##  8:       HR     0       none 20000   0.02390 0.01850 0.00540 0.97235 0.02765       0  0.01850          Hazard Ratio (Cox)
##  9:       HR     0 bonferroni 20000   0.01045 0.00840 0.00205 0.98635 0.01365       0  0.00840          Hazard Ratio (Cox)
## 10:       HR     0       asOF 20000   0.00000 0.00000 0.00000 0.97370 0.02630       0  0.00000          Hazard Ratio (Cox)
## 11:       HR     0        asP 20000   0.00610 0.00470 0.00140 0.98020 0.01980       0  0.00470          Hazard Ratio (Cox)
## 12:       HR     0    delayed 20000   0.00080 0.00035 0.00045 0.97430 0.02570       0  0.00035          Hazard Ratio (Cox)
quantile(sapply(1:1e3, function(iSeed){set.seed(iSeed);mean(rbinom(1e4, size = 1, prob = 0.025))}))
##     0%    25%    50%    75%   100% 
## 0.0199 0.0239 0.0250 0.0262 0.0297 

dt.delayH0.decision$estimand2 <- sapply(dt.delayH0.decision$estimand, switch,
                                          "HR" = "Hazard Ratio (Cox)",
                                          "NTB" = "Net Treatment Benefit (GPC)")
gg.delayH0.type1error <- ggplot(dt.delayH0.decision, aes(x = rho, y = H1, group = design, shape = design, color = design))
gg.delayH0.type1error <- gg.delayH0.type1error + geom_point(size = 3) + geom_line(linewidth = 1)
gg.delayH0.type1error <- gg.delayH0.type1error + facet_wrap(~estimand2)
gg.delayH0.type1error <- gg.delayH0.type1error + geom_hline(yintercept = 0.025, color = "black", linetype = 2)
gg.delayH0.type1error <- gg.delayH0.type1error + labs(y = "Type 1 error")
gg.delayH0.type1error <- gg.delayH0.type1error + theme(text = element_text(size=15), 
                                                           axis.line = element_line(linewidth = 1.25),
                                                           axis.ticks = element_line(linewidth = 2),
                                                           axis.ticks.length=unit(.25, "cm"),
                                                           legend.key.size = unit(2,"line"))
## gg.delayH0.type1error


## ** power 
## *** example 
dfH1.delay <- simTrial(n.E = n.E, n.C = n.C, seed = 1, scale.censoring.C = Inf, scale.C = c(0.75,0.5), scale.E = c(0.8,1))
## plot(dfH1.delay, facet = ~group)
## plot(dfH1.delay, type = "recruitment")

coxH1.delay <- coxph(Surv(timeSurv, statusSurv) ~ group, data = dfH1.delay)
runTrial(coxH1.delay, data = dfH1.delay, interim = 0.75, typeOfDesign = "delayed", maxInfo = 60)
GPCH1.delay <- BuyseTest(group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox), data = dfH1.delay, scoring.rule = "Gehan", trace = FALSE)
runTrial(GPCH1.delay, data = dfH1.delay, interim = 0.75, typeOfDesign = "delayed", maxInfo = 295)

## warper(iSeed, rho = 0, scale.C = c(0.01,0.01), admin.censoring = 0.01, maxInfo = c(GPC = 301, HR = 63))

## *** planned information
ls.planInfoH1.delay <- pblapply(1:min(n.sim,1000), function(iSeed){
    warper(iSeed, rho = 0, scale.C = c(0.75,0.5), scale.E = c(0.8,1), typeOfDesign = "fixed")
}, cl = cl)
dtL.planInfoH1.delay <- as.data.table(do.call(rbind,ls.planInfoH1.delay))
planInfoH1.delay <- dtL.planInfoH1.delay[stage == "final", .(.N, info = mean(1/se^2), info.min = min(1/se^2), info.max = max(1/se^2)), by = "estimand"]
planInfoH1.delay
##    estimand     N      info   info.min info.max
##      <char> <int>     <num>     <num>     <num>
## 1:      NTB  1000 295.16836 281.1074 301.0156
## 2:       HR  1000  59.24351  51.6336  67.3829


## *** simulation
cat(" - power \n")
ls.delayH1 <- pblapply(1:(n.sim/4), function(iSeed){ ## iSeed <- 1
    do.call(rbind,lapply(seqRho, function(iRho){
        warper(iSeed, rho = iRho, scale.C = c(0.75,0.5), scale.E = c(0.8,1), maxInfo = c(GPC = 295, HR = 60))
    }))
}, cl = cl)
dt.delayH1 <- as.data.table(do.call(rbind,ls.delayH1))

## *** information
dtL.delayH1.info <- dt.delayH1[design=="fixed",.(.N, info = mean(1/se^2), info.min = min(1/se^2), info.max = max(1/se^2), info.pc = mean(info.pc)), by = c("estimand","stage","rho")]
dtL.delayH1.info[, info2 := paste0(round(info,2)," (",round(100*info.pc,2),"%) [",round(info.min,2),";",round(info.max,2),"]")]
dtW.delayH1.info <- reshape(dtL.delayH1.info[,.(N, estimand,stage,rho,info2)], direction = "wide", timevar = "stage", idvar = c("estimand","rho"), v.names = "info2")
dtW.delayH1.info[order(dtW.delayH1.info$estimand,),]
dtW.delayH1.info <- reshape(dtL.delayH1.info[,.(N, estimand,stage,rho,info2)], direction = "wide", timevar = "stage", idvar = c("estimand","rho"), v.names = "info2")
dtW.delayH1.info[order(dtW.delayH1.info$estimand,),]
##         N estimand   rho                   info2.interim                  info2.decision                      info2.final
##     <int>   <char> <num>                          <char>                          <char>                           <char>
## 1:  5000       HR -0.50          13.2 (22%) [7.95;18.69]   44.39 (73.99%) [36.59;52.19]      59.26 (98.77%) [49.41;68.2]
## 2: 10000       HR  0.25          13.2 (22%) [8.19;18.18]   44.41 (74.02%) [36.63;52.56]     59.28 (98.79%) [49.91;68.63]
## 3:  5000       HR  0.00          13.2 (22%) [7.63;18.84]       44.4 (74%) [36.75;53.36]     59.26 (98.77%) [49.94;69.92]
## 4:  5000       HR  0.50      13.21 (22.02%) [7.33;18.28]        44.4 (74%) [33.9;52.35]      59.26 (98.76%) [49.6;67.93]
## 5:  5000      NTB -0.50 440.71 (149.39%) [364.65;542.93] 219.9 (74.54%) [202.87;224.71]   293.88 (99.62%) [276.2;299.58]
## 6: 10000      NTB  0.25  415.6 (140.88%) [336.63;507.18] 221.55 (75.1%) [202.11;227.22] 296.01 (100.34%) [275.85;302.37]
## 7:  5000      NTB  0.00  423.35 (143.51%) [345.45;513.4] 220.85 (74.87%) [201.6;226.01] 295.12 (100.04%) [277.03;301.12]
## 8:  5000      NTB  0.50 407.63 (138.18%) [331.22;497.76] 222.6 (75.46%) [199.46;229.38] 297.38 (100.81%) [275.13;305.91]

## *** bound
dtL.delayH1.bound <- dt.delayH1[rho==0,.(rho = rho[1], .N, bound = mean(bound,na.rm=TRUE), bound.min = min(bound,na.rm=TRUE), bound.max = max(bound,na.rm=TRUE)), by = c("estimand","stage","design")]
dtL.delayH1.bound[, bound2 := paste0(round(bound,2)," [",round(bound.min,2),";",round(bound.max,2),"]")]
reshape(dtL.delayH1.bound[,.(rho,.N,estimand,stage,design,bound2)], direction = "wide", timevar = "stage", idvar = c("estimand","design"), v.names = "bound2")
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
## 10:     0    36       HR       asOF 4.66 [3.83;6.18] 2.01 [1.98;2.05] 1.96 [1.96;1.96]
## 11:     0    36       HR        asP  2.41 [2.3;2.58] 2.26 [2.23;2.26] 2.09 [2.04;2.13]
## 12:     0    36       HR    delayed 3.04 [2.81;3.35] 1.54 [1.45;1.64] 1.97 [1.97;1.99]

## *** Decision
dt.delayH1.decision <- dt.delayH1[,.(nSim = (.N/3),
                                       earlyStop = sum(decision %in% c("stop","infoMax"),na.rm=TRUE)/(.N/3),
                                       earlyH0 = sum(!is.na(decision) & decision == "futility" & stage == "decision", na.rm=TRUE)/(.N/3),
                                       earlyH1 = sum(!is.na(decision) & decision == "efficacy" & stage == "decision", na.rm=TRUE)/(.N/3),
                                       H0 = sum(decision == "futility",na.rm=TRUE)/(.N/3),
                                       H1 = sum(decision == "efficacy",na.rm=TRUE)/(.N/3),
                                       maxInfo = sum(decision == "infoMax",na.rm=TRUE)/(.N/3),
                                       reversal = sum(reversal, na.rm=TRUE)/(.N/3)),
                              by = c("estimand","rho","design")]
dt.delayH1.decision[rho==0]
##     estimand   rho     design  nSim earlyStop earlyH0 earlyH1     H0     H1 maxInfo reversal
##       <char> <num>     <char> <num>     <num>   <num>   <num>  <num>  <num>   <num>    <num>
##  1:      NTB     0      fixed  5000    1.0000  0.0000  0.0000 0.0000 0.0000       1   0.0000
##  2:      NTB     0       none  5000    1.0000  0.5740  0.4260 0.5740 0.4260       1   0.0000
##  3:      NTB     0 bonferroni  5000    1.0000  0.6734  0.3266 0.6734 0.3266       1   0.0000
##  4:      NTB     0       asOF  5000    1.0000  0.5740  0.4260 0.5740 0.4260       1   0.0000
##  5:      NTB     0        asP  5000    1.0000  0.5740  0.4260 0.5740 0.4260       1   0.0000
##  6:      NTB     0    delayed  5000    1.0000  0.5740  0.4260 0.5740 0.4260       1   0.0000
##  7:       HR     0      fixed  5000    0.0000  0.0000  0.0000 0.6780 0.3220       0   0.0000
##  8:       HR     0       none  5000    0.1068  0.0374  0.0694 0.6830 0.3170       0   0.0374
##  9:       HR     0 bonferroni  5000    0.0616  0.0250  0.0366 0.7754 0.2246       0   0.0250
## 10:       HR     0       asOF  5000    0.0000  0.0000  0.0000 0.6780 0.3220       0   0.0000
## 11:       HR     0        asP  5000    0.0426  0.0156  0.0270 0.7282 0.2718       0   0.0156
## 12:       HR     0    delayed  5000    0.0068  0.0000  0.0068 0.6830 0.3170       0   0.0000

dt.delayH1.decision$estimand2 <- sapply(dt.delayH1.decision$estimand, switch,
                                          "HR" = "Hazard Ratio (Cox)",
                                          "NTB" = "Net Treatment Benefit (GPC)")
gg.delayH1.power <- ggplot(dt.delayH1.decision, aes(x = rho, y = H1, group = design, color = design))
gg.delayH1.power <- gg.delayH1.power + geom_point(size = 2) + geom_line(linewidth = 1)
gg.delayH1.power <- gg.delayH1.power + facet_wrap(~estimand)
gg.delayH1.power <- gg.delayH1.power + labs(y = "Power")
gg.delayH1.power <- gg.delayH1.power + theme(text = element_text(size=15), 
                                             axis.line = element_line(linewidth = 1.25),
                                             axis.ticks = element_line(linewidth = 2),
                                             axis.ticks.length=unit(.25, "cm"),
                                             legend.key.size = unit(2,"line"))
## gg.delayH1.power
cat("\n")

## * export
cat("Export \n")

saveRDS(dtL.nodelayH0.info, file = file.path("../inst/results","nodelay-H0-info.rds")) ## dtL.nodelayH0.info <- readRDS("~/Github/GSD-GPC/inst/results/nodelay-H0-info.rds")
saveRDS(dtL.nodelayH0.bound, file = file.path("../inst/results","nodelay-H0-bound.rds")) ## dtL.nodelayH0.bound <- readRDS("~/Github/GSD-GPC/inst/results/nodelay-H0-bound.rds")
saveRDS(dt.nodelayH0.decision, file = file.path("../inst/results","nodelay-H0-decision.rds")) ## dt.nodelayH0.decision <- readRDS("~/Github/GSD-GPC/inst/results/nodelay-H0-decision.rds")
ggsave(gg.nodelayH0.type1error, file = file.path("../inst/results","nodelay-H0-gg-type1error.pdf"), width = 12, height = 7)

saveRDS(dtL.nodelayH1.info, file = file.path("../inst/results","nodelay-H1-info.rds")) ## dtL.nodelayH1.info <- readRDS("~/Github/GSD-GPC/inst/results/nodelay-H1-info.rds")
saveRDS(dtL.nodelayH1.bound, file = file.path("../inst/results","nodelay-H1-bound.rds"))  ## dtL.nodelayH1.bound <- readRDS("~/Github/GSD-GPC/inst/results/nodelay-H1-bound.rds")
saveRDS(dt.nodelayH1.decision, file = file.path("../inst/results","nodelay-H1-decision.rds")) ## dt.nodelayH1.decision <- readRDS("~/Github/GSD-GPC/inst/results/nodelay-H1-decision.rds")
ggsave(gg.nodelayH1.power, file = file.path("../inst/results","nodelay-H1-gg-power.pdf"), width = 12, height = 7)

saveRDS(dtL.delayH0.info, file = file.path("../inst/results","delay-H0-info.rds")) ## dtL.delayH0.info <- readRDS("~/Github/GSD-GPC/inst/results/delay-H0-info.rds")
saveRDS(dtL.delayH0.bound, file = file.path("../inst/results","delay-H0-bound.rds")) ## dtL.delayH0.bound <- readRDS("~/Github/GSD-GPC/inst/results/delay-H0-bound.rds")
saveRDS(dt.delayH0.decision, file = file.path("../inst/results","delay-H0-decision.rds"))  ## dt.delayH0.decision <- readRDS("~/Github/GSD-GPC/inst/results/delay-H0-decision.rds")
ggsave(gg.delayH0.type1error, file = file.path("../inst/results","delay-H0-gg-type1error.pdf"), width = 12, height = 7)

saveRDS(dtL.delayH1.info, file = file.path("../inst/results","delay-H1-info.rds")) ## dtL.delayH1.info <- readRDS("~/Github/GSD-GPC/inst/results/delay-H1-info.rds")
saveRDS(dtL.delayH1.bound, file = file.path("../inst/results","delay-H1-bound.rds"))  ## dtL.delayH1.bound <- readRDS("~/Github/GSD-GPC/inst/results/delay-H1-bound.rds")
saveRDS(dt.delayH1.decision, file = file.path("../inst/results","delay-H1-decision.rds")) ## dt.delayH1.decision <- readRDS("~/Github/GSD-GPC/inst/results/delay-H1-decision.rds")
ggsave(gg.delayH1.power, file = file.path("../inst/results","delay-H1-gg-power.pdf"), width = 12, height = 7)
cat("Done \n")

##----------------------------------------------------------------------
### test-GPC.R ends here



