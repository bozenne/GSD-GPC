### test-GPC.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 12 2025 (12:13) 
## Version: 
## Last-Updated: jun 16 2025 (18:27) 
##           By: Brice Ozenne
##     Update #: 62
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

Mplan <- do.call(rbind,pblapply(1:100, function(iSeed){ ## iSeed <- 1
    iDf <- simTrial(n.E = n.E, n.C = n.C, seed = iSeed, scale.censoring.C = Inf)
    iDf$statusTox.bin <- (iDf$statusTox == 1) 
    iBT <- BuyseTest(group ~ tte(timeSurv, statusSurv,threshold = 0.1) + bin(statusTox.bin, operator = "<0"),
                     data = iDf, scoring.rule = "Gehan", trace = FALSE)
    iBT2 <- BuyseTest(group ~ cont(timeSurv, threshold = 0.1) + cont(timeTox, operator = "<0"),
                      data = iDf, scoring.rule = "Gehan", trace = FALSE)
    iOut <- cbind(seed = iSeed,
                  endpoint = c(rownames(confint(iBT)),rownames(confint(iBT2))[2]),
                  rbind(confint(iBT),confint(iBT2)[2,]))
    iOut$endpoint <- factor(iOut$endpoint, levels = unique(iOut$endpoint))
    return(iOut)
}, cl = 10))
rownames(Mplan) <- NULL
Infoplan <- tapply(1/Mplan$se^2, Mplan$endpoint, mean)
## 1/tapply(Mplan$se, Mplan$endpoint, mean)^2

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
warper <- function(iSeed, formula = group ~ tte(timeSurv, statusSurv, threshold = 0.1) + bin(statusTox.bin, operator = "<0"), rho  = 0, ...){ ## iSeed <- 84

    ## ** simulate data
    df.final <- simTrial(n.E = n.E, n.C = n.C, seed = iSeed, scale.censoring.C = Inf, 
                         rho.C = rho, rho.E = rho, ...)
    df.final$statusTox.bin <- (df.final$statusTox == 1) 
    ## plot(df.final)
    df.interim <- subset(df.final, interim =  0.75)
    df.interim$statusTox.bin <- (df.interim$statusTox == 1) 
    ## summary(df.interim)
    ## plot(df.interim)
    df.decision <- df.final[df.final$id %in% df.interim$id,]

    ## ** run GPC
    eBT.interim <- BuyseTest(formula, data = df.interim, scoring.rule = "Gehan", trace = FALSE)
    eBT.decision <- BuyseTest(formula, data = df.decision, scoring.rule = "Gehan", trace = FALSE)
    eBT.final <- BuyseTest(formula, data = df.final, scoring.rule = "Gehan", trace = FALSE)

    name.NTB <- tail(names(coef(eBT.interim)),1) ## name of the last endpoint

    ## ** extact information and test statistic
    eBTinfo.interim <- (1/confint(eBT.interim)[name.NTB,"se"]^2)/Infoplan[name.NTB]
    eBTstat.interim <- confint(eBT.interim)[name.NTB,"estimate"]/confint(eBT.interim)[name.NTB,"se"]

    eBTinfo.decision <- (1/confint(eBT.decision)[name.NTB,"se"]^2)/Infoplan[name.NTB]
    eBTstat.decision <- confint(eBT.decision)[name.NTB,"estimate"]/confint(eBT.decision)[name.NTB,"se"]

    eBTinfo.final <- (1/confint(eBT.final)[name.NTB,"se"]^2)/Infoplan[name.NTB]
    eBTstat.final <- confint(eBT.final)[name.NTB,"estimate"]/confint(eBT.final)[name.NTB,"se"]

    ## ** retrieve boundaries
    bound <- critical.threshold
    ## interim
    if(eBTinfo.interim>1){ ## over-run: go to decision and decision is final
        bound["gsd.of",c("interim","decision","final")] <- c(-Inf,bound["fixed","final"],NA)
        bound["gsd.p",c("interim","decision","final")] <- c(-Inf,bound["fixed","final"],NA)
    }else{
        bound["gsd.of",c("interim","final")] <- getDesignGroupSequential(kMax = 2,
                                                                         alpha = 0.025,
                                                                         sided = 1,
                                                                         informationRates = c(eBTinfo.interim,if(eBTinfo.interim < eBTinfo.final){min(eBTinfo.final,1)}else{1}),
                                                                         typeOfDesign = "asOF")$criticalValues
        bound["gsd.p",c("interim","final")] <- getDesignGroupSequential(kMax = 2,
                                                                        alpha = 0.025,
                                                                        sided = 1,
                                                                        informationRates = c(eBTinfo.interim,if(eBTinfo.interim < eBTinfo.final){min(eBTinfo.final,1)}else{1}),
                                                                        typeOfDesign = "asP")$criticalValues
    
        ## decision (deletion approach)
        if(eBTinfo.decision>1){ ## over-run: decision is final
            bound["gsd.of","decision"] <- bound["fixed","final"]
            bound["gsd.p","decision"] <- bound["fixed","final"]
        }else{
            bound["gsd.of","decision"] <- getDesignGroupSequential(kMax = 2,
                                                                   alpha = 0.025,
                                                                   sided = 1,
                                                                   informationRates = c(eBTinfo.decision,1),
                                                                   typeOfDesign = "asOF")$criticalValues[1]
            bound["gsd.p","decision"] <- getDesignGroupSequential(kMax = 2,
                                                                  alpha = 0.025,
                                                                  sided = 1,
                                                                  informationRates = c(eBTinfo.decision,1),
                                                                  typeOfDesign = "asP")$criticalValues[1]
        }
    }

    ## ** retrieve boundaries
    decision.interim <- ifelse(eBTstat.interim >= bound[,"interim"], eBTstat.decision >= bound[,"decision"], NA)
    decision.final <- ifelse(is.na(decision.interim), eBTstat.final >= bound[,"final"], NA)
    reversal.interim <- (eBTstat.interim >= bound[,"interim"] & eBTinfo.interim < 1) & !is.na(decision.interim) & (decision.interim==FALSE)

    out <- data.frame(seed = iSeed, rho = rho,
                      method = rownames(bound),
                      n.interim = NROW(df.interim),
                      n.pipeline = sum(df.interim$pipeline),
                      estimate.interim = confint(eBT.interim)[name.NTB,"estimate"],
                      se.interim = confint(eBT.interim)[name.NTB,"se"],
                      pvalue.interim = confint(eBT.interim)[name.NTB,"p.value"],
                      info.interim = as.double(eBTinfo.interim),
                      bound.interim = bound[,"interim"],
                      info.decision = as.double(eBTinfo.decision),
                      decision.interim = decision.interim,
                      reversal.interim = reversal.interim,
                      estimate.final = confint(eBT.final)[name.NTB,"estimate"],
                      se.final = confint(eBT.final)[name.NTB,"se"],
                      pvalue.final = confint(eBT.final)[name.NTB,"p.value"],
                      bound.final = bound[,"final"],
                      info.final = as.double(eBTinfo.final),
                      decision.final = decision.final,
                      decision = ifelse(is.na(decision.interim),decision.final,decision.interim))

    return(out)
}

warper(1, formula = group ~ tte(timeSurv, statusSurv, threshold = 0.1) + bin(statusTox.bin, operator = "<0"))
warper(1, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox))

## * simulation study
cl <- makeCluster(50)
parallel::clusterExport(cl, varlist = c("path","n.E","n.C","Infoplan","critical.threshold","warper"))
tempo <- parallel::clusterCall(cl = cl, function(x){
    library(BuyseTest)
    library(rpact)
    sapply(list.files(path,full.names=TRUE), source)
})

## ** type 1 error (normal information flow)
ls.simH0 <- pblapply(1:1000, function(iSeed){
    rbind(warper(iSeed, rho = -0.5, scale.E = NULL, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox)),
          warper(iSeed, rho = -0.25, scale.E = NULL, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox)),
          warper(iSeed, rho = 0, scale.E = NULL, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox)),
          warper(iSeed, rho = 0.25, scale.E = NULL, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox)),
          warper(iSeed, rho = 0.5, scale.E = NULL, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox)))
}, cl = cl)

## one sided test: should match alpha/2=0.025
dt.simH0 <- as.data.table(do.call(rbind,ls.simH0))
dt.simH0[,.(.N,
            info = paste(round(100*mean(info.interim),2),round(100*mean(info.decision),2),round(100*mean(info.final),2),sep=";"),
            stop.interim = 100*mean(!is.na(decision.interim)),
            reversal.interim = 100*sum(reversal.interim, na.rm=TRUE)/.N,
            decision.interim = 100*sum(decision.interim, na.rm=TRUE)/.N,
            decision.final = 100*sum(decision.final, na.rm=TRUE)/.N,
            decision = 100*mean(decision)),
         by = c("rho","method")]

##       rho         method     N              info stop.interim reversal.interim decision.interim decision.final decision
##     <num>         <char> <int>            <char>        <num>            <num>            <num>          <num>    <num>
##  1: -0.50          fixed  1000  71.47;71.55;95.4            0                0                0            2.9      2.9
##  2: -0.50       gsd.none  1000  71.47;71.55;95.4            0                0                0            2.9      2.9
##  3: -0.50         gsd.of  1000  71.47;71.55;95.4            0                0                0            2.1      2.1
##  4: -0.50          gsd.p  1000  71.47;71.55;95.4            0                0                0            1.0      1.0
##  5: -0.50 gsd.bonferroni  1000  71.47;71.55;95.4            0                0                0            1.0      1.0
##  6: -0.25          fixed  1000 71.46;71.57;95.43            0                0                0            2.8      2.8
##  7: -0.25       gsd.none  1000 71.46;71.57;95.43            0                0                0            2.8      2.8
##  8: -0.25         gsd.of  1000 71.46;71.57;95.43            0                0                0            2.4      2.4
##  9: -0.25          gsd.p  1000 71.46;71.57;95.43            0                0                0            1.6      1.6
## 10: -0.25 gsd.bonferroni  1000 71.46;71.57;95.43            0                0                0            1.6      1.6
## 11:  0.00          fixed  1000 71.47;71.62;95.49            0                0                0            3.6      3.6
## 12:  0.00       gsd.none  1000 71.47;71.62;95.49            0                0                0            3.6      3.6
## 13:  0.00         gsd.of  1000 71.47;71.62;95.49            0                0                0            2.6      2.6
## 14:  0.00          gsd.p  1000 71.47;71.62;95.49            0                0                0            1.8      1.8
## 15:  0.00 gsd.bonferroni  1000 71.47;71.62;95.49            0                0                0            1.9      1.9
## 16:  0.25          fixed  1000   71.47;71.7;95.6            0                0                0            3.3      3.3
## 17:  0.25       gsd.none  1000   71.47;71.7;95.6            0                0                0            3.3      3.3
## 18:  0.25         gsd.of  1000   71.47;71.7;95.6            0                0                0            3.1      3.1
## 19:  0.25          gsd.p  1000   71.47;71.7;95.6            0                0                0            1.8      1.8
## 20:  0.25 gsd.bonferroni  1000   71.47;71.7;95.6            0                0                0            1.8      1.8
## 21:  0.50          fixed  1000 71.47;71.87;95.82            0                0                0            3.6      3.6
## 22:  0.50       gsd.none  1000 71.47;71.87;95.82            0                0                0            3.6      3.6
## 23:  0.50         gsd.of  1000 71.47;71.87;95.82            0                0                0            3.1      3.1
## 24:  0.50          gsd.p  1000 71.47;71.87;95.82            0                0                0            2.0      2.0
## 25:  0.50 gsd.bonferroni  1000 71.47;71.87;95.82            0                0                0            2.0      2.0
##       rho         method     N              info stop.interim reversal.interim decision.interim decision.final decision

quantile(sapply(1:1e3, function(iSeed){set.seed(iSeed);mean(rbinom(1e4, size = 1, prob = 0.025))}))
##     0%    25%    50%    75%   100% 
## 0.0199 0.0239 0.0250 0.0262 0.0297 

warper(1, rho=0)


## ** power
ls.simH1 <- pblapply(1:1000, function(iSeed){
    rbind(warper(iSeed, rho = -0.5, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox),
                 scale.C = c(1,0.8), scale.E = c(0.9,5), shape.C = c(1,1), shape.E = c(1,1)),
          warper(iSeed, rho = -0.25, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox),
                 scale.C = c(1,0.8), scale.E = c(0.9,5), shape.C = c(1,1), shape.E = c(1,1)),
          warper(iSeed, rho = 0, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox),
                 scale.C = c(1,0.8), scale.E = c(0.9,5), shape.C = c(1,1), shape.E = c(1,1)),
          warper(iSeed, rho = 0.25, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox),
                 scale.C = c(1,0.8), scale.E = c(0.9,5), shape.C = c(1,1), shape.E = c(1,1)),
          warper(iSeed, rho = 0.5, formula = group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox),
                 scale.C = c(1,0.8), scale.E = c(0.9,5), shape.C = c(1,1), shape.E = c(1,1)))
}, cl = cl)

dt.simH1 <- as.data.table(do.call(rbind,ls.simH1))
dt.simH1[,.(.N,
            info = paste(round(100*mean(info.interim),2),round(100*mean(info.decision),2),round(100*mean(info.final),2),sep=";"),
            stop.interim = 100*mean(!is.na(decision.interim)),
            reversal.interim = 100*sum(reversal.interim, na.rm=TRUE)/.N,
            decision.interim = 100*sum(decision.interim, na.rm=TRUE)/.N,
            decision.final = 100*sum(decision.final, na.rm=TRUE)/.N,
            decision = mean(decision)),
         by = c("rho","method")]

## ** debug power




df.final <- simTrial(n.E = 1e3, n.C = 1e3, seed = 1, scale.censoring.C = Inf,
                     scale.C = c(1,0.8), scale.E = c(0.9,5),
                     shape.C = c(1,1), shape.E = c(1,1), rho.C = -0.5, rho.E = -0.5)
df.interim <- subset(df.final, interim = 0.75)
df.decision <- df.final[df.final$id %in% df.interim$id,]
plot(df.final, facet= ~ group)
plot(df.decision, facet= ~ group)
plot(df.interim, facet= ~ group)
summary(BuyseTest(group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox), scoring.rule = "Gehan", data = df.interim, trace = FALSE))
summary(BuyseTest(group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox), scoring.rule = "Gehan", data = df.decision, trace = FALSE))
summary(BuyseTest(group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox), scoring.rule = "Gehan", data = df.final, trace = FALSE))

confint(BuyseTest(group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox), scoring.rule = "Gehan", data = df.interim, trace = FALSE))
confint(BuyseTest(group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox), scoring.rule = "Gehan", data = df.decision, trace = FALSE))
confint(BuyseTest(group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox), scoring.rule = "Gehan", data = df.final, trace = FALSE))

summary(BuyseTest(group ~ tte(timeSurv, statusSurv) + tte(timeTox, statusTox.bin, operator = "<0"), scoring.rule = "Gehan", data = df.interim))

table(df.interim$statusSurv)
table(df.interim$statusTox)

warper(1, rho = -0.5, scale.E = c(0.75,Inf), formula = group ~ cont(timeSurv) + cont(timeTox, operator = "<0"))



 

## example of reversal
df.final <- simTrial(n.E = 10*n.E, n.C = 10*n.C, seed = 1, scale.censoring.C = Inf, scale.E = c(0.6,1))
df.final$statusTox.bin <- (df.final$statusTox == 1) 

df.interim <- subset(df.final, interim = 0.75)

plot(df.final, facet = ~group)
plot(df.interim, facet = ~group)

eBT.interim <- BuyseTest(group ~ tte(timeSurv, statusSurv, threshold = 0.1) + bin(statusTox.bin, operator = "<0"),
                         data = df.interim, scoring.rule = "Gehan", trace = FALSE)
confint(eBT.interim)

eBT.final <- BuyseTest(group ~ tte(timeSurv, statusSurv, threshold = 0.1) + bin(statusTox.bin, operator = "<0"),
                       data = df.final, scoring.rule = "Gehan", trace = FALSE)
confint(eBT.final)

eBT.final <- BuyseTest(group ~ cont(timeSurv) + cont(timeTox, operator = "<0"),
                       data = df.final, trace = FALSE)
confint(eBT.final)

confint(BuyseTest(group ~ cont(timeSurv, threshold = 0.1) + bin(statusTox.bin, operator = "<0"), data = df.final, trace = FALSE))

## ** double check type 1 error
warperRed <- function(iSeed, rho){
    df.final <- simTrial(n.E = n.E, n.C = n.C, seed = iSeed, scale.censoring.C = Inf,
                         rho.C = rho, rho.E = rho)
    df.final$statusTox.bin <- (df.final$statusTox == 1) 

    eBT.final <- BuyseTest(group ~ tte(timeSurv, statusSurv, threshold = 0.1) + bin(statusTox.bin, operator = "<0"),
                           data = df.final, scoring.rule = "Gehan", trace = FALSE)

    out <- cbind(seed = iSeed, rho = rho, outcome = rownames(confint(eBT.final)), confint(eBT.final))
    return(out)
}

warperRed(1, rho = 0.75)

parallel::clusterExport(cl, varlist = c("warperRed"))
ls.simH0Red <- pblapply(1:20000, function(iSeed){
    rbind(warperRed(iSeed, rho = 0),
          warperRed(iSeed, rho = 0.75))
}, cl = cl)

dt.simH0Red <- as.data.table(do.call(rbind,ls.simH0Red))
dt.simH0Red[, .(.N,
                type1 = mean(p.value <= 0.05),
                type1pos = mean(p.value <= 0.05 & estimate > 0),
                type1neg = mean(p.value <= 0.05 & estimate < 0)), by = c("outcome","rho")]
##          outcome   rho     N   type1 type1pos type1neg
##           <char> <num> <int>   <num>    <num>    <num>
## 1: timeSurv_t0.1  0.00 20000 0.05150  0.02485  0.02665
## 2: statusTox.bin  0.00 20000 0.05145  0.02545  0.02600
## 3: timeSurv_t0.1  0.75 20000 0.05450  0.02775  0.02675
## 4: statusTox.bin  0.75 20000 0.05280  0.02680  0.02600

quantile(sapply(1:1e3, function(iSeed){set.seed(iSeed);mean(rbinom(2e4, size = 1, prob = 0.05))}))
##      0%     25%     50%     75%    100% 
## 0.04460 0.04895 0.05005 0.05100 0.05590 
quantile(sapply(1:1e3, function(iSeed){set.seed(iSeed);mean(rbinom(2e4, size = 1, prob = 0.025))}))
##      0%     25%     50%     75%    100% 
## 0.02155 0.02425 0.02495 0.02575 0.02875 
##----------------------------------------------------------------------
### test-GPC.R ends here



