## helperfunctions for COVID analyses

add.index <- function(dt)  { # why does this return NA values for ANON_ID? 
    dt$case.index <- 1:nrow(dt)
    dt$case.interval <- c(NA, as.integer(diff(dt$SPECIMENDATE)))
    return(dt[, .(ANON_ID, case.index, case.interval)])
}

year.to2020 <- function(x) {
    lubridate::year(x) <- 2020
    return(x)
}

year.to2021 <- function(x) {
    lubridate::year(x) <- 2021
    return(x)
}

cleanID <- function(filename) {
    namestring <- eval(parse(text=filename))
    print("cleaning", namestring)
    if(grepl("\\.csv$", namestring)) {
        x <- fread(namestring)
    } else {
        x <- as.data.table(readRDS(namestring))
    }
    setnames(x, "anon_id", "ANON_ID", skip_absent=TRUE)
    setnames(x, "CC_ANON_ID", "ANON_ID", skip_absent=TRUE)
    x[, ANON_ID := as.integer(gsub(".* ", "", ANON_ID))] # remove date prefix
    namesstring <- gsub("\\.csv$", "\\.rds", namestring) # save with .rds extension
    print(namestring)
    saveRDS(x, file=namestring)
}

checkID <- function(rdsfilename) {
    namestring <- eval(parse(text=rdsfilename))
    if(!is.null(namestring)) {
        x <- readRDS(namestring)
        setnames(x, "anon_id", "ANON_ID", skip_absent=TRUE)
        setnames(x, "CC_ANON_ID", "ANON_ID", skip_absent=TRUE)
        cat(rdsfilename, x[["ANON_ID"]][1], "\n")
    }
}

recode.dmtype <- function(x) {
    r <- as.factor(car::recode(x, 
                        "c(0, 10, 17, 97)='Not diabetic';
                         c(1, 101, 102)='Type 1 diabetes';
                         c(2, 202, 203)='Type 2 diabetes'; 
                         3:9='Other/unknown type';
                         11:16='Other/unknown type';
                         18:96='Other/unknown type';
                         98:100='Other/unknown type'"))
    r <- factor(r, levels=levels(r)[c(1, 3, 4, 2)])
    return(r)
}

inhosp <- function(admission.daysbefore, discharge.daysbefore, lower, upper) {
        as.integer(ifelse(
        (admission.daysbefore >= lower & admission.daysbefore <= upper) | # admitted during interval
        (discharge.daysbefore >= lower & discharge.daysbefore <=upper) | # discharged during interval
        (admission.daysbefore > upper & (is.na(discharge.daysbefore) |
                                         discharge.daysbefore < lower)), # in hosp throughout
        1, 0))
    }

recode.tw <- function(x) {
        car::recode(x,
                    recodes="1='Less recent TW only'; 2='Recent TW only'; 3='Both TWs'",
                    as.factor=TRUE,
                    levels=c("Less recent TW only", "Recent TW only", "Both TWs")) 
    }

group.tw <- function(lessrecent, recent) {
        timewingr <- ifelse(recent==0 & lessrecent==1, 1, # exposed in less recent window only
                     ifelse(recent==1 & lessrecent==0, 2, # exposed in recent window only
                     ifelse(recent==1 & lessrecent==1, 3, # exposed in both windows
                            NA)))
        return(recode.tw(timewingr))
    }

nextdate <- function(dates1, dates2) {
        ## dates1 and dates2 are vectors of equal length
        ## for each nonmissing value in dates1, get the first date in dates2 that is no earlier than dates1
        nextdate <- dates2
        for(i in 1:length(dates1)) {
            if(!is.na(dates1[i])) {
                dates.noearlier <- dates2[dates2 - dates1[i] >= 0]
                nextdate[i] <- as.Date(min(dates.noearlier, na.rm=TRUE))
            }
        }
        return(as.Date(nextdate))
    }

findintercept <- function(alpha, logit1, logit2) {
    ## alpha is intercept
    ## logit1 is correctly calibrated to give expected cases
    ## logit2 is from conditional model
    expected <- sum(invlogit(logit1))
    expected.fullmodel <- sum(invlogit(alpha + logit1 + logit2))
    return(abs(expected.fullmodel - expected))
}

## covid age is the age at  which logit.norm equates to gam.logit
getCOVIDage <- function(x, y) {
    ## left join x (id, sex, logit.norm) on y(COVID.age, sex, gam.logit)
    ## and get covid age from nearest matched row in y
    #setnames(y, "COVID.age", "COVID.age")
    setkeyv(x, c("Sex", "logit.norm"))
    setkeyv(y, c("Sex", "gam.logit"))
    x <- y[x, roll="nearest"]
    return(x)
}

logit <- function(x) log(x / (1 - x))

invlogit <- function(x) {
    return(x / (1 + exp(-x)))
}

icdToInt <- function(x) {
    N <- length(x)
    int3 <- integer(N)
    lastchar <- substr(x, 3, 3)
    lastchar.asc <- DescTools::CharToAsc(lastchar)
    lastchar.digit <- lastchar.asc >= 48 & lastchar.asc <= 57
    int1 <- as.integer(DescTools::CharToAsc(substr(x, 1, 1)))
    int2 <- as.integer(substr(x, 2, 2))
    int3[lastchar.digit] <- as.integer(lastchar[lastchar.digit]) + 10
                                        # integer range 10 to 19
    int3[!lastchar.digit] <- lastchar.asc[!lastchar.digit]
                                        # integer range 65 to 90
    1000 * int1 + 100 * int2 + int3
}

RDStodt <- function(rds.filename, keyname=NULL) {
    ## read RDS file 
    ## convert numeric to int, character with < 20 unique values to factor
    dt <- as.data.table(readRDS(rds.filename))
    numeric.cols <- names(which(unlist(sapply(dt, is.numeric))))
    if(length(numeric.cols) > 0) {
        integer.cols <- names(which(unlist(sapply(dt[, ..numeric.cols],
                                                  function(x) is.numeric(x) &
                                                              isTRUE(all.equal(x, floor(x)))))))
        dt[, (integer.cols) := lapply(.SD, as.integer), .SDcols = integer.cols]
    }

    factor.cols <- names(which(unlist(sapply(dt,
                                             function(x) is.character(x) &
                                                         length(unique(x)) <= 20))))
    if(length(factor.cols) > 0) {
        dt[, (factor.cols) := lapply(.SD, as.factor), .SDcols = factor.cols]
    }
    if(!is.null(keyname)) {
        setkeyv(dt, keyname)
    }
    return(dt)
}

split.strata <- function(data) { # splits data into person-time intervals that map to weeks of calendar time
    ## final split is done by coxph itself. 
    ## (calendar) startweek is used to generate time-updated variable for use as stratifier
    ## data columns must include startday, startweek, tstop, event
    ## loop over startdays
    cutpoints0 <- c(7, 14, 21, 28)
    startdays <- unique(data$startday)
    data.split <- NULL # data frame
    for(day in startdays) {
        data.split.day <- survSplit(formula=Surv(time=tstop, event=event) ~ .,
                                    data=data[startday==day],
                                    cut=cutpoints0 - day %% 7, # cuts the time argument
                                    zero=0, # already set zero at firstdate
                                    end="tstop",
                                    event="event")
        data.split <- rbind(data.split, data.split.day)
    }
    data.split <- as.data.table(data.split)
    data.split[, week:= startweek + floor((1 + tstart) / 7)]
    return(data.split)
}

dateformat <- function(x) { # format date and strip leading zero
     gsub("^0", "", format(x, "%d %B %Y"))
}

format.estcipv <- function(x.ci, x.pval) {
    x <- gsub(" \\(", " \\(95\\% CI ", as.character(x.ci))
    x <- gsub(", ", " to ", x)
    x <- gsub("\\)", paste0(", _p_=\\", x.pval, "\\)"), x)
    x <- gsub("times", "\\\\times", x)
    return(x)
}

`clogistLoglike`  <-  function(n,m,x,beta) {
    ## https://rdrr.io/cran/saws/src/R/clogistLoglike.R
    ## for computational stability should use matrixStats::logSumExp()
    ## m is vector of 0-1 indicators of failure
    ## n is vector of 1s
    M <- sum(m) 
    N <- sum(n)
    if (M==0)  return(0)
    else if (M==N) return(0)
    x <- as.matrix(x)
    eta <-  x %*% beta
    U <- exp(eta)
    if (M==1) return(sum(eta*m) - log(sum(U*n)) )
    if (M > N/2) {
        ## for efficiency, keep loop part of calculation to minimum
        ## by switching m and n-m, beta and -beta
        m <- n - m
        M <- N - M
        U <- 1/U
        eta <-  -eta
    }
    if (M==1) return(sum(eta*m) - log(sum(U*n)) )
    B <- rep(1, N-M+1)
    u <- rep(NA, N)
    count <- 1
    for (a in 1:length(n)){
        u[count:(count + n[a] - 1)] <- U[a]
        count <- count + n[a]
    }
    ## last 2 lines may be written more 
    ## clearly (i.e. more like in Gail et al) BUT LESS EFFICIENTLY as:
                                        #B <- matrix(0,M+1,N+1)
                                        #B[1,] <- 1
                                        #for (i in 1:M){
                                        #for (j in i:(N-M+i)){
                                        #B[i+1,j+1] <-  B[i+1,j]+u[j]*B[i,j]
                                        #}
                                        #} 
                                        #sum(eta*m) - log(B[M+1,N+1])
    for (i in 1:(M-1)) B <-  cumsum(B * u[i:(N-M+i)])
    sum(eta * m) - log(sum(B * u[M:N]))
}

recode.indicator <- function(x) {
    x[is.na(x)] <- 0
    x[x > 1] <- 1
    return(as.factor(x))
}

pnorm.extreme <- function(z, upper=TRUE) {
    ## https://www.johndcook.com/blog/norm-dist-bounds/
    if(upper) { # upper bound on log10 tail prob
        c <- 8 / pi
    } else { # lower bound
        c = 4 
    }
    x <-  2 * sqrt(2 / pi) / (z + sqrt(z^2 + c))
    ln_p = -0.5 * z^2 + log(x)
    log10p <- ln_p / log(10)
    exponent <- floor(log10p)
    coeff <- 10^(log10p - exponent)
    string <- paste0(round(coeff), "E", exponent)
    return(string)
}

logit <- function(x) log(x) - log(1 - x) 

invlogit <- function(x) 1 / (1 + exp(-x))

tolower.exceptfirstchar <- function(x) {
    lowercase.str <- substr(tolower(x), 2, nchar(x))
    x <- paste0(substr(x, 1, 1), lowercase.str)
    return(x)
}

select.union <- function(x, y, covariates, stratum) {
    ## select j as column of x to drop that maximizes loglik of
    ## clogit model using rowSums([x, -j])
    ## returns coeff and loglik of full model without a dropped column,
    loglik.dropcol <- numeric(ncol(x))
    x.u <- rowSums(x)
    covariates <- as.data.frame(covariates)
    names(covariates) <- paste0("c", 1:ncol(covariates)) 
    cl.formula <- as.formula(paste("y ~ x.u", ifelse(is.null(covariates), "", "+"),
                                   paste(names(covariates), collapse=" + "),
                                   "+ strata(stratum)"))
    cl.data <- data.frame(y=y, x.u=x.u, covariates, stratum=stratum)
    
    model.full <- summary(clogit(formula=cl.formula, data=cl.data))
    beta.full <- model.full$coefficients[1, 1]
    loglik.full <- model.full$loglik[2]
    for(j in 1:ncol(x)) {
        x.u <- rowSums(x[, -j, drop=FALSE])
        cl.formula <- as.formula(paste("y ~ x.u", ifelse(is.null(covariates), "", "+"),
                                       paste(names(covariates), collapse=" + "),
                                       "+ strata(stratum)"))
        cl.data <- data.frame(y=y, x.u=x.u, covariates, stratum=stratum)
        
        loglik.dropcol[j] <- summary(clogit(formula=cl.formula, data=cl.data))$loglik[2]
    }

    j.max <- which.max(loglik.dropcol)
    return(data.frame(col.drop=j.max,
                      name=colnames(x)[j.max],
                      beta.full=beta.full,
                      loglik.full=loglik.full))
}

stepwise.union.dropcols <- function(x, y, covariates=NULL, stratum) {
    x.drop <- x
    stepwise.drop <- NULL
    for(j in 1:(ncol(x) - 1)) {
        select.drop <- select.union(x=x.drop, y=y, covariates=covariates, stratum=stratum)
        stepwise.drop <- rbind(stepwise.drop, select.drop)
        x.drop <- x.drop[, -select.drop$col.drop, drop=FALSE]
        ## for each dropped variable, stepwise.drop has the coeff and loglik of the full model before that variable was dropped  
    }
    ## add extra row for last variable retained
    last.name <- colnames(x)[!colnames(x) %in% stepwise.drop$name]
    last.col <- match(last.name, colnames(x))
    last.model <- summary(clogit(y ~ x[, last.col] + strata(stratum)))
    stepwise.drop <- rbind(stepwise.drop,
                           data.frame(col.drop=NA,
                                      name=last.name,
                                      beta.full=last.model$coefficients[1, 1],
                                      loglik.full=last.model$loglik[2]))
    loglik.initial <- stepwise.drop$loglik.full[1]
    ## shift the log likelihood values up one row so that the values are
    ## the log lik after dropping the corresponding variable
    stepwise.drop$loglik.full <- c(stepwise.drop$loglik.full[-1], NA)
    ## append the log-likelihood of the null model
    stepwise.drop$loglik.full[nrow(stepwise.drop)] <- last.model$loglik[1]
    ## scale the log likelihood as difference from loglik of initial full model  
    stepwise.drop$loglik.full <- stepwise.drop$loglik.full - loglik.initial
    return(stepwise.drop[, -1])
}

lookup.names <- data.frame(varname=c("AGE", "sex",
                                     "deathwithin28", "scrip.any", "diag.any", "care.home",
                                     "diag.other", "scripordiag", 
                                     "emerg", "icu.hdu.ccu", "inpat",
                                     "protonpump",
                                     "diabetes.any",
                                     "antihypertensive.any",
                                     "IHD.any",
                                     "CVD",
                                     "heart.other.any",
                                     "circulatory.other",
                                     "oad.any",
                                     "respinf.orTB",
                                     "resp.other",
                                     "cysticfibrosis.any",
                                     "ckd.any",
                                     "neuro.any",
                                     "liver.any",
                                     "listed.any",
                                     "connectivetissue.any",
                                     "mono.poly.neuro",
                                     "epilepsy.any",
                                     "otherneuro.any",
                                     "blood.cancer",
                                     "lung.cancer",
                                     "other.cancer",
                                     "immune.any",
                                     "gabapentinoids",
                                     "antiepileptics.other",
                                     "hh.over18",
                                     "adultsgt2",
                                     "hosp.recent",
                                     "after.letter",
                                     "interval2.shielding",
                                     "interval3.shielding",
                                     "qSIMD.integer",
                                     "preschool.any",
                                     "is.hcw",
                                     "Sgene.dropout",
                                     "av2channels"
                                ),
                           longname=c("Age", "Males",
                                      "Death within 28 days of test",
                                      "Any prescription", "Any admission", "Care home",
                                      "No listed condition, other diagnosis",
                                      "Diagnosis or prescription",
                                      "Emergency admission last year",
                                      "Critical care admission last year",
                                      "Any admission last year",
                                      "Proton pump inhibitor",
                                      "Diabetes (any type)",
                                      "Any antihypertensive",
                                      "Ischaemic heart disease",
                                      "Cerebrovascular disease",
                                      "Other heart disease",
                                      "Other circulatory disease",
                                      "Asthma or chronic airway disease",
                                      "Respiratory infections",
                                      "Other respiratory disease",
                                      "Cystic fibrosis",
                                      "Chronic kidney disease or transplant recipient",
                                      "Neurological (except epilepsy) or dementia",
                                      "Liver disease",
                                      "Any listed condition", 
                                      "Connective tissue disease",
                                      "Neuropathy (mono- or poly-)",
                                      "Epilepsy",
                                      "Other neurological conditions",
                                      "Cancer of blood-forming organs",
                                      "Lung cancer",
                                      "Other cancer",
                                      "Immune deficiency or suppression", 
                                      "Gabapentinoids",
                                      "Other drugs used for epilepsy",
                                      "Number of adults in household",
                                      ">2 adults in household",
                                      "Recent hospital visit/stay",
                                      "Letter sent >= 14 days earlier",
                                      "Shielding-eligible and interval 2",
                                      "Shielding-eligible and interval 3",
                                      "SIMD quintile (integer)",
                                      "At least one child under 5",
                                      "Health-care worker",
                                      "S gene dropout",
                                      "Mean Ct of ORF1ab and N genes"
                                      ))

clean.header <- function(x) {
    x <- gsub("\\.\\.", " (", x)
    x <- gsub("\\.",  ")", x)
    return(x)
}

tabulate.freqs.regressions <- function(varnames, outcome="CASE", data,
                                       model="clogit", stratumvar="stratum") {
    ## table freqs should identify sparse variables to be dropped from varnames
    table.freqs <- univariate.tabulate(varnames=varnames, outcome=outcome, data=data,
                                       drop.reflevel=FALSE, drop.sparserows=FALSE,
                                       drop.emptyrows=TRUE, minrowsum=10)
    ## table freqs will drop rows that do not have at least one non-reference
    ## row with total >= minrowsum
    ## problem is to fix univariate.clogit so that it does not return an error
    ## when clogit returns infinite rate ratio.   
    univariate.table <-
        univariate.clogit(varnames=varnames,
                          outcome=outcome, data=data, add.reflevel=TRUE, model=model,
                          stratumvar=stratumvar)
    multivariate.table <-
        multivariate.clogit(varnames=varnames,
                            outcome=outcome, data=data, add.reflevel=TRUE, model=model,
                            stratumvar=stratumvar)
    table.aug <- combine.tables3(table.freqs, univariate.table, multivariate.table)
    rownames(table.aug) <- replace.names(rownames(table.aug))
    return(table.aug)
}

univariate.clogit <- function(varnames, outcome="CASE", data, add.reflevel=FALSE,
                              model="clogit", stratumvar="stratum") {
    univariate.table <- NULL
    for(i in 1:length(varnames)) {
        xvar <- data[[varnames[i]]]
        numrows.i <- ifelse(is.factor(xvar),
                            length(levels(as.factor(as.integer(xvar)))) - 1, 1)
        ## compute regression only if min count of 5 for one of 2 levels, or > 2 levels
        xtab <- table(xvar)
        if(length(xtab) > 2  | length(xtab) == 2 & min(xtab) >=5) {
            if(model=="clogit") { # conditional logistic regression using stratumvar
                univariate.formula <-
                    as.formula(paste0(outcome, " ~ xvar + strata(", stratumvar, ")"))
                x <- summary(clogit(formula=univariate.formula, data=data))$coefficients
            } else if(model=="logistic") {
                ## unconditional logistic regression: stratumvar must be a factor
                univariate.formula <-
                    as.formula(paste(outcome, "~ xvar + ", stratumvar))
                x <- summary(glm(formula=univariate.formula, data=data,
                                 family="binomial"))$coefficients[2:(1 + numrows.i), , drop=FALSE]
            } else { # Cox regression using data in counting process format stratifying by stratumvar
                ## data must have cols tstart, tstop, event
                univariate.formula <-
                    as.formula(paste0("Surv(time=tstart, time2=tstop, event=event) ~ xvar + strata(", stratumvar, ")"))
                x <- summary(coxph(formula=univariate.formula, data=data))$coefficients
            }
            
            x.colnames <- colnames(x)
        } else { # do not run regression if only 2 levels and min count < 5
            x <- matrix(rep(NA, 5), nrow=1)  ##FIXME: will not work for > 2 levels
        }
        if(!is.factor(xvar)) {
            rownames(x) <- varnames[i]
        } else {
            if(length(levels(xvar)) > 2 & add.reflevel) {
                x <- rbind(rep(NA, ncol(x)), x) # extra line for reference level
                rownames(x)[1] <- with(data, levels(xvar)[1])
                colnames(x) <- x.colnames
            }
        }
        univariate.table <- rbind(univariate.table, x)
    }
    return(univariate.table)
}

multivariate.clogit <- function(varnames, outcome="CASE", data, add.reflevel=FALSE,
                                model="clogit", stratumvar="stratum") {

    colnames <- c(outcome, stratumvar, varnames)
    if(model=="cox") {
        colnames <- c("tstart", "tstop", colnames)
    }
    data.selected <- na.omit(data[, ..colnames])

    ## -stratum ensures that stratum is not included in the model as a covariate
    if(model=="clogit") {
        ## remove stratum from list of covariates
        multivariate.formula <- as.formula(paste0(outcome, "~ . -", stratumvar,
                                                  " + strata(", stratumvar, ")"))
        multivariate.coeffs <- summary(clogit(formula=multivariate.formula,
                                              data=data.selected))$coefficients
    } else if(model=="logistic") { # logistic regression model
        multivariate.formula <- as.formula(paste(outcome, "~ ."))
        multivariate.coeffs <- summary(glm(formula=multivariate.formula,
                                           data=data.selected,
                                           family="binomial"))$coefficients[-1, , drop=FALSE]
        keep.rows <- grep(paste0("^", stratumvar), rownames(multivariate.coeffs), invert=TRUE)
        multivariate.coeffs <- multivariate.coeffs[keep.rows, , drop=FALSE]
    } else { # Cox regression using data in counting process format stratifying by stratumvar
        ## data must have cols tstart, tstop, event
        multivariate.formula <- as.formula(paste0("Surv(time=tstart, time2=tstop, event=event) ~ . -",
                                                  stratumvar,
                                                  " + strata(", stratumvar, ")"))
        print(multivariate.formula)
        multivariate.coeffs <- summary(coxph(formula=multivariate.formula, data=data.selected))$coefficients
    }
 

    multivariate.table <- NULL
    for(i in 1:length(varnames)) {
        ## numrows.i is number of rows in multivariate.coeffs for varnames[i]
        ## numrows.i is 1 for numeric variables
        ## is length(levels) -1 for factor variables
        xvar <- data.selected[[varnames[i]]]
        numrows.i <- ifelse(is.factor(xvar),
                            length(levels(as.factor(as.integer(xvar)))) - 1, 1)
        ## if variable is factor with > 2 levels and add.reflevel,
        ## add extra line to multivariate.table for reference level
        if(is.factor(xvar) & length(levels(xvar)) > 2 &
           add.reflevel) {
            ## add empty line to multivariate.table
            multivariate.table <- rbind(multivariate.table,
                                        rep(NA, ncol(multivariate.coeffs)))
            ## label this empty line as reference level of factor
            rownames(multivariate.table)[nrow(multivariate.table)] <-
                levels(xvar)[1]
        }
        ## label rows of multivariate.coeffs that will be added
        ## this step is run irrespective of add.reflevel
        ## if variable is factor with >2 levels
        if(is.factor(xvar) & length(levels(xvar)) > 2) {
            ## label rows with levels
            ## FIXME: should drop levels with 0 observations
            keep.levels <- which(as.integer(table(xvar)) > 0)
            numrows.kept <- length(keep.levels) - 1
            rownames(multivariate.coeffs)[1:numrows.kept] <- levels(xvar)[keep.levels][-1]
        } else { # if numeric, or factor with 2 levels
            ## label single row with variable name
            rownames(multivariate.coeffs)[1] <- varnames[i]
        }
        ## add rows from multivariate.coeffs    
        multivariate.table <- rbind(multivariate.table,
                                    multivariate.coeffs[1:numrows.i, , drop=FALSE])
        multivariate.coeffs <- multivariate.coeffs[-(1:numrows.i), , drop=FALSE]
    }
    return(multivariate.table)
}

tabulate.bnfsection <- function(chnum, outcome="CASE", data=cc.severe, minrowsum=20) {
    ## get sectioncode in wide format, one col per subchapter
    scrips.section.wide <- data.table::dcast(scrips[scrips$chapternum==chnum, ],
                                         ANON_ID ~ sectioncode, fun.aggregate=length,
                                         value.var="sectioncode")
    colnames(scrips.section.wide)[-1] <- paste0("section.",
                                              as.integer(colnames(scrips.section.wide)[-1]))
  
    ## drop rare sectioncodes
    cols.keep <- c(TRUE, colSums(scrips.section.wide)[-1] >= 20)
    scrips.section.wide <- scrips.section.wide[, ..cols.keep]

    setkey(scrips.section.wide, ANON_ID)
    cc.bnf.section <- merge(data[, c("ANON_ID", "CASE", "stratum")],
                       scrips.section.wide,
                       by="ANON_ID", all.x=TRUE)
    ## now fix colnames and set missing to 0
    bnfcols <- grep("^section.", colnames(cc.bnf.section))
    for(j in bnfcols) {
        cc.bnf.section[[j]][is.na(cc.bnf.section[[j]])] <- 0
        cc.bnf.section[[j]][cc.bnf.section[[j]] > 1] <- 1
        cc.bnf.section[[j]] <- as.factor(cc.bnf.section[[j]])
    }
    
    bnf.section <- colnames(cc.bnf.section)[bnfcols]
    ## FIXME: regressions should be fixed to use only rows retained by univariate.tabulate
    table.bnf.section <- univariate.tabulate(varnames=bnf.section, outcome=outcome, data=cc.bnf.section,
                                        drop.sparserows=TRUE, minrowsum=20)
    if(nrow(table.bnf.section) >= 1) { 
        bnf.section <- rownames(table.bnf.section)  
        subsectioncodes <- as.integer(gsub("section.", "", rownames(table.bnf.section)))
        rownames(table.bnf.section) <-
            bnfcodes$sectionname[match(subsectioncodes, bnfcodes$sectioncode)]
        univariate.bnf.section <- NULL
        for(i in 1:length(bnf.section)) {
            univariate.formula <- as.formula(paste(outcome, "~", bnf.section[i], " + strata(stratum)"))
            x <- summary(clogit(formula=univariate.formula, data=cc.bnf.section))$coefficients
            univariate.bnf.section <- rbind(univariate.bnf.section, x)
        }
        
        table.bnf.section.aug <- combine.tables2(table.bnf.section, univariate.bnf.section)
        return(table.bnf.section.aug) 
    } else {
        return(NULL)
    }
}

tabulate.bnfparas <- function(chnum, outcome="CASE", data=cc.severe, minrowsum=20) {
    ## get sectioncode in wide format, one col per subchapter
    scrips.para.wide <- data.table::dcast(scrips[scrips$chapternum==chnum, ],
                                         ANON_ID ~ paracode, fun.aggregate=length,
                                         value.var="paracode")
    colnames(scrips.para.wide)[-1] <- paste0("para.",
                                              as.integer(colnames(scrips.para.wide)[-1]))

    ## drop rare paragraphcodes
    cols.keep <- c(TRUE, colSums(scrips.para.wide)[-1] >= 20)
    scrips.para.wide[, ..cols.keep]
    
    scrips.para.wide <- as.data.table(scrips.para.wide, key="ANON_ID")
    cc.bnf.para <- scrips.para.wide[data[, c("ANON_ID", "CASE", "stratum")]]
    ## now fix colnames and set missing to 0
    bnfcols <- as.integer(grep("^para.", colnames(cc.bnf.para)))
    cc.bnf.para[, (bnfcols) := lapply(.SD, recode.indicator), .SDcols = bnfcols]

    bnf.para <- colnames(cc.bnf.para)[bnfcols]
    
    table.bnf.para <- univariate.tabulate(varnames=bnf.para, outcome="CASE", data=cc.bnf.para,
                                        drop.sparserows=TRUE, minrowsum=minrowsum)
    ## regressions fixed to use only rows retained by univariate.tabulate
    if(nrow(table.bnf.para) >= 1) { 
        bnf.para <- rownames(table.bnf.para)  
        paracodes <- as.integer(gsub("para.", "", rownames(table.bnf.para)))
        rownames(table.bnf.para) <-
            bnfparacodes$paraname[match(paracodes, bnfparacodes$paracode)]
        univariate.bnf.para <- NULL
        for(i in 1:length(bnf.para)) {
            univariate.formula <- as.formula(paste(outcome, "~", bnf.para[i], " + strata(stratum)"))
            x <- summary(clogit(formula=univariate.formula, data=cc.bnf.para))$coefficients
            univariate.bnf.para <- rbind(univariate.bnf.para, x)
        }

        multivariate.bnf.para <- multivariate.clogit(varnames=bnf.para,
                                                   data=cc.bnf.para, add.reflevel=FALSE)

        table.bnf.para.aug <- combine.tables3(table.bnf.para, univariate.bnf.para,
                                            multivariate.bnf.para)
        return(table.bnf.para.aug) 
    } else {
        return(NULL)
    }
}

merge.bnfsubparas <- function(chnums, data) {
    ## merge with data all drug subparas in chapter numbers in vector chnums
    ## subparacode has 7 digits, leading zeroes stripped in scrips
    scrips.subpara.wide <- data.table::dcast(scrips[scrips$chapternum %in% chnums, ],
                                      ANON_ID ~ bnf_paragraph_code, fun.aggregate=length,
                                      value.var="bnf_paragraph_code")
    names.subparas <-
        bnfsubparacodes$subparaname[match(as.integer(colnames(scrips.subpara.wide)[-1]),
                                          bnfsubparacodes$subparacode)]
    colnames(scrips.subpara.wide)[-1] <- paste("subpara",
                                          as.integer(colnames(scrips.subpara.wide)[-1]),
                                          names.subparas, sep=".")
    ## drop rare subparagraphcodes
    scrips.subpara.wide <- subset(scrips.subpara.wide,
                                  select=colSums(scrips.subpara.wide) > 20)

    setkey(scrips.subpara.wide, ANON_ID)
    cc.bnf.subpara <- scrips.subpara.wide[data]
    ## fix colnames and recode indicator variables
    bnfcols <- as.integer(grep("^subpara.", colnames(cc.bnf.subpara)))
    ## recode indicator variables
    cc.bnf.subpara[, (bnfcols) := lapply(.SD, recode.indicator), .SDcols = bnfcols]
    return(cc.bnf.subpara)
}

tabulate.bnfsubparas <- function(chnum, outcome="CASE", data=cc.severe, minrowsum=20) {
    ## get subparacode in wide format, one col per subchapter
    ## subpara code is misnamed bnf_paragraph_code
    scrips.subpara.wide <- data.table::dcast(scrips[scrips$chapternum==chnum, ],
                                         ANON_ID ~ bnf_paragraph_code, fun.aggregate=length,
                                         value.var="bnf_paragraph_code")
    colnames(scrips.subpara.wide)[-1] <- paste0("subpara.",
                                              as.integer(colnames(scrips.subpara.wide)[-1]))
    ## drop rare subparagraphcodes
    cols.keep <- c(TRUE, colSums(scrips.subpara.wide)[-1] >= 20)
    scrips.subpara.wide <- scrips.subpara.wide[, ..cols.keep]

    scrips.subpara.wide <- as.data.table(scrips.subpara.wide, key="ANON_ID")
    cc.bnf.subpara <- scrips.subpara.wide[data[, c("ANON_ID", "CASE", "stratum")]]
    ## now fix colnames and set missing to 0
    bnfcols <- as.integer(grep("^subpara.", colnames(cc.bnf.subpara)))
    cc.bnf.subpara[, (bnfcols) := lapply(.SD, recode.indicator), .SDcols = bnfcols]

    bnf.subpara <- colnames(cc.bnf.subpara)[bnfcols]
    table.bnf.subpara <- univariate.tabulate(varnames=bnf.subpara, outcome="CASE",
                                             data=cc.bnf.subpara,
                                             drop.sparserows=TRUE, minrowsum=minrowsum)
    ## regressions fixed to use only rows retained by univariate.tabulate
    if(!is.null(table.bnf.subpara) & nrow(table.bnf.subpara) >= 1) { 
        bnf.subpara <- rownames(table.bnf.subpara)
        ## get subpara names
        subparacodes <- as.integer(gsub("subpara.", "", rownames(table.bnf.subpara)))
        rownames(table.bnf.subpara) <-
            bnfsubparacodes$subparaname[match(subparacodes, bnfsubparacodes$subparacode)]
        univariate.bnf.subpara <- NULL
        for(i in 1:length(bnf.subpara)) {
            univariate.formula <- as.formula(paste(outcome, "~", bnf.subpara[i], " + strata(stratum)"))
            x <- summary(clogit(formula=univariate.formula, data=cc.bnf.subpara))$coefficients
            univariate.bnf.subpara <- rbind(univariate.bnf.subpara, x)
        }

        multivariate.bnf.subpara <- multivariate.clogit(varnames=bnf.subpara,
                                                   data=cc.bnf.subpara, add.reflevel=FALSE)

        table.bnf.subpara.aug <- combine.tables3(table.bnf.subpara, univariate.bnf.subpara,
                                            multivariate.bnf.subpara)
        return(table.bnf.subpara.aug) 
    } else {
        return(NULL)
    }
}

tabulate.bnfchemicals <- function(chnum, subpara=NULL, subpara.exclude=NULL, outcome="CASE", data=cc.severe, minrowsum=50) {
    ## get drug names in wide format,one col per approved_name
    chnum.str <- sprintf("%02d", chnum)
    scrips.chem <- scrips[substr(bnf_paragraph_code, 1, 2)==chnum.str]
    if(!is.null(subpara)) {
        scrips.chem <- scrips.chem[bnf_paragraph_code == subpara]
    }
    if(!is.null(subpara.exclude)) {
        scrips.chem <- scrips.chem[bnf_paragraph_code != subpara.exclude]
    }
    scrips.chem.wide <- data.table::dcast(scrips.chem,
                                         ANON_ID ~ approved_name, fun.aggregate=length,
                                         value.var="approved_name")

    subpara.chem <- unique(scrips.chem[, .(approved_name, bnf_paragraph_code)])
    setkey(subpara.chem, bnf_paragraph_code) 
    setcolorder(scrips.chem.wide, subpara.chem[, approved_name])
    
    ## drop rare chemicalgraphcodes
    cols.keep <- c(TRUE, colSums(scrips.chem.wide)[-1] >= 20)
    scrips.chem.wide <- scrips.chem.wide[, ..cols.keep]

    scrips.chem.wide <- as.data.table(scrips.chem.wide, key="ANON_ID")
    cc.bnf.chem <- scrips.chem.wide[data[, c("ANON_ID", "CASE", "stratum")]]
    ## fix colnames and recode indicator variables
    bnfcols <- as.integer(which(colnames(cc.bnf.chem) %in% colnames(scrips.chem.wide)[-1]))
    cc.bnf.chem[, (bnfcols) := lapply(.SD, recode.indicator), .SDcols = bnfcols]

    ## approved names contain spaces that have to be replaced with . for use in formula
    colnames(cc.bnf.chem) <- gsub(" |-", "\\.", colnames(cc.bnf.chem))
    
    bnf.chem <- colnames(cc.bnf.chem)[bnfcols]

    table.bnf.chem <- univariate.tabulate(varnames=bnf.chem, outcome=outcome, data=cc.bnf.chem,
                                        drop.sparserows=TRUE, minrowsum=minrowsum)
    ## regressions fixed to use only rows retained by univariate.tabulate
    if(nrow(table.bnf.chem) >= 1) { 
        bnf.chem <- rownames(table.bnf.chem)
        univariate.bnf.chem <- NULL
        for(i in 1:length(bnf.chem)) {
            univariate.formula <- as.formula(paste(outcome, "~", bnf.chem[i], " + strata(stratum)"))
            x <- summary(clogit(formula=univariate.formula, data=cc.bnf.chem))$coefficients
            univariate.bnf.chem <- rbind(univariate.bnf.chem, x)
        }

        multivariate.bnf.chem <- multivariate.clogit(varnames=bnf.chem,
                                                   data=cc.bnf.chem, add.reflevel=FALSE)

        table.bnf.chem.aug <- combine.tables3(table.bnf.chem, univariate.bnf.chem,
                                            multivariate.bnf.chem)
        return(table.bnf.chem.aug) 
    } else {
        return(NULL)
    }
}

tabulate.icdchapter <- function(chnum, data=cc.severe, minrowsum=20) {
    ## get chapternum in wide format, one col per subchapter
    diagnoses.ch.wide <- data.table::dcast(diagnoses[diagnoses$chapter==chnum, ],
                                         ANON_ID ~ subchapter, fun.aggregate=length,
                                         value.var="subchapter")
    colnames(diagnoses.ch.wide)[-1] <- paste0("subCh.",
                                              as.integer(colnames(diagnoses.ch.wide)[-1]))
    ## drop rare subchapters
    diagnoses.ch.wide <- subset(diagnoses.ch.wide, select=colSums(diagnoses.ch.wide) > 10) # drop=FALSE]

    diagnoses.ch.wide <- as.data.table(diagnoses.ch.wide, key="ANON_ID")
    if(ncol(diagnoses.ch.wide) > 1) {
        cc.icd.ch <- diagnoses.ch.wide[data[, c("ANON_ID", "CASE", "stratum")]]

        ## set missing to 0
        icdcols <- as.integer(grep("^subCh.", colnames(cc.icd.ch)))
        cc.icd.ch[, (icdcols) := lapply(.SD, recode.indicator), .SDcols = icdcols]
        
        icd.ch <- colnames(cc.icd.ch)[icdcols]
        ## FIXME: this function should be fixed to drop rows with small numbers
        ## regressions should be fixed to use only rows retained by univariate.tabulate
        table.icd.ch <- univariate.tabulate(varnames=icd.ch, outcome="CASE", data=cc.icd.ch,
                                            drop.sparserows=TRUE, minrowsum=minrowsum)
        if(!is.null(table.icd.ch)) {
            if(nrow(table.icd.ch) >= 1) { 
                icd.ch <- rownames(table.icd.ch)  ## loop over rows of table.icd.ch 
                subchapternums <- as.integer(gsub("subCh.", "", rownames(table.icd.ch)))
                rownames(table.icd.ch) <- icdsubchapters$name[subchapternums]
                univariate.icd.ch <- NULL
                for(i in 1:length(icd.ch)) {
                    univariate.formula <- as.formula(paste("CASE ~ ", icd.ch[i], " + strata(stratum)"))
                    x <- summary(clogit(formula=univariate.formula, data=cc.icd.ch))$coefficients
                    univariate.icd.ch <- rbind(univariate.icd.ch, x)
                }
                
                ## this should reformat pvalues
                table.icd.ch.aug <- combine.tables2(ftable=table.icd.ch, utable=univariate.icd.ch)
                #browser("tabulateicd")
                # table.icd.ch.aug[, 4] <- pvalue.latex(table.icd.ch.aug[, 4])
                # print(table.icd.ch.aug) # without this line the first two cols are dropped?
            } else {
                return(NULL)
            }
        } else {
            return(NULL)
        }
        } else {
        return(NULL)
    }
}

format.pvalue <- function(z, pvalue) {
    pvalue.na <- is.na(pvalue)
    ## convert to character so that extreme p-values can be represented in scientific notation
    pvalue <- as.character(signif(pvalue, 1))
    ## calculate exact p-values where R outputs 0
    pvalue[!pvalue.na & pvalue == "0"] <- pnorm.extreme(z[!pvalue.na & pvalue == "0"])
    pvalue <- as.character(pvalue.latex(pvalue))
    return(pvalue)
}

combine.tables3 <- function(ftable, utable, mtable)  {# returns single table from freqs, univariate, multivariate 


    ## redo this to use the format.pvalue function
    se.colnum <- ncol(utable) - 2
    z.colnum <- ncol(utable) - 1
    pval.colnum <- ncol(utable)
    
    u.ci <- or.ci(utable[, 1], utable[, se.colnum])

    pvalue.na <- is.na(utable[, pval.colnum])
    ## convert to character so that extreme p-values can be represented in scientific notation
    u.pvalue <- as.character(signif(utable[, pval.colnum], 1))
    
    u.pvalue[!pvalue.na & u.pvalue == "0"] <- pnorm.extreme(utable[, z.colnum][!pvalue.na &
                                                                        u.pvalue == "0"])
    u.pvalue <- as.character(pvalue.latex(u.pvalue))
    
    mult.ci <- or.ci(mtable[, 1], mtable[, se.colnum])
    mult.pvalue <- signif(mtable[, pval.colnum], 1)
    
    pvalue.na <- is.na(mtable[, pval.colnum])
    ## convert to character so that extreme p-values can be represented in scientific notation
    mult.pvalue <- as.character(signif(mtable[, pval.colnum], 1))
    
    mult.pvalue[!pvalue.na & mult.pvalue == "0"] <-
        pnorm.extreme(mtable[, z.colnum][!pvalue.na & mult.pvalue == "0"])
    mult.pvalue <- as.character(pvalue.latex(mult.pvalue))
     
    table.aug <- data.frame(ftable,
                            u.ci, u.pvalue,
                            mult.ci, mult.pvalue)
    return(table.aug)
}

combine.tables2 <- function(ftable, utable)  {# returns single table from freqs, univariate 
    
    u.ci <- or.ci(utable[, 1], utable[, 3]) 
    u.pvalue <- signif(utable[, 5], 1)
    u.pvalue <- as.character(pvalue.latex(u.pvalue))
    
    table.aug <- data.frame(ftable,
                            u.ci, u.pvalue)
    return(table.aug)
}    

## format a vector of pvalues in LaTeX and return a vector of mode character 
pvalue.latex <- function(x, n=1, nexp=1) {
    ## this function has to be able to handle x whether numeric or character 
    pvalue <- sapply(x, function(z) { # sapply returns a vector applying FUN to each element of x

        if (is.na(z) | is.nan(z)) {
            return(NA)
        } else if(as.numeric(z) >= 0.001) {
            ## return character string to one sig fig, not in scientific notation
            return(as.character(signif(as.numeric(z), 1)))
        } else {
            if(is.numeric(z)) {
                ## rounds to 1 sig fig and convert to character string
                z <- sprintf("%.*E", 0, signif(z, n)) # default is 1 sig fig
            } else {
                z <- toupper(z)
            }
            z <- as.numeric(unlist(strsplit(as.character(z), "E"))) # split z at E
            sprintf("\\ensuremath{%.*f\\times 10^{%0*d}}", 0, z[1], nexp, z[2])
        }
    }
    )
    pvalue <- as.character(pvalue)
    #pvalue[grep("\\\\times", pvalue)] <- "<0.001" ## fix for thresholding at 0.001
    return(pvalue)
}


## create and format a 95% confidence interval around an odds/hazard ratio
or.ci <- function(coeff, se, ndigits=2) {
  ci.lo <- exp(coeff - 1.96 * se)
  ci.up <- exp(coeff + 1.96 * se)
  ndigits <- rep(2, length(coeff))
  ndigits[exp(coeff) > 5] <- 1
  x <- sprintf("%.*f (%.*f, %.*f)", ndigits, exp(coeff), ndigits, ci.lo, ndigits, ci.up)
  x[is.na(coeff)] <- ""
  return(x)
}

paste.vectortomatrix <- function(x, y) {
    matrix(paste(x, y), nrow=nrow(x), dimnames=dimnames(x))
}

univariate.tabulate <- function(varnames, outcome="CASE", data, drop.reflevel=TRUE,
                                drop.sparserows=FALSE, drop.emptyrows=FALSE,
                                minrowsum=10, digits=0,
                                colpercent=TRUE) {
    outcomevar <- data[[outcome]]
    table.varnames <- NULL
    for(i in 1:length(varnames)) {
        keep.x <- TRUE
        ## test whether variable is factor or numeric
        z <- data[[varnames[i]]] 
        if(is.numeric(z)) { # median (IQR) for numeric variables 
            x <- tapply(z, outcomevar,
                        function(x) {
                            xq <- quantile(x, probs=c(0.5, 0.25, 0.75), na.rm=TRUE)
                            dec <- max(1, -log10(abs(xq[1])))
                            xq <- round(xq, dec)
                            return(paste0(xq[1], "(", xq[2], "-", xq[3], ")"))
                        }
                        )
            x <- matrix(x, nrow=1)
            rownames(x) <- varnames[i]
        } else { # freqs for factor variables
            x <- table(z, outcomevar)
            ## keep if at least one factor level has row sum >= minrowsum
            keep.x <- !any(rowSums(x) < minrowsum)
            if(drop.emptyrows) {
                x <- x[rowSums(x) > 0, ]
            } 
            if(colpercent) {
                x <- paste.colpercent(x, digits=digits)
            } else {
                x <- paste.rowpercent(x, digits=digits)
            }
            ## rownames are labelled with levels(varname)
            ## if two levels OR drop.reference level, drop reference level
            ## if single row left, label rows with varname
            ## else keep reference level
            if(nrow(x) == 2  | drop.reflevel) {
                x <- x[-1, , drop=FALSE] # drop reference category
                if(nrow(x) == 1)
                    rownames(x) <- varnames[i]
            } 
        }
        if(!drop.sparserows | keep.x) { # rbind this table variable
            table.varnames <- rbind(table.varnames, x)
        }
    }
    if(!is.null(table.varnames)) {
        if(outcome=="CASE") {
            colnames(table.varnames)[1:2] <- c("Controls", "Cases")
        }
        colnames(table.varnames) <- paste0(colnames(table.varnames),
                                           " (", as.integer(table(outcomevar)), ")")
    }
    return(table.varnames)
}

replace.names <- function(varnames) {
    names <- varnames
    found <- names %in% lookup.names$varname
    names[found] <- as.character(lookup.names$longname[match(names[found],
                                                             lookup.names$varname)])
    return(names)
}

cv.predict <- function(nfold, cv.data, lower.formula, upper.formula) {
    cv.predicted <- NULL
    ## maybe more memory-efficient not to pass cv.data as arg
    cv.predicted <-
        foreach(foldnum=1:nfold,
                .combine=rbind,
                .inorder=FALSE) %dopar% traintest.fold(foldnum=foldnum,
                                                       cv.data=cv.data,
                                                       lower.formula=lower.formula,
                                                       upper.formula=upper.formula)
    return(cv.predicted)
}

traintest.fold <- function(foldnum, cv.data,
                           lower.formula, upper.formula) { # stepwise clogit on training fold, predict on test fold
    test <- cv.data[, test.fold] == foldnum
    train.data <- cv.data[!test]
    ## use do.call to refer to the dataset in the calling environment
    ## otherwise step cannot find train.data
    start.model <- do.call(clogit, args=list(formula=upper.formula, data=train.data))
    stepwise.model <- step(start.model,
                           scope=list(lower=lower.formula, upper=upper.formula),
                           direction="both", trace=-1, method="approximate")
    
    test.data  <- cv.data[test]
    x <- t(with(test.data, table(CASE, stratum)))
    numcontrols.stratum <- x[, 1]
    numcases.stratum <- x[, 2]
    a <- table(numcases.stratum, exclude=NULL)
    b <- table(numcontrols.stratum, exclude=NULL)
    if(length(a) > 1) stop("not all strata contain single case")
    unnorm.p <- predict(object=stepwise.model, newdata=test.data,
                        na.action="na.pass", 
                        type="risk", reference="sample")
    ## normalize within each stratum
    norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                            stratum=test.data$stratum,
                                            y=test.data$CASE)
    return(norm.predicted)
}

nonmissing.obs <- function(x, varnames) { ## subset rows nonmissing for varnames
    keep <- rep(TRUE, nrow(x))
    for(j in 1:length(varnames)) {
        keep[is.na(subset(x, select=match(varnames[j], colnames(x))))] <- FALSE
    }
    return(keep)
}

normalize.predictions <- function(unnorm.p, stratum, y) { 
    ## normalize probs so that they sum to 1 within each stratum
    ## returns data frame with 3 columns: stratum, prior.p, posterior.p

    predictions <- data.table(unnorm.p, stratum, y)
    predictions[, unnorm.p := as.numeric(unnorm.p)]
    predictions[, stratum := as.integer(as.character(stratum))]
    predictions <- na.omit(predictions) 

    ## keep only strata that contain a single case 
    numcases.strata <- predictions[y==1, .N, by=stratum]
    numcases.strata <- numcases.strata[N==1]
    setkey(numcases.strata, stratum)
    setkey(predictions, stratum)
    predictions <- predictions[numcases.strata] 
    
    predictions[, stratum.size := .N, by=stratum]
    predictions[, norm.constant := sum(unnorm.p), by=stratum]
    predictions[, prior.p := 1/stratum.size]
    predictions[, posterior.p := unnorm.p /norm.constant]
    predictions <- predictions[, .(prior.p, posterior.p, y)]
    return(predictions)
}

normalize.logrates <- function(unnorm.lograte, stratum, y) { 
    ## normalize logrates so that exp(lograte) sums to 1 within each stratum
    ## returns data frame with 3 columns: stratum, prior.p, posterior.p

    predictions <- data.table(unnorm.lograte, stratum, y)
    predictions[, unnorm.lograte := as.numeric(unnorm.lograte)]
    predictions <- na.omit(predictions) 

    ## keep only strata that contain a single case 
    numcases.strata <- predictions[y==1, .N, by=stratum]
    numcases.strata <- numcases.strata[N==1]
    setkey(numcases.strata, stratum)
    setkey(predictions, stratum)
    predictions <- predictions[numcases.strata] 
    
    predictions[, stratum.size := .N, by=stratum]
    predictions[, log.normconstant := matrixStats::logSumExp(unnorm.lograte), by=stratum]
    predictions[, prior.p := 1/stratum.size]
    predictions[, posterior.p := exp(unnorm.lograte - log.normconstant)]
    predictions <- predictions[, .(stratum, prior.p, posterior.p, y)]
    return(predictions)
}

append.emptyrow <- function(x) {
    empty=matrix(c(rep.int(NA, length(x))), nrow=1, ncol=ncol(x))  
    colnames(empty) = colnames(x)
    return(rbind(x, empty))
}

paste.colpercent <- function(x, digits=0, escape.pct=TRUE) { # paste column percentages into freq table
    x.colpct <- paste0("(", round(100 * prop.table(x, 2), digits))
    if(escape.pct) {
        x.colpct <- paste0(x.colpct, "\\%)")
    } else {
        x.colpct <- paste0(x.colpct, "%)")
    }
    z <- matrix(paste(x, x.colpct), nrow=nrow(x),
                dimnames=dimnames(x))
    return(z)
}

paste.rowpercent <- function(x, digits=0, escape.pct=TRUE) { # paste row percentages into freq table and return a matrix
    x.rowpct <- paste0("(", round(100 * prop.table(x, 1), digits))
    if(escape.pct) {
        x.rowpct <- paste0(x.rowpct, "\\%)")
    } else {
        x.rowpct <- paste0(x.rowpct, "%)")
    }
    z <- matrix(paste(x, x.rowpct), nrow=nrow(x),
                dimnames=dimnames(x))
    return(z)
}

testfolds.bystratum <- function(stratum, y, nfold) {
    ## keep only strata that contain a single case 
    table.strata <- tapply(y, stratum, sum) == 1
    strata.onecase <- as.integer(names(table.strata)[as.integer(table.strata)==1])
    keep <- stratum %in% strata.onecase

    stratum.unique <- unique(stratum[keep])
    N <- length(stratum.unique)
    test.folds <- data.frame(stratum=stratum.unique,
                             test.fold=1 + sample(1:N, size=N) %% nfold)
    return(test.folds)
}

MakeEthnicOnomap <- function(onolytic, geographic){
  ## Take vectors onolytic and geographic to recode ethnicity and return this as a vector
  require(car)
  if(length(onolytic) != length(geographic)) warning("Input vectors should be of same length")
  onolytic <- car::recode(onolytic,
                          "'NOT FOUND'=NA; 'INTERNATIONAL'=NA; 'UNCLASSIFIED'=NA; 'VOID'=NA; 'VOID - FORENAME'=NA; 'VOID INITIAL'=NA")
  eth <- rep("Other", length(onolytic))
  eth[is.na(onolytic)] <- NA
  eth[onolytic=="CHINESE" |
        onolytic=="HONG KONGESE" |
        onolytic=="SINGAPORESE"] <- "Chinese"
  eth[geographic=="EAST ASIA" &
        eth != "Chinese"] <- "Other Asia & Pacific"
  
  eth[geographic=="SOUTH ASIA"] <- "South Asian"
  
  eth[geographic=="BRITISH ISLES"] <- "Britain&Ireland"
  eth[geographic=="CENTRAL EUROPE" |
        geographic=="EASTERN EUROPE" |
        geographic=="NORTHERN EUROPE" |
        geographic=="SOUTHERN EUROPE" |
        onolytic=="AFRIKAANS"] <- "Other Europe"
  
  eth[geographic=="AFRICA" &
        onolytic != "AFRIKAANS" &
        onolytic != "LIBYAN"] <- "Black African"
  eth[geographic=="MIDDLE EAST" &
        onolytic != "MUSLIM"] <- "East Med"
  eth[onolytic == "MUSLIM"] <- "South Asian" # "Muslim, not localized"
  
  eth[onolytic == "BLACK CARIBBEAN"] <- "Black Caribbean"
  
  ## reduce to 5 categories: White, South Asian, Chinese, Black, Other
  ethnic5 <- car::recode(eth, "'Black African'='Black'; 'Black Caribbean'='Black';  'Britain&Ireland'='White';  'Other Europe'='White';  'East Med'='Other'; 'Other Asia & Pacific'='Other'; 'Muslim, not localized'='Other'")
  ethnic5 <- relevel(as.factor(ethnic5), ref="White")
  ethnic5
}
