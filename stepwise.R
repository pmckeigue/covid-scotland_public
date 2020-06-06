
## this doesn't work when wrapped as a function 
cv.predict <- function(nfold, cv.data, lower.formula, upper.formula) {
    
    cv.predicted <- NULL
    for(i in 1:nfold) {
        test <- cv.data$test.fold == i
        test.data  <- cv.data[test, ]
        train.data <- cv.data[!test, ]
        cat(length(which(test)), "observations in test fold", i, "\n")
        
        x <- t(with(test.data, table(CASE, stratum)))
        numcontrols.stratum <- x[, 1]
        numcases.stratum <- x[, 2]
        a <- table(numcases.stratum, exclude=NULL)
        b <- table(numcontrols.stratum, exclude=NULL)
        if(length(a) > 1) stop("not all strata contain single case")
        
        start.model <- clogit(formula=upper.formula, data=train.data)
        stepwise.model <- step(start.model,
                               scope=list(lower=lower.formula, upper=upper.formula),
                               direction="both", trace=-1, method="approximate")
        
        ## normalize within each stratum
        unnorm.p <- predict(object=stepwise.model, newdata=test.data,
                            na.action="na.pass", 
                            type="risk", reference="sample")
        norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                                stratum=test.data$stratum,
                                                y=test.data$CASE)
        cv.predicted <- rbind(cv.predicted, norm.predicted)
    }
    return(cv.predicted)
}

## stepwise regressions for casecontrol.R

lower.varnames <- "care.home"
demog.varnames <- c("care.home", "SIMD.quintile")
listed.varnames <- c(demog.varnames, "dm.type", listed.conditions)
full.varnames <- c(listed.varnames, "diag.any", "scrip.any", drugs)

lower.formula <- as.formula(paste0("CASE ~ ",
                            paste(lower.varnames, collapse="+"), 
                            " + strata(stratum)"))
demog.formula <- as.formula(paste0("CASE ~ ",
                            paste(demog.varnames, collapse="+"), 
                            " + strata(stratum)"))
listed.formula <- as.formula(paste0("CASE ~ ",
                            paste(listed.varnames, collapse="+"), 
                            " + strata(stratum)"))
full.formula <- as.formula(paste0("CASE ~ ",
                            paste(full.varnames, collapse="+"), 
                            " + strata(stratum)"))

nonmissing <- nonmissing.obs(cc.severe, full.varnames)

## fit models
demog.model <- clogit(formula=demog.formula, data=cc.severe[nonmissing, ])
listed.model <- clogit(formula=listed.formula, data=cc.severe[nonmissing, ])
full.model <- clogit(formula=full.formula, data=cc.severe[nonmissing, ])

if(stepwise) {
    ## stepwise for full variable set
    stepwise.full <- step(full.model,
                          scope=list(lower=lower.formula, upper=full.formula),
                          direction="both", method="approximate", trace=-1)
    stepwise.full <- summary(stepwise.full)$coefficients
    rownames(stepwise.full) <- replace.names(rownames(stepwise.full))
    print(stepwise.full)

################ cross-validation of stepwise regression ######
   
    test.folds <- testfolds.bystratum(stratum=cc.severe[nonmissing, ]$stratum,
                                      y=cc.severe[nonmissing, ]$CASE, nfold=nfold)
    ## test.folds has one row per stratum, one column per fold containing indicator vars
    ## merge by stratum adds single column test.fold
    cv.data <- merge(cc.severe[nonmissing, ], test.folds, by="stratum")
    cv.data <- cv.data[, c("test.fold", "CASE", "stratum", full.varnames)]
    cat("cv.data has", nrow(cv.data), "rows\n")

    ## cross-validation procedure below fails when wrapped into a function ?? why
###########################################
    upper.formula <- demog.formula
#####################################
    cv.predicted <- NULL
    for(i in 1:nfold) {
        test <- cv.data$test.fold == i
        test.data  <- cv.data[test, ]
        train.data <- cv.data[!test, ]
        cat(length(which(test)), "observations in test fold", i, "\n")
        
        x <- t(with(test.data, table(CASE, stratum)))
        numcontrols.stratum <- x[, 1]
        numcases.stratum <- x[, 2]
        a <- table(numcases.stratum, exclude=NULL)
        b <- table(numcontrols.stratum, exclude=NULL)
        if(length(a) > 1) stop("not all strata contain single case")
        
        start.model <- clogit(formula=upper.formula, data=train.data)
        stepwise.model <- step(start.model,
                               scope=list(lower=lower.formula, upper=upper.formula),
                               direction="both", trace=-1, method="approximate")
        
        ## normalize within each stratum
        unnorm.p <- predict(object=stepwise.model, newdata=test.data,
                            na.action="na.pass", 
                            type="risk", reference="sample")
        norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                                stratum=test.data$stratum,
                                                y=test.data$CASE)
        cv.predicted <- rbind(cv.predicted, norm.predicted)
    }
###################################
    demog.predicted <- cv.predicted
    demog.predicted <- demog.predicted[demog.predicted$prior.p < 1, ]
    demog.densities <- with(demog.predicted,
                            Wdensities(y, posterior.p, prior.p,
                                       recalibrate=FALSE))
    pander(summary(demog.densities), table.style="multiline",
           split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
           caption="Prediction of severe COVID-19 from demographic variables only")

########################################
    upper.formula <- listed.formula
#####################################
    cv.predicted <- NULL
    for(i in 1:nfold) {
        test <- cv.data$test.fold == i
        test.data  <- cv.data[test, ]
        train.data <- cv.data[!test, ]
        cat(length(which(test)), "observations in test fold", i, "\n")
        
        x <- t(with(test.data, table(CASE, stratum)))
        numcontrols.stratum <- x[, 1]
        numcases.stratum <- x[, 2]
        a <- table(numcases.stratum, exclude=NULL)
        b <- table(numcontrols.stratum, exclude=NULL)
        if(length(a) > 1) stop("not all strata contain single case")
        
        start.model <- clogit(formula=upper.formula, data=train.data)
        stepwise.model <- step(start.model,
                               scope=list(lower=lower.formula, upper=upper.formula),
                               direction="both", trace=-1, method="approximate")
        
        ## normalize within each stratum
        unnorm.p <- predict(object=stepwise.model, newdata=test.data,
                            na.action="na.pass", 
                            type="risk", reference="sample")
        norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                                stratum=test.data$stratum,
                                                y=test.data$CASE)
        cv.predicted <- rbind(cv.predicted, norm.predicted)
    }
###################################
    listed.predicted <- cv.predicted
    listed.predicted <- listed.predicted[listed.predicted$prior.p < 1, ]
    listed.densities <- with(listed.predicted,
                            Wdensities(y, posterior.p, prior.p,
                                       recalibrate=FALSE))
    pander(summary(listed.densities), table.style="multiline",
           split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
           caption="Prediction of severe COVID-19 from demographic + listed variables")

######################################
       upper.formula <- full.formula
#####################################
    ## this fails when wrapped into a function ???
    cv.predicted <- NULL
    for(i in 1:nfold) {
        test <- cv.data$test.fold == i
        test.data  <- cv.data[test, ]
        train.data <- cv.data[!test, ]
        cat(length(which(test)), "observations in test fold", i, "\n")
        
        x <- t(with(test.data, table(CASE, stratum)))
        numcontrols.stratum <- x[, 1]
        numcases.stratum <- x[, 2]
        a <- table(numcases.stratum, exclude=NULL)
        b <- table(numcontrols.stratum, exclude=NULL)
        if(length(a) > 1) stop("not all strata contain single case")
        
        start.model <- clogit(formula=upper.formula, data=train.data)
        stepwise.model <- step(start.model,
                               scope=list(lower=lower.formula, upper=upper.formula),
                               direction="both", trace=-1, method="approximate")
        
        ## normalize within each stratum
        unnorm.p <- predict(object=stepwise.model, newdata=test.data,
                            na.action="na.pass", 
                            type="risk", reference="sample")
        norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                                stratum=test.data$stratum,
                                                y=test.data$CASE)
        cv.predicted <- rbind(cv.predicted, norm.predicted)
    }
###################################

    full.predicted <- cv.predicted
    full.predicted <- full.predicted[full.predicted$prior.p < 1, ]
    full.densities <- with(full.predicted,
                           Wdensities(y, posterior.p, prior.p,
                                      recalibrate=FALSE))
    pander(summary(full.densities), table.style="multiline",
           split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
           caption="Prediction of severe COVID-19 from full variable set")

    if(old) {
    save(stepwise.full,
         demog.densities, listed.densities, full.densities, 
         file="./data/stepwise.RData")
    } else {
            save(stepwise.full,
         demog.densities, listed.densities, full.densities, 
         file="./data/stepwise_15May.RData")

    }

} else {
    if(old) {
        load("./data/stepwise.RData")
    } else {
         load("./data/stepwise_15May.RData")
    }
}

