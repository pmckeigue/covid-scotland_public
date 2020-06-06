
## helperfunctions for COVID analyses

tolower.exceptfirstchar <- function(x) {
    lowercase.str <- substr(tolower(x), 2, nchar(x))
    x <- paste0(substr(x, 1, 1), lowercase.str)
    return(x)
}

select.union <- function(x, y, stratum) {
    ## select j as column of x to drop that maximizes loglik of
    ## clogit model using rowSums([x, -j])
    ## returns coeff and loglik of full model without a dropped column,
    loglik.dropcol <- numeric(ncol(x))
    x.u <- rowSums(x) 
    model.full <- summary(clogit(y ~ x.u + strata(stratum)))
    beta.full <- model.full$coefficients[1, 1]
    loglik.full <- model.full$loglik[2]
    for(j in 1:ncol(x)) {
        x.u <- rowSums(x[, -j, drop=FALSE])
        loglik.dropcol[j] <- summary(clogit(y ~ x.u + strata(stratum)))$loglik[2]
    }

    j.max <- which.max(loglik.dropcol)
    return(data.frame(col.drop=j.max,
                      name=colnames(x)[j.max],
                      beta.full=beta.full,
                      loglik.full=loglik.full))
}

stepwise.union.dropcols <- function(x, y, stratum) {
    x.drop <- x
    stepwise.drop <- NULL
    for(j in 1:(ncol(x) - 1)) {
        select.drop <- select.union(x=x.drop, y=y, stratum=stratum)
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

lookup.names <- data.frame(varname=c("deathwithin28", "scrip.any", "diag.any", "care.home",
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
                                     "alpha_adrenoreceptor6.bnf",
                                     "ace6.bnf",
                                     "angio6.bnf",
                                     "renin_angiotensin6.bnf",
                                     "thiazides6.bnf",
                                     "calcium_channel6.bnf",
                                     "antihypertensive.other",
                                     "anticoagulants6.bnf",
                                     "nsaids6.bnf",
                                     "lipid_regulating6.bnf",
                                     "statins6.bnf",
                                     "hydroxychloroquine6.bnf"
                                     ),
                           longname=c("Death within 28 days of test",
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
                                      "alpha-adrenoreceptor blocker",
                                      "ACE inhibitor", "Angiotensin-II receptor blocker",
                                      "ACE or A-IIR inhibitor",
                                      "Thiazides",
                                      "Calcium channel blocker",
                                      "Other antihypertensive",
                                      "Anticoagulants",
                                      "Non-steroidal anti-inflammatory drugs",
                                      "Lipid-regulating agents",
                                      "Statins",
                                      "Hydroxychloroquine"
                                      ))

clean.header <- function(x) {
    x <- gsub("\\.\\.", " (", x)
    x <- gsub("\\.",  ")", x)
    return(x)
}

tabulate.freqs.regressions <- function(varnames, outcome="CASE", data) {
    ## FIXME: first run should identify sparse variables to be dropped from varnames
    table.freqs <- univariate.tabulate(varnames=varnames, outcome=outcome, data=data,
                                       drop.reflevel=FALSE, drop.sparserows=FALSE)
    univariate.table <-
        univariate.clogit(varnames=varnames, outcome=outcome, data=data, add.reflevel=TRUE)
    multivariate.table <-
        multivariate.clogit(varnames=varnames, outcome=outcome, data=data, add.reflevel=TRUE)
    table.aug <- combine.tables3(table.freqs, univariate.table, multivariate.table)
    rownames(table.aug) <- replace.names(rownames(table.aug))
   return(table.aug)
}   

univariate.clogit <- function(varnames, outcome="CASE", data, add.reflevel=FALSE) {
    univariate.table <- NULL
    for(i in 1:length(varnames)) {
        if(length(table(as.numeric(with(data, eval(str2expression(varnames[i])))))) > 1) {
            univariate.formula <-
                as.formula(paste(outcome, "~", varnames[i], "+ strata(stratum)"))
            x <- summary(clogit(formula=univariate.formula, data=data))$coefficients
            x.colnames <- colnames(x)
        } else {
            x <- rep(NA, 5) ##FIXME: will not work for > 2 levels
            
        }
        if(with(data, is.factor(eval(str2expression(varnames[i])))) &
           with(data, length(levels(eval(str2expression(varnames[i]))))) > 2 &
           add.reflevel) {
           x <- rbind(rep(NA, ncol(x)), x)
           rownames(x)[1] <- with(data, levels(eval(str2expression(varnames[i])))[1])
           colnames(x) <- x.colnames
        }
        univariate.table <- rbind(univariate.table, x)
    }
    return(univariate.table)
}

multivariate.clogit <- function(varnames, outcome="CASE", data, add.reflevel=FALSE) {
    multivariate.formula <- as.formula(paste(outcome, "~",
                                  paste(varnames, collapse=" + "),
                                  "+ strata(stratum)"))
    nonmissing <- nonmissing.obs(data, varnames)
    multivariate.model <- clogit(formula=multivariate.formula, data=data[nonmissing, ])
    multivariate.coeffs <- summary(multivariate.model)$coefficients
    multivariate.table <- NULL
   
    for(i in 1:length(varnames)) {
        ## numrows.i is number of rows in multivariate.coeffs for varnames[i]
        ## numrows.i is 1 for numeric variables
        ## is length(levels) -1 for factor variables
        numrows.i <- 1
        ## set number of rows in multivariate.coeffs for factor variable as num levels -1 
        if(with(data[nonmissing, ], is.factor(eval(str2expression(varnames[i]))))) {
            numrows.i <-
                with(data[nonmissing, ], length(levels(eval(str2expression(varnames[i]))))) - 1
        }
        ## if variable is factor, > 2 levels and add.reflevel
        if(with(data[nonmissing, ], is.factor(eval(str2expression(varnames[i])))) &
           with(data[nonmissing, ], length(levels(eval(str2expression(varnames[i]))))) > 2 &
           add.reflevel) {
            ## add empty line to multivariate.table
            multivariate.table <- rbind(multivariate.table,
                                        rep(NA, ncol(multivariate.coeffs)))
            ## label this empty line as reference level of factor
            rownames(multivariate.table)[nrow(multivariate.table)] <-
                 with(data[nonmissing, ], levels(eval(str2expression(varnames[i])))[1])
        }
        ## label rows of multivariate.coeffs that will be added
        ## this step is run irrespective of add.reflevel
        ## if variable is factor with >2 levels
        if(with(data[nonmissing, ], is.factor(eval(str2expression(varnames[i])))) &
           with(data[nonmissing, ], length(levels(eval(str2expression(varnames[i]))))) > 2) {
            ## label rows with levels
            rownames(multivariate.coeffs)[1:numrows.i] <-
                with(data[nonmissing, ], levels(eval(str2expression(varnames[i])))[-1])
        } else { # if numeric, or factor with 2 levels
            ## label single row with variable name
            rownames(multivariate.coeffs)[1] <-
                with(data[nonmissing, ], varnames[i])
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
    scrips.ch.wide <- reshape2::dcast(scrips[scrips$chapternum==chnum, ],
                                         ANON_ID ~ sectioncode, fun.aggregate=length,
                                         value.var="sectioncode")
    colnames(scrips.ch.wide)[-1] <- paste0("section.",
                                              as.integer(colnames(scrips.ch.wide)[-1]))
    ## drop rare sectioncodes
    scrips.ch.wide <- scrips.ch.wide[, colSums(scrips.ch.wide) > 20]
    
    cc.bnf.ch <- merge(data[, c("ANON_ID", "CASE", "stratum")],
                       scrips.ch.wide,
                       by="ANON_ID", all.x=TRUE)
    ## now fix colnames and set missing to 0
    bnfcols <- grep("^section.", colnames(cc.bnf.ch))
    for(j in bnfcols) {
        cc.bnf.ch[, j][is.na(cc.bnf.ch[, j])] <- 0
        cc.bnf.ch[, j][cc.bnf.ch[, j] > 1] <- 1
        cc.bnf.ch[, j] <- as.factor(cc.bnf.ch[, j])
    }
    
    bnf.ch <- colnames(cc.bnf.ch)[-(1:3)]
    ## FIXME: regressions should be fixed to use only rows retained by univariate.tabulate
    table.bnf.ch <- univariate.tabulate(varnames=bnf.ch, outcome=outcome, data=cc.bnf.ch,
                                        drop.sparserows=TRUE, minrowsum=20)
    if(nrow(table.bnf.ch) >= 1) { 
        bnf.ch <- rownames(table.bnf.ch)  
        subsectioncodes <- as.integer(gsub("section.", "", rownames(table.bnf.ch)))
        rownames(table.bnf.ch) <-
            bnfcodes$sectionname[match(subsectioncodes, bnfcodes$sectioncode)]
        univariate.bnf.ch <- NULL
        for(i in 1:length(bnf.ch)) {
            univariate.formula <- as.formula(paste(outcome, "~", bnf.ch[i], " + strata(stratum)"))
            x <- summary(clogit(formula=univariate.formula, data=cc.bnf.ch))$coefficients
            univariate.bnf.ch <- rbind(univariate.bnf.ch, x)
        }
        
        table.bnf.ch.aug <- combine.tables2(table.bnf.ch, univariate.bnf.ch)
        return(table.bnf.ch.aug) 
    } else {
        return(NULL)
    }
}

tabulate.bnfparas <- function(chnum, outcome="CASE", data=cc.severe, minrowsum=20) {
    ## get sectioncode in wide format, one col per subchapter
    scrips.ch.wide <- reshape2::dcast(scrips[scrips$chapternum==chnum, ],
                                         ANON_ID ~ paracode, fun.aggregate=length,
                                         value.var="paracode")
    colnames(scrips.ch.wide)[-1] <- paste0("para.",
                                              as.integer(colnames(scrips.ch.wide)[-1]))
    ## drop rare paragraphcodes
    scrips.ch.wide <- scrips.ch.wide[, colSums(scrips.ch.wide) > 20]
    
    cc.bnf.ch <- merge(data[, c("ANON_ID", "CASE", "stratum")],
                       scrips.ch.wide,
                       by="ANON_ID", all.x=TRUE)
    ## now fix colnames and set missing to 0
    bnfcols <- grep("^para.", colnames(cc.bnf.ch))
    for(j in bnfcols) {
        cc.bnf.ch[, j][is.na(cc.bnf.ch[, j])] <- 0
        cc.bnf.ch[, j][cc.bnf.ch[, j] > 1] <- 1
        cc.bnf.ch[, j] <- as.factor(cc.bnf.ch[, j])
    }
    
    bnf.ch <- colnames(cc.bnf.ch)[-(1:3)]
    table.bnf.ch <- univariate.tabulate(varnames=bnf.ch, outcome="CASE", data=cc.bnf.ch,
                                        drop.sparserows=TRUE, minrowsum=minrowsum)
    ## regressions fixed to use only rows retained by univariate.tabulate
    if(nrow(table.bnf.ch) >= 1) { 
        bnf.ch <- rownames(table.bnf.ch)  
        paracodes <- as.integer(gsub("para.", "", rownames(table.bnf.ch)))
        rownames(table.bnf.ch) <-
            bnfparacodes$paraname[match(paracodes, bnfparacodes$paracode)]
        univariate.bnf.ch <- NULL
        for(i in 1:length(bnf.ch)) {
            univariate.formula <- as.formula(paste(outcome, "~", bnf.ch[i], " + strata(stratum)"))
            x <- summary(clogit(formula=univariate.formula, data=cc.bnf.ch))$coefficients
            univariate.bnf.ch <- rbind(univariate.bnf.ch, x)
        }

        multivariate.bnf.ch <- multivariate.clogit(varnames=bnf.ch,
                                                   data=cc.bnf.ch, add.reflevel=FALSE)

        table.bnf.ch.aug <- combine.tables3(table.bnf.ch, univariate.bnf.ch,
                                            multivariate.bnf.ch)
        return(table.bnf.ch.aug) 
    } else {
        return(NULL)
    }
}

merge.bnfsubparas <- function(chnums, data) {
    ## merge with data all drug subparas in chapter numbers in vector chnums
    ## subparacode has 7 digits, leading zeroes stripped in scrips
    scrips.ch.wide <- reshape2::dcast(scrips[scrips$chapternum %in% chnums, ],
                                      ANON_ID ~ bnf_paragraph_code, fun.aggregate=length,
                                      value.var="bnf_paragraph_code")
    names.subparas <-
        bnfsubparacodes$subparaname[match(as.integer(colnames(scrips.ch.wide)[-1]),
                                          bnfsubparacodes$subparacode)]
    colnames(scrips.ch.wide)[-1] <- paste("subpara",
                                          as.integer(colnames(scrips.ch.wide)[-1]),
                                          names.subparas, sep=".")
    
    ## drop rare subparagraphcodes
    scrips.ch.wide <- scrips.ch.wide[, colSums(scrips.ch.wide) > 20]
    
    cc.bnf.ch <- merge(data, scrips.ch.wide,
                       by="ANON_ID", all.x=TRUE)
    ## now fix colnames and set missing to 0
    bnfcols <- grep("^subpara.", colnames(cc.bnf.ch))
    for(j in bnfcols) {
        cc.bnf.ch[, j][is.na(cc.bnf.ch[, j])] <- 0
        cc.bnf.ch[, j][cc.bnf.ch[, j] > 1] <- 1
        cc.bnf.ch[, j] <- as.factor(cc.bnf.ch[, j])
    }
    return(cc.bnf.ch)
}

tabulate.bnfsubparas <- function(chnum, outcome="CASE", data=cc.severe, minrowsum=20) {
    ## get subparacode in wide format, one col per subchapter
    ## subpara code is misnamed bnf_paragraph_code
    scrips.ch.wide <- reshape2::dcast(scrips[scrips$chapternum==chnum, ],
                                         ANON_ID ~ bnf_paragraph_code, fun.aggregate=length,
                                         value.var="bnf_paragraph_code")
    colnames(scrips.ch.wide)[-1] <- paste0("subpara.",
                                              as.integer(colnames(scrips.ch.wide)[-1]))
    ## drop rare subparagraphcodes
    scrips.ch.wide <- scrips.ch.wide[, colSums(scrips.ch.wide) > 20]
    
    cc.bnf.ch <- merge(data[, c("ANON_ID", "CASE", "stratum")],
                       scrips.ch.wide,
                       by="ANON_ID", all.x=TRUE)
    ## now fix colnames and set missing to 0
    bnfcols <- grep("^subpara.", colnames(cc.bnf.ch))
    for(j in bnfcols) {
        cc.bnf.ch[, j][is.na(cc.bnf.ch[, j])] <- 0
        cc.bnf.ch[, j][cc.bnf.ch[, j] > 1] <- 1
        cc.bnf.ch[, j] <- as.factor(cc.bnf.ch[, j])
    }
    
    bnf.ch <- colnames(cc.bnf.ch)[-(1:3)]
    table.bnf.ch <- univariate.tabulate(varnames=bnf.ch, outcome="CASE", data=cc.bnf.ch,
                                        drop.sparserows=TRUE, minrowsum=minrowsum)
    ## regressions fixed to use only rows retained by univariate.tabulate
    if(nrow(table.bnf.ch) >= 1) { 
        bnf.ch <- rownames(table.bnf.ch)
        ## get subpara names
        subparacodes <- as.integer(gsub("subpara.", "", rownames(table.bnf.ch)))
        rownames(table.bnf.ch) <-
            bnfsubparacodes$subparaname[match(subparacodes, bnfsubparacodes$subparacode)]
        univariate.bnf.ch <- NULL
        for(i in 1:length(bnf.ch)) {
            univariate.formula <- as.formula(paste(outcome, "~", bnf.ch[i], " + strata(stratum)"))
            x <- summary(clogit(formula=univariate.formula, data=cc.bnf.ch))$coefficients
            univariate.bnf.ch <- rbind(univariate.bnf.ch, x)
        }

        multivariate.bnf.ch <- multivariate.clogit(varnames=bnf.ch,
                                                   data=cc.bnf.ch, add.reflevel=FALSE)

        table.bnf.ch.aug <- combine.tables3(table.bnf.ch, univariate.bnf.ch,
                                            multivariate.bnf.ch)
        return(table.bnf.ch.aug) 
    } else {
        return(NULL)
    }
}

tabulate.bnfchemicals <- function(chnum, outcome="CASE", data=cc.severe, minrowsum=50) {
    ## get drug names in wide format,one col per approved_name
    chnum.str <- sprintf("%02d", chnum)
    scrips.ch <- readRDS(scrips.filename) %>%
        subset(substr(bnf_paragraph_code, 1, 2)==chnum.str)
    scrips.ch.wide <- reshape2::dcast(scrips.ch,
                                         ANON_ID ~ approved_name, fun.aggregate=length,
                                         value.var="approved_name")
    #colnames(scrips.ch.wide)[-1] <- paste0("chemical.",
    #                                          as.integer(colnames(scrips.ch.wide)[-1]))
    ## drop rare chemicalgraphcodes
    scrips.ch.wide <- scrips.ch.wide[, colSums(scrips.ch.wide) > 20]
    
    cc.bnf.ch <- merge(data[, c("ANON_ID", "CASE", "stratum")],
                       scrips.ch.wide,
                       by="ANON_ID", all.x=TRUE)
    ## now fix colnames and set missing to 0
    bnfcols <- 4:ncol(cc.bnf.ch)
    for(j in bnfcols) {
        cc.bnf.ch[, j][is.na(cc.bnf.ch[, j])] <- 0
        cc.bnf.ch[, j][cc.bnf.ch[, j] > 1] <- 1
        cc.bnf.ch[, j] <- as.factor(cc.bnf.ch[, j])
    }

    ## approved names contain spaces that have to be replaced with . for use in formula
    colnames(cc.bnf.ch) <- gsub(" ", "\\.", colnames(cc.bnf.ch))
    
    bnf.ch <- colnames(cc.bnf.ch)[-(1:3)]
    table.bnf.ch <- univariate.tabulate(varnames=bnf.ch, outcome=outcome, data=cc.bnf.ch,
                                        drop.sparserows=TRUE, minrowsum=minrowsum)
    ## regressions fixed to use only rows retained by univariate.tabulate
    if(nrow(table.bnf.ch) >= 1) { 
        bnf.ch <- rownames(table.bnf.ch)
        univariate.bnf.ch <- NULL
        for(i in 1:length(bnf.ch)) {
            univariate.formula <- as.formula(paste(outcome, "~", bnf.ch[i], " + strata(stratum)"))
            x <- summary(clogit(formula=univariate.formula, data=cc.bnf.ch))$coefficients
            univariate.bnf.ch <- rbind(univariate.bnf.ch, x)
        }

        multivariate.bnf.ch <- multivariate.clogit(varnames=bnf.ch,
                                                   data=cc.bnf.ch, add.reflevel=FALSE)

        table.bnf.ch.aug <- combine.tables3(table.bnf.ch, univariate.bnf.ch,
                                            multivariate.bnf.ch)
        return(table.bnf.ch.aug) 
    } else {
        return(NULL)
    }
}

tabulate.icdchapter <- function(chnum, data=cc.severe, minrowsum=20) {
    ## get chapternum in wide format, one col per subchapter
    diagnoses.ch.wide <- reshape2::dcast(diagnoses[diagnoses$chapter==chnum, ],
                                         ANON_ID ~ subchapter, fun.aggregate=length,
                                         value.var="subchapter")
    colnames(diagnoses.ch.wide)[-1] <- paste0("subCh.",
                                              as.integer(colnames(diagnoses.ch.wide)[-1]))
    ## drop rare subchapters
    diagnoses.ch.wide <- diagnoses.ch.wide[, colSums(diagnoses.ch.wide) > 10, drop=FALSE]
    
    if(ncol(diagnoses.ch.wide) > 1) {
        cc.icd.ch <- merge(data[, c("ANON_ID", "CASE", "stratum")],
                           diagnoses.ch.wide,
                           by="ANON_ID", all.x=TRUE)
        ## now fix colnames and set missing to 0
        icdcols <- grep("^subCh.", colnames(cc.icd.ch))
        for(j in icdcols) {
            cc.icd.ch[, j][is.na(cc.icd.ch[, j])] <- 0
            cc.icd.ch[, j][cc.icd.ch[, j] > 1] <- 1
            cc.icd.ch[, j] <- as.factor(cc.icd.ch[, j])
        }
        
        icd.ch <- colnames(cc.icd.ch)[-(1:3)]
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
                print(table.icd.ch.aug) # without this line the first two cols are dropped?
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

combine.tables3 <- function(ftable, utable, mtable)  {# returns single table from freqs, univariate, multivariate 
    
    u.ci <- or.ci(utable[, 1], utable[, 3]) 
    u.pvalue <- signif(utable[, 5], 1)
    u.pvalue <- as.character(pvalue.latex(u.pvalue))
    
    mult.ci <- or.ci(mtable[, 1], mtable[, 3])
    mult.pvalue <- signif(mtable[, 5], 1)
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
    x <- as.numeric(x)
    pvalue <- sapply(x, function(z) { # sapply returns a vector applying FUN to each element of z
        if (is.na(z)) {
            return(NA)
        } else if(z >= 0.001) {
            return(signif(z, 1))
        } else {
            z <- sprintf("%.*E", 0, signif(z, n)) # default is 1 sig fig
            z <- as.numeric(unlist(strsplit(z, "E"))) # split z at E
            sprintf("\\ensuremath{%.*f\\times 10^{%0*d}}", 0, z[1], nexp, z[2])
        }
    })
    return(as.character(pvalue))
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

univariate.tabulate <- function(varnames, outcome="CASE", data, drop.reflevel=TRUE,
                                drop.sparserows=FALSE, minrowsum=10) {
    outcome <- data[, match(outcome, names(data))]
    table.varnames <- NULL
    for(i in 1:length(varnames)) {
        keep.x <- TRUE
        ## test whether variable is factor or numeric
        z <- data[, match(varnames[i], names(data))] 
        if(is.numeric(z)) { # median (IQR) for numeric variables 
            x <- tapply(z, outcome,
                        function(x) return(paste0(median(x, na.rm=TRUE),
                                                  " (", quantile(x, probs=0.25, na.rm=TRUE),
                                                  "-", quantile(x, probs=0.75, na.rm=TRUE), ")")))
            x <- matrix(x, nrow=1)
            rownames(x) <- varnames[i]
        } else { # freqs for factor variables
            x <- table(z, outcome)
            ## keep if at least one factor level has row sum >= minrowsum
            keep.x <- !any(rowSums(x) < minrowsum)
            x <- paste.colpercent(x)
            ## rownames are labelled with levels(varname)
            ## if two levels OR drop.reference level, drop reference level
            ## if single row left, label rows with varname
            ## else keep reference level
            if(nrow(x) == 2  | drop.reflevel) {
                x <- x[-1, , drop=FALSE] # drop reference category
                if(nrow(x) == 1)
                    rownames(x) <- varnames[i]
            } 
            if(!drop.sparserows | keep.x) { # rbind this table variable
                table.varnames <- rbind(table.varnames, x)
            } 
        } 
    }
    
    ## this line breaks code
    if(!is.null(table.varnames)) {
        colnames(table.varnames)[1:2] <- c("Controls", "Cases")
        colnames(table.varnames) <- paste0(colnames(table.varnames),
                                           " (", as.integer(table(outcome)), ")")
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

traintest.fold <- function(i, cv.data,
                           lower.formula, upper.formula) { # stepwise clogit on training fold, predict on test fold
    test.data  <- cv.data[cv.data$test.fold == i, ]
    train.data <- cv.data[cv.data$test.fold != i, ]
    start.model <- clogit(formula=upper.formula, data=train.data)
    cat("train.data dimensions", dim(train.data), "\n")
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
    return(norm.predicted)
}

nonmissing.obs <- function(x, varnames) { ## subset rows nonmissing for varnames
    keep <- rep(TRUE, nrow(x))
    for(j in 1:length(varnames)) {
        keep[is.na(x[, match(varnames[j], colnames(x))])] <- FALSE
    }
    return(keep)
}

normalize.predictions <- function(unnorm.p, stratum, y) { # format a pvalue in latex
    ## normalize probs so that they sum to 1 within each stratum
    ## returns data frame with 4 columns: stratum, normconst, prior.p, posterior.p

    unnorm.p <- as.numeric(unnorm.p)
    stratum <- as.integer(as.character(stratum))
    nonmissing <- !is.na(unnorm.p)

    ## keep only strata that contain a single case 
    table.strata <- tapply(y, stratum, sum) == 1
    strata.onecase <- as.integer(names(table.strata)[as.integer(table.strata)==1])
    keep <- !is.na(unnorm.p) & stratum %in% strata.onecase

    cat(length(which(!keep)), "predictions dropped because stratum does not contain one case\n")

    unnorm.p <- unnorm.p[keep]
    stratum <- as.factor(stratum[keep])
    y <- y[keep]

    ## compute normalizing constant and prior for each stratum, then merge
    norm.constant <- as.numeric(tapply(unnorm.p, stratum, sum))
    stratum.size <- as.numeric(tapply(unnorm.p, stratum, length))
    print(table(stratum.size))
    prior.p <- 1/stratum.size
    norm.constant <- data.frame(stratum=levels(stratum),
                                normconst=norm.constant,
                                prior.p=as.numeric(prior.p))
    
    norm.predicted <- data.frame(unnorm.p, y, stratum)
    norm.predicted <- merge(norm.predicted, norm.constant, by="stratum")
    # normalize probs
    norm.predicted$posterior.p <- norm.predicted$unnorm.p / norm.predicted$normconst
    return(norm.predicted)
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
