
## tables and figures for shielding report


## rollmean.dtN returns a data.table with columns date, avdaily (rolling mean) 
rollmean.dtN <- function(dt, k) {
    ## FIXME: add a by= argument
    N.dates <- dt[, .N, by=SPECIMENDATE]
    setkey(N.dates, SPECIMENDATE)
    all.dates <- data.table(date=seq(from=min(N.dates$SPECIMENDATE),
                                     to=max(N.dates$SPECIMENDATE),
                                     by=1))
    setkey(all.dates, date)

    ## add rows for missing dates by a left join
    all.dates <- N.dates[all.dates]
    all.dates[is.na(N), N := 0]
    return(data.table(date=dt$SPECIMENDATE,
                      avdaily=zoo::rollmean(dt$N, k=k, fill=NA)))
}

rollsum.datewin <- function(dt, k, datevar) {
    ## returns a table of rolling sums of width k centred on date
    ## dt is a dataset of totals (N) by date, which may have no rows for some dates in the
    ## range of interest
    ## add rows for missing dates by a left join of all.dates with dt
    all.dates <- data.table(date=seq(from=min(dt[[datevar]]),
                                     to=max(dt[[datevar]]),
                                     by=1))
    setkey(all.dates, date)
    # setkeyv(dt, datevar) can't set physical key here
    all.dates <- dt[all.dates]

    return(data.table(date=all.dates[[datevar]],
                      rsum=zoo::rollsum(all.dates[, N], k=k, fill=0, align="center")))
}

freqs.slidingwindow <- function(dt, winsize=7, datevar=NULL, categoricvar=NULL,
                                groupvar=NULL) {
    ## returns a table of freqs for categoric var by groupvar over sliding windows of
    ## width winsize centred at date
    freqs.bydate <- dt[, .N, by=c(datevar, categoricvar, groupvar)] 
    setkeyv(freqs.bydate, datevar)
    freqs.tw <- freqs.bydate[, rollsum.datewin(dt=.SD, k=winsize, datevar=datevar),
                         by=c(categoricvar, groupvar),
                         .SDcols=c(datevar, "N")]
    ## compute sum over levels of categoricvar in each window, by groupvar
    freqs.N <- freqs.tw[, .(rsum.allcategories = sum(rsum)), by=c("date", groupvar)]
    ## left join freqs.tw with freqs.N
    setkeyv(freqs.N, c("date", groupvar))
    setkeyv(freqs.tw, c("date", groupvar))
    freqs.tw <- freqs.N[freqs.tw]
    freqs.tw[, p := rsum / rsum.allcategories]
    freqs.tw[, se.p := sqrt(p * (1 - p) / rsum.allcategories)]
    return(freqs.tw) # one row for each level of categoricvar
}

datadir <- "data/2021-02-18/"
load(paste0(datadir, "cc.all.RData"))

theme_set(theme_gray(base_size = 14))

lastdate <- as.Date("2021-02-01")

## Table 1 - shielding cohort by eligibility category x age

shielded.full[, shield.group6 := car::recode(shield.group,
                                             recode="'Pregnant with heart disease'=NA",
                                             levels=levels(shield.group)[-6])]
                                             
bygroup <- with(shielded.full, table(shield.group6))
bygroup <- as.integer(c(bygroup, sum(bygroup)))

shieldgroups.byage <- with(shielded.full,
                           as.matrix(as.data.frame.matrix(table(shield.group6, agegr20))))
shieldgroups.byage <- rbind(shieldgroups.byage, colSums(shieldgroups.byage))
rownames(shieldgroups.byage)[nrow(shieldgroups.byage)] <- "All shielding categories"

shieldgroups.byage.rowpct <- round(100 * shieldgroups.byage /
                                   rowSums(shieldgroups.byage))
shieldgroups.byage.rowpct <- paste0("(", shieldgroups.byage.rowpct, "\\%)")
table.shieldgroups.byage <- paste.vectortomatrix(shieldgroups.byage,
                                                 shieldgroups.byage.rowpct)
table.shieldgroups.byage <- data.frame(table.shieldgroups.byage,
                                       All=bygroup)
colnames(table.shieldgroups.byage) <- c(levels(shielded.full$agegr20), "All")


## Table 2 - shielding cohort by shielding category x outcome
table.shield.casegr <- univariate.tabulate(outcome="casegr3",
                                           varnames=c("shield.group6"),
                                           data=shielded.full, drop.reflevel=FALSE,
                                           digits=1, colpercent=FALSE)

table.shield.casegr.rowsums <- rowSums(with(shielded.full, table(shield.group6, casegr3)))

table.shield.casegr.rowsums <- c(sum(table.shield.casegr.rowsums),
                                     table.shield.casegr.rowsums)

## add row for all shielding categories with row percentages
table.casegr <- with(shielded.full, table(casegr3))
table.casegr <- paste0(table.casegr, " (",
                       round(100 * table.casegr / sum(table.casegr), 1),
                       "\\%)")
table.shield.casegr <- rbind(table.casegr, table.shield.casegr)
rownames(table.shield.casegr) <- replace.names(rownames(table.shield.casegr))
rownames(table.shield.casegr)[1] <- "All shielding categories"
colnames(table.shield.casegr) <- gsub("\\(.+\\)", "", colnames(table.shield.casegr))
table.shield.casegr <- data.frame(table.shield.casegr, All=table.shield.casegr.rowsums)
colnames(table.shield.casegr) <- gsub("\\.+", " ", colnames(table.shield.casegr))
colnames(table.shield.casegr) <- gsub(" non ", ", non-", colnames(table.shield.casegr))


###################################

## Table 3 - shielding cohort by date of letter
table.date.letter <- as.matrix(as.data.frame.matrix(
    with(shielded.full, table(shield.group6, Date.Sent))))
col5name <- colnames(table.date.letter)[5]
table.date.letter <- rbind(table.date.letter, colSums(table.date.letter))
rownames(table.date.letter)[nrow(table.date.letter)] <- "All shielding categories"
                           
totals <-  rowSums(table.date.letter)
lastcol <- ncol(table.date.letter)

table.date.letter <- cbind(table.date.letter[, 1:4], rowSums(table.date.letter[, 5:lastcol]))
date.letter.rowpct <- round(100 * table.date.letter / totals)
date.letter.rowpct <- paste0("(", date.letter.rowpct, "\\%)")

table.date.letter <- paste.vectortomatrix(table.date.letter, date.letter.rowpct)
table.date.letter <- cbind(table.date.letter, totals)
colnames(table.date.letter)[5:6] <- c(col5name, "All")

table.date.letter <- as.data.frame(table.date.letter)
colnames(table.date.letter)[1:5] <- format(as.Date(colnames(table.date.letter)[1:5]),
                                           "%d %b")
colnames(table.date.letter) <- gsub("^0", "", colnames(table.date.letter))
colnames(table.date.letter)[5] <- paste(colnames(table.date.letter)[5], "or later")

#####################################
## Table 4: rate ratios asssociated with eligibility category by advice interval

table.shielding.byinterval <-
    cbind(
        tabulate.freqs.regressions(varnames=c("shield.group", "shield.any"),
                                   outcome="CASE",
                                   data=cc.severe[care.home=="Independent" &
                                               SPECIMENDATE <= as.Date("2020-04-17")])[, 1:4], 
        tabulate.freqs.regressions(varnames=c("shield.group", "shield.any"), 
                                   outcome="CASE",
                                   data=cc.severe[care.home=="Independent" &
                                               SPECIMENDATE > as.Date("2020-04-17")])[, 1:4]) 
rownames(table.shielding.byinterval)[nrow(table.shielding.byinterval)] <- "All shielding categories"

with(cc.all[shield.group != "No shielding"], table(hh.over18gr, agegr20))

## calculate estimated proportion of all cases infected after arrival of each batch of letters
severe.infected.date <- table(cc.severe[care.home=="Independent" & vaxstatus=="No vaccine by presentation date"][CASE==1, SPECIMENDATE - 7])
dateby.props <- cumsum(severe.infected.date) / sum(severe.infected.date)
dateby.props <- dateby.props[match(levels(as.factor(cc.all$Date.Sent)), names(dateby.props))]


#####################################

## calculate estimated proportion of all cases infected after arrival of each batch of letters
severe.infected.date <- table(cc.severe[care.home=="Independent" & vaxstatus=="No vaccine by presentation date"][CASE==1, SPECIMENDATE - 7])
dateby.props <- cumsum(severe.infected.date) / sum(severe.infected.date)
dateby.props <- dateby.props[match(levels(as.factor(cc.all$Date.Sent)), names(dateby.props))]

#########################################
## exclude care home residents

## Table 5 - regression of severe case status on shielding group and covariates
table.shielded.severecases.nocare <-
    tabulate.freqs.regressions(varnames=c("shield.group", 
                                           "preschool.any", "hh.schoolagegr",
                                          "hh.over18gr", 
                                          "SIMD.quintile", "occup", "hosp.recent"),
                               outcome="CASE",
                               data=cc.severe[care.home=="Independent" & vaxstatus=="No vaccine by presentation date"])

## regression of severe case status on listedgr3 and covariates to get coeff for adultsgt1
table.listedgr3.severecases.nocare <-
    tabulate.freqs.regressions(varnames=c("listedgr3",
                                           "preschool.any", "hh.schoolagegr",
                                          "adultsgt1", 
                                          "SIMD.quintile", "occup", "hosp.recent"),
                               outcome="CASE",
                               data=cc.severe[care.home=="Independent" & vaxstatus=="No vaccine by presentation date"])

table.listedgr3.wave1 <-
    tabulate.freqs.regressions(varnames=c("listedgr3",
                                           "preschool.any", "hh.schoolagegr",
                                          "adultsgt1", 
                                          "SIMD.quintile", "occup", "hosp.recent"),
                               outcome="CASE",
                               data=cc.severe[care.home=="Independent" & vaxstatus=="No vaccine by presentation date" & SPECIMENDATE < as.Date("2020-09-01")])

table.listedgr3.wave2 <-
    tabulate.freqs.regressions(varnames=c("listedgr3",
                                           "preschool.any", "hh.schoolagegr",
                                          "adultsgt1", 
                                          "SIMD.quintile", "occup", "hosp.recent"),
                               outcome="CASE",
                               data=cc.severe[care.home=="Independent" & vaxstatus=="No vaccine by presentation date" & SPECIMENDATE >= as.Date("2020-09-01")])

                               
## regression of severe case status on listedgr3 and covariates to get coeff for hh.schoolage.any
table.listedgr3.severecases.nocare.anyschoolage <-
    tabulate.freqs.regressions(varnames=c("listedgr3",
                                           "preschool.any", "hh.schoolage.any",
                                          "hh.over18gr", 
                                          "SIMD.quintile", "occup", "hosp.recent"),
                               outcome="CASE",
                               data=cc.severe[care.home=="Independent" & vaxstatus=="No vaccine by presentation date"])

### PARF for severe cases associated with recent hospital exposure and with probable HCAI

r.rapid <- summary(clogit(data=cc.severe[care.home=="Independent" & vaxstatus=="No vaccine by presentation date"],
               CASE ~ rapid.recent + strata(stratum)))$coefficients[2]
p.rapid <- with(cc.severe[care.home=="Independent" & CASE==1], mean(as.integer(rapid.recent)))

parf.rapid <- p.rapid * (r.rapid - 1) / r.rapid

r.hcai <- summary(clogit(data=cc.severe[care.home=="Independent" & vaxstatus=="No vaccine by presentation date"],
               CASE ~ prob.hcai + strata(stratum)))$coefficients[2]
p.hcai <- with(cc.severe[care.home=="Independent" & CASE==1], mean(as.integer(prob.hcai)))

parf.hcai <- p.hcai * (r.hcai - 1) / r.hcai

## Table 6 - regression in shielded eligible only
summary(cc.severe[listedgr3 == "Eligible for shielding" & care.home=="Independent",
                  .(shieldedonly.group, hh.over18gr, hosp.recent)])
table.shieldedonly.severecases.nocare <-
    tabulate.freqs.regressions(varnames=c("shieldedonly.group", 
                                           "preschool.any", "hh.schoolagegr",
                                          "hh.over18gr3", 
                                          "qSIMD.integer", "hosp.recent"),
                               outcome="CASE",
                               data=cc.severe[listedgr3 == "Eligible for shielding" &
                                              care.home=="Independent"])


###########################################
## Figure 1: severe cases by date of presentation
winsize.casedates <- 3

casedates.allgr <-  cc.severe[CASE==1 & care.home=="Independent", .N, by=SPECIMENDATE]
setkey(casedates.allgr, SPECIMENDATE)
casedates.allgr <- casedates.allgr[, rollmean.dtN(dt=.SD, k=winsize.casedates),
                                   .SDcols=c("SPECIMENDATE", "N")]

## use by= to get rolling means by category
casedates.byelig <- cc.severe[CASE==1 & care.home=="Independent", .N, by=.(SPECIMENDATE, listedgr3)]
setkey(casedates.byelig, SPECIMENDATE)
casedates.byelig <- casedates.byelig[, rollmean.dtN(dt=.SD, k=winsize.casedates), by=listedgr3,
                       .SDcols=c("SPECIMENDATE", "N")]
                                 
p.byelig <- ggplot(data=casedates.byelig,
                aes(x=date,
                    y=avdaily, fill=listedgr3)) +
    geom_area() +
    scale_y_continuous(expand=c(0, 0)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = as.Date("2020-03-01"),
                                                  to = lastdate, by = "month"),
                                         "%d %b")
                 ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
    theme(legend.position = c(0.5, 0.7)) +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values=c("grey", "blue", "red")) +
    xlab(paste0("Presentation date: mid-point of ", winsize.casedates, "-day window")) +
         ylab("Daily severe cases")
p.byelig

#######################################################################
## Figure 2: plot rate ratios associated with shielding eligibility, recent hospital exposure and number of children

firstdate <- as.Date("2020-03-01")

## for coeffs we have to specify explicitly the midpoint of the time window 
winsize <- 21 # should be an odd number for date.midpoint to be exactly centred
startdates <- with(cc.severe, firstdate:max(SPECIMENDATE) - winsize)
enddates <- startdates + winsize

## 2 (a) plot rate ratio associated with shielding eligibility by time window

coeffs.timewindow <- NULL
for(timewin in 1:length(startdates)) {
    tdata <- cc.severe[care.home=="Independent" &
                       as.integer(SPECIMENDATE) >= startdates[timewin] &
                       as.integer(SPECIMENDATE) <= enddates[timewin],
                       .(CASE, listedgr3, stratum)]
    if(length(with(tdata, table(CASE, listedgr3))) > 2) {
        date.midpoint <- startdates[timewin] + floor(winsize / 2)
        coeffs <- c(date.midpoint,
                    tryCatch(summary(clogit(formula=CASE ~ listedgr3 + strata(stratum), data=tdata))$coefficients,
                             error=function(cond) return(matrix(rep(NA, 10), nrow=2))
                             )
                    )
        coeffs.timewindow <- rbind(coeffs.timewindow, coeffs)
    }
}

coeffs.timewindow <- data.table(coeffs.timewindow)
coeffs.long <- melt(coeffs.timewindow, idvars=1, measure=list(2:3, 4:5, 6:7, 8:9, 10:11),
                    value.name=c("coeff", "rateratio", "se.coeff", "Z", "pvalue"))
colnames(coeffs.long)[1:2] <- c("date.midpoint", "riskgroup")
coeffs.long[, date.midpoint := as.Date(date.midpoint, origin="1970-01-01")]
coeffs.long[, riskgroup := car::recode(riskgroup,
                                       "1='Moderate risk condition'; 2='Eligible for shielding';",
                                       as.factor=TRUE)]
 
coeffs.long[date.midpoint >= as.Date("2020-06-01") &
                  date.midpoint < as.Date("2020-09-01"), coeff := NA]

num.letters <- with(shielded.full, table(Date.Sent))
num.letters <- num.letters / num.letters[1]
dates.letters <- as.Date(names(num.letters))
coeffs.long[date.midpoint >= as.Date("2020-06-01") & date.midpoint <= as.Date("2020-09-30"), coeff := NA]
                                  
p.rateratio <-
    ggplot(data=coeffs.long, aes(x=date.midpoint, y=coeff, color=riskgroup)) +
    geom_line(size=0.01 * coeffs.long[, se.coeff]^-2) +
    scale_color_manual(name="Risk category",
                       values=c("red", "blue")) + 
    theme(legend.position = c(0.47, 0.82)) +
    theme(legend.title = element_text(size=10), legend.text=element_text(size=10)) + 
    scale_y_continuous(breaks=log(c(2, 3, 4, 5, 6, 7, 8, 10, 12)),
                       labels=c(2, 3, 4, 5, 6, 7, 8, 10, 12),
                       limits=log(c(2, 12)), expand=c(0, 0)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                     format.Date(seq.Date(from = as.Date("2020-03-01"),
                              to = lastdate, by = "month"),
                     "%d %b")
                 ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
    labs(x=paste0("Presentation date: mid-point of ", winsize, "-day window"),
         y="Rate ratio (log scale)",
         title="Association of severe COVID-19 with risk category",
         tag="(a)",
         caption="Arrows indicate dates that shielding advice letters were sent, with thickness proportional to number of letters") + 
    annotate(geom="segment",
             x = dates.letters,
             xend = dates.letters,
             y=log(2),
             yend=log(2.3),
             size=as.numeric(num.letters), 
             arrow=arrow(ends="first", length=unit(0.1, "inches")))
p.rateratio

#####################################

## 2 (b) plot rate ratio associated with recent hospital exposure by time window
## coeffs.hosp calculated from unadjusted coefficient
coeffs.hosp.timewindow <- NULL
for(timewin in 1:length(startdates)) {
    tdata <- cc.severe[care.home=="Independent" &
                       as.integer(SPECIMENDATE) >= startdates[timewin] &
                       as.integer(SPECIMENDATE) <= enddates[timewin], .(CASE, hosp.recent,
                                                                        stratum)]
    if(
        length(with(tdata, table(CASE))) > 1 & 
        length(with(tdata, which(hosp.recent==TRUE))) > 0
    ) {
        date.midpoint <- startdates[timewin] + floor(winsize / 2)
        coeffs.hosp <- c(date.midpoint,
                         summary(clogit(formula=CASE ~ hosp.recent + strata(stratum),
                                        data=tdata))$coefficients)
        coeffs.hosp.timewindow <- rbind(coeffs.hosp.timewindow, coeffs.hosp)
    }
}

colnames(coeffs.hosp.timewindow) <- c("date.midpoint", "coeff", "rateratio", "se.coeff", "Z", "pvalue")

coeffs.hosp.timewindow <- data.table(coeffs.hosp.timewindow)
coeffs.hosp.timewindow[, date.midpoint := as.Date(date.midpoint, origin="1970-01-01")]
coeffs.hosp.timewindow <- coeffs.hosp.timewindow[se.coeff < 0.5]
coeffs.hosp.timewindow[date.midpoint >= as.Date("2020-06-01") &
                  date.midpoint <= as.Date("2020-09-30"), coeff := NA]

p.hosp.rateratio <-
    ggplot(data=coeffs.hosp.timewindow, aes(x=date.midpoint, y=coeff)) +
    geom_line(size=0.01 * coeffs.hosp.timewindow[, se.coeff]^-2) +
    scale_y_continuous(breaks=log(c(5, 10, 20, 30, 50)),
                           labels=c(5, 10, 20, 30, 50)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                     format.Date(seq.Date(from = as.Date("2020-03-01"),
                              to = lastdate, by = "month"),
                     "%d %b")
                 ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
    labs(x=paste0("Presentation date: mid-point of ", winsize, "-day window"),
         y="Rate ratio (log scale)",
         title="Association of severe COVID-19 with recent hospital exposure",
         tag="(b)")

#################################################################################
## 2 (c) plot rate ratio associated with number of children by time window

coeffs.household.timewindow <- NULL
for(timewin in 1:length(startdates)) {
    tdata <- cc.severe[care.home=="Independent" &
                       as.integer(SPECIMENDATE) >= startdates[timewin] &
                       as.integer(SPECIMENDATE) <= enddates[timewin], .(CASE, hh.schoolage, hh.over18,
                                                                        qSIMD.integer,
                                                                        stratum)]
    if(
        length(with(tdata, table(CASE))) > 1 & 
        length(with(tdata, which(hh.schoolage > 0))) > 0
    ) {
        date.midpoint <- startdates[timewin] + floor(winsize / 2)
        coeffs.household <-
            data.frame(date.midpoint=rep(date.midpoint, 2),
                       hhvar=c("per child", "per adult"),
                       summary(clogit(formula=CASE ~ hh.schoolage + hh.over18 + qSIMD.integer 
                                      + strata(stratum),
                                      data=tdata))$coefficients[1:2, ])
        coeffs.household.timewindow <- rbind(coeffs.household.timewindow, coeffs.household)
    }
}

coeffs.household.timewindow <- data.table(coeffs.household.timewindow)

colnames(coeffs.household.timewindow) <- 
    c("date.midpoint", "RR", "coeff", "rateratio", "se.coeff", "Z", "pvalue")
coeffs.household.timewindow[, date.midpoint := as.Date(date.midpoint, origin="1970-01-01")]
coeffs.household.timewindow <- coeffs.household.timewindow[se.coeff < 0.5]
coeffs.household.timewindow[date.midpoint >= as.Date("2020-06-01") &
                  date.midpoint < as.Date("2020-09-01"), coeff := NA]

coeffs.adults.timewindow <- coeffs.household.timewindow[RR=="adults"]
coeffs.children.timewindow <- coeffs.household.timewindow[RR=="children"]
coeffs.household.timewindow[date.midpoint >= as.Date("2020-06-01") & date.midpoint <= as.Date("2020-09-30"), coeff := NA]

p.household.rateratio <-
    ggplot(data=coeffs.household.timewindow, aes(x=date.midpoint, y=coeff, color=RR)) +
    geom_line(size=0.05 * coeffs.household.timewindow[, se.coeff]^-1) +
    xlim(as.Date(c("2020-03-01", "2020-11-30"))) +   
    scale_y_continuous(breaks=log(c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8)),
                       labels=c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8), 
                       limits=log(c(0.6, 1.85)), expand=c(0, 0)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                     format.Date(seq.Date(from = as.Date("2020-03-01"),
                              to = lastdate, by = "month"),
                     "%d %b")
                 ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
    theme(legend.position = c(0.5, 0.5)) + 
    geom_hline(size=0.2, yintercept=0) +
    labs(x=paste0("Presentation date: mid-point of ", winsize, "-day window"),
         y="Rate ratio (log scale)", 
         color='Rate ratio',
         title="Association with adults and school-age children in household",
         caption="Data from 1 June to 30 September 2020 are omitted because the numbers are small", 
         tag="(b)") 

#################################################################################

## Figure 3: frequency of recent hospital exposure in controls(excluding care home residents, population attributable risk fraction excluding care home residents, and fraction of cases in care home residents

## 3 (a) frequency of recent hospital exposure
winsize.hosp <- 7
freqs.hosp.ctrls.tw <- freqs.slidingwindow(dt=cc.severe[CASE==0 & care.home=="Independent"],
                                winsize=winsize.hosp,
                                datevar="SPECIMENDATE",
                                categoricvar="hosp.recent",
                                groupvar="listedgr3")
freqs.hosp.ctrls.tw <- freqs.hosp.ctrls.tw[hosp.recent==TRUE]
freqs.hosp.ctrls.tw[rsum.allcategories <5, `:=`(p=NA, se.p=NA)]
freqs.hosp.ctrls.tw[date >= as.Date("2020-06-01") & date <= as.Date("2020-09-30"), p := NA]

p.hosp <-
    ggplot(data=freqs.hosp.ctrls.tw,
           aes(x=date, y=p, color=listedgr3, group=listedgr3)) +
    geom_line() + # size=0.002 * (freqs.hosp.ctrls.tw$rsum.allcategories)^0.5) + 
    scale_y_continuous(limits=c(0, max(freqs.hosp.ctrls.tw$p, na.rm=TRUE)), expand=c(0, 0)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = as.Date("2020-03-01"),
                                                  to = lastdate, by = "month"),
                                         "%d %b")
                             ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
    theme(legend.position = c(0.4, 0.5)) +
    theme(legend.title = element_blank()) +
    scale_color_manual(values=c("black", "blue", "red")) +
    #guides(fill = guide_legend(reverse = TRUE)) + 
    labs(x=paste0("Presentation date: mid-point of ", winsize, "-day window"),
         y="Frequency", 
         title="Recent hospital exposure in controls matched to severe cases",
         tag="(a)")

###########################################################

## 3 (b) fraction attributable to recent hospital exposure excluding care home residents

## compute daily counts of cases and controls with recent hospital exposure
cc.hosp.recent <- cc.severe[CASE==0 & care.home=="Independent" & hosp.recent==TRUE,
                        .N, by=.(SPECIMENDATE, CASE)]
setkeyv(cc.hosp.recent, c("SPECIMENDATE", "CASE"))

## compute daily counts of all cases and controls
cc.counts.all <- cc.severe[care.home=="Independent", .N, by=.(SPECIMENDATE, CASE)]
setkeyv(cc.counts.all, c("SPECIMENDATE", "CASE"))

#########################

## compute rolling mean of avdaily in controls with recent hospital exposure
controls.hosp.rollmean <- cc.hosp.recent[CASE==0, rollmean.dtN(dt=.SD, k=winsize),
                                  .SDcols=c("SPECIMENDATE", "N")]
setnames(controls.hosp.rollmean, "avdaily", "avdaily.exposed")
setkey(controls.hosp.rollmean, date)
## compute rolling mean of avdaily in all controls
controls.all.rollmean <- cc.counts.all[CASE==0, rollmean.dtN(dt=.SD, k=winsize),
                               .SDcols=c("SPECIMENDATE", "N")]
setnames(controls.all.rollmean, "avdaily", "avdaily.total")
setkey(controls.all.rollmean, date)
## left join av counts in all controls with av counts in hosp.recent
controls.all.rollmean <- controls.hosp.rollmean[controls.all.rollmean]
## compute frequency of exposure in controls
controls.all.rollmean[, expfreq := avdaily.exposed / avdaily.total]
setkey(controls.all.rollmean, date)

######################################

## compute rolling mean of avdaily in cases with recent hospital exposure
cases.hosp.rollmean <- cc.hosp.recent[CASE==0, rollmean.dtN(dt=.SD, k=winsize),
                                  .SDcols=c("SPECIMENDATE", "N")]
setnames(cases.hosp.rollmean, "avdaily", "avdaily.exposed")
setkey(cases.hosp.rollmean, date)
## compute rolling mean of avdaily in all cases
cases.all.rollmean <- cc.counts.all[CASE==1, rollmean.dtN(dt=.SD, k=winsize),
                               .SDcols=c("SPECIMENDATE", "N")]
setnames(cases.all.rollmean, "avdaily", "avdaily.total")
setkey(cases.all.rollmean, date)
## left join av counts in all cases with av counts in hosp.recent
cases.all.rollmean <- cases.hosp.rollmean[cases.all.rollmean]
## compute frequency of exposure in cases
cases.all.rollmean[, case.expfreq := avdaily.exposed / avdaily.total]
setnames(cases.all.rollmean, "avdaily.exposed", "avdaily.exposed.case")
setnames(cases.all.rollmean, "avdaily.total", "avdaily.total.case")
setkey(cases.all.rollmean, date)

######################################

## merge with coeffs
setkey(coeffs.hosp.timewindow, date.midpoint)
coeffs.hosp.timewindow <- controls.all.rollmean[coeffs.hosp.timewindow]
coeffs.hosp.timewindow <- cases.all.rollmean[coeffs.hosp.timewindow]

## Miettinen's definition
coeffs.hosp.timewindow[, popattr.frac := case.expfreq * (rateratio - 1) / rateratio]

coeffs.hosp.timewindow[date.midpoint >= as.Date("2020-06-01") &
                      date.midpoint <= as.Date("2020-09-30"), popattr.frac := NA]

coeffs.hosp.timewindow[date >= as.Date("2020-06-01") & date <= as.Date("2020-09-30"), popattr.frac := NA]

p.fracattr <-
    ggplot(data=coeffs.hosp.timewindow, aes(x=date, y=popattr.frac)) +
    geom_line(size=0.01 * coeffs.hosp.timewindow[, se.coeff]^-2) +
    scale_y_continuous(limits=c(0.2, 0.7), expand=c(0, 0)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = as.Date("2020-03-01"),
                                                  to = lastdate, by = "month"),
                                         "%d %b")
                             ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
    labs(x=paste0("Presentation date: mid-point of ", winsize, "-day window"),
         y="Population attributable risk fraction",
         title="Fraction of severe cases attributable to hospital exposure",
         caption="Data from 1 June to 30 September 2020 are omitted because the numbers are small.",
         tag="(c)")

#############################################################################
## Supplementary Figure 1: case fatality rate

casedates.all <- cc.all[CASE==1, .N, by=.(SPECIMENDATE, fatalcase, shield.any)]
## cast to wide format
casedates.all <- dcast(casedates.all, SPECIMENDATE ~ fatalcase + shield.any, value.var="N")
## set NA to 0
for(j in 2:5) {
    set(casedates.all, which(is.na(casedates.all[[j]])), j, 0)
}
## compute sliding window
winsize.fatality <- 7
for(j in names(casedates.all)[2:5]) {
    casedates.all[, (j) := zoo::rollmean(casedates.all[[j]], k=winsize.fatality, fill=NA)]
}

colnames(casedates.all)[2:5] <- c("nonfatal.inelig", "nonfatal.elig",
                                  "fatal.inelig", "fatal.elig")

## compute N and  casefatality
casedates.all[, N.inelig := fatal.inelig + nonfatal.inelig]
casedates.all[, N.elig := fatal.elig + nonfatal.elig]

casedates.all[, casefatality.inelig := fatal.inelig / N.inelig]
casedates.all[, casefatality.elig := fatal.elig / N.elig]

## concatenate inelig and elig
casedates.inelig <- casedates.all[, .(SPECIMENDATE, N.inelig, casefatality.inelig)]
casedates.elig <- casedates.all[, .(SPECIMENDATE, N.elig, casefatality.elig)]
colnames(casedates.inelig) <- colnames(casedates.elig) <- c("SPECIMENDATE", "N", "casefatality")
casedates.inelig[, shield.elig := FALSE]
casedates.elig[, shield.elig := TRUE]
casedates.all <- rbind(casedates.inelig, casedates.elig)
setkey(casedates.all, SPECIMENDATE)

p.casefatality <- ggplot(data=casedates.all,
       aes(x=SPECIMENDATE, y=casefatality, group=shield.elig, color=shield.elig)) +
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = as.Date("2020-05-01"), by = "week"),
                 labels=gsub("^0", "", 
                     format.Date(seq.Date(from = as.Date("2020-03-01"),
                                          to = as.Date("2020-05-01"),
                                          by = "week"),
                                 "%d %b")
                     ),
                 limits=c(as.Date("2020-03-01"), as.Date("2020-05-01"))) +
    scale_y_continuous(expand=c(0, 0)) + 
    geom_line(size=ifelse(casedates.all[["shield.elig"]], 0.4, 0.02) *
                  casedates.all[["N"]]^0.5 
              ) +
    scale_color_manual(name="Shielding status",
                       labels=c("Ineligible","Eligible"),
                       values=c("blue", "red")) +
    theme(legend.position = c(0.6, 0.8)) + 
    xlab(paste0("Presentation date: midpoint of ", winsize.fatality, "-day window")) +
    ylab("Proportion of cases that were fatal") 

############################################################################
winsize.casedates <- 3
## use by= to get rolling means by category

casedates.bysource <- cc.severe[CASE==1, .N, by=.(SPECIMENDATE, exp.group)]
setkey(casedates.bysource, SPECIMENDATE)
casedates.bysource <- casedates.bysource[, rollmean.dtN(dt=.SD, k=winsize.casedates),
                                         by=exp.group,
                                         .SDcols=c("SPECIMENDATE", "N")]
                                 
p.bysource <- ggplot(data=casedates.bysource[date < as.Date("2021-01-20")],
                     aes(x=date,
                         y=avdaily, fill=exp.group)) +
    geom_area() +
    scale_x_date(breaks = seq.Date(from = firstdate,
                                   to = lastdate, by = "month"),
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = firstdate,
                                                  to = lastdate, by = "month"),
                                         "%d %b")
                 ),
                 limits=c(firstdate, lastdate)) +
    theme(legend.position = c(0.1, 0.8)) +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values=c("grey", "orange", "green")) +
    xlab(paste0("Specimen date: mid-point of ", winsize.casedates, "-day window")) +
         ylab("Daily severe cases")

###########################################################

## case-crossover analysis

varnames.twgr <- c("rapid.timewingr",
                   "disch.timewingr",
                   "daycase.timewingr",
                   "opd.timewingr",
                   "psych.timewingr")
varnames.keep <- c("CASE", "stratum", varnames.twgr)
cc.severe.testpos <- cc.severe[stratum %in% cc.severe[testpositive.case==TRUE, .N,
                                                      by=stratum][["stratum"]],
                               ..varnames.keep]  

freqs.tw <- univariate.tabulate(varnames=varnames.twgr,
                                outcome="CASE",
                                data=cc.severe.testpos,
                                drop.reflevel=FALSE)
coeffs.tw <- univariate.clogit(varnames=varnames.twgr,
                               outcome="CASE",
                               data=cc.severe.testpos,
                               add.reflevel=TRUE)

table.tw <- combine.tables2(freqs.tw, coeffs.tw)
table.tw <- data.frame(interval=rep(c("Less recent interval only", "Recent interval only",
                                   "Both intervals"), 5),
                       table.tw)

###########################################################

table.hcai.rapid.recent <- with(cc.severe[testpositive.case==TRUE],
                                                 table(hcai, rapid.recent))
colnames(table.hcai.rapid.recent) <- c("No recent hospital exposure",
                                "Recent hospital exposure")

table.hcai.rapid.tw <- with(cc.severe[testpositive.case==TRUE],
                                      table(hcai, rapid.timewingr, exclude=NULL))[, 1:2]
colnames(table.hcai.rapid.tw) <- c("Less recent time window only",
                                "Recent time window only")

table.hcai.rapid <- data.frame(as.data.frame.matrix(table.hcai.rapid.recent),
                               as.data.frame.matrix(table.hcai.rapid.tw))
table.hcai.rapid <- rbind(table.hcai.rapid, colSums(table.hcai.rapid))
rownames(table.hcai.rapid)[6] <- "All ECDC categories"

table.hcai.rapid.recent <- rbind(paste.colpercent(table.hcai.rapid.recent),
                                 colSums(table.hcai.rapid.recent))
rownames(table.hcai.rapid.recent)[6] <- "All ECDPC categories"

#######################################################
table.casegr <- paste.colpercent(with(cc.all, table(casegr, hosp.recent)))
table.hcai.hosp.recent <- paste.colpercent(with(cc.all, table(hcai, hosp.recent)))
##################################################

## restriction to those with underlying cause or critical care for ARHAI query
strata.testpos <- cc.severe[testpositive.case==TRUE &
                            (covid_ucod==1 | criticalcare==TRUE), .N, by=stratum]

testpos.coeffs <- tabulate.freqs.regressions(varnames=c("listedgr3", "hosp.recent", "occup"),
                           outcome="CASE",
                           data=cc.severe[care.home=="Independent" &
                                          stratum %in% strata.testpos$stratum==TRUE])

table.rapid.coeffs <- tabulate.freqs.regressions(varnames=c("listedgr3", "rapid.recent", "occup"),
                           outcome="CASE",
                           data=cc.severe[care.home=="Independent" & vaxstatus=="No vaccine by presentation date"])

table.daycase.coeffs <- tabulate.freqs.regressions(varnames=c("listedgr3", "daycase.recent", "occup"),
                           outcome="CASE",
                           data=cc.severe[care.home=="Independent" & vaxstatus=="No vaccine by presentation date"])

table.opd.coeffs <- tabulate.freqs.regressions(varnames=c("listedgr3", "opd.recent", "occup"),
                           outcome="CASE",
                           data=cc.severe[care.home=="Independent" & vaxstatus=="No vaccine by presentation date"])

#############################################################

winsize.nrs <- 15
freqs.nrs.deaths.tw <- freqs.slidingwindow(dt=cc.all[CASE==1 & fatalcase],
                                winsize=winsize.nrs,
                                datevar="SPECIMENDATE",
                                categoricvar="cod.case",
                                groupvar=NULL)
freqs.nrs.deaths.tw <- freqs.nrs.deaths.tw[cod.case==1]
freqs.nrs.deaths.tw[rsum.allcategories <20, `:=`(p=NA, se.p=NA)]
freqs.nrs.deaths.tw[date >= as.Date("2020-06-01") & date <= as.Date("2020-09-30"), p := NA]

p.nrs <-
    ggplot(data=freqs.nrs.deaths.tw,
           aes(x=date, y=p)) +
    geom_line() + # size=0.002 * (freqs.nrs.deaths.tw$rsum.allcategories)^0.5) + 
    scale_y_continuous(limits=c(0, max(freqs.nrs.deaths.tw$p, na.rm=TRUE)), expand=c(0, 0)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = as.Date("2020-03-01"),
                                                  to = lastdate, by = "month"),
                                         "%d %b")
                             ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
       xlab(paste0("Presentation date: mid-point of ", winsize.nrs, "-day window")) +
    ylab("Frequency") +
    ggtitle("Proportion of fatal cases ascertained only through death certificate")
p.nrs

###########################################################


rmarkdown::render("transmissionBMED-D-21-00640hc.Rmd")
#rmarkdown::render("transmission.Rmd", output_file=paste0("transmission_", linkdate, ".pdf"))
