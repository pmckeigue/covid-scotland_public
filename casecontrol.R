## analysis script for case-control study
rm(list=ls())

library(car)
library(survival)
library(MASS)
library(wevid)
library(rmarkdown)
library(pander)
library(ggplot2)
library(doParallel)
library(reshape2)
library(readxl)
library(DescTools)
library(icd.data)
library(gam)
library(dplyr)

registerDoParallel(cores=2)

bsource("helperfunctions.R")

old <- TRUE
old <- FALSE # uncomment to use 15 May linkage

nrs <- TRUE
nrs <- FALSE # uncomment to exclude NRS-only deaths

stepwise <- TRUE   
#stepwise <- FALSE ## uncomment to save time if old version still valid



if(old) {
    cc.all <- readRDS("./data/CC_linked_ANON_20200501 (2).rds")
    diagnoses <- readRDS("./data/CC_SMR01_ICD10_x25_ANON_20200501.rds") 
    procedures <- readRDS("./data/CC_SMR01_OPCS4_MAIN.x25_ANON_20200501.rds")
    scrips <- readRDS("./data/CC_PIS_x15_ANON_20200501.rds")
    controls.status <- readRDS("./data/CC_CHI_CHECK_2020-05-05.rds")
} else {
    cc.all <- readRDS("./data/CC_linked_ANON_20200515.rds")
    diagnoses <- readRDS("./data/CC_SMR01_ICD10_x25_ANON_20200515.rds")
    procedures <- readRDS("./data/CC_SMR01_OPCS4_MAIN.x25_ANON_20200515.rds")
    scrips.filename <- "./data/CC_PIS_x15_ANON_20200515.rds" 
    scrips <- readRDS(scrips.filename)[, c("ANON_ID", "bnf_paragraph_code",
                                           "bnf_paragraph_description")] 
    controls.status <- readRDS("./data/CC_CHI_CHECK_2020-05-15.rds")

    ## scrips should be 2 tables to save space
    ## one record per scrip
    scripvars <- c("ANON_ID", "dispensed_date", "bnf_paragraph_code",
                   "formulation_code",
                   "item_strength", "item_strength_uom", "item_code",
                   "num_items", "quantity")
    ## lookup table for subpara code and item code
    drugvars <- c("bnf_paragraph_code", "bnf_paragraph_description",
                  "item_code", "approved_name")
    ## other files for diagnoses, procedures and scrips within the time limit
    ## data/CC_SMR01_ICD10_25_ANON_20200515.rds
    ## data/CC_SMR01_OPCS4_MAIN.25_ANON_20200501.rds
    ## data/CC_PIS_15_ANON_20200515.rds

    scrips.protonpump <-
        readRDS(scrips.filename) %>%
        subset(bnf_paragraph_code=="0103050")
}

#################### add fields to scrips ########################################

## for now, just keep id, paragraph code, paragraph_description

scrips$bnf_paragraph_description <- as.factor(scrips$bnf_paragraph_description)
cat("scrips object uses", object.size(scrips) * 1E-6, "MB\n")

## we have 7 digits on scrips, giving resolution to subpara level only
length(table(substr(scrips$bnf_paragraph_code, 1, 2))) # chapter
length(table(substr(scrips$bnf_paragraph_code, 1, 4))) # chapter, section
length(table(substr(scrips$bnf_paragraph_code, 1, 6)))  # chapter, section, paragraph
length(table(scrips$bnf_paragraph_code)) # 537 groups

## we need integer variables chapter, sectioncode, paracode for use with reshape2::dcast
scrips$chapternum <- as.integer(substr(scrips$bnf_paragraph_code, 1, 2))
scrips$sectioncode <- as.integer(substr(scrips$bnf_paragraph_code, 1, 4))
scrips$paracode <- as.integer(substr(scrips$bnf_paragraph_code, 1, 6))

## recode scrips$bnf.chapter values > 14 or NA to 14
scrips$chapternum[is.na(scrips$chapternum)] <- 14
scrips$chapternum[scrips$chapternum > 14] <- 14

################### short names for ICD chapters ########################

icdchapters <- data.frame(names(icd10_chapters),
                             t(matrix(as.character(unlist(icd10_chapters)), nrow=2)))
colnames(icdchapters) <- c("name", "start", "end")
icdchapters$shortname <- gsub("Diseases of the ", "", icdchapters$name)
icdchapters$shortname <- gsub("Certain ", "", icdchapters$shortname)
icdchapters$shortname <- gsub("conditions originating in the ", "",
                              icdchapters$shortname)
icdchapters$shortname <- gsub("Factors influencing ", "", icdchapters$shortname)
truncate.at <- StrPos(icdchapters$shortname, " |,|-") - 1
truncate.at[is.na(truncate.at)] <- nchar(icdchapters$shortname[is.na(truncate.at)])
icdchapters$shortname <- substr(icdchapters$shortname, 1, truncate.at) 
icdchapters$shortname <- gsub("^health$", "Health_factors", icdchapters$shortname)
icdchapters$start <- as.character(icdchapters$start)
icdchapters$end <- as.character(icdchapters$end)
icdsubchapters <- data.frame(names(icd10_sub_chapters),
                             t(matrix(as.character(unlist(icd10_sub_chapters)), nrow=2)))
colnames(icdsubchapters) <- c("name", "start", "end")
icdsubchapters$start <- as.character(icdsubchapters$start)
icdsubchapters$end <- as.character(icdsubchapters$end)

######################################################################

source("bnfcodes.R")

###############################################

names(cc.all) <- gsub("CASE_NO", "stratum", names(cc.all))
names(cc.all) <- gsub("^SEX$", "sex", names(cc.all))
names(cc.all) <- gsub("imumune", "immune", names(cc.all))
names(cc.all) <- gsub("^ethnic$", "ethnic.old", names(cc.all))
names(cc.all) <- gsub("^CAREHOME$", "care.home", names(cc.all))
names(cc.all) <- gsub("^simd$", "SIMD.quintile", names(cc.all))
names(cc.all) <- gsub("DATE_OF_DEATH", "Date.Death", names(cc.all))
names(cc.all) <- gsub("^age$", "AGE", names(cc.all))

## HAI is based on the ECDC definition of nosocomial infection

cc.all$scrip.any <- as.factor(as.integer(cc.all$ANON_ID %in% scrips$ANON_ID))
cc.all$diag.any <- as.factor(as.integer(cc.all$ANON_ID %in% diagnoses$ANON_ID))

## exclude controls not observable

controls.status <- controls.status[, c("ANON_ID", "CHI.EXTENDED.STATUS", "CHI.DATE.OF.DEATH",
                                   "CHI.Explanation")]

cc.all <- merge(cc.all, controls.status, by="ANON_ID", all.x=TRUE)

print(with(cc.all, table(CHI.Explanation, scrip.any)))
print(with(cc.all, table(CHI.Explanation, diag.any)))

## exclude controls classified on CHI database as no longer current
cc.all <- subset(cc.all, subset=is.na(CHI.Explanation) |
                             CHI.Explanation=="Current and no history" |
                             CHI.Explanation=="Current with history")


## exclude controls already dead on date of test of case they were matched to
controls.deceased <- with(cc.all, CASE==0 &
                                  !is.na(Date.Death) &
                                  Date.Death <= SPECDATE)
cc.all <- cc.all[!controls.deceased, ]                         

cc.all$stratum <- as.integer(cc.all$stratum)
## check that each stratum contains a single case
cat("checking that each stratum contains a single case ...")
table.strata <- tapply(cc.all$CASE, cc.all$stratum, sum) == 1
strata.onecase <- as.integer(names(table.strata)[as.integer(table.strata)==1])
keep <- cc.all$stratum %in% strata.onecase
cc.all <- cc.all[keep, ]
cat("done:", length(which(!keep)), "observations dropped\n")

cc.all$SIMD.quintile <- car::recode(cc.all$SIMD.quintile, "'Unknown'=NA")


cc.all$scripordiag <- as.integer(with(cc.all, as.integer(diag.any)==2 |
                                              as.integer(scrip.any)==2))

cc.all$sex <- car::recode(as.factor(cc.all$sex), "1='Male'; 2='Female'")
cc.all <- within(cc.all, sex <- relevel(sex, ref="Female"))

cc.all$agegr20 <- as.factor(car::recode(as.integer(cc.all$AGE),
                              "0:39='0-39'; 40:59='40-59';
                                  60:74='60-74'; 75:hi='75 or more'"))

cc.all$agegr3 <-
    as.factor(car::recode(cc.all$AGE,
                          "0:59='0-60 years'; 60:74='60-74 years'; 75:hi='75+ years'"))

cc.all$care.home[cc.all$NURSINGHOME==1] <- 1
cc.all$care.home <- as.factor(car::recode(cc.all$care.home, "0='Independent'; 1='Care/nursing home'"))
cc.all <- within(cc.all, care.home <- relevel(care.home, ref="Independent"))

## all cases have nonmissing SPECDATE
## controls are recoded to have same SPECDATE as the case they were matched to
cc.all$deathwithin28 <- 0
cc.all$deathwithin28[with(cc.all,
                          CASE==1 &
                          !is.na(Date.Death) & 
                          Date.Death - SPECDATE >= 0 &
                          Date.Death - SPECDATE <= 28)] <- 1

## fatality rates by type of unit
## 34% in icu, 38% in icu.hdu.ccu,
## 24% in those in hdu but never in icu or icu.hdu.ccu (a small group)
cc.all$unitcategory <- numeric(nrow(cc.all))
cc.all$unitcategory[cc.all$icu==0 & cc.all$hdu==0 & (cc.all$inhosp | cc.all$adm28==1)] <- 0
cc.all$unitcategory[cc.all$icu==0 & cc.all$hdu==1] <- 1
cc.all$unitcategory[cc.all$icu==1 & cc.all$hdu==0] <- 2
cc.all$unitcategory[cc.all$icu==1 & cc.all$hdu==1] <- 3
cc.all$unitcategory[cc.all$CASE==0] <- NA
cc.all$unitcategory <- car::recode(cc.all$unitcategory,
                                   "0='Hospitalized, no HDU or ICU'; 1='HDU only'; 2='ICU only'; 3='HDU and ICU'")

if(!old) {
    print(with(cc.all[cc.all$CASE==1, ], paste.colpercent(table(deathwithin28, unitcategory))))
}

## check this is correct: all those with icu==1 or icu.hdu.ccu==1 should be coded as severe

## integer values > 1 for icu and inhosp may represent days from test to entry
## values of 0 must be for those not admitted, as there are no missing values
with(cc.all[cc.all$CASE==1, ], table(adm28, exclude=NULL))
with(cc.all[cc.all$CASE==1, ], table(icu, exclude=NULL))
if(!old) {
    with(cc.all[cc.all$CASE==1, ], table(hdu, exclude=NULL))
    with(cc.all[cc.all$CASE==1, ], table(icu, hdu, exclude=NULL))
}

## coding of case groups -- check this is correct
if(old) {
    cc.all$nrs_covid_case <- 0
    cc.all$hdu <- 0
} # old dataset did not include test-neg cases ascertained through death registration

cc.all$group <- NA
## Cases ascertained only through NRS coded as D
cc.all$group[cc.all$CASE==1 & cc.all$nrs_covid_case==1] <- "D"

## Cases not tested in hospital and not admitted within 29 days coded as C  
cc.all$group[cc.all$CASE==1 & cc.all$nrs_covid_case==0 &
             cc.all$adm28==0 & cc.all$inhosp==0] <- "C"

## Cases tested in hospit
cc.all$group[cc.all$CASE==1 & cc.all$nrs_covid_case==0 &
             (cc.all$adm28 > 0 | cc.all$inhosp > 0)] <- "B"
table(cc.all$deathwithin28, cc.all$group, exclude=NULL)

cc.all$group[cc.all$CASE==1 & cc.all$nrs_covid_case==0 &
                 (cc.all$icu==1 | cc.all$hdu==1 | cc.all$deathwithin28==1)]  <- "A"

print(table(cc.all$CASE, cc.all$group, exclude=NULL))

## assign controls to same group as matched case i.e. A, B, C and create a new variable named casegroup
casegroups <- cc.all[cc.all$CASE==1, ][, c("stratum", "group")]
colnames(casegroups)[2] <- "casegroup"
cc.all <- merge(cc.all, casegroups, by=c("stratum"), all.x=T)
table(cc.all$CASE, cc.all$casegroup)
with(cc.all[cc.all$CASE==1, ], table(casegroup, deathwithin28, exclude=NULL))

cc.all$fatalcase <- as.integer(cc.all$CASE==1 & cc.all$deathwithin28==1)

if(old) {
    cc.all$casegroup <- car::recode(cc.all$casegroup,
                                    "'A'='Critical care or fatal'; 'B'='Hospitalised, not severe'; 'C'='Test-positive, not hospitalised'")
} else {
    cc.all$casegroup <- car::recode(cc.all$casegroup,
                                    "'A'='Critical care or fatal'; 'B'='Hospitalised, not severe'; 'C'='Test-positive, not hospitalised'; 'D'='COVID mention on death cert, no test'")
}   

cc.all$casegroup <- as.factor(cc.all$casegroup)
cc.all <- within(cc.all, casegroup <- relevel(casegroup, ref="Critical care or fatal"))

print(paste.colpercent(with(cc.all[cc.all$CASE==1, ], table(care.home, casegroup))))

######## coding ethnicity ##############################

OnolyticsType <- cc.all$OnolyticsType
GeographicalArea <- cc.all$GeographicalArea
ethnic.smr <- as.character(cc.all$ETHNIC_SMR_LAST)

source("ethnic_assign.R")

cc.all$ethnic5.smr <- ethnic5.smr

## recode SMR ethnicity to 4 categories: White, Black, South Asian, Other
cc.all$ethnic4.smr <- as.factor(car::recode(cc.all$ethnic5.smr, "'Chinese'='Other'"))
cc.all <- within(cc.all, ethnic4.smr <- factor(ethnic4.smr,
                                               levels=levels(ethnic4.smr)[c(4, 3, 1, 2)]))   
#if(length(OnolyticsType) > 0) {
    cc.all$ethnic5 <- ethnic5
   
    ## tabulate ONOMAP ethnicity against SMR ethnicity
    table.ethnic <- table(cc.all$ethnic5, cc.all$ethnic5.smr, exclude=NULL)
    
    tn <- table(cc.all$ethnic5, cc.all$ethnic5.smr)
    SouthAsian.sensitivity <- 100 * tn[5, 2] / sum(tn[, 2])
    SouthAsian.specificity <- 100 * (sum(tn[, -2]) - sum(tn[5, ]) + tn[5, 2]) / sum(tn[, -2])
    sum.xtabulate <- sum(tn)
                                                       
    table.ethnic <- paste.colpercent(table.ethnic)

    ## recode ONOMAP ethnicity to 4 categories: White, South Asian, Chinese, Other
    cc.all$ethnic4 <- car::recode(cc.all$ethnic5, "'Black'='Other'")
    cc.all <- within(cc.all, ethnic4 <- relevel(as.factor(ethnic4), ref="White"))
    cc.all$ethnic4 <- factor(cc.all$ethnic4, levels=levels(cc.all$ethnic4)[c(1, 4, 2, 3)])
   
    ## recode ONOMAP ethnicity to 3 categories: White, South Asian, Other
    cc.all$ethnic3 <- car::recode(cc.all$ethnic4, "'Chinese'='Other'")
    cc.all <- within(cc.all, ethnic3 <- relevel(as.factor(ethnic3), ref="White"))
    cc.all$ethnic3 <- factor(cc.all$ethnic3, levels=levels(cc.all$ethnic3)[c(1, 3, 2)])
#}


########################################################################################
###  diabetes based on combining Sci-Diabetes records with ICD codes and drugs 
## add in BNF codes 6.1 for diabetes drugs and
## E10 to E14 for diabetes

## immune.any includes primary immunodeficiency and secondary immunosuppression

########## coding listed conditions ####################
ids.icd.diabetes <- unique(diagnoses$ANON_ID[grep("^E1[0-4]", diagnoses$ICD10)])

## 802 other immunomodulating drugs
## Methotrexate and chloroquine appear in musculoskeletal chapter 
ids.bnf.diabetes <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 601])

ids.diabetes.extra <- unique(c(ids.icd.diabetes, ids.bnf.diabetes))

# recode diabetes type
cc.all$dm.type <- as.integer(cc.all$dm.type)
## missing recoded as zero
cc.all$dm.type[is.na(cc.all$dm.type)] <- 0

## add in extra cases notified directly from SCI-Diabetes register, without assignment
## of diabetes type from SDRN database
if(!old) {
    cc.all$dm.type[cc.all$dm.type==0 & cc.all$diab.reg==1] <- 3
}

print(table(cc.all$dm.type, cc.all$ANON_ID %in% ids.diabetes.extra))

browser("diabetes.extra")

## code diagnoses detected from discharges or BNF codes as unknown type
## we could classify those not on insulin as definite Type 2 but Helen says no

## REVISION: for consistency with the diabetes paper, we will not include the extra cases identified through diagnostic codes or drug codes
# cc.all$dm.type[cc.all$dm.type==0 & cc.all$ANON_ID %in% ids.diabetes.extra] <- 3

cc.all$dm.type <- 
  as.factor(car::recode(cc.all$dm.type, 
                        "c(0, 10, 17, 97)='Not diabetic';
                         c(1, 101, 102)='Type 1 diabetes';
                         c(2, 202, 203)='Type 2 diabetes'; 
                         3:9='Other/unknown type';
                         11:16='Other/unknown type';
                         18:96='Other/unknown type';
                         98:100='Other/unknown type'"))

cc.all <- within(cc.all, dm.type <- relevel(dm.type, ref="Not diabetic"))
cc.all$dm.type <- factor(cc.all$dm.type, levels=levels(cc.all$dm.type)[c(1, 3, 4, 2)])

## define an indicator variable for any diabetes
cc.all$diabetes.any <- as.integer(cc.all$dm.type != "Not diabetic")
cc.all$diabetes.any <- as.factor(car::recode(cc.all$diabetes.any,
                                        "0='Not diabetic'; 1='Diabetic'"))
cc.all <- within(cc.all, diabetes.any <- relevel(diabetes.any, ref="Not diabetic"))

####################################################################

## code care home residents as 1, other as 0

cc.all$cats3 <-  as.integer(cc.all$care.home == "Care/nursing home")
cc.all$cats3[cc.all$cats3==0 & cc.all$listed.any==1] <- 2 
cc.all$cats3[cc.all$cats3==0 & cc.all$listed.any==0] <- 3
cc.all$cats3 <- as.factor(car::recode(cc.all$cats3, "1='Care/nursing home';
                               2='Independent living, listed condition';
                               3='Independent living, no listed condition'"))     

####### drug exposure specifically relevant to proton pump ##################

ids.protonpump <- unique(scrips$ANON_ID[scrips$bnf_paragraph_code == "0103050"])
cc.all$protonpump <- as.factor(as.integer(cc.all$ANON_ID %in% ids.protonpump))
cc.all$y.protonpump <- as.integer(cc.all$protonpump =="1")

TW <- 120 # time window in days

if(!old) {
    ## merge SPECDATE into scrips.protonpump
    scrips.protonpump <- merge(scrips.protonpump, cc.all[, c("ANON_ID", "SPECDATE")],
                               by="ANON_ID", all.x=TRUE) 
    scrips.protonpump$daysbefore <- as.integer(scrips.protonpump$SPECDATE -
                                               scrips.protonpump$dispensed_date)
    ## minimum value is 16 days -- the 15-day cutoff has been applied to the scrips table

    scrips.protonpump$dispensing.days <- as.integer(scrips.protonpump$SPECDATE - 15 - 
                                                    as.Date("2019-06-01"))
    scrips.protonpump$intervalTWday <- ceiling((scrips.protonpump$daysbefore - 15) / TW)
    scrips.protonpump$intervalTWday[scrips.protonpump$intervalTWday > 3] <- 3

## https://www.whocc.no/atc_ddd_index/?code=A02BC&showdescription=yes gives defined daily doses
##
##ATC code  	Name  	DDD 	 U 	Adm.R	 Note
##A02BC05 	esomeprazole 	30 	mg 	O
##A02BC03 	lansoprazole 	30 	mg 	O 
##A02BC01 	omeprazole 	20 	mg 	O 
##A02BC02 	pantoprazole 	40 	mg 	O 	
##A02BC04 	rabeprazole 	20 	mg 	O 	
    DDD5 <- c(30, 30, 20, 40, 20)

    scrips.protonpump$dose <- scrips.protonpump$item_strength * scrips.protonpump$quantity 

    ## calculate dose of each drug over entire period
    dose.protonpump <-  reshape2::dcast(data=scrips.protonpump,
                                        formula=ANON_ID + dispensing.days ~ approved_name,
                                        value.var="dose",
                                        fun.aggregate=sum)

    for(j in 1:5) { # loop over approved names to divide by DDD x dispensing.days 

        dose.protonpump[, j + 2] <- dose.protonpump[, j + 2] /
            (DDD5[j] * dose.protonpump$dispensing.days) 
    }
    dose.protonpump$DDDs.all <- rowSums(dose.protonpump[, -(1:2)])

    ###############  average dose in each TW-day interval ######################## 
    doseTWday.protonpump <-  reshape2::dcast(data=scrips.protonpump,
                                    formula=ANON_ID + intervalTWday ~ approved_name,
                                    value.var="dose",
                                    fun.aggregate=sum)
    ## strictly, each interval should be divided by the number of days in that interval for that individual
    ## for now, use TW days
    for(j in 1:5) { # loop over approved names to divide by DDD x TW 
        doseTWday.protonpump[, j + 2] <- doseTWday.protonpump[, j + 2] /
            (DDD5[j] * TW)
    }
    doseTWday.protonpump$DDDsTWday.all <- rowSums(doseTWday.protonpump[, -(1:2)])
    # print(summary(doseTWday.protonpump))

    ## drop columns for individual drugs, and cast again to get one column for each interval
    doseTWday.protonpump <- doseTWday.protonpump[, c("ANON_ID", "intervalTWday", "DDDsTWday.all")]
    doseTWday.protonpump <- doseTWday.protonpump[!is.na(doseTWday.protonpump$intervalTWday), ]

    doseTWday.protonpump.wide <-  reshape2::dcast(data=doseTWday.protonpump,
                                    formula=ANON_ID ~ intervalTWday,
                                    value.var="DDDsTWday.all",
                                    fun.aggregate=sum)
    colnames(doseTWday.protonpump.wide)[-1] <- c("DDD.interval1", "DDD.interval2", "DDD.interval3") 

    cc.all <- merge(cc.all, dose.protonpump, by="ANON_ID", all.x=TRUE)
    cc.all$DDDs.all[is.na(cc.all$DDDs.all)] <- 0
    cc.all$ESOMEPRAZOLE[is.na(cc.all$ESOMEPRAZOLE)] <- 0
    cc.all$LANSOPRAZOLE[is.na(cc.all$LANSOPRAZOLE)] <- 0
    cc.all$OMEPRAZOLE[is.na(cc.all$OMEPRAZOLE)] <- 0
    cc.all$PANTOPRAZOLE[is.na(cc.all$PANTOPRAZOLE)] <- 0
    cc.all$RABEPRAZOLE[is.na(cc.all$RABEPRAZOLE)] <- 0
    cc.all$dispensing.days <- as.integer(cc.all$SPECDATE - as.Date("2019-06-01"))
                                        #cc.all$DDDs.average <- cc.all$DDDs.all / cc.all$dispensing.days
    cc.all$DDDsgr <- 0.5 * ceiling(2 * cc.all$DDDs.all)
    cc.all$DDDsgr <- as.factor(car::recode(cc.all$DDDsgr, "2:hi='2 or more'"))
    
    cc.all <- merge(cc.all, doseTWday.protonpump.wide, by="ANON_ID", all.x=TRUE)
    cc.all$DDD.interval1[is.na(cc.all$DDD.interval1)] <- 0
    cc.all$DDD.interval2[is.na(cc.all$DDD.interval2)] <- 0
    cc.all$DDD.interval3[is.na(cc.all$DDD.interval3)] <- 0
}

ids.nonopioid.analgesic <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 40701])
cc.all$nonopioid.analgesic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                   ids.nonopioid.analgesic))
ids.antiplatelet <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 209])
cc.all$antiplatelet <- as.factor(as.integer(cc.all$ANON_ID %in% ids.antiplatelet))

ids.nsaid <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 100101])
cc.all$nsaid <- as.factor(as.integer(cc.all$ANON_ID %in% ids.nsaid))

ids.opioid.analgesic <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 40702])
cc.all$opioid.analgesic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                      ids.opioid.analgesic))

ids.antipsychotic <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 40201])
cc.all$antipsychotic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                      ids.antipsychotic))

ids.osmotic.laxative <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 10604])
cc.all$osmotic.laxative <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                ids.osmotic.laxative))

ids.anticoagulant.any <-
    unique(scrips$ANON_ID[scrips$bnf_paragraph_code == "0208010" |
                          scrips$bnf_paragraph_code == "0208020"])
cc.all$anticoagulant.any <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                ids.anticoagulant.any))
if(!old) {
ids.hydroxychloroquine <- readRDS(scrips.filename) %>%
    subset(subset=approved_name=="HYDROXYCHLOROQUINE SULFATE", select=ANON_ID)
ids.hydroxychloroquine <- unique(ids.hydroxychloroquine$ANON_ID)
cc.all$hydroxychloroquine <- as.factor(as.integer(cc.all$ANON_ID %in% ids.hydroxychloroquine))
}
############# some tables for ethnicity report ##########################
########################################################################

testpositives.ethnic.smr <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                  table(ethnic5.smr, casegroup)), 1)

table.testpositives.demog.ethnicsmr <-
    tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home", "SIMD.quintile"),
                               data=cc.all[!is.na(cc.all$ethnic5.smr), ])


table.hospitalized.demog.ethnicsmr <-
    tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home", "SIMD.quintile"),
                               data=cc.all[(cc.all$casegroup=="Hospitalised, not severe" |
                                            cc.all$casegroup=="Critical care or fatal") &
                                           !is.na(cc.all$ethnic5.smr), ])

if(FALSE) {
## tabulate ethnicity by case group
testpositives.ethnic <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                              table(ethnic4, casegroup)), 1)

testpositives.carehome <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                table(ethnic4, care.home)), 0)

testpositives.healthboard <- t(paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                table(ethnic4, HBRES_NAME)), 0))

table.testpositives.demog <-
    tabulate.freqs.regressions(varnames=c("ethnic4", "care.home", "SIMD.quintile"),
                               data=cc.all)

table.hospitalized.demog <-
    tabulate.freqs.regressions(varnames=c("ethnic4", "care.home", "SIMD.quintile"),
                               data=cc.all[cc.all$casegroup=="Hospitalised, not severe" |
                                           cc.all$casegroup=="Critical care or fatal", ])
}

############ save hospitalized non-severe for later use as training dataset #############

hosp <- cc.all$casegroup=="Hospitalised, not severe"
## merge BNF chapters, one variable per subpara
chnums = 1:13
cc.hosp <- merge.bnfsubparas(chnums=chnums, data=cc.all[hosp, ])

saveRDS(cc.hosp, file="cchosp.rds")
rm(cc.hosp)

########## restrict to severe cases and matched controls ###################### 
cat("Restricting to severe cases and matched controls\n")
if(!nrs) {
    cc.severe <- cc.all[cc.all$casegroup=="Critical care or fatal", ]
} else {
    cc.severe <- cc.all[cc.all$casegroup=="Critical care or fatal" |
                      cc.all$casegroup=="COVID mention on death cert, no test", ]
}

rm(cc.all)

## merge drugs, one variable per chapter
scrips.wide <- reshape2::dcast(scrips, ANON_ID ~ chapternum, fun.aggregate=length, 
                               value.var="chapternum")
shortnames.cols <-  bnfchapters$shortname[match(as.integer(colnames(scrips.wide)[-1]),
                                                as.integer(bnfchapters$chapternum))]
colnames(scrips.wide)[-1] <- paste("BNF", colnames(scrips.wide)[-1], shortnames.cols,
                                   sep="_")
cc.severe <- merge(cc.severe, scrips.wide, by="ANON_ID", all.x=TRUE)

bnfcols <- grep("^BNF", colnames(cc.severe))
for(j in bnfcols) {
    cc.severe[, j][is.na(cc.severe[, j])] <- 0
    cc.severe[, j][cc.severe[, j] > 1] <- 1
    cc.severe[, j] <- as.factor(cc.severe[, j])
}

## merge BNF chapters, one variable per subpara
chnums = 1:13
cc.severe <- merge.bnfsubparas(chnums=chnums, data=cc.severe)
subparacols <- grep("^subpara\\.", names(cc.severe))
x <- apply(cc.severe[, subparacols], 2, as.character)
x <- apply(x, 2, as.integer)
cc.severe$numdrugs.subpara <- rowSums(x)
cc.severe$numdrugsgr <- 3 * ceiling(cc.severe$numdrugs.subpara / 3)
cc.severe$numdrugsgr <- as.factor(car::recode(cc.severe$numdrugsgr,
                                              "3='1 to 3'; 6='4 to 6'; 9='7 to 9';
                                              12='10 to 12'; 15:hi='>12'"))
cc.severe$numdrugsgr <- factor(cc.severe$numdrugsgr,
                               levels=levels(cc.severe$numdrugsgr)[c(2, 3, 5, 6, 4, 1)])
table(cc.severe$numdrugsgr)

cc.severe$numdrugs.notppi <- cc.severe$numdrugs.subpara - cc.severe$y.protonpump
cc.severe$numdrugs.notppi.gr <- 3 * ceiling(cc.severe$numdrugs.notppi / 3)
cc.severe$numdrugs.notppi.gr <- as.factor(car::recode(cc.severe$numdrugs.notppi.gr, "15:hi='>12'"))
cc.severe$numdrugs.notppi.gr <- factor(cc.severe$numdrugs.notppi.gr,
                               levels=levels(cc.severe$numdrugs.notppi.gr)[c(2, 4:6, 3, 1)])
table(cc.severe$numdrugs.notppi.gr)

cardiovasc.subparacols <- grep("^subpara\\.2", names(cc.severe))
x <- apply(cc.severe[, cardiovasc.subparacols], 2, as.character)
x <- apply(x, 2, as.integer)
cc.severe$numdrugs.cardiovasc <- rowSums(x)
cc.severe$numdrugs.notcardiovasc <- cc.severe$numdrugs.subpara - cc.severe$numdrugs.cardiovasc

## merge ICD diagnoses

## chapter is assigned as the position of the first element in icdchapters$start that x is greater than or equal to
unique.diagnoses <- as.character(unique(diagnoses$ICD10))
chapter <- integer(length(unique.diagnoses))
subchapter <- integer(length(unique.diagnoses))

for(i in 1:length(unique.diagnoses)) {
    chapter[i] <- min(which(substr(unique.diagnoses[i], 1, 3) <= icdchapters$end))
    ## subchapter is row in icdsubchapters table 
    subchapter[i] <- min(which(substr(unique.diagnoses[i], 1, 3) <= icdsubchapters$end))
}
unique.diagnoses <- data.frame(ICD10=unique.diagnoses, chapter=chapter, subchapter=subchapter)
diagnoses <- merge(diagnoses, unique.diagnoses, by="ICD10", all.x=TRUE)

diagnoses.wide <- reshape2::dcast(diagnoses, ANON_ID ~ chapter, fun.aggregate=length,
                                  value.var="chapter")
colnames(diagnoses.wide)[-1] <-
    paste0("Ch.", as.integer(colnames(diagnoses.wide)[-1]), "_", 
           icdchapters$shortname[as.integer(colnames(diagnoses.wide)[-1])])
## drop rare chapters
diagnoses.wide <- diagnoses.wide[, colSums(diagnoses.wide) > 20]

cc.severe <- merge(cc.severe, diagnoses.wide, by="ANON_ID", all.x=TRUE)
icdcols <- grep("^Ch.", colnames(cc.severe))
for(j in icdcols) {
    cc.severe[, j][is.na(cc.severe[, j])] <- 0
    cc.severe[, j][cc.severe[, j] > 1] <- 1
    cc.severe[, j] <- as.factor(cc.severe[, j])
}

cc.severe$num.icdchapters <- rowSums(matrix(as.integer(as.matrix(cc.severe[, icdcols])),
                                            nrow=nrow(cc.severe)))
cc.severe$num.icdchapters.gr <-
    as.factor(car::recode(cc.severe$num.icdchapters,
                          "0='No discharge records'; 1:2='1-2 ICD-10 chapters'; 3:hi='3 or more chapters'")
                         )
cc.severe <- within(cc.severe, num.icdchapters.gr <- relevel(num.icdchapters.gr, ref="No discharge records"))

###########################################

source("comorbidity.R")

## 8 listed conditions designated by NHS
listed.conditions <- c("dm.type", "IHD.any", "heart.other.any", "oad.any",
                       "ckd.any", "neuro.any", "liver.any", "immune.any")
#if(!old) {
#    listed.conditions <- c(listed.conditions, "circulatory.other")
#}

############ extract predefined disease categories #################
## as these are coded as factors, lowest level will be 1

cc.severe$listed.any <-
    as.factor(as.integer(with(cc.severe,
                              diabetes.any=="Diabetic" | IHD.any==1 |
                              heart.other.any==1 |
                              ckd.any==1 | oad.any==1 |
                              neuro.any==1 | liver.any==1 | immune.any==1)))

####################################################################

## code care home residents as 1, other as 0
cc.severe$cats3 <-  as.integer(cc.severe$care.home == "Care/nursing home")
cc.severe$cats3[cc.severe$cats3==0 & cc.severe$listed.any==1] <- 2 
cc.severe$cats3[cc.severe$cats3==0 & cc.severe$listed.any==0] <- 3
cc.severe$cats3 <- as.factor(car::recode(cc.severe$cats3, "1='Care/nursing home';
                               2='Independent living, listed condition';
                               3='Independent living, no listed condition'"))     

########### variable lists for tabulating

demog <- c("ethnic3", "SIMD.quintile", "care.home")

if(length(OnolyticsType)==0) {
      demog <- c("SIMD.quintile", "care.home")
}

demog.smr <- c("ethnic4.smr", "SIMD.quintile", "care.home")

bnf.chapternames <- colnames(cc.severe)[bnfcols]
drugs <- bnf.chapternames

icd.chapternames <- colnames(cc.severe)[icdcols]
conditions <- icd.chapternames

### incidence and mortality using national population estimates #####

case.freqs <- with(cc.severe[cc.severe$CASE==1, ], table(AGE, sex, exclude=NULL))
death.freqs <- with(cc.severe[cc.severe$fatalcase==1, ], table(AGE, sex, exclude=NULL))

save(case.freqs, death.freqs, file="casefreqs.agesex.RData")

source("incidencemortality.R")

###############################################################

if(FALSE) { # coding antihypertensives
antihypertensive.classes <- c("vasodilator_ah6.bnf", 
                              "centrally_acting_ah6.bnf",  
                              "adrenergic_neurone_block6.bnf",  
                              "alpha_adrenoceptor6.bnf", 
                              "ace6.bnf",  
                              "angio6.bnf",  
                              "renin_angiotensin6.bnf",  
                              "thiazides6.bnf",  
                              "calcium_channel6.bnf")
antihypertensives <- c(antihypertensive.classes[4:9], "antihypertensive.other")
}

########################################################

table.severe.demog <-
    tabulate.freqs.regressions(varnames=demog,
                               data=cc.severe)

cc.severe$diag.other <- as.integer(cc.severe$listed.any==0 & cc.severe$diag.any==1)
cc.severe$diag.other <- as.factor(car::recode(cc.severe$diag.other,
                                  "0='Listed condition or no admission';
                                   1='No listed condition, but other admission diagnosis'"))

table.scripordiag <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- with(cc.severe[cc.severe$agegr20==agegr, ],
              table(scripordiag, CASE))
    colnames(x) <- c("Controls", "Cases")
    x <- paste.colpercent(x)
    y <- with(cc.severe[cc.severe$agegr20==agegr, ],
              table(listed.any, CASE))
    y <- paste.colpercent(y)
    y <- y[2, , drop=FALSE]
    rownames(y) <- "Any listed condition"
    z <- with(cc.severe[cc.severe$agegr20==agegr, ],
              table(listed.any==0 & diag.any==1, CASE))
    z <- paste.colpercent(z)
    z <- z[2, , drop=FALSE]
    rownames(z) <- "No listed condition, but other diagnosis"
    
    x <- rbind(x[1, , drop=FALSE], y, z, x[2, , drop=FALSE])
    table.scripordiag <- cbind(table.scripordiag, x)
}


table.scripordiag.fatal <- NULL
for(agegr in levels(cc.severe$agegr3)) {
    x <- with(cc.severe[cc.severe$agegr3==agegr, ],
              table(scripordiag, fatalcase))
    colnames(x) <- c("Controls", "Fatal cases")
    x <- paste.colpercent(x)
    # x <- x[1, , drop=FALSE]
    rownames(x) <- c("No scrip or diagnosis", "Scrip or diagnosis")
    table.scripordiag.fatal <- rbind(table.scripordiag.fatal, x)
}

cc.severe$scripordiag <- as.factor(cc.severe$scripordiag)
varnames.listed <- c("care.home", "scrip.any", "diag.any", "listed.any", "diag.other",
                     "scripordiag", listed.conditions)

table.agegr <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- univariate.tabulate(varnames=varnames.listed, 
                             outcome="CASE",
                             data=cc.severe[cc.severe$agegr20==agegr, ],
                             drop.reflevel=FALSE)
    table.agegr <- cbind(table.agegr, x)
}
                                        #

freqs.all <- univariate.tabulate(varnames=c("deathwithin28", varnames.listed, "listed.any"),
                             outcome="CASE",
                             data=cc.severe,
                             drop.reflevel=FALSE)

keep.varnames <- logical(length(varnames.listed))
for(i in 1:length(varnames.listed)) {
    x <- cc.severe[, match(varnames.listed[i], colnames(cc.severe))]
    exposed <- as.integer(x) > 1
    a <- with(cc.severe[exposed, ], table(agegr3, CASE))
    keep.varnames[i] <- !any(as.integer(a)==0)
}

tables.agegr <- vector("list", length(levels(cc.severe$agegr3)))
for(i in 1:length(levels(cc.severe$agegr3))) {
    agegr <- levels(cc.severe$agegr3)[i]
    tables.agegr[[i]] <-
        tabulate.freqs.regressions(varnames=c("care.home", "scrip.any",
                                              "diag.any", listed.conditions),
                                   data=cc.severe[cc.severe$agegr3==agegr, ])
}

table.agegr.all <- tabulate.freqs.regressions(varnames=c("care.home", "scrip.any",
                                              "diag.any", listed.conditions),
                                               data=cc.severe)

## demographic vars
table.demog.aug <- tabulate.freqs.regressions(varnames=demog, data=cc.severe)

## separate analysis using SMR ethnicity 
table.ethnicsmr <- univariate.tabulate(varnames="ethnic4.smr", outcome="CASE",
                                       data=cc.severe[!is.na(cc.severe$ethnic4.smr), ],
                                       drop.reflevel=FALSE)
univariate.ethnicsmr <-
    univariate.clogit(varnames="ethnic4.smr",
                      data=cc.severe[!is.na(cc.severe$ethnic4.smr), ],
                      add.reflevel=TRUE)
table.ethnicsmr.aug <- combine.tables2(table.ethnicsmr, univariate.ethnicsmr)
rownames(table.ethnicsmr.aug) <- replace.names(rownames(table.ethnicsmr.aug))

#### 5 ethnic groups for HPS report
table.ethnic5smr <-
    tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home",
                                          "SIMD.quintile"),
                               outcome="CASE",
                               data=cc.severe[!is.na(cc.severe$ethnic5.smr), ])
rownames(table.ethnic5smr) <- replace.names(rownames(table.ethnic5smr))

## listed conditions
table.listed.conditions.lt60 <-
    tabulate.freqs.regressions(varnames=listed.conditions,
                               data=cc.severe[cc.severe$AGE < 60, ])
table.listed.conditions.ge60 <-
    tabulate.freqs.regressions(varnames=listed.conditions,
                               data=cc.severe[cc.severe$AGE >= 60, ])

## full multivariate model variables -- for this use dm.type rather than diabetes.any 
multivariate.all <-
    multivariate.clogit(varnames=c(demog, "dm.type", listed.conditions[-1],
                                   "diag.any", conditions, "scrip.any", drugs,
                                   "protonpump"),
                        data=cc.severe, add.reflevel=TRUE)

################# restrict to those without listed conditions #############

nocare <- cc.severe$care.home=="Independent"
notlisted <- cc.severe$listed.any == 0
nocare.notlisted <- nocare & notlisted

## conditions
table.conditions.aug <- tabulate.freqs.regressions(varnames=conditions, 
                                                   data=cc.severe[notlisted, ])
cat("Tabulating ICD subchapter diagnoses ...")
## tabulate subchapters in ICD chapters of interest
table.icdchapter2 <- tabulate.icdchapter(chnum=2, data=cc.severe[notlisted, ])
table.icdchapter7 <-  tabulate.icdchapter(chnum=7, data=cc.severe[notlisted, ])
table.icdchapter11 <- tabulate.icdchapter(chnum=11, data=cc.severe[notlisted, ])

table.icdsubchapters <- NULL
for(i in 1:20) {
    table.icdsubchapters <-
        rbind(table.icdsubchapters,
              tabulate.icdchapter(chnum=i, data=cc.severe[notlisted, ], minrowsum=50))
}
#table.icdsubchapters <- table.icdsubchapters[grep("ensuremath",
#                                                  table.icdsubchapters$u.pvalue), ]
cat("done\n")

#########################################################################

## drugs 
table.drugs.aug <- tabulate.freqs.regressions(varnames=drugs, 
                                              data=cc.severe[notlisted, ])

#############################################################################

## tabulate scrip.any effect by carehome and listed.any

## code care home residents as 1, other as 0
cc.severe$cats3 <-  as.integer(cc.severe$care.home == "Care/nursing home")
cc.severe$cats3[cc.severe$cats3==0 & cc.severe$listed.any==1] <- 2 
cc.severe$cats3[cc.severe$cats3==0 & cc.severe$listed.any==0] <- 3
cc.severe$cats3 <- as.factor(car::recode(cc.severe$cats3, "1='Care/nursing home';
                               2='Independent living, listed condition';
                               3='Independent living, no listed condition'"))     

################################################################

## backwards selection of smallest subset of BNF chapters that explains most of the scrip.any effect
    
## tabulate associations with drug chapters in those not in care homes and without listed conditions 
table.drugs.nocare.notlisted <- tabulate.freqs.regressions(varnames=drugs, 
                                                           data=cc.severe[nocare.notlisted, ])

######## stepwise regressions use saved version #####################

nfold <- 4
source("stepwise.R")

#####################################################################

if(!old) source("pharmaco.R")

#####################################################
if(old) {
    rmarkdown::render("casecontrol.Rmd", output_file="casecontrol5May.pdf")
} else  {
    if(!nrs) {
        rmarkdown::render("casecontrol.Rmd", output_file="casecontrol15May.pdf")
        rmarkdown::render("pharmaco.Rmd", output_file="pharmaco.pdf")
    } else {
        rmarkdown::render("casecontrol.Rmd", output_file="casecontrol15May_withNRS.pdf")
    } 
}

## remove large objects from memory

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

if(old) {
    save.image(file="casecontrol5May.RData")
} else {
    if(!nrs) {
        save.image(file="casecontrol15May.RData")
    } else {
        saveRDS(table.agegr, file="table.agegr.withNRS.rds")
    }
}
