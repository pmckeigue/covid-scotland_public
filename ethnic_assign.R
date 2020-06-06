######## coding ethnicity ##############################

## dataset is named cc.all

## this script should be edited as a function that will work with any dataset, with user-specified names for source variables and output variables

## Source variables: 

## ethnic.smr - raw SMR categories
## OnolyticsType - raw Onomap types based on name classification
## Geographical Area - regional groupings of OnolyticsType

## Output variables: ethnic5.smr, eth, eth5

## this script

## (1) recodes ETHNIC_smr to specify Chinese as separate category, then collapses to 5 categories: White, South Asian, Chinese, Black, Other

## (2) derives a new variable "eth" from Onomap types

## (3) collapses eth to eth5, which has five categories as above

## We cannot combine SMR coding with Onomap coding because this will introduce non-differential missclassification of ethnicity between cases and controls, especially for those categorised as Black.

## The Onolytics type "MUSLIM" is coded as South Asian - this is reasonable in Scotland where most people with Muslim names are South Asian.


## label the codes used in the raw SMR variable ethnic.smr
##  https://www.ndc.scot.nhs.uk/Dictionary-A-Z/Definitions/index.asp?Search=E&ID=243&Title=Ethnic%20Group
#'1A'='Scottish';
#'1B'='Other British';
#'1C'='Irish';
#'1K'='Gypsy/ Traveller';
#'1L'='Polish';
#'1Z'='Any other white ethnic group';
#'2A'='Any mixed or multiple ethnic groups';
#'3F'='Pakistani, Pakistani Scottish or Pakistani British';
#'3G'='Indian, Indian Scottish or Indian British';
#'3H'='Bangladeshi, Bangladeshi Scottish or Bangladeshi British';
#'3J'='Chinese, Chinese Scottish or Chinese British';
#'3Z'='Other Asian, Asian Scottish or Asian British';
#'4D'='African, African Scottish or African British';
#'4Y'='Other African';
#'5C'='Caribbean, Caribbean Scottish or Caribbean British';
#'5D'='Black, Black Scottish or Black British';
#'5Y'='Other Caribbean or Black';
#'6A'='Arab, Arab Scottish or Arab British';
#'6Z'='Other ethnic group';
#'98'='Refused/Not provided by patient';
#'99'='Not Known'")

## extra codes: 1[DEFGHJ] White, 3[ABC] South Asian, 3[E] Chinese, 3D Other Asian, 4[ABCE]  African

ethnic5.smr <- character(length(ethnic.smr))
ethnic5.smr[grep("^1[A-Z]$", ethnic.smr)] <- "White"
ethnic5.smr[grep("^2[A-Z]$", ethnic.smr)] <- "Other"
ethnic5.smr[grep("^3[ABCFGH]$", ethnic.smr)] <- "South Asian"
ethnic5.smr[grep("^3[EJ]$", ethnic.smr)] <- "Chinese"
ethnic5.smr[grep("^3[DZ]$", ethnic.smr)] <- "Other"
ethnic5.smr[grep("^4[ABCDEY]$", ethnic.smr)] <- "Black"
ethnic5.smr[grep("^5[ABCDY]$", ethnic.smr)] <- "Black"
ethnic5.smr[grep("^6[AZ]$", ethnic.smr)] <- "Other"
ethnic5.smr[grep("^9", ethnic.smr)] <- NA
ethnic5.smr[ethnic5.smr==""] <- NA
ethnic5.smr <- as.factor(ethnic5.smr)
ethnic5.smr <- factor(ethnic5.smr, levels=levels(ethnic5.smr)[c(5, 4, 2, 1, 3)])

if(length(OnolyticsType) > 0) {
    OnolyticsType <- car::recode(OnolyticsType,
                            "'NOT FOUND'=NA; 'INTERNATIONAL'=NA; 'UNCLASSIFIED'=NA; 'VOID'=NA; 'VOID - FORENAME'=NA; 'VOID INITIAL'=NA")
    
    table(OnolyticsType[GeographicalArea=="SOUTH ASIA"])
    table(OnolyticsType[GeographicalArea=="AFRICA"])
    table(OnolyticsType[GeographicalArea=="BRITISH ISLES"])
    table(OnolyticsType[GeographicalArea=="EAST ASIA"])
    table(OnolyticsType[GeographicalArea=="MIDDLE EAST"])
    
    eth <- rep("Other", nrow(cc.all))
    eth[is.na(OnolyticsType)] <- NA
    eth[OnolyticsType=="CHINESE" |
        OnolyticsType=="HONG KONGESE" |
        OnolyticsType=="SINGAPORESE"] <- "Chinese"
    eth[GeographicalArea=="EAST ASIA" &
        eth != "Chinese"] <- "Other Asia & Pacific"
    
    eth[GeographicalArea=="SOUTH ASIA"] <- "South Asian"
    
    eth[GeographicalArea=="BRITISH ISLES"] <- "Britain&Ireland"
    eth[GeographicalArea=="CENTRAL EUROPE" |
        GeographicalArea=="EASTERN EUROPE" |
        GeographicalArea=="NORTHERN EUROPE" |
        GeographicalArea=="SOUTHERN EUROPE" |
        OnolyticsType=="AFRIKAANS"] <- "Other Europe"
    
    eth[GeographicalArea=="AFRICA" &
        OnolyticsType != "AFRIKAANS" &
        OnolyticsType != "LIBYAN"] <- "Black African"
    eth[GeographicalArea=="MIDDLE EAST" &
        OnolyticsType != "MUSLIM"] <- "East Med"
    eth[OnolyticsType == "MUSLIM"] <- "South Asian" # "Muslim, not localized"
    
    eth[OnolyticsType == "BLACK CARIBBEAN"] <- "Black Caribbean"
    
    ## reduce to 5 categories: White, South Asian, Chinese, Black, Other
    ethnic5 <- car::recode(eth, "'Black African'='Black'; 'Black Caribbean'='Black';  'Britain&Ireland'='White';  'Other Europe'='White';  'East Med'='Other'; 'Other Asia & Pacific'='Other'; 'Muslim, not localized'='Other'") 
    ethnic5 <- relevel(as.factor(ethnic5), ref="White")
}
