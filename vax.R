---
title: "Briefing note: efficacy of vaccination in those eligible for shielding in Scotland"
header-includes:
  \usepackage{newunicodechar}
  \let\origquote\quote
  \def\quote{\origquote\itshape}
  \usepackage{graphicx}
  \usepackage{geometry}
  \usepackage{longtable}
  \usepackage{booktabs}
  \usepackage{float}
  \floatplacement{figure}{H}
  \usepackage{array}
  \usepackage{threeparttable}
  \usepackage{longtable}
  \usepackage{lscape}
  \usepackage{pdfpages}
output: pdf_document
always_allow_html: true
urlcolor: blue
linenumbers: false
linkcolor: cyan
---

```{r vax, echo=FALSE, warning=FALSE, message=FALSE}
## vaccination analysis

datadir <- "data/2021-03-16/"

vax.firstdoses <- vacc[, .N, by=c("weekdose1", "vax_type_1")]
ggplot(data=vax.firstdoses, aes(x=weekdose1, y=N, color=vax_type_1)) +
    geom_line() 

keep.vars <- c("CASE", "care.home", "hh.over18", "listedgr3", "shield.any", "shield.group",
               listed.conditions, "shield.group", "vaxstatus", "stratum")

cc.vax <- na.omit(cc.severe[SPECIMENDATE > as.Date("2020-12-01"), ..keep.vars])
cc.vax[, vax14 := as.integer(vaxstatus) > 1]

summary(clogit(data=cc.vax, formula = CASE ~ care.home + hh.over18 + listedgr3 / vax14 + strata(stratum)))$coefficients

table.shieldgroups <- summary(clogit(data=cc.vax, formula = CASE ~ care.home + hh.over18 + shield.group / vax14 + strata(stratum)))$coefficients
table.shieldgroups <- data.table(effect=rownames(table.shieldgroups), table.shieldgroups)
table.shieldgroups[, rateratio := or.ci(coef, `se(coef)`)]
table.shieldgroups[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
table.shieldgroups <- table.shieldgroups[, .(effect, rateratio, pvalue)]

cc.vax[`dm.typeOther/unknown type` == 1, `dm.typeType 2 diabetes` := 1]
cc.vax[, `dm.typeOther/unknown type` := NULL]

controls.vax <- cc.kept[CASE==0 & care.home=="Independent", .(AGE, COVID.age, vax_dose_1)]
controls.vax[, vaxnow := !is.na(vax_dose_1)]
controls.vax[, vaxnow := car::recode(vaxnow,
                                   recodes="FALSE='Not vaccinated';
                                            TRUE='At least one dose'",
                                   as.factor=TRUE,
                                   levels=c("Not vaccinated", "At least one dose"))]

theme_set(theme_gray(base_size = 16))
p.covidage <- ggplot(data=controls.vax, aes(x=AGE, y=COVID.age, color=vaxnow)) +
    geom_point(size=0.5, position="jitter") +
    scale_x_continuous(breaks=seq(20, 80, by=10), limits=c(15, 90)) + 
    scale_y_continuous(breaks=seq(20, 80, by=10), limits=c(15, 90)) +
    theme(legend.title = element_blank()) +
    theme(legend.position=c(0.15, 0.85)) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(title=paste("COVID age versus calendar age by vaccination status on",
                        format.Date(max(vacc$vax_dose_1), '%d %b %Y')),
         caption="Controls matched to severe cases, not resident in care homes. Limits of COVID age are 16 and 90 years",
         x="Calendar age (years)",
         y="COVID age (years)")
p.covidage

png("COVIDage_vaxstatus.png")
p.covidage
dev.off()
rm(p.covidage)

knitr::kable(table.reinfection[, -(3:4)], 
             escape=FALSE, 
             booktabs=TRUE,
			 label="reinfection",
             row.names=TRUE,
             align=rep("r", 5),
             caption="Non-fatal cases in first wave, by reinfection in second wave",
             col.names=c(clean.header(colnames(table.reinfection)[1:2]),
                          "Rate ratio (95\\% CI)", "\\ensuremath{p}-value")) %>%
    column_spec(column=1, width="3cm") %>%
    column_spec(column=2, width="2.5cm") %>%
    column_spec(column=3, width="2.5cm") %>%
   kable_styling(latex_options="HOLD_position", "scale_down")

```


```{r shieldgroups, echo=FALSE, warning=FALSE, message=FALSE}

knitr::kable(table.shieldgroups, 
             escape=FALSE, 
             booktabs=TRUE,
			 label="fatal",
             row.names=FALSE,
             align=c("l", "r, "r),
             caption="Rate ratios associated with at least one dose of any vaccine at least 14 days before",
             col.names=c("Effect", 
                         "Rate ratio (95\\% CI)", "\\ensuremath{p}-value")) %>%
    column_spec(column=1, width="3cm") %>%
    kable_styling(latex_options="HOLD_position", "scale_down")

```
