
scotpop <- read_excel("./Scotland_midyearpop_est2019.xlsx")

# load("casefreqs.RData") ## national cumulative cases and deaths by sex and one year age group

####### incidence and mortality using national population estimates ######################

case.freqs <- data.frame(Age=as.integer(rownames(case.freqs)),
                         Females=as.integer(case.freqs[, 1]),
                         Males=as.integer(case.freqs[, 2]))
case.long <- reshape2::melt(case.freqs, id="Age")
colnames(case.long) <- c("Age", "Sex", "Cases")

death.freqs <- data.frame(Age=as.integer(rownames(death.freqs)),
                         Females=as.integer(death.freqs[, 1]),
                         Males=as.integer(death.freqs[, 2]))
death.long <- reshape2::melt(death.freqs, id="Age")
colnames(death.long) <- c("Age", "Sex", "Deaths")

scotpop.long <- reshape2::melt(scotpop[, -2], id="Age")
colnames(scotpop.long) <- c("Age", "Sex", "Population")

discrim <- merge(scotpop.long, case.long, by=c("Age", "Sex"), all.x=TRUE)
discrim <- merge(discrim, death.long, by=c("Age", "Sex"), all.x=TRUE)
discrim$Cases[is.na(discrim$Cases)] <- 0
discrim$Deaths[is.na(discrim$Deaths)] <- 0

discrim$Sex <- as.factor(discrim$Sex)

discrim$Noncases <- discrim$Population - discrim$Cases
y.cases <- cbind(as.integer(discrim$Cases), as.integer(discrim$Noncases))

discrim$Survivors <- discrim$Population - discrim$Deaths
y.deaths <- cbind(as.integer(discrim$Deaths), as.integer(discrim$Survivors))

cases.model <- glm(formula=y.cases ~ Sex + Age, family="binomial", data=discrim)
deaths.model <- glm(formula=y.deaths ~ Sex + Age, family="binomial", data=discrim)

cases.model.coeffs <- summary(cases.model)$coefficients
deaths.model.coeffs <- summary(deaths.model)$coefficients
logistic.coeffs <- data.frame(severecase=cases.model.coeffs[, 1],
                              death=deaths.model.coeffs[, 1])

male <- discrim$Sex=="Males"
female <- discrim$Sex=="Females"
gam.model.MaleDeaths <- gam::gam(formula=y.deaths[male, ] ~ s(Age), family=binomial("logit"),
                                 data=discrim[male, ])
gam.model.FemaleDeaths <- gam::gam(formula=y.deaths[female, ] ~ s(Age), family=binomial("logit"),
                                   data=discrim[female, ])
gam.model.MaleCases<- gam::gam(formula=y.cases[male, ] ~ s(Age), family=binomial("logit"),
                               data=discrim[male, ])
gam.model.FemaleCases <- gam::gam(formula=y.cases[female, ] ~ s(Age), family=binomial("logit"),
                                  data=discrim[female, ])

gam.male <- data.frame(Cases=car::logit(gam.model.MaleCases$fitted.values),
                       Deaths=car::logit(gam.model.MaleDeaths$fitted.values),
                       Age=discrim$Age[male])
gam.male.long <- reshape2::melt(data=gam.male, id="Age")
colnames(gam.male.long)[2] <- "Status"
gam.male.long$Sex <- "Males"

gam.female <- data.frame(Cases=car::logit(gam.model.FemaleCases$fitted.values),
                       Deaths=car::logit(gam.model.FemaleDeaths$fitted.values),
                       Age=discrim$Age[female])
gam.female.long <- reshape2::melt(data=gam.female, id="Age")
colnames(gam.female.long)[2] <- "Status"
gam.female.long$Sex <- "Females"
gam <- rbind(gam.male.long, gam.female.long)
     
###############################################################

logodds.posterior <- predict(object=cases.model, newdata=discrim, type="link")
logodds.prior <- log(sum(discrim$Cases) / sum(discrim$Noncases))
log.likratio <- logodds.posterior - logodds.prior
discrim$W <- log.likratio / log(2)
lambda1 <- sum(discrim$W * discrim$Cases) / sum(discrim$Cases)
lambda0 <- sum(-discrim$W * discrim$Noncases) / sum(discrim$Noncases)
cases.Lambda.agesex <- 0.5 * (lambda0 +  lambda1)


logodds.posterior <- predict(object=deaths.model, newdata=discrim, type="link")
logodds.prior <- log(sum(discrim$Deaths) / sum(discrim$Survivors))
log.likratio <- logodds.posterior - logodds.prior
discrim$W <- log.likratio / log(2)
lambda1 <- sum(discrim$W * discrim$Deaths) / sum(discrim$Deaths)
lambda0 <- sum(-discrim$W * discrim$Survivors) / sum(discrim$Survivors)
deaths.Lambda.agesex <- 0.5 * (lambda0 +  lambda1)

