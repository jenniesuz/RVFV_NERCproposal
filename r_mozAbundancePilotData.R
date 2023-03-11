library(ggplot2)
library(lme4)
library(plyr)
library(MASS)


#************ read in data *************
readxl::excel_sheets("SEP22_PILOT_CATCH_DATA.xlsx")


# combine spreadsheets
village <- data.frame(readxl::read_excel("SEP22_PILOT_CATCH_DATA.xlsx","Village Table"))
trapSet <- data.frame(readxl::read_excel("SEP22_PILOT_CATCH_DATA.xlsx","Trap Set - September"))
catch <- data.frame(readxl::read_excel("SEP22_PILOT_CATCH_DATA.xlsx","Trap Catch"))

dat <- merge(village,trapSet,by="V_ID",all.y=T)
dat <- merge(dat,catch,by="TS_ID",all.y=T)

# focus only on culex - need to include zero for traps where no culex caught

culex <- dat[dat$Genus %in% "Culex",] # 61 of 120 rows
length(unique(culex$TS_ID)) # 28 traps
length(unique(dat$TS_ID))   # of 54 traps
notCulex <- dat[!dat$TS_ID %in% unique(culex$TS_ID),]
# confirm result
length(unique(notCulex$TS_ID))
notCulex$Genus <- "Culex"
notCulex$Total.Females <- 0

datC <- rbind.data.frame(culex,notCulex)

dat <- ddply(datC,.(District,Village,Date.deployed,TS_ID),summarise,moz=sum(Total.Females))
#******************************************

dat$type <- NA
dat$wardType[dat$Village %in% "Mapea"] <- "mixed"
dat$wardType[dat$Village %in% "Mtakuja"] <- "rural"
dat$wardType[dat$Village %in% "Mheza"] <- "mixed"
dat$wardType[dat$Village %in% "Ngarenaro St"] <- "mixed"
dat$wardType[dat$Village %in% "Esilalaei"] <- "rural"
dat$wardType[dat$Village %in% "Ruvu darajani"] <- "rural"
dat$wardType[dat$Village %in% "Kikavu chini"] <- "rural"
dat$wardType[dat$Village %in% "Kisutu"] <- "mixed"

hist(dat$moz)
hist(log10(dat$moz)+1)
qqnorm(dat$moz)
qqnorm(log10(dat$moz+1))

ggplot(dat) +
  geom_boxplot(aes(x=District,y=moz))
ggplot(dat) +
  geom_boxplot(aes(x=Village,y=moz))
ggplot(dat) +
  geom_boxplot(aes(x=wardType,y=moz))

length(dat$moz[dat$moz==0]) #26

median(dat$moz[dat$wardType %in% "rural"])
#[1] 0
median(dat$moz[dat$wardType %in% "mixed"])
#[1] 1.5
mean(dat$moz[dat$wardType %in% "mixed"])
#[1] 7.866667
mean(dat$moz[dat$wardType %in% "rural"])
#[1] 4.583333

# probably going to need zero-inflated negative binomial mixed effects model
# unless more villages have presence this time around

# analyse presence absence

dat$present <- NA
dat$present[dat$moz>0] <- 1
dat$present[dat$moz==0] <- 0

xtabs(present~wardType,data=dat)
length(dat$moz[dat$wardType %in% "mixed"]) # 30
length(dat$moz[dat$wardType %in% "rural"]) # 24
20/30
8/24

presAbsMod1 <- glmer(present~wardType + (1|Village),family="binomial",data=dat)
presAbsMod2 <- glmer(present~(1|Village),family="binomial",data=dat)
presAbsMod3 <- glm(present~wardType,family="binomial",data=dat)
AIC(presAbsMod1,presAbsMod2,presAbsMod3)

# take away absences and look at abundance where present
abund <- dat[dat$present==1,]
abundMod1 <- glmer.nb(moz~wardType + (1|Village),data=abund)
abundMod2 <- glmer.nb(moz~ (1|Village),data=abund)
abundMod3 <- glm.nb(moz~wardType,data=abund)
AIC(abundMod1,abundMod2,abundMod3)

