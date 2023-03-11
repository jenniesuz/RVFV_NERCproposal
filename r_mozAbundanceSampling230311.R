# Simulations to assess numbers of traps required to detect differences 
# in mosquito abundance between rural 
# and peri-urban areas 

library("ggplot2")
library("GLMMmisc")
library("parallel")

#*****************simulate data***************
simAbundanceDatFunc <- function(traps=1:24               # how many traps per day per landtype
                                ,days=1:3               # how many days per area per landtype
                                ,area=LETTERS[1:4]      # how many areas 
                                ,daysRand=1.8           # days random effect
                                ,areaRand=1.8           # area random effect
                                ,interc=50){          # area random effect

  pr <- c("periurban","rural") # removed grassland rural

  moz.data <- expand.grid(traps=traps,pr=pr,area=area,days=days)
   lengthMD <- length(moz.data[,1])
   numReps <- lengthMD/max(traps)
    days <- sapply(1:numReps,function(x){
      return(c(rep(x,max(traps))))
    })
    moz.data$days <- c(days)
    
  moz.data <- sim.glmm(
      design.data = moz.data,
      fixed.eff =
        list(
          intercept = log(interc),      
          pr =
            log(c(periurban = 1,               
              rural = 0.5))
              ),
         rand.V =
            inv.mor(
               c(days = daysRand, # 1.8
                area = areaRand)), #1.8             
         distribution = "poisson")      
 return(moz.data)
}



dat1 <- simAbundanceDatFunc()

ggplot(dat1) +
  geom_boxplot(aes(x=pr,y=log(response)))


#*********************************************************************









#************************Simulations to determine number of traps************
pwrFunc <- function(...){
  library("GLMMmisc")
days <- 1:4
traps <- 1:5
pr <- c("periurban","rural") # removed grassland rural
area <- LETTERS[1:4]

moz.data <-
expand.grid(traps=traps,pr=pr,area=area,days=days)

moz.data$days[(moz.data$area %in% "B")&(moz.data$days %in% 1)] <- 4
moz.data$days[(moz.data$area %in% "B")&(moz.data$days %in% 2)] <- 5
moz.data$days[(moz.data$area %in% "B")&(moz.data$days %in% 3)] <- 6

moz.data$days[(moz.data$area %in% "C")&(moz.data$days %in% 1)] <- 7
moz.data$days[(moz.data$area %in% "C")&(moz.data$days %in% 2)] <- 8
moz.data$days[(moz.data$area %in% "C")&(moz.data$days %in% 3)] <- 9

moz.data$days[(moz.data$area %in% "D")&(moz.data$days %in% 1)] <- 10
moz.data$days[(moz.data$area %in% "D")&(moz.data$days %in% 2)] <- 11
moz.data$days[(moz.data$area %in% "D")&(moz.data$days %in% 3)] <- 12

moz.data<-
  sim.glmm(
    design.data = moz.data,
    fixed.eff =
      list(
        intercept = log(10),      
        pr =
          log(                     
            c(periurban = 1,               
              rural = 0.5))
        ),
    rand.V =
      inv.mor(
        c(days = 1.8,
          area = 1.8)),            
    distribution = "poisson")      

moz.pois <-
    lme4::glmer(response ~ pr + (1 | area) + (1| days) ,
                family = "poisson", data = moz.data)
output <- summary(moz.pois)
cfs<-output$coefficients
prrural<-cfs[2,4]
pvalrural <- 0


if(length(output$optinfo$conv$lme4$messages)==0){
if(prrural<=0.05){pvalrural<-1}
}
return(c(pvalrural))
}


sim.res <- mclapply(1:1000, pwrFunc, mc.cores =1)

l<-do.call(rbind,sim.res)
l2<-rowSums(l)
length(l2[l2==1])

# 785/1000 = 5 traps over three nights per area
# 803/1000 = 6 traps over three nights per area per land use

# 4 days 5 traps per area per land use type # 160 # 991

# 872/ 1000