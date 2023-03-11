# Simulations to assess numbers of traps required to detect differences 
# in mosquito abundance between rural 
# and peri-urban areas 

library("ggplot2")
library("GLMMmisc")
library("parallel")

#*****************simulate data***************
simAbundanceDatFunc <- function(traps=1:6               # how many traps per day per landtype
                                #,days=1                 # how many days per area per landtype
                                ,area=LETTERS[1:12]    # how many 'areas' - village 
                                ,daysRand=1.8           # traps random effect
                                ,areaRand=4             # area random effect
                                ,interc=5){            # area random effect

  pr <- c("periurban","rural") # removed grassland rural

  moz.data <- expand.grid(traps=traps,pr=pr,area=area)#,days=days)
   #lengthMD <- length(moz.data[,1])
   #numReps <- lengthMD/max(traps)
  #  days <- sapply(1:numReps,function(x){
  #    return(c(rep(x,max(traps))))
  #  })
   # moz.data$days <- c(days)
    
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
               c(#days = daysRand, # 1.8
                area = areaRand)), #1.8             
         distribution = "negbinomial",
      theta=0.5)      
 return(moz.data)
}



dat1 <- simAbundanceDatFunc(areaRand=1.5)

hist(dat1$response)

ggplot(dat1) +
  geom_boxplot(aes(x=area,y=response))

ggplot(dat1) +
  geom_boxplot(aes(x=pr,y=response))



#*********************************************************************









#************************Simulations to determine number of traps************
pwrFunc <- function(...){
  library("GLMMmisc")
#days <- 1
traps <- 1:8
pr <- c("periurban","rural") # removed grassland rural
area <- LETTERS[1:12]

moz.data <-
expand.grid(traps=traps,pr=pr,area=area)#,days=days)
#lengthMD <- length(moz.data[,1])
#numReps <- lengthMD/max(traps)
#days <- sapply(1:numReps,function(x){
#  return(c(rep(x,max(traps))))
#})
#moz.data$days <- c(days)

moz.data<-
  sim.glmm(
    design.data = moz.data,
    fixed.eff =
      list(
        intercept = log(5),      
        pr =
          log(                     
            c(periurban = 1,               
              rural = 0.5))
        ),
    rand.V =
      inv.mor(
        c(#days = 1.3,
          area = 1.5)),            
    distribution = "negbinomial",
    theta=0.5)      

moz.nb <-
    lme4::glmer.nb(response ~ pr + (1 | area) 
                   , data = moz.data
                  ,control = glmerControl(optimizer ="Nelder_Mead"))
output <- summary(moz.nb)
cfs<-output$coefficients
prrural<-cfs[2,4]
pvalrural <- 0


if(length(output$optinfo$conv$lme4$messages)==0){
if(prrural<=0.05){pvalrural<-1}
}
return(c(pvalrural))
}


sim.res <- mclapply(1:100, pwrFunc, mc.cores =1)

l<-do.call(rbind,sim.res)
l2<-rowSums(l)
length(l2[l2==1])

# c. 85 with 10 traps
