# Simulations to assess numbers of individual blood fed mosquitoes that would be required
# to detect a difference between peri-urban and and rural areas
library("ggplot2")
library("GLMMmisc")
library("parallel")

#***************************************************************8
#*###BLOODMEALS#***********

pwrFunc <- function(...){
  library("GLMMmisc")
  
  inds <- 20
  pr <- c("periurban","rural")
  area <- LETTERS[1:4]
  numTraps<-1:20
  
  moz.data <-
    expand.grid(inds=inds,pr=pr,area=area,numTraps=numTraps)
  
  moz.data$response <- inds
  moz.data$n <- moz.data$response
  
  moz.data<-
    sim.glmm(
      design.data = moz.data,
      fixed.eff =
        list(
          intercept = qlogis(0.50),      
          pr =
            log(                     
              c(rural= 1,               
                periurban = 1.5))
        ),
      rand.V =
        inv.mor(
          c(area = 1.3
            ,numTraps=1.3)),            
      distribution = "binomial")      
  
  moz.bin <-
    lme4::glmer(cbind(response, n - response) ~ pr + (1 | area) + (1 | numTraps) ,
                family = "binomial", data = moz.data)
  output <- summary(moz.bin)
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
sum(l2)
