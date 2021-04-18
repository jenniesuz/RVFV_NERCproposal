# Simulations to assess numbers of traps required to detect differences 
# in mosquito abundance between rural 
# and peri-urban areas and between grassland and cultivated crop habitats

library("ggplot2")
library("GLMMmisc")
library("parallel")


# what would a biologically significant difference in mosquito abundance?
# how does R0 change with mosquito: host ratio - if assume a trap catch relates

# R0 that keeps moz:host ratio constant, but varies proportion of bm on humans, x= rho.human
R0<-function(x,num.human=50,num.ruminant=900,rho.ruminant,human.m.trans=0,bite=1/3,ruminant.m.trans,m.h.trans=0.05,EIP=14,gamma=1/5,surv=0.9){ 
  rho.human<-1-rho.ruminant
  R.0.human<-(x/num.human)*(((rho.human*bite)^2)*m.h.trans*human.m.trans*(surv^EIP)/(-(log(surv))*gamma))
  R.0.ruminant<-(x/num.ruminant)*(((rho.ruminant*bite)^2)*m.h.trans*ruminant.m.trans*(surv^EIP)/(-(log(surv))*gamma))
  R.0<-sum(R.0.ruminant,R.0.human)
  #prop.inf<-1-(1/R.0)
  return(R.0)
}

num.moz.t<-seq(from=10000,to=2000000,by=100) # from=1000000,to=10000000,by=10000) to get above 1
moz2hostt<-num.moz.t/1200

x<-seq(from=0.8,to=1,0.01)

R01<-sapply(num.moz.t,function(x){R0(x,rho.ruminant=0.9,bite=1/3,EIP=14,ruminant.m.trans=0.5)})

R02<-sapply(num.moz.t,function(x){R0(x,rho.ruminant=0.8,bite=1/3,EIP=14,ruminant.m.trans=0.5)})

R03<-sapply(num.moz.t,function(x){R0(x,rho.ruminant=0.7,bite=1/3,EIP=14,ruminant.m.trans=0.5)})

R04<-sapply(num.moz.t,function(x){R0(x,rho.ruminant=0.6,bite=1/3,EIP=14,ruminant.m.trans=0.5)})

R05<-sapply(num.moz.t,function(x){R0(x,rho.ruminant=0.5,bite=1/3,EIP=14,ruminant.m.trans=0.5)})



dat <- cbind.data.frame(moz2host=rep(moz2hostt,3)
                        ,Ro=c(R01,R02,R03)
                        ,rho=c(rep("1",length(R01))
                               ,rep("5",length(R02))
                               ,rep("10",length(R03))
                               
                        )
)

dat$rho <- as.character(dat$rho)
dat$rho <- factor(dat$rho, levels=c("1","5","10"))

# with very relaxed assumptions (e.g. short bite interval short EIP high moz:hot ratio)
ggplot(dat) +
  geom_line(aes(x=moz2host,y=Ro,linetype=rho)) +
  scale_x_continuous(trans='log',breaks=c(10,20,50,100,200,500,1000),limits=c(10,1000)) +
  #ylim(0,20) +
  geom_hline(yintercept=1) +
  labs(x="Mosquitoes per host"
       ,y="R0"
       ,linetype="% of bloodmeals \n on dead-end hosts") +
  theme_set(theme_bw()) +
  theme(
    axis.line = element_line(color = 'black')
    ,text=element_text(size=12)
    ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
    ,axis.text=element_text(size=12)
    ,legend.key.size = unit(0.8,"line")
    ,legend.background = element_blank()
    ,legend.text=element_text(size=12)
    ,legend.position =c(0.3,0.8)
    #,legend.title = element_blank()
    ,strip.background = element_rect(colour="white", fill="white")
    ,panel.border = element_blank()
  )

#************************Simulations to determine number of traps************
pwrFunc <- function(...){
  library("GLMMmisc")
days <- 1:3
traps <- 1:10
pr <- c("periurban","rural")
habitat <- c("grassland","crop")
area <- LETTERS[1:4]

moz.data <-
expand.grid(traps=traps,habitat=habitat,pr=pr,area=area,days=days)

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
              rural = 0.5)),
        habitat=log(
          c(grassland=1,
            crop=2)
        )),
    rand.V =
      inv.mor(
        c(days = 1.8,
          area = 1.8)),            
    distribution = "poisson")      

moz.pois <-
    lme4::glmer(response ~ pr + habitat + (1 | area) + (1| days) ,
                family = "poisson", data = moz.data)
output <- summary(moz.pois)
cfs<-output$coefficients
prrural<-cfs[2,4]
prhabitat<-cfs[3,4]
pvalrural <- 0
pvalhabitat <- 0

if(length(output$optinfo$conv$lme4$messages)==0){
if(prrural<=0.05){pvalrural<-1}
if(prhabitat<=0.05){pvalhabitat<-1}
}
return(c(pvalrural,pvalhabitat))
}


sim.res <- mclapply(1:1000, pwrFunc, mc.cores =1)

l<-do.call(rbind,sim.res)
l2<-rowSums(l)
length(l2[l2==2])


