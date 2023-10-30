################################################################################/
#Sentinel Prey Data Analysis and Visualisation
#2023 TN
#modified from SR and TN 2022
################################################################################/

#Set Up-------------------------------------------------------------------------
library(tidyverse)
library(gamlss)
library(gratia)
library(mgcv)
library(googlesheets4)
library(DHARMa)

#load data
sentPrey<-read.csv("sentinelPreyData.csv")
sentPrey<-read_sheet("https://docs.google.com/spreadsheets/d/1nPHKuQKl2J0DRYot5sx8afq6LNKcX68v7SOxvBCtRac/edit?usp=sharing")
write_csv(sentPrey, "sentinelPreyData.csv")
#center the lat/lon
sentPrey<-sentPrey %>%
  group_by(BLID) %>% mutate(cLon=mean(lon_dup),cLat=mean(lat_dup)) %>%
  ungroup()%>% mutate(lon_dup=lon_dup-cLon,lat_dup=lat_dup-cLat)

#make sure vectors are represented properly
sentPrey <- sentPrey %>%
  mutate(BLID=as.factor(BLID)) %>%
  mutate(year=as.factor(year)) %>%
  mutate(type=as.factor(type))

#make NX give negative distances
sentPrey <- sentPrey %>% 
  mutate(dist = case_when(site == "N1" ~ -dist,
                          site == "N2" ~ -dist,
                          site == "N3" ~ -dist,
                          site == "N3" ~ -dist,
                          site == "N4" ~ -dist,
                          site == "N5" ~ -dist,
                          site == "N6" ~ -dist,
                          TRUE ~ dist))
#remove columns that wont be needed at all
sentPrey<-sentPrey %>%
  dplyr::select(-...1, #sometimes -X instead
                -trapID,
                -BTID,
                -trapPassID,
                -stnDate,
                -midDate)

#remove some rows of sentinel prey that were raised off the ground
sentPrey<-sentPrey %>% filter(type=="GD") %>%
  droplevels()


#make some new sub-sets of data ----------------------------------------------
#only within crop data
sentPreyCrop<-filter(sentPrey, dist>=0)

#crop vs non crop data
sentPreyN<-sentPrey %>% 
  filter(site != "C1",
         site != "C2",
         site != "C3",
         site != "C3",
         site != "C4",
         site != "C5",
         site != "C6",
         site != "C7",
         site != "C8",
         site != "C9") %>%
  mutate(site="nonCrop")

sentPreyC<-sentPrey %>%
  filter(site != "N1",
         site != "N2",
         site != "N3",
         site != "N3",
         site != "N4",
         site != "N5",
         site != "N6") %>%
  mutate(site="Crop")

sentPreyNC<-bind_rows(sentPreyC,sentPreyN) %>%
  dplyr::select(-dist)%>%
  mutate(site=as.factor(site))%>%
  drop_na()

rm(sentPreyC)
rm(sentPreyN)


#sentinelPrey by distance ----------------------------------------

#Base 
formSP1 <- as.formula(chewingInsect ~ s(dist) + #Distance from edge
                        s(GDD) + #Growing degree day
                        ti(dist,GDD) + #Distance:GDD interaction
                        year + #Year 
                        s(lon_dup,lat_dup,by=BLID) + #Within-field distances
                        s(BLID,bs='re')) #Between-field

sentPreyDistGAM1 <- gam(formula = formSP1,
                        family = nb,
                        data=sentPreyCrop)
write_rds(sentPreyDistGAM1, "sentPreyDistGAM1.rds")

#With abundance (beetCount)
formSP2 <- as.formula(chewingInsect ~ s(dist, bs="ts") + #Distance from edge
                        s(GDD, bs="ts") + #Growing degree day
                        ti(dist,GDD, bs="ts") + #Distance:GDD interaction
                        year + #Year 
                        s(lon_dup,lat_dup,by=BLID, bs="ts") + #Within-field distances
                        s(beetCount, bs="ts") + #abundance of beetles
                        ti(beetCount, dist, bs="ts") + #Distance:abundance interaction
                        s(BLID,bs='re')) #Between-field

sentPreyDistGAM2 <- gam(formula = formSP2,
                        family = nb,
                        data=sentPreyCrop)
write_rds(sentPreyDistGAM2, "sentPreyDistGAM_final.rds")

#sentinelPrey by crop vs non-crop---------------------------------------------

#Base
formSP3 <- as.formula(chewingInsect ~ site*GDD + #site:GDD interaction and main
                        year + #Year 
                        s(lon_dup,lat_dup,by=BLID) + #Within-field distances
                        s(BLID,bs='re')) #Between-field 

sentPreyNCGAM1 <- gam(formula = formSP3,
                        family = nb,
                        data=sentPreyNC)
write_rds(sentPreyNCGAM1, "sentPreyNCGAM1.rds")

#Smooth abundance (beetCount)
formSP4 <- as.formula(chewingInsect ~ site + #crop or non-crop
                        s(GDD, bs="ts") + #growing degree day
                        s(GDD, by=site, bs="ts") + #crop non crop:GDD interaction
                        year + #Year 
                        s(lon_dup,lat_dup,by=BLID, bs="ts") + #Within-field distances
                        s(beetCount, bs="ts") + #abundance
                        s(beetCount, by=site, bs="ts") + #crop non crop:abundance interaction
                        s(BLID,bs='re')) #Between-field

sentPreyNCGAM2 <- gam(formula = formSP4,
                      family = nb,
                      data=sentPreyNC)
write_rds(sentPreyNCGAM2, "sentPreyNCGAM_final.rds")


#visualize sentinelPrey by distance ------------------------------------------
m1gam<-read_rds("sentPreyDistGAM_final.rds")
#test the model
simulateResiduals(m1gam, plot=T)

summary(m1gam)
newdat <- expand.grid(dist=seq(0,200,by=5),GDD=c(300,500,700),
                      beetCount=10, year='2021',BLID='41007',lon_dup=0,lat_dup=0) 
newdat <- predict.gam(m1gam,newdata=newdat,se.fit = TRUE,
                      exclude = c('s(BLID)',paste0('s(lon_dup,lat_dup):BLID',levels(sentPreyCrop$BLID)))) %>% 
  do.call('data.frame',.) %>% 
  mutate(upr=fit+se.fit*1.96,lwr=fit-se.fit*1.96) %>% 
  mutate(across(c(fit,upr,lwr),exp)) %>% 
  bind_cols(dplyr::select(newdat,dist,GDD),.)

cols <- c('blue','purple','red')

(p <- newdat %>% mutate(GDD=factor(GDD,labels = c('Early (300)','Mid (500)','Late (700)'))) %>% 
    ggplot(aes(x=dist))+geom_ribbon(aes(ymax=upr,ymin=lwr,fill=GDD),alpha=0.2)+
    geom_line(aes(y=fit,col=GDD))+
    geom_text(data=sentPreyCrop,aes(x=dist,y=0.05),label='|',position=position_jitter(width = 1, height=0),alpha=0.4,size=2)+
    labs(x='Distance from nearest non-crop vegetation area',y='Bite marks per caterpillar')+
    xlim(0,200)+
    scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
    theme(legend.position = c(0.85,0.85),legend.background = element_rect(colour='grey'))
)

ggsave('./figures/biteMarks.png',p,width = 10,height=6)
ggsave('./figures/biteMarks2.png',p,width = 6,height=6)

(p <- newdat %>% mutate(GDD=factor(GDD,labels = c('Early (300)','Mid (500)','Late (700)'))) %>% 
    ggplot()+geom_ribbon(aes(x=dist,ymax=upr,ymin=lwr),alpha=0.2)+
    geom_line(aes(x=dist,y=fit))+
    geom_text(data=dplyr::select(sentPreyCrop,dist),aes(x=dist,y=0.05),label='|',position=position_jitter(width = 1, height=0),alpha=0.4,size=2)+
    facet_wrap(~GDD)+
    labs(x='Distance from nearest non-crop vegetation area (m)',y='Bite marks per sentinel')+
    xlim(0,200)+
    scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
    theme_bw()
)
p
ggsave('./figures/biteMarks1.png',p,width = 10,height=6, scale=0.7)
ggsave('./figures/biteMarks3.svg',p,width = 10,height=6)

#visualize the interaction effect of beetle abundance and dist on bite marks
p<-draw(smooth_estimates(m1gam, "ti(beetCount,dist)"))+
  labs(y="Distance from nearest non-crop vegetation area (m)", 
       x="Number of beetles",
       title = "")+
  ylim(0.4,200)+ #clean up the plot and only use reliable data
  theme_minimal()+
  coord_flip()
p
ggsave('./figures/biteMarksByAbund.png',p,width = 10,height=6, scale=0.7)

# Visualize sentinelPrey by crop vs. non-crop ----------------------------------
m2gam<-read_rds("sentPreyNCGAM_final.rds")
summary(m2gam)
#test the model
simulateResiduals(m2gam, plot=T)

newdat <- expand.grid(site=c('nonCrop', 'Crop'), GDD=c(300,500,700),
                      beetCount=10, year='2021',BLID='41007',lon_dup=0,lat_dup=0) 
newdat <- predict.gam(m2gam,newdata=newdat,se.fit = TRUE,
                      exclude = c('s(BLID)',paste0('s(lon_dup,lat_dup):BLID',levels(sentPreyNC$BLID)))) %>% 
  do.call('data.frame',.) %>% 
  mutate(upr=fit+se.fit*1.96,lwr=fit-se.fit*1.96) %>% 
  mutate(across(c(fit,upr,lwr),exp)) %>% 
  bind_cols(dplyr::select(newdat,site,GDD),.)

cols <- c('blue','purple','red')

(p <- newdat %>% mutate(GDD=factor(GDD,labels = c('Early (300)','Mid (500)','Late (700)'))) %>% 
    ggplot(aes(x=site))+
    geom_errorbar(aes(ymax=upr,ymin=lwr,fill=GDD),alpha=0.2)+
    geom_point(aes(y=fit,col=GDD))+
    labs(x='Non Crop vs Crop',y='Bite marks per caterpillar')+
    scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
    theme(legend.position = c(0.85,0.85),legend.background = element_rect(colour='grey'))
)

ggsave('./figures/biteMarks.png',p,width = 10,height=6)
ggsave('./figures/biteMarks2.png',p,width = 6,height=6)

(p <- newdat %>% mutate(GDD=factor(GDD,labels = c('Early (300)','Mid (500)','Late (700)'))) %>% 
    ggplot()+
    geom_errorbar(aes(x=site,ymax=upr,ymin=lwr),alpha=0.9, width=0.1, position = position_dodge(0.9), linewidth=1)+
    geom_point(aes(x=site,y=fit), pch=19, size=4)+
    facet_wrap(~GDD)+
    labs(x='Location of Sampling',y='Bite marks per sentinel')+
    scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
    theme_bw()
)

ggsave('./figures/biteMarks2.png',p,width = 10,height=6, scale=0.7)
ggsave('./figures/biteMarks3.svg',p,width = 10,height=6)



# New code (Oct 30) ----------------------------------------------------------------


sentPreyAb2 <- sentPrey %>% 
  # mutate(isNonCrop=grepl('N',as.character(site))) %>%
  mutate(isNonCrop=factor(ifelse(grepl('N',site),'NonCrop','Crop'))) %>% 
  mutate(cDist = ifelse(grepl('N',site),0,dist)) %>% 
  na.omit()

isComplete <- complete.cases(sentPreyAb2)

mForm <- as.formula(chewingInsect ~ 
                      #"Main effects"
                      s(beetCount) + #Count of beetles
                      s(cDist) + #Distance from edge
                      isNonCrop + #Difference b/w crop and non crop (at dist=0)
                      s(GDD,by=isNonCrop) + #median growing degree day 
                      s(lon_dup,lat_dup,by=BLID) + #Within-field distances
                      #"Interactions"
                      ti(GDD,cDist) + #GDD:Distance
                      ti(beetCount,cDist)+ #Count and distance
                      s(BLID,bs='re')) #Between-field re

library(parallel)
cl <- makeCluster(5)
gamSent <- bam(formula = mForm,family='nb',data=sentPreyAb2,cluster = cl,samfrac=0.1)
stopCluster(cl)

summary(gamSent)
plot(gamSent)

par(mfrow=c(2,2));gam.check(gamN);abline(0,1,col='red');par(mfrow=c(1,1))

library(DHARMa)
gamSent_simRes <- simulateResiduals(gamSent,n=1000)
plot(gamSent_simRes)
testZeroInflation(gamSent_simRes) 
testOutliers(gamSent_simRes, type = "binomial")
testCategorical(gamSent_simRes, sentPreyAb2$isNonCrop[isComplete])
testCategorical(gamSent_simRes, sentPreyAb2$BLID[isComplete])

#Check for site-to-site spatial AC
reCoefs <- gamSent$coefficients[grepl('s(BLID).',names(gamSent$coefficients),fixed = TRUE)]

sentPreyAb2 %>% dplyr::select(BLID,cLon,cLat) %>% distinct() %>% 
  mutate(ranef=reCoefs) %>% 
  ggplot(aes(x=cLon,y=cLat))+geom_point(aes(size=abs(ranef),col=ranef))+
  scale_color_gradient2(low='red',high='blue')

#Calculate Moran's I
blid_weight <- sentPreyAb2 %>% dplyr::select(cLon,cLat) %>% distinct() %>% 
  dist() %>% as.matrix() 
blid_weight <- 1/blid_weight
diag(blid_weight) <- 0
ape::Moran.I(reCoefs,blid_weight) #No problem here (cite in paper)

#Check for temporal AC
sentPreyAb2 %>% na.omit() %>% mutate(res=residuals(gamSent,type='deviance')) %>% 
  ggplot(aes(x=GDD,y=res))+geom_point()+
  geom_smooth(method='loess',se=TRUE)

sentPreyAb2 %>% na.omit() %>%  mutate(res=residuals(gamSent,type='deviance')) %>% 
  ggplot(aes(x=GDD,y=res))+geom_point()+
  geom_smooth(method='loess',se=TRUE)+
  facet_wrap(~BLID)

#Make some partial effects plots
library(ggeffects)
theme_set(theme_bw())

ggpredict(gamSent,terms=c('cDist','isNonCrop')) %>% 
  data.frame() %>% 
  filter(x==0|group=='Crop')

#Predictions for noncrop
dontUse <- c('BLID','lon_dup','lat_dup')
useGDD <- c(300,500,700)
gddLabs <- c('Early (300)','Mid (500)','Late (700)')
predDat1 <- data.frame(cDist=0,isNonCrop=sentPreyAb2$isNonCrop[1], #Noncrop
                       GDD=useGDD,lon_dup=0,lat_dup=0,beetCount=14,
                       BLID=sentPreyAb2$BLID[1]) %>% 
  bind_cols(.,predict(gamSent,newdata=.,se.fit=TRUE,exclude=dontUse)) %>% 
  mutate(upr=fit+se.fit*1.96,lwr=fit-se.fit*1.96) %>% 
  mutate(across(c(fit,upr,lwr),~exp(.x)/2)) %>% 
  mutate(cDist=-5) %>% mutate(GDD=factor(GDD,labels=gddLabs))

predDat2 <- data.frame(cDist=0:200,isNonCrop=sentPreyAb2$isNonCrop[605], #Crop
                       GDD=useGDD,lon_dup=0,lat_dup=0,beetCount=14,
                       BLID=sentPreyAb2$BLID[1]) %>% 
  bind_cols(.,predict(gamSent,newdata=.,se.fit=TRUE,exclude=dontUse)) %>% 
  mutate(upr=fit+se.fit*1.96,lwr=fit-se.fit*1.96) %>% 
  mutate(across(c(fit,upr,lwr),~exp(.x)/2)) %>% 
  mutate(GDD=factor(GDD,labels=gddLabs))

ggplot(data=predDat2,aes(x=cDist,y=fit))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.2)+
  geom_line(linewidth=1)+
  geom_pointrange(data=predDat1,aes(ymax=upr,ymin=lwr),size=1)+
  facet_wrap(~GDD,ncol=3)+
  labs(x='Distance',y='Bites per week')

ggpredict(gamSent,terms=c('cDist[0:200]','beetCount[0,100]')) %>% 
  data.frame() %>% mutate(group=factor(group,labels=c('Low (0)','High (50)'))) %>% 
  ggplot(aes(x=x,y=predicted))+
  geom_ribbon(aes(ymax=conf.high,ymin=conf.low,fill=group),alpha=0.2)+
  geom_line(aes(col=group),linewidth=2)+
  scale_colour_manual(values=c('red','blue'),aesthetics = c('colour','fill'))+
  labs(x='Distance',y='Bites per week',colour='Beetles\nper\nweek',
       fill='Beetles\nper\nweek')+
  theme(legend.position = c(0.8,0.8),legend.background = element_rect(color='black'))

save(list = c('sentPrey','gamSent'),file = './sentData_SR.Rdata')
