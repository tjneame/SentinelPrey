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

#remove some rows that have data input issues
sentPrey<-sentPrey %>% filter(type=="GD"|type=="RA") %>%
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
                        type+ #RA or GD
                        s(BLID,bs='re')) #Between-field

sentPreyDistGAM1 <- gam(formula = formSP1,
                        family = nb,
                        data=sentPreyCrop)
write_rds(sentPreyDistGAM1, "sentPreyDistGAM1.rds")

#Smooth abundance (beetCount)
formSP2 <- as.formula(chewingInsect ~ s(dist) + #Distance from edge
                        s(GDD) + #Growing degree day
                        ti(dist,GDD) + #Distance:GDD interaction
                        year + #Year 
                        s(lon_dup,lat_dup,by=BLID) + #Within-field distances
                        type+ #RA or GD
                        s(beetCount)+ #abundance
                        s(BLID,bs='re')) #Between-field

sentPreyDistGAM2 <- gam(formula = formSP2,
                        family = nb,
                        data=sentPreyCrop)
write_rds(sentPreyDistGAM2, "sentPreyDistGAM2.rds")

#sentinelPrey by crop vs non-crop---------------------------------------------

#Base
formSP3 <- as.formula(chewingInsect ~ site*GDD + #site:GDD interaction and main
                        year + #Year 
                        type + #GD or RA
                        s(lon_dup,lat_dup,by=BLID) + #Within-field distances
                        s(BLID,bs='re')) #Between-field 

sentPreyNCGAM1 <- gam(formula = formSP3,
                        family = nb,
                        data=sentPreyNC)
write_rds(sentPreyNCGAM1, "sentPreyNCGAM1.rds")

#Smooth abundance (beetCount)
formSP4 <- as.formula(chewingInsect ~ site*GDD + #Distance:GDD interaction and main
                        year + #Year 
                        s(lon_dup,lat_dup,by=BLID) + #Within-field distances
                        type + #GD or RA
                        s(beetCount) + #abundance
                        s(BLID,bs='re')) #Between-field

sentPreyNCGAM2 <- gam(formula = formSP4,
                      family = nb,
                      data=sentPreyNC)
write_rds(sentPreyNCGAM2, "sentPreyNCGAM2.rds")


#visualize sentinelPrey by distance ------------------------------------------
m1gam<-read_rds("sentPreyDistGAM2.rds")

newdat <- expand.grid(dist=seq(0,200,by=5),GDD=c(300,500,700), type='RA',
                      beetCount=0, year='2021',BLID='41007',lon_dup=0,lat_dup=0) 
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
    labs(x='Distance from field edge',y='Bite marks per caterpillar')+
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
    labs(x='Distance from field edge',y='Bite marks per caterpillar')+
    xlim(0,200)+
    scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
    theme(legend.position = c(0.15,0.85),legend.background = element_rect(colour='grey'))
)

ggsave('./figures/biteMarks3.png',p,width = 10,height=6)
ggsave('./figures/biteMarks3.svg',p,width = 10,height=6)

# Visualize sentinelPrey by crop vs. non-crop
m2gam<-read_rds("sentPreyNCGAM2.rds")
newdat <- expand.grid(site=c('Crop', 'nonCrop'), GDD=c(300,500,700), type='RA',
                      beetCount=0, year='2021',BLID='41007',lon_dup=0,lat_dup=0) 
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
    labs(x='Crop vs Non Crop',y='Bite marks per caterpillar')+
    scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
    theme(legend.position = c(0.85,0.85),legend.background = element_rect(colour='grey'))
)

ggsave('./figures/biteMarks.png',p,width = 10,height=6)
ggsave('./figures/biteMarks2.png',p,width = 6,height=6)

(p <- newdat %>% mutate(GDD=factor(GDD,labels = c('Early (300)','Mid (500)','Late (700)'))) %>% 
    ggplot()+
    geom_errorbar(aes(x=site,ymax=upr,ymin=lwr),alpha=0.2)+
    geom_point(aes(x=site,y=fit))+
    facet_wrap(~GDD)+
    labs(x='Crop vs Non Crop',y='Bite marks per caterpillar')+
    scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
    theme(legend.position = c(0.15,0.85),legend.background = element_rect(colour='grey'))
)

ggsave('./figures/biteMarks3.png',p,width = 10,height=6)
ggsave('./figures/biteMarks3.svg',p,width = 10,height=6)

