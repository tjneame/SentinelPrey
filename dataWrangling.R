                        #########################################/
                              #Sentinel Prey Data Wrangling
                                       #2023 TN
                                #modified from SR and TN 2022
                        #########################################/
#The code for sentinel prey written in September 2022 by SR and TN needs to be 
#updated and new models need to be run
#therefore new data wrangling needs to happen

#Set Up-------------------------------------------
#Load the libraries
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
theme_set(theme_bw())

#Functions
geom2cols <- function(d,x=lon,y=lat,removeGeom=TRUE,epsg=NA){
  require(sf); require(dplyr)
  if(!('sf' %in% class(d))) stop('Data must be an sf object')
  if(is_grouped_df(d)){ #If data are grouped
    warning('Data are grouped. Will ungroup data before converting 
geometry column.')
    d <- ungroup(d)
  }
  if(!is.na(epsg)) d <- st_transform(d,epsg) #Transform to new CRS
  d <- d %>% 
    mutate({{x}}:=st_coordinates(.)[,1],{{y}}:=st_coordinates(.)[,2]) #Make 
  #new columns from coordinates
  if(removeGeom) d <- st_drop_geometry(d) #Drop geometry
  return(d)
}

#Load things ------------------------------------

setwd("~/Documents/School/MastersData/SentinelPrey") #Set working directory
# setwd("/media/rsamuel/ROBINSON8GB/SentinelPrey") 
#setwd("~/Documents/Stats for Tobyn 2022/SentinelPrey")

#Read in shapefiles
fieldWhole <- read_sf("FieldWhole2.shp") %>% st_transform(3402)
grass <- read_sf("Grass2.shp") %>% st_transform(3402)
infra <- read_sf("Infrastructure2.shp") %>% st_transform(3402)
trees <- read_sf("Trees2.shp") %>% st_transform(3402)
wet <- read_sf("Wetlands2.shp") %>% st_transform(3402)

#Read in csvs
sampleLocs <- read.csv('SampleSiteLocations.csv',stringsAsFactors = FALSE,strip.white = TRUE) %>% 
  mutate(lat_dup=lat,lon_dup=lon) %>%
  st_as_sf(coords=c('lon','lat')) %>% #Add spatial feature info
  st_set_crs(4326) %>%
  st_transform(3402) %>%
  geom2cols(.,x=lon_dup, y=lat_dup, removeGeom = FALSE, epsg=3402)

sentPrey <- read.csv('sentinelPreyCaterpillars.csv',stringsAsFactors = FALSE,strip.white = TRUE)

ACISDat <- read.csv('ACISDailyData-20210601-20220831-PID131006430.csv')

# #Looks OK
# ggplot(fieldWhole)+geom_sf()
# ggplot(grass)+geom_sf()
# ggplot(infra)+geom_sf()
# ggplot(trees)+geom_sf()
# 
# #Take a look at BLID 41121
# ggplot()+
#   geom_sf(data=filter(fieldWhole,BLID==41121))+ #Field
#   geom_sf(data=filter(grass,BLID==41121),fill='green')+ #Grass
#   geom_sf(data=filter(trees,BLID==41121),fill='darkgreen')+ #Trees
#   geom_sf(data=filter(sampleLocs,BLID==41121))

#Distance from each point to boundary lines of grass/tree polygons--------------

# #Quicker example at a single field
# treeLine <- st_cast(filter(trees,BLID==41021),'MULTILINESTRING')
# grassLine <- st_cast(filter(grass,BLID==41021),'MULTILINESTRING')
# edgeLine <- st_union(treeLine,grassLine)
# 
# ggplot()+
#   geom_sf(data=filter(fieldWhole,BLID==41021))+ #Field
#   geom_sf(data=filter(grass,BLID==41021),fill='green')+ #Grass
#   geom_sf(data=filter(trees,BLID==41021),fill='darkgreen')+ #Trees
#   geom_sf(data=filter(sampleLocs,BLID==41021))+
#   geom_sf(data=edgeLine,col='red')
# 
# filter(sampleLocs,BLID==41021) %>% 
#   mutate(dist=apply(st_distance(filter(sampleLocs,BLID==41021),edgeLine),1,min))

#Turns polygons into lines (multi-line strings)
treeLine <- st_cast(trees,'MULTILINESTRING') %>% st_geometry()
grassLine <- st_cast(grass,'MULTILINESTRING')%>% st_geometry()
wetLine <- st_cast(wet,'MULTILINESTRING')%>% st_geometry()
infraLine <- st_cast(infra,'MULTILINESTRING')%>% st_geometry()

#Gets distances from noncrop lines
sampleLocs$dist <- sapply(list(treeLine,grassLine,wetLine,infraLine),function(x){
  apply(st_distance(sampleLocs,x),1,min)}) %>% 
  apply(.,1,min)

#Conglomerate data and clean----------------------------------------------------
#Merges sample distances with prey records
sampleLocs <- sampleLocs %>% unite(col = trapID, c(BLID,stationID),sep='-',remove = FALSE) #Creates a "Trap ID" column
sentPrey <- sentPrey %>% unite(col = trapID, c(BLID,site),sep='-',remove = FALSE) #Creates a "Trap ID" column

#Merges the tables together by "trapID" 
sentPreyDist <- st_drop_geometry(sampleLocs) %>%
  dplyr::select(trapID, lon_dup, lat_dup, dist) %>% 
  full_join(sentPrey,by = 'trapID')

#Merges prey records and sample distances with GDD data
sentPreyDist <- sentPreyDist %>% 
  unite(col = stnDate, c(weatherStn,midDate),sep='-',remove = FALSE) #Creates a "station date" column called stnDate
ACISDat <- ACISDat %>% 
  unite(col=stnDate, c(Station_Name,Date), sep = "-", remove = FALSE) #Creates a "station date" column called stnDate
sentPreyDistGDD <- left_join(sentPreyDist, ACISDat, by = 'stnDate') #Merges the tables together by "stnDate"

#turns dataframe in to a Tibble
sentPreyDistGDD <- as_tibble(sentPreyDistGDD)
is_tibble(sentPreyDistGDD)

#remove columns that are duplicates or don't mean anything
sentPreyF <- dplyr::select(sentPreyDistGDD, -GDDSourceFlag, 
                           -Notes, -Station_Name, -Date, 
                           -notRecovered, -peircingInsect, -arachnid,
                           -aves, -non.rodentMammal, -rodent, -weatherStn) %>%
  rename(pass=return) %>%
  drop_na()

#add abundance data!

abun<-read.csv("/Users/tobynneame/Documents/School/MastersData/beetleDataAnalysis/beetleDataAbundance.csv")%>%
  dplyr::select(trapPassID,beetCount)
sentPreyF<-sentPreyF %>% unite(col = trapPassID, c(BLID,site,pass),sep='-',remove = FALSE) #make trapPassID column
sentPreyF<-left_join(sentPreyF,abun, by="trapPassID")

write.csv(sentPreyF, "sentinelPreyData.csv")

