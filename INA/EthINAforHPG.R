library(readxl)
library(rgeoboundaries)
library(sf)
library(INA)
library(raster)

#------------------------------------------------------------------------------

#lets map out Ethiopia potato using MapSPAM data

potatoSPAM <- raster('/home/betherton/EthiopiaINA/spam2017V1r1_SSA_gr_H_POTA_A.tif')
west_ext <- extent(33, 47, 3, 15)
#west_ext <- extent(37, 44, 8, 12)
potatoSPAM <- crop(potatoSPAM, west_ext)
eth_region <- gb_adm1("ethiopia", type = "sscu",version="3_0_0") #regional boundaries
eth_zones <- gb_adm2("ethiopia", type = "sscu",version="3_0_0") #zonal boundaries

#Lets map out the pH tiff

topsoil_pH <- raster('/home/betherton/EthiopiaINA/T_PH_H2O.tif')
topsoil_Eth <- crop(topsoil_pH, west_ext)
topsoil_pH <- resample(topsoil_Eth,potatoSPAM,method = "bilinear")
ph_values <- getValues(topsoil_pH)
ph_coords<-cbind(coordinates(topsoil_pH),ph_values);colnames(ph_coords)<-c("x","y","ph")

#lets set up our data matrix with geo pts and region/zonal boundaries

potato_values <- getValues(potatoSPAM)
coords<-cbind(coordinates(potatoSPAM),potato_values)
coords<-cbind(coords,ph_coords[,3]);colnames(coords)<-c("x","y","density","ph")
coords_geo<-coords[,c(1,2)]
pnts<-data.frame(name=1:dim(coords_geo)[1],latitude=coords_geo[,2],longitude=coords_geo[,1])
sp_pnts<-st_as_sf(pnts,coords=c('longitude','latitude'))
st_crs(sp_pnts)=4326
sp_points=st_transform(sp_pnts,crs=st_crs(eth_zones))
sp_points$local_zone<-apply(st_intersects(eth_zones,sp_points,sparse=FALSE),
                            2,function(col){eth_zones[which(col),]$shapeName})  #assign each point its zone
sp_points$region<-apply(st_intersects(eth_region,sp_points,sparse=FALSE),
                        2,function(col){eth_region[which(col),]$shapeName})  #assign each point its region
sp_points$pH<-cbind(coords[,4]);sp_points$density<-cbind(coords[,3])
sp_points<-sp_points[which(sp_points$local_zone!="character(0)"),] #filtering out bodies of water
sp_points<-sp_points[which(is.na(sp_points$density)==FALSE),] #filtering out regions with NA cropland density 
new_coords <- do.call(rbind, st_geometry(sp_points$geometry)) %>% 
  as.matrix() %>% setNames(c("lon","lat")) #a new coord set excluding the no data areas
rm(sp_pnts,pnts,coords,coords_geo,topsoil_Eth)

#----------------------------------------------------------------------------

#Bacterial Wilt Incidence, rate of new occurrence

sp_points$incidence <- ifelse(sp_points$pH>6.550966, 0, sp_points$pH*(-3.5033)+22.95)

#-----------------------------------------------------------------------------

#BW Prevalence, individuals in population 

{
  sp_points$disease_prev<-0#sample(0:100,dim(sp_points)[1],replace=TRUE)
  #probs taken from https://peerj.com/articles/14661/
  sp_points$disease_prev[which(sp_points$local_zone=="Awi/Agew")]<-25
  sp_points$disease_prev[which(sp_points$local_zone=="Arsi")]<-sample(0.8:91.6,length(which(sp_points$local_zone=="Arsi")),replace=TRUE)
  sp_points$disease_prev[which(sp_points$local_zone=="East Shewa")]<-sample(0.8:63,length(which(sp_points$local_zone=="East Shewa")),replace=TRUE)
  sp_points$disease_prev[which(sp_points$local_zone=="Gamo Gofa")]<-sample(10:80,length(which(sp_points$local_zone=="Gamo Gofa")),replace=TRUE)
  sp_points$disease_prev[which(sp_points$local_zone=="Gurage")]<-sample(0:23,length(which(sp_points$local_zone=="Gurage")),replace=TRUE)
  sp_points$disease_prev[which(sp_points$local_zone=="Hadiya")]<-sample(0:89,length(which(sp_points$local_zone=="Hadiya")),replace=TRUE)
  sp_points$disease_prev[which(sp_points$local_zone=="West Shewa")]<-sample(10.87:90,length(which(sp_points$local_zone=="West Shewa")),replace=TRUE)
  sp_points$disease_prev[which(sp_points$local_zone=="Jimma")]<-sample(0:100,length(which(sp_points$local_zone=="Jimma")),replace=TRUE)
  sp_points$disease_prev[which(sp_points$local_zone=="Keffa")]<-sample(21:78,length(which(sp_points$local_zone=="Keffa")),replace=TRUE)
  sp_points$disease_prev[which(sp_points$local_zone=="Alaba")]<-42.5
  sp_points$disease_prev[which(sp_points$local_zone=="North Gonder")]<-sample(50:100,length(which(sp_points$local_zone=="North Gonder")),replace=TRUE)
  sp_points$disease_prev[which(sp_points$local_zone=="Sidama")]<-62.5
  sp_points$disease_prev[which(sp_points$local_zone=="South Gonder")]<-20
  sp_points$disease_prev[which(sp_points$local_zone=="South West Shewa")]<-60
  sp_points$disease_prev[which(sp_points$local_zone=="West Arsi")]<-sample(25:75,length(which(sp_points$local_zone=="West Arsi")),replace=TRUE)
  sp_points$disease_prev[which(sp_points$local_zone=="West Gojam")]<-sample(66.7:100,length(which(sp_points$local_zone=="West Gojam")),replace=TRUE)
  sp_points$disease_prev[which(sp_points$local_zone=="West Shewa")]<-sample(1.5:82.5,length(which(sp_points$local_zone=="West Shewa")),replace=TRUE)
} #using the data from peerj, assign a relative intital disease prevalence across each pixel
sp_points$incidence[which(is.na(sp_points$incidence)==TRUE)]<-0 #if no pH data is present, assume no incidence 
#if the pH is high enough, then BW can be present at a node as a functions of it's disease prevalence 
for(i in 1:dim(sp_points)[1]){
  if(sp_points$incidence[i]!=0){
    sp_points$init_loc[i]<-sample(c(0,1),1,prob=c(((100-sp_points$disease_prev[i])/100),(sp_points$disease_prev[i]/100)))
  }
  else{
    sp_points$init_loc[i]<-0
  }
  print(i)
}

#---------------------------------------------------------------------------

#prob of disease estab = disease incidence / no of suceptible individuals 
#roughly 6873 potatoes per pixel, from FAO, and given each pixel is roughly
#0.09818 ha in size

#sp_points$estab<-ifelse(sp_points$density==0,0,sp_points$incidence/((sp_points$density/100)*6873))
#if there is incidence or the pH is high enough, then disease establishment is possible
sp_points$estab<-ifelse(sp_points$incidence!=0,1,0)
print("Estab!")
#----------------------------------------------------------------------------

EKE<- read_excel("/home/betherton/EthiopiaINA/Expert Knowledge Elicitation (EKE)_ Potato Ethiopia (Responses).xlsx")
EKE_6<-EKE[,227:275]
EKE_7<-EKE[,276:343]
colnames(EKE_7)<-c("Years","Confidence","perc_Tigray","perc_Afar","perc_Amhara",
                   "perc_Gumuz","perc_Somalia","perc_Oromia","perc_Gambela","perc_SN",
                   "perc_Sidama","perc_SW","perc_Harari","perc_Addis","perc_Diredawa",
                   "season_Tigray","season_Afar","season_Amhara",
                   "season_Gumuz","season_Somalia","season_Oromia","season_Gambela","season_SN",
                   "season_Sidama","season_SW","season_Harari","season_Addis","season_Diredawa",
                   "incomplete_report",
                   "cost_Tigray","cost_Afar","cost_Amhara",
                   "cost_Gumuz","cost_Somalia","cost_Oromia","cost_Gambela","cost_SN",
                   "cost_Sidama","cost_SW","cost_Harari","cost_Addis","cost_Diredawa",
                   "estab_Tigray","estab_Afar","estab_Amhara",
                   "estab_Gumuz","estab_Somalia","estab_Oromia","estab_Gambela","estab_SN",
                   "estab_Sidama","estab_SW","estab_Harari","estab_Addis","estab_Diredawa",
                   "comm_Tigray","comm_Afar","comm_Amhara",
                   "comm_Gumuz","comm_Somalia","comm_Oromia","comm_Gambela","comm_SN",
                   "comm_Sidama","comm_SW","comm_Harari","comm_Addis","comm_Diredawa")#replace colnames
subst<-EKE_7[,56:68]
comm_mat<-matrix(0,nrow=13,ncol=13)
colnames(comm_mat)<-c("Tigray","Afar","Amhara",
                      "Gumuz","Somalia","Oromia","Gambela","SN",
                      "Sidama","SW","Harari","Addis","Diredawa")
rownames(comm_mat)<-c("Tigray","Afar","Amhara",
                      "Gumuz","Somalia","Oromia","Gambela","SN",
                      "Sidama","SW","Harari","Addis","Diredawa")
ls<-c("Tigray","Afar","Amhara",
      "Benishangul-Gumuz","Somalia","Oromia","Gambela","Southern Nations, Nationalities and Peoples'",
      "Sidama","Southwest Ethiopia Peoples'","Harari","Addis Ababa","Diredawa")
for(i in 1:13){
  comm_mat[i,1]<-comm_mat[i,1]+length(which(grepl(ls[i],subst$comm_Tigray)==TRUE))
  comm_mat[i,2]<-comm_mat[i,2]+length(which(grepl(ls[i],subst$comm_Afar)==TRUE))
  comm_mat[i,3]<-comm_mat[i,3]+length(which(grepl(ls[i],subst$comm_Amhara)==TRUE))
  comm_mat[i,4]<-comm_mat[i,4]+length(which(grepl(ls[i],subst$comm_Gumuz)==TRUE))
  comm_mat[i,5]<-comm_mat[i,5]+length(which(grepl(ls[i],subst$comm_Somalia)==TRUE))
  comm_mat[i,6]<-comm_mat[i,6]+length(which(grepl(ls[i],subst$comm_Oromia)==TRUE))
  comm_mat[i,7]<-comm_mat[i,7]+length(which(grepl(ls[i],subst$comm_Gambela)==TRUE))
  comm_mat[i,8]<-comm_mat[i,8]+length(which(grepl(ls[i],subst$comm_SN)==TRUE))
  comm_mat[i,9]<-comm_mat[i,9]+length(which(grepl(ls[i],subst$comm_Sidama)==TRUE))
  comm_mat[i,10]<-comm_mat[i,10]+length(which(grepl(ls[i],subst$comm_SW)==TRUE))
  comm_mat[i,11]<-comm_mat[i,11]+length(which(grepl(ls[i],subst$comm_Harari)==TRUE))
  comm_mat[i,12]<-comm_mat[i,12]+length(which(grepl(ls[i],subst$comm_Addis)==TRUE))
  comm_mat[i,13]<-comm_mat[i,13]+length(which(grepl(ls[i],subst$comm_Diredawa)==TRUE))
}

#lets make our bio adj matrix now

subst<-EKE_6[,1:13]
biomat<-matrix(0,nrow=13,ncol=13)
colnames(subst)<-c("bioTigray","bioAfar","bioAmhara",
                   "bioGumuz","bioSomalia","bioOromia","bioGambela","bioSN",
                   "bioSidama","bioSW","bioHarari","bioAddis","bioDiredawa")
colnames(biomat)<-c("Tigray","Afar","Amhara",
                    "Gumuz","Somalia","Oromia","Gambela","SN",
                    "Sidama","SW","Harari","Addis","Diredawa")
rownames(biomat)<-c("Tigray","Afar","Amhara",
                    "Gumuz","Somalia","Oromia","Gambela","SN",
                    "Sidama","SW","Harari","Addis","Diredawa")
ls<-c("Tigray","Afar","Amhara",
      "Benishangul-Gumuz","Somalia","Oromia","Gambela","Southern Nations, Nationalities and Peoples'",
      "Sidama","Southwest Ethiopia Peoples'","Harari","Addis Ababa","Diredawa")
for(i in 1:13){
  biomat[i,1]<-(biomat[i,1]+length(which(grepl(ls[i],subst$bioTigray)==TRUE)))
  biomat[i,2]<-(biomat[i,2]+length(which(grepl(ls[i],subst$bioAfar)==TRUE)))
  biomat[i,3]<-(biomat[i,3]+length(which(grepl(ls[i],subst$bioAmhara)==TRUE)))
  biomat[i,4]<-(biomat[i,4]+length(which(grepl(ls[i],subst$bioGumuz)==TRUE)))
  biomat[i,5]<-(biomat[i,5]+length(which(grepl(ls[i],subst$bioSomalia)==TRUE)))
  biomat[i,6]<-(biomat[i,6]+length(which(grepl(ls[i],subst$bioOromia)==TRUE)))
  biomat[i,7]<-(biomat[i,7]+length(which(grepl(ls[i],subst$bioGambela)==TRUE)))
  biomat[i,8]<-(biomat[i,8]+length(which(grepl(ls[i],subst$bioSN)==TRUE)))
  biomat[i,9]<-(biomat[i,9]+length(which(grepl(ls[i],subst$bioSidama)==TRUE)))
  biomat[i,10]<-(biomat[i,10]+length(which(grepl(ls[i],subst$bioSW)==TRUE)))
  biomat[i,11]<-(biomat[i,11]+length(which(grepl(ls[i],subst$bioHarari)==TRUE)))
  biomat[i,12]<-(biomat[i,12]+length(which(grepl(ls[i],subst$bioAddis)==TRUE)))
  biomat[i,13]<-(biomat[i,13]+length(which(grepl(ls[i],subst$bioDiredawa)==TRUE)))
}

rm(subst,EKE,EKE_7,EKE_6)

colnames(biomat)<-c("Tigray","Afar","Amhara",
                    "Beneshangul Gumu","Somali","Oromia","Gambela","SNNPR",
                    "Sidama","SW","Hareri","Addis Ababa","Dire Dawa")
rownames(biomat)<-c("Tigray","Afar","Amhara",
                    "Beneshangul Gumu","Somali","Oromia","Gambela","SNNPR",
                    "Sidama","SW","Hareri","Addis Ababa","Dire Dawa")

big_comm_mat<-matrix(0,nrow=dim(sp_points)[1],ncol=dim(sp_points)[1])
big_bio_mat<-matrix(0,nrow=dim(sp_points)[1],ncol=dim(sp_points)[1])

for(i in 1:dim(sp_points)[1]){
  for(j in 1:dim(sp_points)[1]){
    big_bio_mat[i,j]<-biomat[which(colnames(biomat)==sp_points$region[i]),which(rownames(biomat)==sp_points$region[j])]
  }
  print(i)
} #may take a minute. It should count to 3462 (with the current crop)
colnames(comm_mat)<-c("Tigray","Afar","Amhara",
                      "Beneshangul Gumu","Somali","Oromia","Gambela","SNNPR",
                      "Sidama","SW","Hareri","Addis Ababa","Dire Dawa")
rownames(comm_mat)<-c("Tigray","Afar","Amhara",
                      "Beneshangul Gumu","Somali","Oromia","Gambela","SNNPR",
                      "Sidama","SW","Hareri","Addis Ababa","Dire Dawa")
for(i in 1:dim(sp_points)[1]){
  for(j in 1:dim(sp_points)[1]){
    big_comm_mat[i,j]<-comm_mat[which(colnames(comm_mat)==sp_points$region[i]),which(rownames(comm_mat)==sp_points$region[j])]
  }
  print(i)
} #may take a minute. Counts to 3462 with current Ethiopia crop

big_bio_mat<-big_bio_mat/max(big_bio_mat)
big_comm_mat<-big_comm_mat/max(big_comm_mat)
big_bio_mat[which(big_bio_mat>0.6)]<-1
big_bio_mat[which(big_bio_mat<0.6)]<-0
big_comm_mat[which(big_comm_mat>0.6)]<-1
big_comm_mat[which(big_comm_mat<0.6)]<-0

#------------------------------------------------------------------------------

#Lets run INA

ethiopia.ina <- INAscene(nreals=10, #lets start small
                         ntimesteps=10, #lets just do 10 times steps for simplicity
                         doplot=FALSE, #and were plotting the outputs
                         readgeocoords=TRUE,
                         geocoords = new_coords, #from mapSPAM (a smaller subset representing 6 regions)
                         readinitinfo=F, 
                         initinfo=NA, 
                         initinfo.norp='num', #initial info was essentially everywhere
                         initinfo.n=5, 
                         initinfo.p=NA, 
                         initinfo.dist='random', 
                         readinitbio=TRUE, 
                         initbio=as.vector(sp_points$init_loc), #bio presence based on EKE survey
                         readseam=TRUE, 
                         seam=big_comm_mat, #comm adj mat based on EKE
                         readbpam=TRUE, 
                         bpam=big_bio_mat, #bio mat based on EKE
                         readprobadoptvec=F, 
                         probadoptvec=NA, 
                         probadoptmean=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), #probability of adoption is high
                         probadoptsd=0.1,
                         readprobestabvec=TRUE, 
                         probestabvec=as.vector(sp_points$estab), 
                         #probestabmean=c(0.5,0.8), #probability of estab is high w/out mgmt
                         #probestabsd=0.1, #we can modify this estab based off pH?
                         maneffdir='decrease_estab', #mgmt decreases estab
                         maneffmean=0.5, 
                         maneffsd=0.1, 
                         usethreshman=F, 
                         maneffthresh=NA, 
                         sampeffort=NA,
                         outputvol="less"                   
)
eth_out<-ethiopia.ina$multdetails
eth_multout<-ethiopia.ina$multout

save.image("/home/betherton/EthiopiaINA/EthOUTImg2.RData")
saveRDS(eth_out,file="/home/betherton/EthiopiaINA/EthOUT2.RDS")

