#setwd("~/Documents/projects/TB_Bedmap/")

library(ggplot2)
library(reshape2)
library(scales)
library(rpart)
library(caret)
library(broom)
library(foreach)
library(dplyr)


source("./Code/HelperFunctions.R")

set.seed(1)

###############################################
# FUNCTIONS
##############################################

getResults<-function(fitDat = NULL,roundVar = 2)
{
  sumDat  <- summary(fitDat)
  fitTab  <- sumDat$coefficients
  
  pvals <- round(fitTab[2:nrow(fitTab),"Pr(>|z|)"],(roundVar+1))
  OR <- round(exp(fitTab[2:nrow(fitTab),"Estimate"]),(roundVar+1))
  
  confTab <- round(exp(confint(fitDat)),roundVar)
  
  UCI  <- confTab[2:nrow(confTab),2]
  LCI <- confTab[2:nrow(confTab),1]
  
  sumRes  <- cbind(OR,LCI,UCI,pvals)
  
  return(sumRes)
}

convertORtoRR<-function(OR = NULL,p = 0.1)
{
  return(OR/(1 -0.1 + (0.1 * OR)))
}

###############################################
# LOAD THE DATA DO SOME CLEAN UP !
##############################################


roomLayout<-read.csv(file="./Data/Room_Coordinates_2008.csv",header=T)
bedCoord<-read.csv(file="./Data/Bed_Coordinates_2008.csv",header=T)
dat<-read.csv(file="./Data/Data_for_analysis_deidentified_R_ready_3.csv",header=T,stringsAsFactors=F)
dat[dat$KEY==1,]$CaseStatus_shrt<-"Index case"

#wrangle the data  a bit
dat$Bed<-dat$value
dat$DateOfVisit<-as.Date(gsub("X","",dat$variable),format="%m.%d.%y")

#Sort out client outcome levels
dat$CaseStatus_shrt<-factor(dat$CaseStatus_shrt,levels=c("Index case","Active","Latent","Uninfected","PriorInfection","Unknown-TSTnotRead","Unknown-noDB","Unknown-noScreen"))
dat<-dat[!is.na(dat$Bed),]

#### adding the bed co-ordinate data to the demographic data
dat<-merge(dat,bedCoord,by="Bed")



###############################################
# CALCULATE SOME DISTANCES!
##############################################


#calculating the distance between the index case
#and a given client on given day of index case visit
#get information about the index case
#note, the clients need to be in the same
#zones to be assigned a distance measure

#idxData is the data for the Index case
idxData = dat[dat$KEY==1,]
distances = c()

for(i in 1:nrow(dat))
{
  #compare find where the index case was on the day of
  #that particular client was in the shelter
  idx<-which(idxData$DateOfVisit==dat[i,]$DateOfVisit)
  
  #
  if (!(is.na(dat[i,]$X) || is.na(dat[i,]$Y)))
  {
    #check to see if in the same zone as index
    if(dat[i,]$Zone == idxData[idx,]$Zone){
      D =sqrt((dat[i,]$X-idxData[idx,]$X)^2 + (dat[i,]$Y-idxData[idx,]$Y)^2 + (dat[i,]$Z-idxData[idx,]$Z)^2)
    }
    else
    {
      #if they are not in the same zone assign some maximial distance (about 24 m is the largest distance)
      #between two points in the shetler
      D = 24
    }
    distances = c(distances,D)
  }
  else{
    #in the weird possibility that someone isn't in the shelter in the
    # master data (shouldn't be the case) just give them an NA
    distances = c(distances,NA)
  }
}
dat$distFromIndex <- distances


# add exposure time to index case (will just be how often a case shows up)
# and average distance,minimum distance and the cumulative distance over the says

temp<-unique(dat$KEY)

dat$ExpoDats<-rep(0,nrow(dat))
dat$avgDist<-rep(0,nrow(dat))
dat$sumDist<-rep(0,nrow(dat))
dat$minDist<- rep(0,nrow(dat))

#iterate through each of the clients and
# get their exposure and distance data
for(i in 1:length(temp))
{
  tempDat<-which(dat$KEY == temp[i])
  
  dat[tempDat,]$ExpoDats<-length(tempDat)
  
  dat[tempDat,]$avgDist<-rep(mean(dat[tempDat,]$distFromIndex,na.rm=T),length(tempDat))
  dat[tempDat,]$minDist<-rep(min(dat[tempDat,]$distFromIndex,na.rm=T),length(tempDat))
  dat[tempDat,]$sumDist<-sum(dat[tempDat,]$distFromIndex,na.rm=T)
}



#######################################
#     ASSESS SOME DEMOGRAPHICS
#######################################

dat.first<-dat[dat$IndexCaseVisit=="First",]
dat.second<-dat[dat$IndexCaseVisit=="Second",]

#Find which clients were there both stays
both = intersect(dat.first$KEY,dat.second$KEY)
both.tab<-table(dat[match(both,dat$KEY),]$CaseStatus_shrt)

#Find which clients were there only for the first stay
firstOnly = setdiff(dat.first$KEY,dat.second$KEY)
first.tab<-table(dat[match(firstOnly,dat$KEY),]$CaseStatus_shrt)

#Find which clients wer there only for the second stay
secondOnly=setdiff(dat.second$KEY,dat.first$KEY)
second.tab<-table(dat[match(secondOnly,dat$KEY),]$CaseStatus_shrt)

#summarize the information, and sum the rows as well


caseDemo <-merge(data.frame(both.tab),merge(data.frame(first.tab),data.frame(second.tab),by = "Var1",all = T),by="Var1",all=T)
colnames(caseDemo)<-c("Outcome","BothVisits","First only","Second only")

caseDemo$Total<-apply(caseDemo,1,function(x){sum(as.numeric(x[2:4]),na.rm=T)})
caseDemo<-caseDemo[,c(1,5,2,3,4)]

write.csv(caseDemo,file="./Tables/CaseDemoByVisit.csv",quote=T)

#Add a label indicating whether a client was there
# only for the first visit, the second or the both
dat$Visits = rep(NA,nrow(dat))

for(key in unique(dat$KEY))
{
  visitTimes = "Second only"
  if(key %in% both ){visitTimes = "Both"}
  else if (key %in% firstOnly){visitTimes = "First only"}
  dat[which(dat$KEY == key),]$Visits <- visitTimes
}
dat$Visits<-factor(dat$Visits,levels=c("Both","First only","Second only"))



#######################################
#     VISUALIZE SOME CLIENT MOVEMENT
#######################################

#According to the date determine the following:
# 1) How many new clients entered the shelter
# 2) How many clients were not yet in the shelter
# 3) How many clients left the shelter:
# 4) How many clients changed location
# 5) How many clients stayed in the same location

#set up temporary holding variable
caseID<-sort(unique(dat$KEY))
prevDatStatus<-rep(0,length(caseID))
prevBedNum<-rep(0,length(caseID))
names(prevDatStatus)<-caseID

temp<-dat[match(unique(dat$KEY),dat$KEY),c("KEY","CaseStatus_shrt")]
caseStat<-temp[order(temp$KEY),]$CaseStatus_shrt

#set up the shelter matrix
shelterMove<-c()

prevDate<-0

for(date in as.character(sort(unique(dat$DateOfVisit))))
{
  #date <- sort(unique(dat$DateOfVisit))[1]
  # 1) How many new clients entered the shelter
  subData<-subset(dat,dat$DateOfVisit == date)
  subData<-subData[order(subData$KEY),]
  
  inShelter<-as.numeric(caseID %in% subData$KEY)
  
  # If it's the first day, there can be no inner movement
  # because we don't have that data
  if(sum(abs(prevDatStatus)) > 0)
  {
    moveStatus<-subData$Bed == prevBedNum[inShelter == 1] 
    
    #move status = FALSE means, the client moved to a new bed
    #Code 3 means : already in the shelter and moved
    #Code 2 means : already in shelter and stayed in the same place
    #Code 1 means : newly arrived to the shetler on that day
    idx.noMoved<-prevDatStatus[inShelter==1]>0 & moveStatus 
    idx.moved<-prevDatStatus[inShelter==1]>0 & !(moveStatus) 
  
    #now update the inshelter codes
    #it's possible that people left and came to the same bed
    #even if the dates weren't consequctive; but I am going
    #to assume that they were still in the shelter and the index
    #case did the moving
    if(sum(idx.noMoved)>0){inShelter[which(inShelter >0)[idx.noMoved]]<-2}
   
    # look at those cases with label 3, if the days weren't consecutive
    # then it means the client probably left the shelter.
    if((prevDate == (as.Date(date)-1)) & sum(idx.moved)>0)
    {
      #so only assume they've actually moved around
      #if it's consecutive days, otherwise
      #they've probably left.
      inShelter[which(inShelter >0)[idx.moved]]<-3
    }
    #now update the previous bed information
    prevBedNum[inShelter>0]<-subData$Bed
  
    #for those clients that have *moved out* of the shetler, set
    #their bed status to zero
    idx.moveOut<-(inShelter==0 & prevDatStatus>=1)
    prevBedNum[idx.moveOut]<-0
    inShelter[idx.moveOut]<- -1
  
  }
  else{
    #Record bed number from the first visit
    prevBedNum[subData$KEY]<-subData$Bed
  }
  #Tabulate the status
  
  prevDatStatus<-inShelter
  prevDate<-as.Date(date)
  
  inShelter<-factor(inShelter,levels=c(-1,0,1,2,3))
  tabDat<-table(caseStat,inShelter)
  
  #collapse the unknowns
  tabDat<-prop.table(rbind(tabDat[1:4,],colSums(tabDat[5:8,])),1)
  rownames(tabDat)<-c(rownames(tabDat)[1:4],"Unknown")  
  
  print(tabDat)
  shelterMove<-rbind(shelterMove,cbind(rep(date,times=nrow(tabDat)),tabDat))

} 

temp<-rownames(shelterMove)
shelterMove<-data.frame(shelterMove,stringsAsFactors=FALSE)
shelterMove$CaseStatus<-temp
colnames(shelterMove)<-c("Date","-1","0","1","2","3","CaseStatus")
shelterMove<-melt(shelterMove,id.vars=c("Date","CaseStatus"))

shelterMove$value<-as.numeric(as.character(shelterMove$value))
shelterMove$Date<-factor(shelterMove$Date)
shelterMove$CaseStatus<-factor(shelterMove$CaseStatus,levels=c("Index Case","Active","Latent","Uninfected","Unknown"))
#tile plot to show the movement over time

movement<-ggplot(data=shelterMove,aes(x=Date,y=variable))+
  geom_tile(aes(fill=value),colour="white")+
  scale_fill_gradient(low="white",high = "black",name = "% of all outcome group",labels = c("0%","25%","50%","75%","100%"))+
  facet_grid(CaseStatus~.)+
  theme_bw()+
  xlab("")+
  ylab("Client Movement Status")+
  theme(axis.text.x=element_text(angle=90))

ggsave(movement,file="Client_Movement_Status_by_Date.pdf")




### VISUALIZATION MOVEMENT WITHING THE ROOM GEOGRAPHY

#a. Visualizing the movement by co-ordinate location in the dorm
pallete = c("white","#d73027","#fc8d59","#4575b4","#e9a3c9","#fee090","#fee090","#fee090")#

moveByRoom <- ggplot(data=dat,aes(x=X,y=Y)) +
  geom_polygon(data=roomLayout,aes(x=X,y=Y,group=Zone),fill="black",colour="white",alpha=0.7)+
  geom_point(aes(colour=CaseStatus_shrt,size=CaseStatus_shrt),alpha=0.75,position=position_jitter(height=0.3,width=0.3))+
  geom_line(data=dat[dat$KEY==1,],aes(x=X,y=Y),colour="white",size=2,alpha=0.4)+
  scale_size_manual(values=c(5,4,4,4,1,1,1,1),name = "Client outcome")+
  scale_colour_manual(values=pallete,name="Client outcome")+
  facet_grid(IndexCaseVisit~.)+
  xlab("Room width (m)")+
  ylab("Room height (m)")+
  theme_bw()+
  theme(legend.key=element_rect(fill="#00000070"),panel.background=element_rect(fill="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank())

moveByRoom <- moveByRoom + theme(legend.text=element_text(size=12,family="Helvetica",face="plain"),
                   legend.title=element_text(size=12,face="plain",family="Helvetica"),
                   axis.text = element_text(family="Helvetica",size=12),
                   axis.title=element_text(family="Helvetica",size=12))

ggsave(moveByRoom,file="./revisions/Rev 1.3//Client_Movement_By_Room.jpeg",units="in",height=4,width=7,dpi=300)



#b. Visualizing the movement of the index case according to bed
pallete = c("#000000","#d73027","#fc8d59","#4575b4","#e9a3c9","#fee090","#fee090","#fee090")

moveByBed <- ggplot(data= dat,aes(x=as.character(DateOfVisit),y=Bed,group=KEY)) +
  geom_line(aes(colour = CaseStatus_shrt,alpha=CaseStatus_shrt))+
  geom_point(aes(colour = CaseStatus_shrt,size=CaseStatus_shrt,alpha=CaseStatus_shrt),pch=21,fill= "white")+
  scale_colour_manual(values=pallete,name="Client Outcome") +
  scale_size_manual(values = c(3,2,2,2,1,1,1,1),name= "Client Outcome")+
  scale_alpha_manual(values=c(1.0,0.9,0.9,0.9,0.5,0.5,0.5,0.5),name= "Client Outcome")+
  theme_bw()+
  ylab("Bed Number")+
  xlab("Date of Visit")+
  theme(axis.text.x= element_text(angle = 90),axis.text=element_text(size=14),axis.title=element_text(face="bold",size=16))
  #theme(axis.text.x = element_blank())
ggsave(moveByBed,file="./figures/Client_Movement_By_BedNum.svg",units="cm",height=14,width=26.5)




#################################################################
#    EVALULATE THE CURRENT TIME AND DISTANCE DATA
#################################################################
#using a multinom regression I will check how the time and distance
#variables stack up

temp<-dat[match(unique(dat$KEY),dat$KEY),]

#remove the index case
temp<-temp[!(temp$KEY == 1),]

#some quick checks:
temp$caseLabel <- as.numeric(temp$CaseStatus_shrt %in% c("Active","Latent","Uninfected"))

# 1. Total case in Shetler, and having a case label:
summary(glm(caseLabel~ExpoDats,data = temp,family = binomial(link="logit")))


# 2. Whether you are more likle to have a case label in visit 1 vs visit 2
chisq.test(with(temp,table(caseLabel,Visits)))

#Ok - now do analysis with clients that have known outcomes
datEval<-dat[dat$CaseStatus_shrt %in% c("Active","Latent","Uninfected"),]
datEval.Uni<-datEval[match(unique(datEval$KEY),datEval$KEY),]


#c. Visualizing movement according to proximity from the index case

moveByDist <-ggplot(data= datEval,aes(x=as.factor(DateOfVisit),y=distFromIndex,group=KEY)) +
  geom_line(aes(colour = CaseStatus_shrt,alpha=CaseStatus_shrt))+
  geom_point(aes(colour = CaseStatus_shrt,fill=CaseStatus_shrt,alpha=CaseStatus_shrt),pch=21)+
  scale_colour_manual(values=pallete[2:9],name="Client outcome") +
  scale_fill_manual(values=pallete[2:9],name="Client outcome")+
  scale_size_manual(values = c(3,2,2,2,1,1,1,1),name= "Client outcome")+
  scale_alpha_manual(values=c(1.0,0.9,0.9,0.9,0.5,0.5,0.5,0.5),name= "Client outcome")+
  scale_x_discrete(labels=c("V1-D1","V1-D6","V1-D7","V1-D8","V1-D9","V1-D10","V1-D11","V2-D1","V2-D5","V2-D6","V2-D9","V2-D10"))+
  #scale_x_date(labels = date_format("%b-%d-%y"),minor_breaks=NULL,breaks=unique(dat$DateOfVisit))+
  coord_polar()+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.text.y = element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14,face="plain"),
        axis.text.x = element_text(face="bold",size=14,colour = c(rep("black",7),rep("azure4",5))),
        axis.ticks.y=element_blank(),panel.border=element_rect(fill=NA,colour=NA))
  

moveByDist<-moveByDist + theme(legend.text=element_text(size=12,family="Helvetica",face="plain"),
                   legend.title=element_text(size=12,face="plain",family="Helvetica"),
                   axis.text = element_text(family="Helvetica",size=12),
                   axis.title=element_text(family="Helvetica",size=12))

ggsave(moveByDist,file="./revisions/Rev 1.3/Client_Movement_By_Distance.jpeg",units="in",height=5,width=7)







regressionResults<-c()

# ..... EXPOSURE DAYS EVALUATION

#Assuming that the index case infected *all* of the latent cases, calculated the OR
# Infected vs. Uninfected
fit.glm<-glm(as.numeric(CaseStatus_shrt != "Uninfected")~ExpoDats,data = datEval.Uni,family = binomial(link="logit"))
regressionResults<-rbind(regressionResults, cbind(rep("Univariable,Exposure Days",times=1),getResults(fit.glm)))


# ..... SUM OF DISTANCE EVALUATION
fit.glm<-glm(as.numeric(CaseStatus_shrt != "Uninfected")~sumDist,data = datEval.Uni,family = binomial(link="logit"))
regressionResults<-rbind(regressionResults, cbind(rep("Univariable,Sum Distance",times=1),getResults(fit.glm)))

# ..... AVERAGE DISTANCE EVALUATION
fit.glm<-glm(as.numeric(CaseStatus_shrt != "Uninfected")~avgDist,data = datEval.Uni,family = binomial(link="logit"))
regressionResults<-rbind(regressionResults, cbind(rep("Univariable,Average Distance",times=1),getResults(fit.glm)))


# ..... MINIMUM DISTANCE EVALUATION
fit.glm<-glm(as.numeric(CaseStatus_shrt != "Uninfected")~minDist,data = datEval.Uni,family = binomial(link="logit"))
regressionResults<-rbind(regressionResults, cbind(rep("Univariable,Min Distance",times=1),getResults(fit.glm)))

# ..... EXPOSURE DAYS / SUM RATIO EVALUATION
datEval.Uni$expRatio <-with(datEval.Uni,sumDist/ExpoDats)
fit.glm<-glm(as.numeric(CaseStatus_shrt != "Uninfected")~expRatio,data = datEval.Uni,family = binomial(link="logit"))
regressionResults<-rbind(regressionResults, cbind(rep("Univariable,expoRatio",times=1),getResults(fit.glm)))


write.csv(file = "./tables/RegressionResults.csv",printPretty(regressionResults),quote=FALSE)

#re-evalute the above under differen scenarios [respoding to review request]

#1. Infected = 36, uninfected = 22
#gc()
#2. Infected = 34, uninfected = 24

maxInfect = table(datEval.Uni$CaseStatus_shrt)["Latent"]

datEval.Uni$binClass <- as.character(mapvalues(datEval.Uni$CaseStatus_shrt,from=as.character(unique(datEval.Uni$CaseStatus_shrt)),to=c(0,NA,1)))


totalInfectedPior<-c()
for(i in 1:30){
    totalChoices<-choose(31,i)
    
    #sometimes total choices is a huge number (> 1 million), when that happens, just try 10,000 combos
    #the max is 31 chooose 15, which has a total 300,540,195 combinations..that would take a long time to run
    
    if(totalChoices > 10000){totalChoices = 10000}
    
    #storing the temporary results
    tempResults <-  list(expo=list(OR=c(),pval=c()),
                         distAvg=list(OR=c(),pval=c()),
                         distMin=list(OR=c(),pval=c()),
                         distexpRatio=list(OR=c(),pval=c()))
    
    for(trail in 1:totalChoices){
      
      idxLatent<-which(datEval.Uni$CaseStatus_shrt == "Latent")
      
      #notinfect implies that the latent cases had TB when they came to the shelter
      #infected implies that the latent cases were infected by the index case
      notinfected <- sample(idxLatent,size=i,replace=FALSE)
      infected <- setdiff(idxLatent,notinfected)
      
      datEval.Uni$binClass[infected]<-1
      datEval.Uni$binClass[notinfected]<-0
      
      
      # ..... SUM OF DISTANCE EVALUATION
      fit.glm<-tidy(glm(as.numeric(binClass)~ExpoDats,data = datEval.Uni,family = binomial(link="logit")))
      tempResults$expo$OR  <-  c(tempResults$expo$OR,exp(fit.glm$estimate[2]))
      tempResults$expo$pval  <-  c(tempResults$expo$pval,fit.glm$p.value[2])
      
      # ..... AVERAGE DISTANCE EVALUATION
      fit.glm<-tidy(glm(as.numeric(CaseStatus_shrt != "Uninfected")~avgDist,data = datEval.Uni,family = binomial(link="logit")))
      tempResults$distAvg$OR  <-  c(tempResults$distAvg$OR,exp(fit.glm$estimate[2]))
      tempResults$distAvg$pval <-  c(tempResults$distAvg$pval,fit.glm$p.value[2])
      
      # ..... MINIMUM DISTANCE EVALUATION
      fit.glm<-tidy(glm(as.numeric(CaseStatus_shrt != "Uninfected")~minDist,data = datEval.Uni,family = binomial(link="logit")))
      tempResults$distMin$OR  <-  c(tempResults$distMin$OR,exp(fit.glm$estimate[2]))
      tempResults$distMin$pval <-  c(tempResults$distMin$pval,fit.glm$p.value[2])
      
      
      # ..... EXPOSURE DAYS / SUM RATIO EVALUATION
      datEval.Uni$expRatio <-with(datEval.Uni,sumDist/ExpoDats)
      fit.glm<-tidy(glm(as.numeric(CaseStatus_shrt != "Uninfected")~expRatio,data = datEval.Uni,family = binomial(link="logit")))
      
      tempResults$distexpRatio$OR  <-  c(tempResults$distexpRatio$OR ,exp(fit.glm$estimate[2]))
      tempResults$distexpRatio$pval  <-  c(tempResults$distexpRatio$pval,fit.glm$p.value[2])
      
    }
    
    #now summarize the results for each
    #first with with OR
    ranges<-sapply(tempResults,function(x){
      tmp  <-  summary(x[[1]])
      list(OR = list(min=tmp[1],max=tmp[6],Median=tmp[3]),
           pVal = list(totalSig=sum(x[[2]]<0.05),totalNotSig=sum(x[[2]]>0.05)))
    })
    
    totalInfectedPior<-rbind(totalInfectedPior,cbind(rep(i,nrow(ranges)),rep(totalChoices,nrow(ranges)),ranges))
}

load(file="TotalInfected.Rda")

##### first, summarize the p-value results
pVal <- totalInfectedPior[which(rownames(totalInfectedPior) == "pVal"),]

pVal<-foreach(i=1:nrow(pVal),.combine=rbind)%do%{
  x  <- unlist(pVal[i,])
  tmp <- x[3:length(x)]
  tmpName<-t(sapply(names(tmp),function(y){unlist(strsplit(y,"\\."))}))
  
  cbind(rep(x[1],nrow(tmpName)),rep(x[2],nrow(tmpName)),tmpName[,1],tmpName[,2],tmp)
}

colnames(pVal)<-c("totalPre","totalPerm","varofInt","sigStatus","numSig")
pVal <- data.frame(pVal)

#there was no circumstance it which the distance measures were significant
#only exposre distance was significant in the scenarios, so I will just filter with that

pVal <- pVal %>% filter(varofInt == "expo")
pVal$perSig <- as.numeric(as.character(pVal$numSig))/as.numeric(as.character(pVal$totalPerm))
  
ggplot(pVal,aes(x=as.numeric(as.character(totalPre)),y=perSig))+
  geom_bar(stat="identity",colour="black",aes(fill=sigStatus))+
  scale_fill_manual(values=c("white","black"),name=expression(paste(italic("P"), "< 0.05")),labels=c("False","True"))+
  scale_x_continuous(breaks=seq(from=1,to=maxInfect,by=1))+
  scale_y_continuous(labels = percent)+
  ylab("Percentage of runs with significant outcomes")+
  xlab("Total number of clients with prior infection")+
  theme_bw() + 
  theme(legend.text=element_text(size=12,family="Helvetica",face="plain"),
        legend.title=element_text(size=12,face="plain",family="Helvetica"),
        axis.text = element_text(family="Helvetica",size=12),
        axis.title=element_text(family="Helvetica",size=12),axis.text.x=element_text(angle=90,vjust=0.5))

ggsave("./revisions/Rev 1.3//SensitivityAnalysis_pvalue.jpeg",unit="in",height=4,width=7,dpi=300)

##### next, summarize the OR results
OR <- totalInfectedPior[which(rownames(totalInfectedPior) == "OR"),]

OR<-foreach(i=1:nrow(OR),.combine=rbind)%do%{
  x  <- unlist(OR[i,])
  tmp <- x[3:length(x)]
  tmpName<-t(sapply(names(tmp),function(y){unlist(strsplit(y,"\\."))}))
  
  cbind(rep(x[1],nrow(tmpName)),rep(x[2],nrow(tmpName)),tmpName[,1],tmpName[,2],tmp)
}

colnames(OR)<-c("totalPre","totalPerm","varofInt","varVal","value")
OR <- data.frame(OR)

#seperate the  min and max and median
tmpMedian <- filter(OR,varofInt == "expo") %>% filter(varVal == "Median")
tmpMin <- filter(OR,varofInt =="expo") %>% filter(varVal == "min")
tmpMax <- filter(OR,varofInt =="expo") %>% filter(varVal == "max")

#now make it all one happy thing
tmp<-data.frame(totalPre=factor(tmpMedian$totalPre,levels=sort(as.numeric(as.character(tmpMedian$totalPre)))),
                totalPerm = tmpMedian$totalPerm,
                Median = as.numeric(as.character(tmpMedian$value)),
                Min = as.numeric(as.character(tmpMin$value)),
                Max = as.numeric(as.character(tmpMax$value)))


ggplot(tmp,aes(x=totalPre,y=Median,group=totalPre))+
  geom_point(size=5)+
  geom_errorbar(aes(ymax=Max,ymin=Min,group=totalPre))+
  geom_hline(yintercept=1,colour="blue",linetype="dashed")+
  ylab("Odds ratio")+
  xlab("Total number of pre-infected clients")+
  theme_bw()+
  theme(legend.text=element_text(size=12,family="Helvetica",face="plain"),
          legend.title=element_text(size=12,face="plain",family="Helvetica"),
          axis.text = element_text(family="Helvetica",size=12),
          axis.title=element_text(family="Helvetica",size=12),
          axis.text.x = element_text(angle=90,vjust=0.5))

ggsave("revisions/Rev 1.3//SensitivityAnalysis_oddsRatio.jpeg",units="in",dpi=300,height=4,width=7)

#################################################################
#    Deriving and Evaluating cutoffs  for exisiting variables
#################################################################

# Since exposure time did appear to be significant, let's look into whether we could
# also derive a  cutpoint for this

#Now see how in this study

FitOrg <- rpart(CaseStatus_shrt~ExpoDats,data=datEval.Uni,method="class",control=rpart.control(maxdepth=1,minsplit=5))

#checking rouhgly how mucha difference this would make

fit.glm<-glm(as.numeric(CaseStatus_shrt != "Uninfected")~(ExpoDats >= FitOrg$splits[,"index"]),data = datEval.Uni,family=binomial(link = "logit"))
getResults(fit.glm)

#A binary response
datEval.Uni$binResponse<-as.numeric(datEval.Uni$CaseStatus_shrt != "Uninfected")
results<-bootOptimism(modelFrame=model.frame(CaseStatus_shrt~ExpoDats,data=datEval.Uni),binResponse=datEval.Uni$binResponse)

write.csv(file="./tables/PerformanceEval_optimismPenalized.csv",results)



#plotting the relationship between exposure and total distance per client
ggplot(data=datEval.Uni,aes(x=ExpoDats,y=sumDist/ExpoDats))+
  geom_point(aes(colour = CaseStatus_shrt,shape=Visits),size=5,alpha=0.6,position=position_jitter(height=0.25,width=0.25))+
  scale_x_discrete(breaks=(1:12))+
  scale_colour_manual(values = c("#d73027","#fc8d59","#4575b4"),name="Client outcome")+
  guides(shape = guide_legend("Presence during index case stays"))+
  geom_vline(xintercept=4.5,linetype="dashed")+
  ylab("Cumulative distance to the index case over time (m/day)")+
  xlab("Total days exposed to index case")+
  theme_bw()+
  theme(legend.text=element_text(size=12,family="Helvetica",face="plain"),
                     legend.title=element_text(size=12,face="plain",family="Helvetica"),
                     axis.text = element_text(family="Helvetica",size=10),
                     axis.title=element_text(family="Helvetica",size=10))
  

ggsave(file="./revisions/Rev 1.3/TotalDaysExposed_cumulativeDistance.jpeg",dpi=300,units="in",height=4, width=7)




