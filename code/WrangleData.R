# The file required for this piece of code has not been made available as only
# santizied data was released
rename_group<-function(someDat=NULL,visit=NULL){
  #die if some data is not provided
  if(is.null(someDat) || length(someDat)<2){stop("!needs more than one date to work")}
  if(is.null(visit)){stop("! Visit information is not provided")}
  
  #just in case recode visit
  if(visit=="First" || visit==1){
    visit=1
  }else{visit=2}
  
  #get min date
  startDate<-min(as.Date(someDat,format="%Y-%m-%d"))
  endDate<-max(as.Date(someDat,format="%Y-%m-%d"))
  
  indexDates<-seq(from=startDate, to=endDate, by=1)
  
  #safety check
  if(length(indexDates)>100){stop("Something is very wrong")}
  
  renamed<-sapply(someDat,function(x){
    day<-which(x==indexDates)
    sprintf("V%d-D%d",visit,day)
  })
  
  return(renamed)
}

############
# MAIN
#############
library(dplyr)
dat<-read.csv(file="Data/Data_for_analysis_deidentified_R_ready_3.csv",header=T,stringsAsFactors=F)

#rename the dates as V#-D# where V= index case stay (1 or 2), and D= day during stay period


temp1<-subset(dat,IndexCaseVisit=="First")
temp1$DateOfVisit2<-rename_group(temp1$DateOfVisit,1)

temp2<-subset(dat,IndexCaseVisit=="Second")
temp2$DateOfVisit2<-rename_group(temp2$DateOfVisit,2)

dat<-rbind(temp1,temp2)
dat$Bed<-dat$value
dat$DateOfVisit<-dat$DateOfVisit2



write.csv(file="./Data/Data_for_analysis_deidentified_R_ready_sanitized.csv",dat[,c("KEY","CaseStatus_shrt","Bed","DateOfVisit","IndexCaseVisit")],quote=F)

