# The file required for this piece of code has not been made available as only
# santizied data was released

library(reshape)
library()
#read in the de-identified data
dat<-read.csv(file="./Data/Data_for_analysis_deidentified.csv",header=T)

temp<-melt(dat,id.var=c("KEY","DiagDate","Case.Status","CaseStatus_shrt"))
temp$DateOfVisit<- as.Date(gsub('X','',temp$variable),format="%m.%d.%y")
temp$IndexCaseVisit<-c("First","Second")[factor(temp$DateOfVisit <= "2008-01-11",levels=c("TRUE","FALSE"))]

write.csv(file="./Data/Data_for_analysis_deidentified_R_ready.csv",temp,quote=T)

