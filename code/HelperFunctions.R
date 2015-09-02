accuracy <- function(truth=NULL,pred=NULL)
{
  #create a 2 by 2 table of the predicted vs truth
  twoBYtwo<-table(truth,pred)
  return((twoBYtwo[1,1] + twoBYtwo[2,2])/sum(twoBYtwo))
}


getmultiTable<-function(test=NULL,roundFact=2)
{
  ncol<-ncol(coef(test))
  pTable<-c()
  
  #odds Ratio
  OR<-round(exp(coef(test)[,2:ncol]),roundFact)
  if(ncol>2){
    nameLevel=rownames(OR)
    measures = colnames(coef(test))[2:ncol]
  }else{
    nameLevel=names(OR)
    measures= colnames(coef(test))[2]
  }
  
  nRepL<-length(nameLevel)
  
  pTable<-rbind(pTable,cbind(nameLevel,rep("OR",nRepL),OR))

  #lower confidence interval
  LCI<-round(exp(rbind(confint(test)[,,"Active"][2:ncol,1],confint(test)[,,"Latent"][2:ncol,1])),roundFact)
  pTable<-rbind(pTable,cbind(nameLevel,rep("LCI",nRepL),LCI))

  
  #upper confidence interval
  UCI<-round(exp(rbind(confint(test)[,,"Active"][2:ncol,2],confint(test)[,,"Latent"][2:ncol,2])),roundFact)
  pTable<-rbind(pTable,cbind(nameLevel,rep("UCI",nRepL),UCI))

  #pvalue
  z <- summary(test)$coefficients/summary(test)$standard.errors
  p <- round((1 - pnorm(abs(z), 0, 1)) * 2,roundFact)
  
  pTable<-rbind(pTable,cbind(nameLevel,rep("p-val",nRepL),p[,2:ncol]))

  
  #putting it together and printing it out
  colnames(pTable)<-c("Level","Measure",measures)
  pTable<-melt(data.frame(pTable,stringsAsFactors=FALSE),id.vars=c("Level","Measure"),measure.vars=measures)

  temp<-c()
  #ugly formatting to reshape this thing
  for(uLevel in unique(pTable$Level))
  {
    temp2<-subset(pTable,pTable$Level==uLevel)
    for(var in unique(temp2$variable))
    {
      temp3<-subset(temp2,temp2$variable==var)
      temp<-rbind(temp,c(uLevel,var,as.character(temp3$value)))
    }
  }
  
  return(temp)
}


printPretty<-function(pTable)
{
  #prints a regression table all pretty like
  pretty<-apply(pTable,1,function(x){
    #sprintf("%s,%s,%s (%s - %s),%s",x[2],x[3],x[4],x[5],x[6],x[7])
    sprintf("%s,%s (%s - %s),%s",x[1],x[2],x[3],x[4],x[5])
  })
  
  return(pretty)
}

modelFitting<-function(ModelOrg=NULL,ModelSamp=NULL,binResponseOrg=NULL,binResponseSamp=NULL)
{
    if(is.null(ModelSamp))
    {
      FitOrg <- rpart(ModelOrg,method="class",control=rpart.control(maxdepth=1,minsplit=5))
     
      if(is.null(FitOrg$splits)){return(c(NA,NA))}
      else{
        sens<-specificity(data=factor(as.numeric(ModelOrg[,2]>FitOrg$splits[4])),reference=factor(binResponseOrg),positive="1")
        spec<-specificity(data=factor(as.numeric(ModelOrg[,2]>FitOrg$splits[4])),reference=factor(binResponseOrg),negative="0")
        
        return(c(sens,spec))
      }
    }
    else{
      #Fit a model using the subsampled data
      FitOrg <- rpart(ModelSamp,method="class",control=rpart.control(maxdepth=1,minsplit=5))
      
      #If an optimal cutpoint can't be found, then reutrn NA
      if(is.null(FitOrg$splits)){return(c(NA,NA))}
      else{
        #Evalaute the sensitivity and speficity on the bootstrapped sample
        sens<-sensitivity(data=factor(as.numeric(ModelSamp[,2]>FitOrg$splits[4])),reference=factor(binResponseSamp),positive="1")
        spec<-specificity(data=factor(as.numeric(ModelSamp[,2]>FitOrg$splits[4])),reference=factor(binResponseSamp),negative="0")
        
        #Evaluate the sensitivity and specifcity on the original sample
        sensORG<-sensitivity(data=factor(as.numeric(ModelOrg[,2]>FitOrg$splits[4])),reference=factor(binResponseOrg),positive="1")
        specORG<-specificity(data=factor(as.numeric(ModelOrg[,2]>FitOrg$splits[4])),reference=factor(binResponseOrg),negative="0")
        
        #return the difference (i.e optimism) in sensitivity and specificity between the bootsrapped
        #and original sample
        return(c(sens-sensORG,spec-specORG))
      }
    }
}


bootOptimism<-function(modelFrame=NULL,binResponse=NULL)
{  
  originalPerf<-modelFitting(ModelOrg=modelFrame,binResponseOrg=binResponse)

  optimism<-c()
  
  #run 1000 boostrapped iterations
  for(i in 1:10000)
  {
    #resample the model with replacement
    sampIdx<-sample(x=1:nrow(modelFrame),replace=TRUE,size=nrow(modelFrame))
    optimism<-rbind(optimism,modelFitting(ModelOrg=modelFrame,ModelSamp=modelFrame[sampIdx,],binResponseOrg=binResponse,binResponseSamp=binResponse[sampIdx]))
  }
  
  sens.OPTADJUST<-originalPerf[1]-mean(optimism[,1],na.rm=T)
  spec.OTPADJUST<-originalPerf[2]-mean(optimism[,2],na.rm=T)
  
  #return a vector with : overfitting sensitivity, optimism adjusted sensitivity, overfitted spec, and optimism adjusted spec
  #also return the percentage of iterations that were not NA
  return(c(originalPerf[1],sens.OPTADJUST,originalPerf[2],spec.OTPADJUST,sum(!is.na(optimism[,1]))/10000))
}

