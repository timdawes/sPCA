install.packages("survival")
install.packages("survivalROC")
install.packages("superpc")
install.packages("maxstat")

library(survival)
library(survivalROC)
library(maxstat)
library(superpc)

 
    optimal.point.by.sens.and.spec<- function (x) {m<- which.min(sqrt((x$FP)^2 + (1-x$TP)^2)); return(list(x.coord=x$FP[m], y.coord= x$TP[m], cutoff=x$cut.values[m]))}
 
     median.survival.times<- function(fit) {
       if (length(fit$strata)>1) {
         c<- c(0,cumsum(fit$strata))
         n<- length(fit$n)
         f<- rep(0,n)
             for (i in 1:n) {
               times<- (c[i]+1):c[i+1]
               y<- fit$surv[times]
               t1<- fit$time[times]
               times<- times[which(y>0 & t1>0)]
               if (length(times)>0) {
                                     y<- fit$surv[times]
                                     t1<- fit$time[times]
                                     if (min(y)!=max(y)) {model<- lm(log(y) ~ t1); f[i]<- (log(2) + model$coef[1]) / -model$coef[2]} else {f[i]<- NA}
               }
             }
      return(range(na.omit(f)))} else {return(rep(NA,2))}}


      
    first<- 38
    last<- 29335
    load("Data.RData")
    loop.data<- matrix(0,nrow=nrow(data), ncol=ncol(data))
         
    no.perms<- 1000
    no.subjects<- 256
    k.fold<- 8

    score<- matrix(0,nrow=(no.perms+1),ncol=6)
    coxph.features.survival<- matrix(0, nrow=8, ncol=29298)

    binary<- rep(0,no.subjects)

    fit.groups.discrete<- rep(0,no.subjects)
    fit.groups.continuous<- rep(0,no.subjects)

    optimal.split.point<- matrix(0,nrow=no.perms+1, ncol=8)
    cox.model.2d.AT<- list()
    cox.model.3d.AT<- list()

    cv.threshold<- matrix(0,nrow=no.perms+1, ncol=8)
  
    AUC<- matrix(0,nrow=(no.perms+1), ncol=2, dimnames=list(1:(no.perms+1),c("1D","3D")))


# Change this according to what 1D predictor you want to use
    OneDPredictor<- c("RVEF")
    if (OneDPredictor=="SMWTDIST") {data<- data[which(data$SMWTDIST!=999),]}
    Predictor<- data$RVEF
    if (match(OneDPredictor,c("CICATH","PVRI","RAP","RVEF","SMWTDIST","SVoverESV","RVEDVI","RVESVI")) %in% c(1,4,5,6)) {Predictor<- max(Predictor) - Predictor}
    options(warn=-1)



  for (loops in 1:no.perms) {
   
    # Original data on loop 1, otherwise permute it
    if (loops==1) {loop.data<- data.frame(data)}
    if (loops>1) {cat("\n \n randomising...");
                  loop.data<- data.frame(sapply(1:last, function(x) {data[sample(1:no.subjects), x]}));
                  colnames(loop.data)<- colnames(data)} 
       
        for (i in 1:8) {
        
            test.subjects<- ((i-1)*(no.subjects/8)+1):(i*(no.subjects/8))
            train.subjects<- setdiff(1:no.subjects,test.subjects)
            cv.subjects<- train.subjects
            
            train.set<- loop.data[train.subjects,]
            cv.set<- loop.data[cv.subjects,]
            test.set<- loop.data[test.subjects,]
          
            # 3D Prediction: SPCA training and testing
                    
                          data.train<- list(x=t(train.set[,first:last]),y=train.set[,5],censoring.status=train.set[,3],featurenames=paste("v",1:29298,sep=""))
                          data.cv<- list(x=t(cv.set[,first:last]),y=cv.set[,5],censoring.status=cv.set[,3],featurenames=paste("v",1:29298,sep=""))
                          data.test<- list(x=t(test.set[,first:last]),y=test.set[,5],censoring.status=test.set[,3],featurenames=paste("v",1:29298,sep=""))
          
                      # Train the PCA
                          cat(", training...")
                          train.obj<- superpc.train(data.train, type="survival")
                          if (loops==1) {coxph.features.survival[i,]<- train.obj$feature}
                          
                      # Cross-validate the model to optimize the threshold
                          cat("cross validating...")
                          cv.obj<- superpc.cv(train.obj, n.threshold=10, folds=1, n.components=1, data.cv)
                          cv.threshold[loops, i]<- cv.obj$thresholds[which.max(cv.obj$scor[1,])]
        
                      # Test the SPCA
                          # Discrete outcomes (to compare with 1D measures later)
                                fit.groups.discrete[test.subjects]<- superpc.predict(train.obj, data.train, data.test,
                                    threshold=cv.threshold[loops,i], n.components=1, prediction.type="discrete")$v.pred.1df
              
            # Single-dimension prediction: (such as from RVEF, CO, PVRi etc)
                    
                                cat("\n Optimising ",OneDPredictor,"'s cutoff...",sep="")
                         
                                # Find the cutoff which maximises sensitivity and specificity in the training group
                                    AUC.1D.train<- survivalROC(Stime=train.set[,5], status=train.set[,3], marker=Predictor[train.subjects], predict.time=365, span=0.2)
                                    optimised<- optimal.point.by.sens.and.spec(AUC.1D.train)
                                    optimal.split.point[loops,i]<- optimised$cutoff
                                    binary[which(Predictor[test.subjects]>optimal.split.point[loops,i])]<- 1
                                    binary[which(Predictor[test.subjects]<=optimal.split.point[loops,i])]<- 0
                                    KM.function.1d<- survfit(Surv(data.test$y,data.test$censoring.status)~binary[test.subjects])
                                    
              
        } # End of loop
    
    
    # Now calculate the strength of prediction from all the data
    
            # Collect the data needed for 2D and 3D PH models
            cat("\n Forming CoxPH model...",sep="")
            survival.data.test<- data.frame(cbind(loop.data$Days, loop.data$Outcome, fit.groups.discrete, binary,
                                                  loop.data$AGEATCMR, loop.data$SEX, loop.data$TYPEOFPH,
                                                  loop.data$SMWTDIST, loop.data$DIAGFC, loop.data$RVEDP, loop.data$CI,
                                                  loop.data$MPAP, loop.data$PVR, loop.data$RAP,
                                                  loop.data$RVEDVI, loop.data$RVESVI, loop.data$SVoverESV))
            
            colnames(survival.data.test)<- c("Days","Outcome","ThreeDPredictorPC1","BinaryPredictor","AGEATCMR","SEX","TYPEOFPH","SMWTDIST","DIAGFC",
                                             "RVEDP","CI","MPAP","PVR","RAP","RVEDVI","RVESVI","SVoverESV")
            
            # ROC AUCs for For 1D and 3D data (respectively)
                AUC.1D.test<- survivalROC(Stime=survival.data.test$Days, status=survival.data.test$Outcome, marker=survival.data.test$BinaryPredictor, predict.time=365, method="KM")
                AUC[loops,1]<- AUC.1D.test$AUC
                AUC.3D.test<- survivalROC(Stime=survival.data.test$Days, status=survival.data.test$Outcome, marker=survival.data.test$ThreeDPredictorPC1, predict.time=365, method="KM")
                AUC[loops,2]<- AUC.3D.test$AUC
                            
            # Calculate Cox PH from binary status    
                cox.model.2d.AT[[loops]]<- coxph(Surv(Days,Outcome) ~ AGEATCMR + TYPEOFPH + BinaryPredictor, data=survival.data.test)
                
            # Calculate Cox PH from the 3D status
                cox.model.3d.AT[[loops]]<- coxph(Surv(Days,Outcome) ~ AGEATCMR + TYPEOFPH + ThreeDPredictorPC1, data=survival.data.test)
                            
            score[loops,1:2]<- c(cox.model.2d.AT[[loops]]$scor, cox.model.3d.AT[[loops]]$score)
            score[loops,3:4]<- median.survival.times(survfit(Surv(Days,Outcome) ~ BinaryPredictor, data=survival.data.test))
            score[loops,5:6]<- median.survival.times(survfit(Surv(Days,Outcome) ~ BinaryPredictor + ThreeDPredictorPC1, data=survival.data.test))
               
}



