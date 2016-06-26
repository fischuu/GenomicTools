mdr <- function(X, status, fold=2, t=NULL, cv=0, cvr=0.75, top=3, NAasValues=TRUE, fix=NULL, verbose=FALSE){

    if(cv > 0) stop("Cross-validation is currently under development!")
  
    if(!is.null(fix) & cv>0) stop("Cross-validation in combination with fix SNPs is currently not supported.")

    if(is.data.frame(X)){
      X <- as.matrix(X)
      if(verbose) message("X casted from data.frame to matrix.")
    } 
    if(!is.matrix(X)) stop("Matrix required as input for X")
    
    N <- nrow(X)
    cvRes <- list()
    
    if(is.character(fix)){
      fix <- which((colnames(X)==fix)==TRUE)
    }

    ifelse(is.null(fix), fix <- -1, fix <- fix - 1)
    
    if(is.null(t)) t <- table(status)[2]/table(status)[1]
    if(nrow(X)!=length(status)) stop ("nrow(X) / length(status) mismatch!\n")

    res <- mdr.C(X=X, fold=fold, status=status, t=t, cv=0, cvp=cvr, top=top, na=as.numeric(NAasValues), fix)  

    if(cv>0){
      indices <- 1:nrow(X)
      for(i in 1:cv){ 
          cvSub <- list()
         	trainSet <- sample(indices,floor(length(status)*cvr))
	        testSet <- indices[-trainSet]
          trainModel <- mdr.C(X=X[trainSet,],fold=fold,status=status[trainSet],t=t,cv=0,cvp=cvr,top=top,na=as.numeric(NAasValues),fix)
          tempModel <- list(mdr=trainModel,fold=fold,t=t,cv=0,top=top,fix=fix,X=X,status=status,cv=cvRes)
	        for(foldRun in 1:fold){
          	  trainSet <- mdrEnsemble(tempModel,data=X[trainSet,],new.status=status[trainSet],fold=foldRun)$cv[[foldRun]]
              testSet <- mdrEnsemble(tempModel,data=X[testSet,],new.status=status[testSet],fold=foldRun)$cv[[foldRun]]
	            cvSub[[foldRun]] <- list(train=trainSet,test=testSet)
	        }
        cvRes[[i]] <- cvSub
      }
    }

    result <- list(mdr=res,fold=fold,t=t,cv=cv,top=top,fix=fix,X=X,status=status,cvRes=cvRes)

    class(result) <- "mdr"
    result
} 

