mdr <- function(X, status, fold=2, t=NULL, top=3, NAasValues=TRUE, fix=NULL, verbose=FALSE){

    if(is.data.frame(X)){
      X <- as.matrix(X)
      if(verbose) message("X casted from data.frame to matrix.")
    } 
    if(!is.matrix(X)) stop("Matrix required as input for X")
    
    N <- nrow(X)
    
    if(is.character(fix)){
      fix <- which((colnames(X)==fix)==TRUE)
    }

    ifelse(is.null(fix), fix <- -1, fix <- fix - 1)
    
    if(is.null(t)) t <- table(status)[2]/table(status)[1]
    if(nrow(X)!=length(status)) stop ("nrow(X) / length(status) mismatch!\n")

    res <- mdr.C(X=X, fold=fold, status=status, t=t, cv=0, cvp=1, top=top, na=as.numeric(NAasValues), fix)  


    result <- list(mdr=res,fold=fold,t=t,top=top,fix=fix,X=X,status=status)

    class(result) <- "mdr"
    result
} 

