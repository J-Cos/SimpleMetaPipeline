IdentifyLastInputPresent<-function(vector){

    x<-list()

    for (i in 1:length(vector)) {
        object<-vector[i]
        x[[i]]<-try(assign(value=readRDS(file=file.path(path, "IntermediateOutputs", paste0(dataname,"_",object ,".RDS"))), x=object, envir=parent.frame()), silent=TRUE)
        tryresults<-lapply(x, class)  
    }

    tryresults<-lapply(x, class)
    startingstep<-min(which(tryresults=="try-error"))

    return(startingstep)
}
