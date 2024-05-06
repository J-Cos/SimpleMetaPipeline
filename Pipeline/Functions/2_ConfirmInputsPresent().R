ConfirmInputsPresent<-function(object){
    if (!exists(object)) {
        assign(value=readRDS(file=file.path(path, "IntermediateOutputs", paste0(dataname,"_",object ,".RDS"))), x=object, envir=parent.frame())
    } else { 
        print(paste0(object, " already present")) 
    }
}