CacheOutput<-function(object){
    if (!is.null(object)) {
        saveRDS(object, file=file.path(path, "IntermediateOutputs", paste0(dataname,"_", deparse(substitute(object)),".RDS")))
    }
}