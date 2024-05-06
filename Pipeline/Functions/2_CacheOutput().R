CacheOutput<-function(object){
    # saves even if null to ensure that result printing is able to confirm which steps were not chosen to be completed
    saveRDS(object, file=file.path(path, "IntermediateOutputs", paste0(dataname,"_", deparse(substitute(object)),".RDS")))
}