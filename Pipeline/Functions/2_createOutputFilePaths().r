#N.b. parent.frame is required here as otherwise function searches for variable in global env (where it is defined), instead of the RunDADA2 env.
createOutputFilePaths<- function(SampleNames=parent.frame()$SampleNames, Run=parent.frame()$Run, dataname=parent.frame()$dataname, suffix, outputDirectoryPrefix){
    output <- file.path(path, "IntermediateOutputs", paste0(dataname, outputDirectoryPrefix) ,paste0("Run",Run), paste0(dataname, SampleNames, suffix))
    names(output) <- SampleNames
    return(output)
}