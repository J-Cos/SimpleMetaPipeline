        createOutputFilePaths<- function(SampleNames=parent.frame()$SampleNames, dataname=parent.frame()$dataname, suffix, outputDirectoryPrefix){
            output <- file.path(path, "IntermediateOutputs", paste0(outputDirectoryPrefix) ,paste0(dataname, SampleNames, suffix)) #parent.frame is required here as otherwise function searches for variable in global env (where it is defined), instead of the RunDADA2 env.
            names(output) <- SampleNames
            return(output)
        }