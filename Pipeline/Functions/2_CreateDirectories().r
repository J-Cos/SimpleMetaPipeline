CreateDirectories<-function() { 
    dir.create((file.path(path, "FASTQs")))
    dir.create((file.path(path, "IntermediateOutputs")))
    dir.create((file.path(path, "Results")))
    dir.create((file.path(path, "Data")))
}
