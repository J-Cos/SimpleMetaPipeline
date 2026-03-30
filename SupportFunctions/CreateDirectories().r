#func to creae directories

CreateDirectories<-function(path) { 
    dir.create((file.path(path, "FASTQs")))
    dir.create((file.path(path, "IntermediateOutputs")))
    dir.create((file.path(path, "Results")))
    dir.create((file.path(path, "Data")))
    dir.create((file.path(path, "Data", "BlastDBs")))
    dir.create((file.path(path, "Data", "Classifiers")))
    dir.create((file.path(path, "Data", "Raw")))


}
