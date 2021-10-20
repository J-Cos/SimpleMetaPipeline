GetTrainingSet<-function(Type, RefLibrary){

    if (Type== "Load") {

        trainingSet<-get(load(file.path(path, "Data", RefLibrary)))

        print("Training Set loaded")
    } else if (Type=="Create"){

        #do flexible training of classifier based on range of reference library inputs

    }
    return(trainingSet)
}