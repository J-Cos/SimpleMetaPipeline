GetTrainingSet<-function(Type, RefLibrary){

    if (Type== "Load") {

        trainingSet<-get(load(file.path(path, "Data", RefLibrary)))
        print("Training Set loaded")
        return(trainingSet)

    } 
    if (Type=="Create"){

        print("Creating Training Set")
        #do flexible training of classifier based on range of reference library inputs
        print("Training Set created")
        return(trainingSet)

    } else {
        print("No training set created or loaded")
        return(NULL)
    } 
}