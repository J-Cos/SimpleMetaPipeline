GetTrainingSet<-function(Type, RefLibrary){

    if (Type== "Load") {

        trainingSet<-get(load(file.path(path, "Data", RefLibrary)))
        print("Training Set loaded")
        return(trainingSet)

    } else {
        print("No training set created or loaded")
        return(NULL)
    } 
}