        loadFastq<-function(pattern){
            output<-sort(list.files(file.path(path, "FASTQs"), pattern=pattern, full.names = TRUE))
            return(output)
        }