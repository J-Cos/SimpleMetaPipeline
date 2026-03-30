# SCRIPT
    # Train IDtaxa classifier from BIOCODE fasta


#1 set seed for replication
    set.seed(0.1)

#2 dependencies    
    #installed software
        #vsearch
        #swarm v2
    #R packages
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load( tidyverse,
                    seqinr,
                    DECIPHER,
                    Biostrings
                    )


#3. parameters
    path<-"../../BioinformaticPipeline_Env"
    BiocodeFasta<-"BIOCODE_JAN21_ER_edit.fasta"
    libraryname<-"Biocode_ER"
    #Functions
        function_files<-list.files(file.path(path, "BioinformaticPipeline", "SupportFunctions"))
        sapply(file.path(path, "BioinformaticPipeline",  "SupportFunctions" ,function_files),source)

#4. format inputs
    seqs <- readDNAStringSet(file.path(path, "Data", "Raw", BiocodeFasta) )

    #taxonomy formatting
        Taxa_df<-CreateTaxaDfFromBiocodeFasta(path=file.path(path, "Data", "Raw"), BiocodeFasta=BiocodeFasta)
        Taxonomy_df<-Taxa_df[,-(1:2)]
        df_args <- c(Taxonomy_df, sep=";")
        Taxonomy_noroot<-do.call(paste, df_args)
        Taxonomy<- paste0("Root; Eukaryota;",Taxonomy_noroot)
        names(seqs)<-Taxonomy

    #sequence formatting
    seqs <- RemoveGaps(seqs)
    seqs <- OrientNucleotides(seqs)

#5. train classifier 
    IdtaxaClassifierOutputs<-TrainIdtaxaClassifier(seqs=seqs, maxGroupSize=10, maxIterations=100)



#6. save 
    trainingSet<-IdtaxaClassifierOutputs[[1]]
    save(trainingSet, file=file.path(path, "Data", "Classifiers", paste0(libraryname, "_IdtaxaClassifier.Rdata")))

    pdf(file = file.path(path, "Data", "Classifiers", paste0(libraryname, "_IdtaxaClassifierTrainingSet.pdf")))   # The directory you want to save the file in
    plot(IdtaxaClassifierOutputs[[1]])
    dev.off()