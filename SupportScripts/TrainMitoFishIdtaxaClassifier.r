# SCRIPT
    # Train IDtaxa classifier from MitoFish fasta

#0 cluster settings
        #on HPC?
        HPC<-FALSE
        if (HPC==TRUE)   {setwd("/rds/general/user/jcw120/home/BioinformaticPipeline_Env")} #necessary as it appears different job classes have different WDs.

        #CRAN mirror
            r = getOption("repos")
            r["CRAN"] = "http://cran.us.r-project.org"
            options(repos = r)



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
 
# Functions
        #function_files<-list.files(file.path(path, "BioinformaticPipeline", "SupportFunctions"))
        #sapply(file.path(path, "BioinformaticPipeline",  "SupportFunctions" ,function_files),source)

    TrainIdtaxaClassifier<-function(seqs, taxid=NULL, maxGroupSize, maxIterations, allowGroupRemoval=FALSE) {
        
        #get some group objects to manipulate
        groups<-names(seqs)
        groupCounts <- table(names(seqs))
        u_groups <- names(groupCounts) # unique groups

        #set maximum number of sequences assigned to each taxa 
        maxGroupSize <- maxGroupSize # max sequences per label (>= 1)
        
        #identify seqs for removal based on this max group size
        remove <- logical(length(seqs))
        for (i in which(groupCounts > maxGroupSize)) {
            index <- which(groups==u_groups[i])
            keep <- sample(length(index),
                            maxGroupSize)
            remove[index[-keep]] <- TRUE   
        }

        #classifier training loop
        probSeqsPrev <- integer() # suspected problem sequences from prior iteration
        for (i in seq_len(maxIterations)) {
            cat("Training iteration: ", i, "\n", sep="")
            
            # train the classifier
            trainingSet <- LearnTaxa(train=seqs[!remove],
                                    taxonomy=names(seqs)[!remove],
                                    verbose=TRUE)
            
            #something prevents this line working on the HPC - it does work locally
            #save(trainingSet, file=file.path(parent.frame()$path, "Data", "Classifiers", paste0(parent.frame()$libraryname, "_IdtaxaClassifier_Iteration:", i, ".Rdata")))
                                    
            # look for problem sequences
            probSeqs <- trainingSet$problemSequences$Index
            if (length(probSeqs)==0) {
                cat("No problem sequences remaining.\n")
                break
            } else if (length(probSeqs)==length(probSeqsPrev) && all(probSeqsPrev==probSeqs)) {
                cat("Iterations converged.\n")
                break
            }
            if (i==maxIterations)
                break
            probSeqsPrev <- probSeqs
            
            # remove any problem sequences
            index <- which(!remove)[probSeqs]
            remove[index] <- TRUE # remove all problem sequences
            if (!allowGroupRemoval) {
                # replace any removed groups
                missing <- !(u_groups %in% groups[!remove])
                missing <- u_groups[missing]
                if (length(missing) > 0) {
                    index <- index[groups[index] %in% missing]
                    remove[index] <- FALSE # don't remove
                }
            }
        }

        return(list(trainingSet, probSeqs))
    }
    ReformatMitoFishTaxa <- function(list_item) {
        string<-lapply(strsplit(list_item, split= "\t"), '[', 2)
        string<-str_replace(string, "cellularOrganisms", "Root")
        return(string)
    }

    MakeTaxDf <- function(seqs){
        splits<-strsplit(names(seqs), ";")
        Tax_df<-plyr::ldply(splits, rbind )
        return(Tax_df)
    }

#4. format inputs
    #MitoFish
        seqs <- readDNAStringSet(file.path(path, "Data", "Raw", "MitoFish.fasta") )
        tax<-as.list(names(seqs))
        names(seqs)<-unlist(lapply(tax, ReformatMitoFishTaxa))
            print("MitoFish seqs prepared")

    #format  seqs
        #sequence formatting
        seqs <- RemoveGaps(seqs)
            print("gaps removed")

        seqs <- OrientNucleotides(seqs)
            print("nucleotides reoriented")

#6. train classifier 
    IdtaxaClassifierOutputs<-TrainIdtaxaClassifier(seqs=seqs, maxGroupSize=Inf, maxIterations=3)


#7. save 
    trainingSet<-IdtaxaClassifierOutputs[[1]]
    save(trainingSet, file=file.path(path, "Data", "Classifiers", paste0("MitoFish_IdtaxaClassifier.Rdata")))

    pdf(file = file.path(path, "Data", "Classifiers", paste0("MitoFish_IdtaxaClassifierTrainingSet.pdf")))   # The directory you want to save the file in
    plot(IdtaxaClassifierOutputs[[1]])
    dev.off()