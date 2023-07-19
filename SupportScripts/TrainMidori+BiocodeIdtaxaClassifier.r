# SCRIPT
    # Train IDtaxa classifier from BIOCODE fasta

#0 cluster settings
        #on HPC?
        HPC<-TRUE
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
    path<-"../BioinformaticPipeline_Env"
    biocodeVersion<-"JAN21"
    midoriVersion="256"
    libraryname<-"ARMS_classifier"
    maxGroupSize<-3
    maxIterations<-3
    LargeTerrestrialOrders<-c("Araneae", "Hemiptera", "Coleoptera", "Hymenoptera", "Lepidoptera", "Diptera")
    LargestMarineOrder<-"Decapoda"
    subsampleLargeTerrestrialOrders<-"SubsampleTerrestrial" # either 1) "SubsampleTerrestrial" or 2) ""

    BiocodeFasta<-paste0("BIOCODE_",biocodeVersion,"_ER_edit.fasta")
    MidoriFasta=paste0("MIDORI2_UNIQ_NUC_SP_GB", midoriVersion,"_CO1_RDP.fasta")

# Functions
        #function_files<-list.files(file.path(path, "BioinformaticPipeline", "SupportFunctions"))
        #sapply(file.path(path, "BioinformaticPipeline",  "SupportFunctions" ,function_files),source)

    TrainIdtaxaClassifier<-function(seqs, taxid=NULL, maxGroupSize, maxIterations=1, allowGroupRemoval=TRUE, probSeqsPrev) {
        
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
        print(paste0(sum(remove), " sequences removed based on max group size"))

        #classifier training loop
        for (i in seq_len(maxIterations)) {
            
            # train the classifier
            trainingSet <- LearnTaxa(train=seqs[!remove],
                                    taxonomy=names(seqs)[!remove],
                                    verbose=TRUE)
            
            # look for problem sequences
            probSeqs <- trainingSet$problemSequences$Index
            if (length(probSeqs)==0) {
                cat("No problem sequences remaining.\n")
                break
            } else if (length(probSeqs)==length(probSeqsPrev) && all(probSeqsPrev==probSeqs)) {
                cat("Iterations converged.\n")
                break
            }

            
            # remove any problem sequences
            index <- which(!remove)[probSeqs]
            remove[index] <- TRUE # remove all problem sequences
            
            #test if any rgoups entirely removed
            missing <- !(u_groups %in% groups[!remove])
            missing <- u_groups[missing]
            if (!allowGroupRemoval) {
                # replace any removed groups
                if (length(missing) > 0) {
                    index <- index[groups[index] %in% missing]
                    remove[index] <- FALSE # don't remove
                }
            } else if (allowGroupRemoval){
                print(paste0(length(missing), " groups removed"))
            }
        }
                
        return(list("trainingSet"=trainingSet, "seqs_filtered"=seqs[!remove], "probSeqs"=probSeqs))
    }
    ReformatMidoriTaxa <- function(list_item) {
        string<-lapply(strsplit(list_item, split= "\t"), '[', 2)
        ranks<-strsplit(unlist(string), split= ";")
        ranks_split<-strsplit(unlist(ranks), split= "_")
        ranks_split_nonum<-list()
        for (i in 1:length(ranks_split)) {
            ranks_split_nonum[[i]]<-ranks_split[[i]][-length(ranks_split[[i]])]
            ranks_split_nonum[[i]]<-paste(ranks_split_nonum[[i]], collapse="_")
        }
        new_list_item<-paste(unlist(ranks_split_nonum), collapse=";")
        #capitalise root!
        new_list_item<-sub("^.", "R", new_list_item)

        return(new_list_item)
    }
    ReformatBiocodeTaxa <- function(list_item) {
        taxstring<-lapply(strsplit(list_item, split= "\\|"), '[', 2)
        #remove semicolons from names so they don't interfere downstream (replace with underscore)
        string_no_semicolons<-gsub('[;]', '_', taxstring)
        ranks<-strsplit(unlist(string_no_semicolons), split= ",")
        ranks_split<-lapply(strsplit(unlist(ranks), split= ": "), '[',2)
        if (length(ranks_split)==6) { #6 because the 6th rank of biocode is the species
            ranks_split[length(ranks_split)]<-paste0( ranks_split[(length(ranks_split)-1)]," ", ranks_split[length(ranks_split)])
        }
        list_item_noroot<-paste0(unlist(ranks_split), collapse=";")
        new_list_item<-paste0("Root;Eukaryota;",list_item_noroot)
        return(new_list_item)
    }
    MakeTaxDf <- function(seqs){
        splits<-strsplit(names(seqs), ";")
        Tax_df<-plyr::ldply(splits, rbind )
        return(Tax_df)
    }
    CheckIfPreviousIncompleteRunExists<-function() {
        ClassifierOutput<-NULL
        for (iter in (maxIterations-1):1){
            PreviousIteration<-try(readRDS(file.path(path, "Data", "Classifiers", paste0(libraryname, "_IdtaxaClassifier_Iteration_", iter, ".RDS"))), silent=TRUE)
            if(! "try-error" %in% class(PreviousIteration)){
                ClassifierOutput<-PreviousIteration
                break 
            }
            iter<-0
        }
        return(list("ClassifierOutput" = ClassifierOutput, "iter"=iter))
    }

# CHECK IF A PREVIOUS INCOMPLETE RUN EXISTS

    ClassifierOutput<-CheckIfPreviousIncompleteRunExists()
    print(paste0("Starting at iteration ", ClassifierOutput[["iter"]]))

    if(is.null(ClassifierOutput[["ClassifierOutput"]])){
        #4. format inputs
            #specify that no problem sequences are known
                probSeqs<-integer()
            #biocode
                Bseqs <- readDNAStringSet(file.path(path, "Data", "Raw", BiocodeFasta) )
                Btax<-as.list(names(Bseqs))
                names(Bseqs)<-unlist(lapply(Btax, ReformatBiocodeTaxa))
                    print("Biocode seqs prepared")

            #midori
                Mseqs <- Biostrings::readDNAStringSet(file.path(path, "Data", "Raw", MidoriFasta))
                    #Mseqs <- Mseqs[sample(1:length(Mseqs), length(Bseqs))] #line for testing purposes
                Mtax<-as.list(names(Mseqs))
                names(Mseqs)<-unlist(lapply(Mtax, ReformatMidoriTaxa))
                        print("Midori names reformatted")

            #combine midori and biocode
                seqs<-c(Bseqs, Mseqs)

            #format combined seqs
                #sequence formatting
                seqs <- RemoveGaps(seqs)
                    print("gaps removed")

                seqs <- OrientNucleotides(seqs)
                    print("nucleotides reoriented")

            # Checkpoint - save seqs to output for inspection
                writeXStringSet(seqs, file=file.path(path, "Data", "Classifiers", paste0(libraryname, "_m", midoriVersion, "_b", biocodeVersion, "_combinedSeqs.fasta")), format="fasta", width=10000)

        #5. subsample large terrestrial orders
            if (subsampleLargeTerrestrialOrders=="SubsampleTerrestrial"){

                Tax_df<-MakeTaxDf(seqs=seqs)       

                SizeOfLargestMarineOrder<-sum(Tax_df[,5] == LargestMarineOrder, na.rm=TRUE) # na removal discounts seqs of unknown order (e.g. a sequence only labelled "arthropoda")

                TerrestrialSeqsToRemove<-c()
                for (i in 1: length(LargeTerrestrialOrders)){
                    TerrestrialOrderSeqs<-which(Tax_df[,5] == LargeTerrestrialOrders[i])
                    TerrestrialOrderSeqsToRemove<-sample(TerrestrialOrderSeqs, length(TerrestrialOrderSeqs)-SizeOfLargestMarineOrder, replace=FALSE)
                    TerrestrialSeqsToRemove<-c(TerrestrialSeqsToRemove, TerrestrialOrderSeqsToRemove)
                }
                seqs<-seqs[-TerrestrialSeqsToRemove,]

            # Checkpoint - save seqs to output for inspection
                writeXStringSet(seqs, file=file.path(path, "Data", "Classifiers", paste0(libraryname,"_",subsampleLargeTerrestrialOrders, "_m", midoriVersion, "_b", biocodeVersion,  "_combinedSeqs.fasta")), format="fasta", width=10000)

            }
    } else {
        seqs<-ClassifierOutput[["ClassifierOutput"]][["seqs_filtered"]]
        probSeqs<-ClassifierOutput[["ClassifierOutput"]][["probSeqs"]]
    }

#6. train classifier
    for (iter in (ClassifierOutput[["iter"]]+1):maxIterations){
        cat("Training iteration: ", iter, "\n", sep="")
        IdtaxaClassifierOutputs<-TrainIdtaxaClassifier(seqs=seqs, maxGroupSize=maxGroupSize, probSeqsPrev=probSeqs)
        saveRDS(IdtaxaClassifierOutputs, file=file.path(parent.frame()$path, "Data", "Classifiers", paste0(parent.frame()$libraryname, "_IdtaxaClassifier_Iteration_", iter, ".RDS")))
        print(paste0(length(IdtaxaClassifierOutputs[["probSeqs"]]) , " problem sequences removed"))
        print(paste0(length(IdtaxaClassifierOutputs[["seqs_filtered"]]), " sequences remaining"))
        if (length(IdtaxaClassifierOutputs[["probSeqs"]])==0){break}
        seqs<-IdtaxaClassifierOutputs[["seqs_filtered"]]
        probSeqs<-IdtaxaClassifierOutputs[["probSeqs"]]
    }   

#7. save 
    trainingSet<-IdtaxaClassifierOutputs[[1]]
    save(trainingSet, file=file.path(path, "Data", "Classifiers", paste0(libraryname,"_",subsampleLargeTerrestrialOrders, "_m", midoriVersion, "_b", biocodeVersion, "_Final_IdtaxaClassifier.Rdata")))

    pdf(file = file.path(path, "Data", "Classifiers", paste0(libraryname,"_",subsampleLargeTerrestrialOrders, "_m", midoriVersion, "_b", biocodeVersion, "_IdtaxaClassifierTrainingSet.pdf")))   # The directory you want to save the file in
    plot(IdtaxaClassifierOutputs[[1]])
    dev.off()