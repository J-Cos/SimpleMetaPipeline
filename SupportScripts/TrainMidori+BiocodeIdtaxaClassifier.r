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
    BiocodeFasta<-"BIOCODE_JAN21_ER_edit.fasta"
    MidoriFasta="MIDORI_UNIQ_NUC_SP_GB247_CO1_RDP.fasta"
    libraryname<-"ARMS_classifier"

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


#4. format inputs
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
    writeXStringSet(seqs, file=file.path(path, "Data", "Classifiers", paste0(libraryname, "_combinedSeqs.fasta")), format="fasta", width=10000)


#5. train classifier 
    IdtaxaClassifierOutputs<-TrainIdtaxaClassifier(seqs=seqs, maxGroupSize=10, maxIterations=3)


#6. save 
    trainingSet<-IdtaxaClassifierOutputs[[1]]
    save(trainingSet, file=file.path(path, "Data", "Classifiers", paste0(libraryname, "_Final_IdtaxaClassifier.Rdata")))

    pdf(file = file.path(path, "Data", "Classifiers", paste0(libraryname, "_IdtaxaClassifierTrainingSet.pdf")))   # The directory you want to save the file in
    plot(IdtaxaClassifierOutputs[[1]])
    dev.off()