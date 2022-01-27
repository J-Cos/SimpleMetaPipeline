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
            
            save(trainingSet, file=file.path(path, "Data", "Classifiers", paste0(libraryname, "_IdtaxaClassifier_Iteration:", i, ".Rdata")))

                                    
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
    CreateTaxaDfFromBiocodeFasta<-function(path, BiocodeFasta) {
        
        #readfasta as strings
        BiocodeFasta <- readLines(file.path(path,BiocodeFasta))

        # if odd number of rows then the last is a fragement and should be removed
        if (length(BiocodeFasta)%%2==1) {
            BiocodeFasta<-BiocodeFasta[1:length(BiocodeFasta)-1]
        }


        BiocodeFastaNames<-BiocodeFasta[seq(from=1, to=length(BiocodeFasta), by=2)]
        
        Tax_df<-as.data.frame( matrix(NA, ncol =8 , nrow = length(BiocodeFastaNames)) )
        names(Tax_df)<-c(    "NumID", "ID", "Phylum", "Class", "Order", "Family", "Genus", "Species")
        
        for (i in 1:length(BiocodeFastaNames)) {

            string<-sub('.', '', BiocodeFastaNames[i])
            string<-strsplit(string, split= "\\|")
            string<-unlist(string)

            Tax_df$NumID[i]<-string[[3]]
            Tax_df$ID[i]<-string[1]

            x<-unlist(strsplit(string[2], split=","))
            x<-unlist(lapply(FUN=strsplit, x, ":", 2))
            ranks<-rep(NA, 6)
            for (j in 1: length(ranks)) {
                ranks[j]<-x[ (j*2) ]
            }

            Tax_df[i,3:8]<-ranks
        }

        return(Tax_df)
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
        new_list_item<-paste(unlist(ranks_split_nonum), collapse="; ")
        #capitalise root!
        new_list_item<-sub("^.", "R", new_list_item)

        return(new_list_item)
    }



#4. format inputs
    #biocode
        Bseqs <- readDNAStringSet(file.path(path, "Data", "Raw", BiocodeFasta) )

        #taxonomy formatting
            Taxa_df<-CreateTaxaDfFromBiocodeFasta(path=file.path(path, "Data", "Raw"), BiocodeFasta=BiocodeFasta)
            Taxonomy_df<-Taxa_df[,-(1:2)]
            df_args <- c(Taxonomy_df, sep=";")
            Taxonomy_noroot<-do.call(paste, df_args)
            Taxonomy<- paste0("Root; Eukaryota:",Taxonomy_noroot)
            names(Bseqs)<-Taxonomy
                print("Biocode seqs prepared")

    #midori
        Mseqs <- Biostrings::readDNAStringSet(file.path(path, "Data", "Raw", MidoriFasta))
        
        #taxonomy formatting
            tax<-as.list(names(Mseqs))
            names(Mseqs)<-unlist(lapply(tax, ReformatMidoriTaxa))
                print("Midori names reformatted")

    #combine midori and biocode
        seqs<-c(Bseqs, Mseqs)

    #format combined seqs
        #sequence formatting
        seqs <- RemoveGaps(seqs)
            print("gaps removed")

        seqs <- OrientNucleotides(seqs)
            print("nucleotides reoriented")

    
#5. train classifier 
    IdtaxaClassifierOutputs<-TrainIdtaxaClassifier(seqs=seqs, maxGroupSize=10, maxIterations=3)


#6. save 
    trainingSet<-IdtaxaClassifierOutputs[[1]]
    save(trainingSet, file=file.path(path, "Data", "Classifiers", paste0(libraryname, "_Final_IdtaxaClassifier.Rdata")))

    pdf(file = file.path(path, "Data", "Classifiers", paste0(libraryname, "_IdtaxaClassifierTrainingSet.pdf")))   # The directory you want to save the file in
    plot(IdtaxaClassifierOutputs[[1]])
    dev.off()