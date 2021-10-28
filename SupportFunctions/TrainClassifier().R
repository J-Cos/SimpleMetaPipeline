
    seqs_path<-file.path(path, "Data", "microgreen_algaebase_biocomPipe.fasta") #fasta downloaded from algaebase
    seqs <- readDNAStringSet(seqs_path)
    rank_path<-file.path(path, "Data","taxid.txt")
    taxid <- read.table(rank_path, 
                        header=FALSE, 
                        col.names=c('Index','Name','Parent','Level','Rank'),
                        sep="*", # asterisks delimited
                        quote="", # preserve quotes
                        stringsAsFactors=FALSE)

    taxid<-taxid[!(taxid$Rank=="species"),]#remove species level

    #renaming unnamed ranks, should go in make taxid script not here
    #taxid[[2]][str_detect( as.vector(taxid[[2]]), "Unnamed_rank")]<-"Unnamed_rank"

    seqs <- RemoveGaps(seqs)
    seqs <- OrientNucleotides(seqs)

    # obtain the taxonomic assignments
    groups <- names(seqs) # sequence names
    groups<-sapply(regmatches(groups, regexpr(";", groups), invert = TRUE), getElement, 2)

    #change naming of groups to what IDTAXA expects
    groups<-paste0("Root;", groups)
    groups<-gsub("\\s*\\([^\\)]+\\)","",groups) #string manipulation 1
    groups<-gsub("(.*);;.*","\\1",groups) #string manipulation 2
    names(seqs)<-groups

    #get some group objects to manipulate
    groupCounts <- table(groups)
    u_groups <- names(groupCounts) # unique groups

    #set maximum number of sequences assigned to each taxa (in this case infinite as database is small)
    maxGroupSize <- Inf # max sequences per label (>= 1)
    
    #identify seqs for removal based on this max group size
    remove <- logical(length(seqs))
    for (i in which(groupCounts > maxGroupSize)) {
        index <- which(groups==u_groups[i])
        keep <- sample(length(index),
                        maxGroupSize)
        remove[index[-keep]] <- TRUE   
    }

    #classifier training loop
    maxIterations <- 3 # must be >= 1
    allowGroupRemoval <- FALSE
    probSeqsPrev <- integer() # suspected problem sequences from prior iteration
    for (i in seq_len(maxIterations)) {
        cat("Training iteration: ", i, "\n", sep="")
        
        # train the classifier
        trainingSet <- LearnTaxa(train=seqs[!remove],
                                taxonomy=names(seqs)[!remove],
                                rank=taxid)
                                
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




save(trainingSet, file="trainingSet.Rdata")


    #save tables
    IdtaxaTables<-list()
    IdtaxaTables[[1]]<-as.data.frame(length(probSeqs))
    rownames(IdtaxaTables[[1]])<- "Number of Problem Sequences"
    t<-tableGrob(IdtaxaTables)
    ggsave( plot=t, file=file.path(path, "Figs",paste0(dataname,"ClassifierTrainingTables.pdf")) )
    
    pdf(file = file.path(path, "Figs", paste0(dataname, "_TrainingSet.pdf")))   # The directory you want to save the file in
    plot(trainingSet)
    dev.off()
