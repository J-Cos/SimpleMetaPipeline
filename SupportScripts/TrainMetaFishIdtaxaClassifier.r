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
                    ape,
                    Biostrings,
                    DECIPHER
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
    ReformatMetaFishTaxa <- function(list_item) {
        #string<-lapply(strsplit(list_item, split= ";"), '[', 2)
        string<-str_replace(list_item, "Chordata", "Root")
        return(string)
    }

    MakeTaxDf <- function(seqs){
        splits<-strsplit(names(seqs), ";")
        Tax_df<-plyr::ldply(splits, rbind )
        return(Tax_df)
    }

#4. get metfish input
        # load REMOTE references and cleaning scripts (requires internet connection)
        source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-load-remote.R")
        source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-clean.R")

        # change 'metabarcode' argument as appropriate:
        reflib.sub <- subset_references(df=reflib.cleaned, metabarcode="12s.miya")

        # [OPTIONAL] taxonomically dereplicate and filter on sequence length
        # 'proplen=0.5' removes sequences shorter than 50% of median sequence length
        # 'proplen=0' retains all sequences
        reflib.sub <- derep_filter(df=reflib.sub, derep=FALSE, proplen=0.5)

        # write out reference library in FASTA format
            # create labels
                df.labs <- reflib.sub %>% 
                    dplyr::mutate(kingdom="Animalia") %>%
                    dplyr::arrange(kingdom,phylum,class,order,family,genus,sciNameValid,dbid) %>% 
                    dplyr::mutate(labelDadaTaxonomy=paste0(phylum,";",class,";",order,";",family,";",genus,";",sciNameValid))
                    # get version
                    gbv <- paste0("v",unique(pull(df.labs,genbankVersion))[1])
                    mbc <- unique(pull(df.labs,metabarcode))
                
            # make filename
                filename.dada.taxonomy <- paste("MetaFish",mbc,"dada.taxonomy",gbv,"fasta",sep=".")

            # convert to fasta file and write out
                ape::write.FASTA(tab2fas(df=df.labs,seqcol="nucleotides",namecol="labelDadaTaxonomy"), file=file.path("../Data/Raw", filename.dada.taxonomy))
                


#5. format inputs
        seqs <- readDNAStringSet(file.path(path, "Data", "Raw", filename.dada.taxonomy) )
        tax<-as.list(names(seqs))
        names(seqs)<-unlist(lapply(tax, ReformatMetaFishTaxa))
            print("MitoFish seqs prepared")

    #format  seqs
        #sequence formatting
        seqs <- RemoveGaps(seqs)
            print("gaps removed")

        seqs <- OrientNucleotides(seqs)
            print("nucleotides reoriented")

#6. train classifier 
    IdtaxaClassifierOutputs<-TrainIdtaxaClassifier(seqs=seqs, maxGroupSize=10, maxIterations=3)


#7. save 
    trainingSet<-IdtaxaClassifierOutputs[[1]]
    save(trainingSet, file=file.path(path, "Data", "Classifiers", paste0("MetaFish_IdtaxaClassifier.Rdata")))

    pdf(file = file.path(path, "Data", "Classifiers", paste0("MetaFish_IdtaxaClassifierTrainingSet.pdf")))   # The directory you want to save the file in
    plot(IdtaxaClassifierOutputs[[1]])
    dev.off()