#Control Script Template
    # general settings
        path <-""
            # If running the pipeline as standard (starting in the BioinformaticPipeline_Env directory) leave this blank. This designates the path to the working directory within which you have the following folders:
                # BioinformaticPipeline - get this from github and then create your own control file from this template - don't modify anything else
                # FASTQs -directory with your unmerged multiplexed raw FASTQ files, in a subdirectory named to match dataname argument
                # ReferenceLibraries - directory with your taxonomic reference library formatted as i) IDTAXA classsifier and/or ii) BLAST database
                # IntermediateOutputs - this directory will be populated by the pipeline as it runs, it will enable the pipeline to be run over multiple sessions as the output from each module is saved here.
                # Results - this is where final results will be saved
        dataPath<-NULL
            # if NULL (default) then the fastqs will be found in the FASTQs directory of the BioinformaticPipeline_Env directory. Normally there is no reason to modify this.
            # In the case of large datasets in external storage or research data stores on HPCs you can specify the path to the fastqs here. Note that all other components of the 
            # pipeline will work as normal and all outputs will be written to the BioinformaticPipeline_Env/IntermediateOutputs and BioinformaticPipeline_Env/Results.
        dataname="Example"
            # this should be the name you associate with this set of fastqs,  the inout fastq folder should be labelled with this
            # and all output and result files will be labelled with this name.
        multithread=4
            # number of threads to use on parallelised steps within the pipeline

    # Cutadapt settings https://cutadapt.readthedocs.io/en/stable/
        UseCutadapt=FALSE
            # either TRUE or FALSE: determines whether cutadapt is run
        FWD ="ACTG"
            ## your forward primer sequence
        REV = "ACTG"
            ## your reverse primer sequence

    # DADA2 settings https://benjjneb.github.io/dada2/index.html
        NumberOfRuns<-2
            #Number of sequencing runs used to generate the data for this experiment, e.g. 1 or 2
        truncLen=list(c(150,150), c(150,150))
            #the lengths at which the sequences should be truncated (forward and reverse) for each run. If single-end or pre-merged 
            # sequences are input only a single number is required for each entry in the list.
        trimLeft=list(c(0,0), c(0,0))
            # the # of bases to remove from the start of the sequences (forward and reverse) for each run. If single-end or pre-merged 
            # sequences are input only a single number is required for each entry in the list.
        maxN=0
            # after truncation seqs with more than this number of Ns are discarded (DADA2 does not allow Ns so this must be kept to 0 if DADA2 outputs are desired)
        maxEE=c(2,2)
            #  After truncation, reads with higher than ‘maxEE’ "expected errors"  
            # will be discarded (forward and reverse). Expected errors are calculated from the nominal definition 
            # of the quality score: EE = sum(10^(-Q/10)). If single-end or pre-merged 
            # sequences are input only a single number is required.
        truncQ=2
            # Truncate reads at the first instance of a quality score less than or equal to ‘truncQ’. 
        DesiredSequenceLengthRange=NULL
            #sequence length range to keep enter as e.g. 360:270, if NULL all sequence lengths are kept. 
            # Sequences that are much longer or shorter than expected
            # may be the result of non-specific priming. This removal is analogous to “cutting a band” 
            # in-silico to get amplicons of the targeted length. 
        pool="pseudo"
            #TRUE, FALSE, or pseudo. pseudo pooling approximates the effect of denoising with pooled samples, but with
            # linearly increasing computational time (ca. doubled compared to no pooling)
        MixedOrientation="FALSE"
            # TRUE or FALSE. FALSE (default) means dada2 proceeds in the normal fashion. TRUE means that _RO and _FO 
            # suffixes will be looked for in the sample names to enable samples containing reverse orientation (RO) sequences 
            # to be reverse complemented and then merged into the matching sample containing the forward orientated sequences.
        ReadType="Paired-end"
            # default "Paired-end". Other options: "Single-read" and "Premerged". Premerged refers to paired-end data which has been merged prior to input into the pipeline.
            # This is not recommended but in some historic datasets it is unavoidable. In this case DADA2's inflateErr(inflation=3) function is used to correct distortions
            # to error rates introduced by premerging the paired-end reads.

    # lulu settings1 https://github.com/tobiasgf/lulu
        MatchRate1=90 #as a %, default 84
            # % matching bases to consider clustering OTUs if co-occurence seen. 
        MinRelativeCo1 = 0.95 #as a decimal, default 0.95
        #minimum_relative_cooccurence: minimum co-occurrence rate – i.e. the
         # lower rate of occurrence of the potential error explained by
         # co-occurrence with the potential parent for considering error
         # state.
        RatioType1 = "min" # options: "min" and "avg"
        #sets whether a potential error must have lower
        #  abundance than the parent in all samples ‘min’ (default), or
        #  if an error just needs to have lower abundance on average
        #  ‘avg’. Choosing lower abundance on average over globally
        #  lower abundance will greatly increase the number of
        #  designated errors. This option was introduced to make it
        #  possible to account for non-sufficiently clustered
        #  intraspecific variation, but is not generally recommended, as
        #  it will also increase the potential of clustering
        #  well-separated, but co-occuring, sequence similar species. 

     # cluster settings https://github.com/torognes/vsearch https://github.com/torognes/swarm
        linkage = "complete"
             #either complete (vsearch) or single (swarm)
        differences=1
            #number of base differences at which swarm clustering will be performed (1=default)                               
        SimilarityThreshold = 0.97
            #%age similarity at which to cluster sequences as a decimal
        threads=20
            # number of threads available 

    # lulu settings2 https://github.com/tobiasgf/lulu
        MatchRate2=90 #as a %, default 84
            # % matching bases to consider clustering OTUs if co-occurence seen. 
        MinRelativeCo2 = 0.95 #as a decimal, default 0.95
        #minimum_relative_cooccurence: minimum co-occurrence rate – i.e. the
         # lower rate of occurrence of the potential error explained by
         # co-occurrence with the potential parent for considering error
         # state.
        RatioType2 = "min" # options: "min" and "avg"
        #sets whether a potential error must have lower
        #  abundance than the parent in all samples ‘min’ (default), or
        #  if an error just needs to have lower abundance on average
        #  ‘avg’. Choosing lower abundance on average over globally
        #  lower abundance will greatly increase the number of
        #  designated errors. This option was introduced to make it
        #  possible to account for non-sufficiently clustered
        #  intraspecific variation, but is not generally recommended, as
        #  it will also increase the potential of clustering
        #  well-separated, but co-occuring, sequence similar species. 

    # IDTAXA settings https://www.bioconductor.org/packages/release/bioc/vignettes/DECIPHER/inst/doc/ClassifySequences.pdf
        IDTAXA =TRUE
            #whether to run IDTAXA assignment, TRUE or FALSE
        desiredranks<-c("rootrank", "domain", "phylum", "class", "order", "family", "genus")
            #determine this based on the training set you are using (accessible with trainingSet[[3]] once the trainingSet is loaded into R). If this trainingSet
            #does not contain ranks then this is ignored, and it is assumed all reference sequences have an assignment at 
            #all ranks. If this is not the case assignment will be unusable. If the trainingSet does contain ranks then only 
            #the ranks specified here are reported, this enables easy comparison between assignments as all sequences receive the same number of taxxonomic ranks
            # . e.g. looking at class assignments for all sequences.
            #some examples for widely used reference libraries:
                #SILVA 18s
                    #desiredranks<-c("rootrank", "domain", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "infraphylum", "superclass", "class", "subclass", "infraclass", "superorder",  "order", "suborder", "superfamily", "family", "subfamily", "genus")
                #PR2 18s
                    #desiredranks<-c("rootrank", "kingdom", "division", "phylum", "class", "order", "family", "genus")
                #GTDB 16s
                    #desiredranks<-c("rootrank", "domain", "phylum", "class", "order", "family", "genus")
        trainingSet= "GTDB_r207-mod_April2022.RData"
            #IDTAXA classifier to use, should be the file name for your pretrained classifier. Note these commonly used pre-trained classifiers can be downloaded here:
            # http://www2.decipher.codes/Downloads.html
        SeqsToAssign ="ESVs"
            #whether to make assignments to "ESVs", "OTUs", or "cOTUs". default of assigning to ESVs is recommended to make full use of multi-algorithm agreement.
        threshold=40 
            # %age confidence of assignment required to record assignment
            # 30=low confidence, 40=moderate, 50 = high, 60= very high
        parallel=FALSE
            #whether idtaxa should be run in parallel using the number of cores specififed with the multithread argument.
            #each core will classify 500 sequences at a time using parallel::mclapply. Setting this to TRUE is recommended 
            # on an an HC with more than 10 cores. Note that this is a workaround as Idtaxa's inbuilt parallelisation seems
            # unable to use more than around 10 cores.
    #BLAST https://www.ncbi.nlm.nih.gov/books/NBK1734/
        Blast= FALSE
            #run blast or not, TRUE or FALSE
        dbname= "MidoriBiocodeBlastDB" 
            # name of the blast db you wish to use, should be the directory name of your prebuilt blast database
        assignmentThresholds=c(0, 0, 80, 85, 90, 95, 97)
            #similarity thresholds at which blast assignments should be made. 
            # read from highest rank to lowest rank. Requires a value for each rank in dataset.
        

# run pipeline
    source(file.path(path, "BioinformaticPipeline", "Pipeline", "Pipeline.R"))





