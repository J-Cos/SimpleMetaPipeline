#settings
    # general settings
        path <-"../../BioinformaticPipeline_Env"
            # this should be the path to the working directory within which you have the following folders:
                # BioinformaticPipeline - get this from github and then create your own control file from this template - don't modify anything else
                # FASTQs - fill this file with your unmerged multiplexed raw FASTQ files
                # ReferenceLibraries - fill this file with your taxonomic reference library formatted either as i) ... or ii) ...
                # IntermediateOutputs - this will be populated by the pipeline as it runs, it will enable the pipeline to be run over multiple sessions as the output from each module is saved here.
                # Results - this is where final results will be saved
        dataname="23s_test"
            # this should be the name you associate with this set of fastqs,  the inout fastq folder should be labelled with this
            # and all output and result files will be labelled with this name.
        multithread=TRUE
            #on windows set multithread=FALSE

    # Cutadapt settings
        UseCutadapt=TRUE
            # either TRUE or FALSE: determines whether cutadapt is run
        FWD ="ACCTGCGGARGGATCA"
            ## CHANGE ME to your forward primer sequence
        REV = "GAGATCCRTTGYTRAAAGTT",  
            ## CHANGE ME...


    # DADA2 settings
        truncLen=c(249,212)
            #the length at which the sequences should be truncated (forward and reverse)
        trimLeft=c(20,20)
            # the # of bases to remove from the start of the sequences
        maxN=0
            # after truncation seqs with more than this number of Ns are discarded (DADA2 does not allow Ns)
        maxEE=c(2,2)
            #  After truncation, reads with higher than ‘maxEE’ "expected errors"  
            # will be discarded. Expected errors are calculated from the nominal definition 
            # of the quality score: EE = sum(10^(-Q/10))
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

     # cluster settings   
        linkage = "complete"
             #either complete (vsearch) or single (swarm)
        differences=1
            #number of base differences at which swarm clustering will be performed (1=default)                               
        SimilarityThreshold = 0.97
            #%age similarity at which to cluster sequences as a decimal
        threads=4
            # number of threads available    
            
    # lulu settings
        MatchRate=84 #as a %, default 84
            # % matching bases to consider clustering OTUs if co-occurence seen. 
        MinRelativeCo = 0.95 #as a decimal, default 0.95
        #minimum_relative_cooccurence: minimum co-occurrence rate – i.e. the
         # lower rate of occurrence of the potential error explained by
         # co-occurrence with the potential parent for considering error
         # state.
        RatioType = "min" # options: "min" and "avg"
        #sets whether a potential error must have lower
        #  abundance than the parent in all samples ‘min’ (default), or
        #  if an error just needs to have lower abundance on average
        #  ‘avg’. Choosing lower abundance on average over globally
        #  lower abundance will greatly increase the number of
        #  designated errors. This option was introduced to make it
        #  possible to account for non-sufficiently clustered
        #  intraspecific variation, but is not generally recommended, as
        #  it will also increase the potential of cluster
        #  well-separated, but co-occuring, sequence similar species.

    # IDTAXA settings
        Type ="No"  
            #whether to "Create" or "Load" a training set, or perform "No Assignment"
        RefLibrary= "microgreen_trainingSet.Rdata" 
            #ref library to load if loading
        SeqsToAssign ="cOTUs"
            #whether to assign to "ESVs", "OTUs", or "cOTUs"
        threshold=60
            # %age confidence of assignment required to record assignment

# run pipeline
    #source(file.path(path, "BioinformaticPipeline", "Pipeline", "DSLI_SQL_Pipeline.R"))
    #source(file.path(path, "BioinformaticPipeline", "Pipeline", "DLSL_Pipeline.R"))
    #source(file.path(path, "BioinformaticPipeline", "Pipeline", "CDLSL_Pipeline.R"))
    source(file.path(path, "BioinformaticPipeline", "Pipeline", "CDLCL_Pipeline.R"))


