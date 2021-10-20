#settings
    # general settings
        path <-"../../BioinformaticPipeline_Env"
            # this should be the path to the working directory within which you have the following folders:
                # BioinformaticPipeline - get this from github and then create your own control file from this template - don't modify anything else
                # FASTQs - fill this file with your unmerged multiplexed raw FASTQ files
                # ReferenceLibraries - fill this file with your taxonomic reference library formatted either as i) ... or ii) ...
                # IntermediateOutputs - this will be populated by the pipeline as it runs, it will enable the pipeline to be run over multiple sessions as the output from each module is saved here.
                # Results - this is where final results will be saved
        dataname="test"
            # this should be the name you associate with this set of fastqs, all output and result files will be labelled with this name.


    # DADA2 settings
        truncLen=c(150,150)
            #the length at which the sequences should be truncated (forward and reverse)
        trimLeft=c(0,0)
            # the # of bases to remove from the start of the sequences
        maxN=0
            # after truncation seqs with more than this number of Ns are discarded (DADA2 does not allow Ns)
        maxEE=c(2,2)
            #  After truncation, reads with higher than ‘maxEE’ "expected errors"  
            # will be discarded. Expected errors are calculated from the nominal definition 
            # of the quality score: EE = sum(10^(-Q/10))
        truncQ=2
            # Truncate reads at the first instance of a quality score less than or equal to ‘truncQ’.  
        pool="pseudo"
            #TRUE, FALSE, or pseudo. pseudo pooling approximates the effect of denoising with pooled samples, but with
            # linearly increasing computational time (ca. doubled compared to no pooling)
        multithread=TRUE
            #on windows set multithread=FALSE

    # swarmv2 settings
        differences=1
            #number of base differences at which swarm clustering will be performed (1=default)                               
            

    # lulu settings
        MatchlistRate=90 #as a %
            # % matching bases to consider clustering OTUs if co-occurence seen. 

    # IDTAXA settings
        Type ="Load"  
            #whether to "Create" or "Load" a training set
        RefLibrary= "SILVA_SSU_r138_2019.RData" 
            #ref library to load if loading
        SeqsToAssign ="ESVs"
            #whether to assign to "ESVs", "OTUs", or "cOTUs"
        threshold=60
            # %age confidence of assignment required to record assignment

# run pipeline
    source(file.path(path, "BioinformaticPipeline", "Pipeline", "DSLI_SQL_Pipeline.R"))
