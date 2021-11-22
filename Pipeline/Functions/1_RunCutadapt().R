
#RunCutadapt(FWD ="ACCTGCGGARGGATCA",  ## CHANGE ME to your forward primer sequence
#            REV = "GAGATCCRTTGYTRAAAGTT",  ## CHANGE ME...
#            multithread=TRUE) #make FALSE for windows

RunCutadapt<-function(multithread, FWD, REV, dataname=NULL, UseCutadapt=FALSE) {
    
    if (UseCutadapt == FALSE) return(NULL)
    else if (UseCutadapt == TRUE){

        # Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
        fnFs<-loadFastq(FASTQ_folder=dataname, pattern="_R1_001.fastq")
        fnRs<-loadFastq(FASTQ_folder=dataname, pattern="_R2_001.fastq")

        FWD.orients <- allOrients(FWD)
        REV.orients <- allOrients(REV)


        # Extract sample names, assuming filenames have format: SAMPLENAME_RN_001.fastq.gz
        SampleNames <- sapply(strsplit(basename(fnFs), "_R1_001.fastq"), `[`, 1)
        
        # Create paths for filtered outputs
        filtFs<-createOutputFilePaths(suffix="_F_filt.fastq.gz", outputDirectoryPrefix="_prefilteredsequences")
        filtRs<-createOutputFilePaths(suffix="_R_filt.fastq.gz", outputDirectoryPrefix="_prefilteredsequences")

        #The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. Next we are going to “pre-filter” the sequences just to remove those with Ns, but perform no other filtering.
        out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, rm.phix=TRUE, compress=TRUE, multithread=multithread, verbose=TRUE) 


        # We are now ready to count the number of times the primers appear in the forward and reverse read, 
        #while considering all possible primer orientations. Identifying and counting the primers on one set of 
        #paired end FASTQ files is sufficient, 
        #assuming all the files were created using the same library preparation, so we’ll just process the first sample.
    
        PrimerCountSamplesList<-list()
        for (i in seq_along(filtFs)) {
            PrimerCountSamplesList[[i]]<-rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtFs[[i]]), 
                FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = filtRs[[i]]), 
                REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtFs[[i]]), 
                REV.ReverseReads = sapply(REV.orients, primerHits, fn = filtRs[[i]]))
        }

        path.cut <- file.path(path, "IntermediateOutputs", paste0(dataname,"_CutadaptedSeqs"))
        if(!dir.exists(path.cut)) dir.create(path.cut)
        fnFs.cut <- file.path(path.cut, basename(fnFs))
        fnRs.cut <- file.path(path.cut, basename(fnRs))

        FWD.RC <- dada2:::rc(FWD)
        REV.RC <- dada2:::rc(REV)

        # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
        R1.flags <- paste("-g", FWD, "-a", REV.RC) 

        # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
        R2.flags <- paste("-G", REV, "-A", FWD.RC) 

        # Run Cutadapt
        for(i in seq_along(fnFs)) {
        system2("cutadapt", args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                                    "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                                    filtFs[i], filtRs[i])) # input files
        }

        PrimerCountCutadaptedSamplesList<-list()
        for (i in seq_along(filtFs)) {
            PrimerCountCutadaptedSamplesList[[i]]<-rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[i]]), 
            FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[i]]), 
            REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[i]]), 
            REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[i]]))
        }


        return(list("PrimerCountSamplesList"=PrimerCountSamplesList, "PrimerCountCutadaptedSamplesList"=PrimerCountCutadaptedSamplesList))
    }

}
