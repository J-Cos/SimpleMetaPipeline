        RunDADA2<-function(truncLen=NULL, trimLeft=NULL, maxN=0, maxEE=c(2,2), 
                            truncQ=2, DesiredSequenceLengthRange=NULL, dataname=NULL, multithread, pool,
                            UseCutadapt=FALSE, MixedOrientation=FALSE, NumberOfRuns=1) {
            #stop function if necessary arguments blank                    
            if (is.null(truncLen) | is.null(trimLeft)) stop("You must specify truncLen and trimLeft")

            #make lists to populate with outputs from each Run
            RunESVtables<-list()
            RunDadaPlots<-list()
            RunDadaTables<-list()

            #loop over all runs
            for (Run in 1:NumberOfRuns) {

                #get filenames of fastas from right folder depending on if cutadapt run
                if (UseCutadapt==FALSE) {
                    # Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
                    fnFs<-sort(list.files(file.path(path, "FASTQs", dataname, paste0("Run", Run)), pattern="_R1_001.fastq", full.names = TRUE))
                    fnRs<-sort(list.files(file.path(path, "FASTQs", dataname, paste0("Run", Run)), pattern="_R2_001.fastq", full.names = TRUE))
                } else if (UseCutadapt==TRUE) {
                    fnFs<-sort(list.files(file.path(path, "IntermediateOutputs", paste0(dataname,"_CutadaptedSeqs")), pattern="_R1_001.fastq", full.names = TRUE))
                    fnRs<-sort(list.files(file.path(path, "IntermediateOutputs", paste0(dataname,"_CutadaptedSeqs")), pattern="_R2_001.fastq", full.names = TRUE))
                }


                # Extract sample names, assuming filenames have format: SAMPLENAME_RN_001.fastq.gz
                SampleNames <- sapply(strsplit(basename(fnFs), "_R1_001.fastq"), `[`, 1)
                
                # Create paths for filtered outputs
                filtFs<-createOutputFilePaths(suffix="_F_filt.fastq.gz", outputDirectoryPrefix="_filteredsequences")
                filtRs<-createOutputFilePaths(suffix="_R_filt.fastq.gz", outputDirectoryPrefix="_filteredsequences")

                #filter and trim
                out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=truncLen[[Run]],
                                    trimLeft=trimLeft[[Run]], maxN=maxN, maxEE=maxEE, truncQ=truncQ, 
                                    rm.phix=TRUE, compress=TRUE, multithread=multithread, verbose=TRUE) 
                errF <- learnErrors(filtFs, multithread=multithread)
                errR <- learnErrors(filtRs, multithread=multithread)

                #denoise
                dadaFs <- dada(filtFs, err=errF, multithread=multithread,pool=pool)
                dadaRs <- dada(filtRs, err=errR, multithread=multithread,pool=pool)

                #merge
                mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

                #make sequence table
                seqtab <- makeSequenceTable(mergers)

                #cut to specific length if needed
                if (!is.null(DesiredSequenceLengthRange)) {
                    seqtab <- seqtab[,nchar(colnames(seqtab)) %in% DesiredSequenceLengthRange]
                }
            

                #remove chimeras
                seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                                    multithread=multithread, verbose=TRUE)

                #make sequence length distribution table
                SeqLengthDist<-table(nchar(getSequences(seqtab.nochim )))

                #create read tracking table 
                track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)) # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
                colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
                rownames(track) <- SampleNames
                RunDadaTables[[paste0("Run", Run)]]<-track

                #create dada plots
                RunDadaPlots[[Run]]<-list()
                RunDadaPlots[[Run]][[1]]<-plotQualityProfile(fnFs) #example
                RunDadaPlots[[Run]][[2]]<-plotQualityProfile(fnRs) #example
                RunDadaPlots[[Run]][[3]]<-plotQualityProfile(filtFs) #example
                RunDadaPlots[[Run]][[4]]<-plotQualityProfile(filtRs) #example
                RunDadaPlots[[Run]][[5]]<-plotErrors(errF, nominalQ=TRUE)
                RunDadaPlots[[Run]][[6]]<-plotErrors(errR, nominalQ=TRUE)

                #create main output of ESVtable    
                ESVtable<-seqtab.nochim

                # standardise orientations if each sample consisted of FO_R1, FO_R2, RO_R1, RO_R2 (i.e read were in mixed orientations)
                if (MixedOrientation) {
                    #get Reverse oriented subsamples, filter out sequences only appearing in forward orientation, reverse complement remaining, and reinsert them into seq table
                    RevSamples<-ESVtable[!grepl('_FO',  rownames(ESVtable), fixed=T),]
                    RevESVtable<-RevSamples[,colSums(RevSamples)!=0]            
                    ReorientedSeqs<- reverseComplement(DNAStringSet(colnames(RevESVtable)))
                    ReorientedESVtable<-RevESVtable
                    colnames(ReorientedESVtable)<- as.character(ReorientedSeqs)

                    # get forward oriented reads, filter out sequences only appearing in reverse orientation
                    ForSamples<-ESVtable[!grepl('_RO',  rownames(ESVtable), fixed=T),]
                    ForESVtable<-ForSamples[,colSums(ForSamples)!=0]            

                    # remove orientation identifiers in sample names as all seqs are now in the same orientation
                    rownames(ReorientedESVtable)<- gsub('_RO', '', rownames(ReorientedESVtable))
                    rownames(ForESVtable)<- gsub('_FO', '', rownames(ForESVtable))

                    #merge forward and reverse complemented reverses, this merges both across sequences (standard) and across samples ( achieved through repeats=sum)
                    ESVtable<-mergeSequenceTables(ForESVtable, ReorientedESVtable, repeats="sum")
                }

                RunESVtables[[Run]]<-ESVtable
            } #finish looping over all runs

            #merge runs
            ESVtable <- mergeSequenceTables(tables=RunESVtables, repeats="Error") %>% collapseNoMismatch # merge all sequence tables and combine sequences of different lengths with all matching bases
            DadaPlots<-unlist(RunDadaPlots, recursive=FALSE)
            DadaTables<-do.call("rbind", RunDadaTables)

            #format ESV data 
            ESVtable<-t(ESVtable) # make samples columns
            ESVtable<-cbind(paste0("ESV_",1:dim(ESVtable)[1]), rownames_to_column(as.data.frame(ESVtable))) #add esv id column
            names(ESVtable)[1:2]<-c("ESV", "Sequence") 
            names(ESVtable)[3:length(names(ESVtable))]<-paste0("Sample_", names(ESVtable)[3:length(names(ESVtable))]) # giving all samples "Sample_" prefix

            

            return(list("SeqDataTable"=ESVtable, "SecondaryOutputs"=list("DadaPlots"=DadaPlots, "DadaTables"=DadaTables, "SeqLengthDist"=SeqLengthDist)))
        }