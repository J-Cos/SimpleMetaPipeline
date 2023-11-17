        RunDADA2<-function(truncLen=NULL, trimLeft=NULL, maxN=0, maxEE=c(2,2), 
                            truncQ=2, DesiredSequenceLengthRange=NULL, dataname=NULL, multithread, pool,
                            UseCutadapt=FALSE, MixedOrientation=FALSE, NumberOfRuns=1, ReadType="Paired-end") {
            #stop function if necessary arguments blank                    
            if (is.null(truncLen) | is.null(trimLeft)) stop("You must specify truncLen and trimLeft")
            if(!ReadType %in% c("Paired-end", "Single-read", "Premerged")) stop("ReadType not set correctly, must be one of: 'Paired-end', 'Single-read', 'Premerged' ")

            #make lists to populate with outputs from each Run
            RunESVtables<-list()
            RunDadaPlots<-list()
            RunDadaTables<-list()

            #loop over all runs
                #get all run names
                RunNames<-list.files(file.path(dataPath, "FASTQs", dataname))

                #first check if any runs completed in a previous exited run, if so load these and append uncompleted runs to this object
                if(file.exists(file=file.path(path, "IntermediateOutputs", paste0(dataname,"_RunESVtables.RDS")))) {
                    RunESVtables<-readRDS(file=file.path(path, "IntermediateOutputs", paste0(dataname,"_RunESVtables.RDS") ))
                    RunDadaPlots<-readRDS(file=file.path(path, "IntermediateOutputs", paste0(dataname,"_RunDadaPlots.RDS") ))
                    RunDadaTables<-readRDS(file=file.path(path, "IntermediateOutputs", paste0(dataname,"_RunDadaTables.RDS") ))
                }
                #adjust number of runs so loop over uncompleted runs only
                if(exists("RunESVtables")) {NumberCompletedRuns<-length(RunESVtables)}
            
            if (NumberCompletedRuns<NumberOfRuns) {
                for (Run in (NumberCompletedRuns+1):NumberOfRuns) {

                    print( paste0( "RUN ", Run))

                    if(ReadType=="Paired-end"){

                        #get filenames of fastas from right folder depending on if cutadapt run
                        if (UseCutadapt==FALSE) {
                            # Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
                            fnFs<-sort(list.files(file.path(dataPath, "FASTQs", dataname, RunNames[Run]), pattern="_R1_001.fastq", full.names = TRUE))
                            fnRs<-sort(list.files(file.path(dataPath, "FASTQs", dataname, RunNames[Run]), pattern="_R2_001.fastq", full.names = TRUE))
                        } else if (UseCutadapt==TRUE) {
                            fnFs<-sort(list.files(file.path(path, "IntermediateOutputs", paste0(dataname,"_CutadaptedSeqs")), pattern="_R1_001.fastq", full.names = TRUE))
                            fnRs<-sort(list.files(file.path(path, "IntermediateOutputs", paste0(dataname,"_CutadaptedSeqs")), pattern="_R2_001.fastq", full.names = TRUE))
                        }

                        # Extract sample names, assuming filenames have format: SAMPLENAME_RN_001.fastq.gz
                        SampleNames <- sapply(strsplit(basename(fnFs), "_R1_001.fastq"), `[`, 1) %>%
                                        paste0(.,"__Run", Run)
                        
                        # Create paths for filtered outputs
                        filtFs<-createOutputFilePaths(suffix="_F_filt.fastq.gz", outputDirectoryPrefix="_filteredsequences")
                        filtRs<-createOutputFilePaths(suffix="_R_filt.fastq.gz", outputDirectoryPrefix="_filteredsequences")

                        #filter and trim
                        out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=truncLen[[Run]],
                                            trimLeft=trimLeft[[Run]], maxN=maxN, maxEE=maxEE, truncQ=truncQ, 
                                            rm.phix=TRUE, compress=TRUE, multithread=multithread, verbose=TRUE) 
                        print("filterAndTrim Complete")
                        errF <- learnErrors(filtFs, multithread=multithread)
                        print("learnErrors Forwards Complete")
                        errR <- learnErrors(filtRs, multithread=multithread)
                        print("learnErrors Reverses Complete")

                        #denoise
                        dadaFs <- dada(filtFs, err=errF, multithread=multithread,pool=pool)
                        print("denoise Forwards Complete")
                        dadaRs <- dada(filtRs, err=errR, multithread=multithread,pool=pool)
                        print("denoise Reverses Complete")

                        #merge
                        mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap=5)
                        print("merging Complete")

                        #create dada plots
                            RunDadaPlots[[Run]]<-list()
                            RunDadaPlots[[Run]][[1]]<-plotQualityProfile(fnFs) #example
                            RunDadaPlots[[Run]][[2]]<-plotQualityProfile(fnRs) #example
                            RunDadaPlots[[Run]][[3]]<-plotQualityProfile(filtFs) #example
                            RunDadaPlots[[Run]][[4]]<-plotQualityProfile(filtRs) #example
                            RunDadaPlots[[Run]][[5]]<-plotErrors(errF, nominalQ=TRUE)
                            RunDadaPlots[[Run]][[6]]<-plotErrors(errR, nominalQ=TRUE)

                        #create read tracking table 
                        track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN)) # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
                        colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")

                    } else if (ReadType!="Paired-end"){
                        #get merged filenames
                        fns<-sort(list.files(file.path(dataPath, "FASTQs", dataname, RunNames[Run]), pattern=".fastq", full.names = TRUE))

                        # Extract sample names, assuming filenames have format: SAMPLENAME_RN_001.fastq.gz
                        SampleNames <- sapply(strsplit(basename(fns), ".fastq"), `[`, 1) %>%
                                        paste0(.,"__Run", Run)

                        # Create paths for filtered outputs
                        filts<-createOutputFilePaths(suffix="_filt.fastq.gz", outputDirectoryPrefix="_filteredsequences")

                        #filter and trim
                        out <- filterAndTrim(fns, filts, truncLen=truncLen[[Run]],
                                            trimLeft=trimLeft[[Run]], maxN=maxN, maxEE=maxEE, truncQ=truncQ, 
                                            rm.phix=TRUE, compress=TRUE, multithread=multithread, verbose=TRUE) 
                        print("filterAndTrim Complete")

                        err <- learnErrors(filts, multithread=multithread)
                        print("learnErrors Complete")

                        if (ReadType=="Premerged"){
                            #inflate errors to deal with problems introduced to error scores by merging forward and reverse befoer running dada2
                            inflatedErr<-inflateErr(err, 3)
                        }

                        #denoise
                        dada <- dada(filts, err=inflatedErr, multithread=multithread,pool=pool)
                        print("denoise Complete")

                        #rename to mergers to feed into next step making sequence table
                        mergers <- dada    
                                                
                        #create dada plots
                            RunDadaPlots[[Run]]<-list()
                            RunDadaPlots[[Run]][[1]]<-plotQualityProfile(fns) #example
                            RunDadaPlots[[Run]][[2]]<-plotQualityProfile(filts) #example
                            RunDadaPlots[[Run]][[3]]<-plotErrors(inflatedErr, nominalQ=TRUE)

                        #create read tracking table 
                        track <- cbind(out, sapply(dada, getN)) # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
                        colnames(track) <- c("input", "filtered", "denoised")
                    }

                    #finish read tracking table (these steps not specific to premerged or non-premerged)
                    rownames(track) <- SampleNames
                    RunDadaTables[[paste0("Run", Run)]]<-track

                    #make sequence table
                    seqtab <- makeSequenceTable(mergers)

                    #append run number to sample names to enable differentiation of samples included in multiple runs
                    row.names(seqtab)<-SampleNames

                    #cut to specific length if needed
                    if (!is.null(DesiredSequenceLengthRange)) {
                        seqtab <- seqtab[,nchar(colnames(seqtab)) %in% DesiredSequenceLengthRange]
                    }
                
                    #make sequence length distribution table
                    SeqLengthDist<-table(nchar(getSequences(seqtab )))

                    #create main output of ESVtable    
                    ESVtable<-seqtab

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
                    CacheOutput(RunESVtables)
                    CacheOutput(RunDadaPlots)
                    CacheOutput(RunDadaTables)
                } #finish looping over all runs
            }

            #merge runs
            if (NumberOfRuns>1) {
                ESVtable <- mergeSequenceTables(tables=RunESVtables, repeats="sum") #%>% collapseNoMismatch # merge all sequence tables and combine sequences of different lengths with all matching bases
                DadaPlots<-unlist(RunDadaPlots, recursive=FALSE)
                DadaTables<-do.call("rbind", RunDadaTables)
                SeqLengthDist<-table(nchar(getSequences(ESVtable)))

            } else if (NumberOfRuns==1) {
                ESVtable <- RunESVtables[[1]]
                DadaPlots<-RunDadaPlots[[1]]
                DadaTables<-RunDadaTables[[1]]
                SeqLengthDist<-table(nchar(getSequences(ESVtable)))
            }

            #post merging chimera removal
                ESVtable<- removeBimeraDenovo(ESVtable, method="consensus", 
                                                multithread=multithread, verbose=TRUE)
                #add chimera removal tracking colunm to table

                    chimericData<-rowSums(ESVtable) %>% 
                        as.data.frame %>% 
                        rownames_to_column("rowname_for_merging")%>%
                        rename("."= "non-chimeric, applied post reorientation and merging of FO and RO") 

                    DadaTables<-DadaTables%>%  
                        as.data.frame %>% 
                        rownames_to_column %>%
                        as_tibble %>%
                        mutate(rowname_for_merging=str_replace(rowname, "_FO","")) %>%
                        mutate(rowname_for_merging=str_replace(rowname_for_merging, "_RO","")) %>% 
                        left_join(.,chimericData) %>%
                        select(!rowname_for_merging) %>%
                        arrange(rowname)

            #format ESV data 
            ESVtable<-t(ESVtable) # make samples columns
            ESVtable<-cbind(paste0("ESV_",1:dim(ESVtable)[1]), rownames_to_column(as.data.frame(ESVtable))) #add esv id column
            names(ESVtable)[1:2]<-c("ESV", "Sequence") 
            names(ESVtable)[3:length(names(ESVtable))]<-paste0("Sample_", names(ESVtable)[3:length(names(ESVtable))]) # giving all samples "Sample_" prefix

            return(list("SeqDataTable"=ESVtable, "SecondaryOutputs"=list("DadaPlots"=DadaPlots, "DadaTables"=DadaTables, "SeqLengthDist"=SeqLengthDist)))
        }