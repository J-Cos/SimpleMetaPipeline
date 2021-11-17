        RunDADA2<-function(truncLen=NULL, trimLeft=NULL, maxN=0, maxEE=c(2,2), truncQ=2, DesiredSequenceLengthRange=NULL, dataname=NULL, multithread, pool) {
            #stop function if necessary arguments blank                    
            if (is.null(truncLen) | is.null(trimLeft)) stop("You must specify truncLen and trimLeft")

            # Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
            fnFs<-loadFastq(FASTQ_folder=dataname, pattern="_R1_001.fastq")
            fnRs<-loadFastq(FASTQ_folder=dataname, pattern="_R2_001.fastq")

            # Extract sample names, assuming filenames have format: SAMPLENAME_RN_001.fastq.gz
            SampleNames <- sapply(strsplit(basename(fnFs), "_R1_001.fastq"), `[`, 1)
            
            # Create paths for filtered outputs
            filtFs<-createOutputFilePaths(suffix="_F_filt.fastq.gz", outputDirectoryPrefix="filteredsequences")
            filtRs<-createOutputFilePaths(suffix="_R_filt.fastq.gz", outputDirectoryPrefix="filteredsequences")

            #filter and trim
            out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=truncLen,
                                trimLeft=trimLeft, maxN=maxN, maxEE=maxEE, truncQ=truncQ, 
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

            #make sequence length distribution table
            SeqLengthDist<-table(nchar(getSequences(seqtab)))

            #remove chimeras
            seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                                multithread=multithread, verbose=TRUE)

            #create read tracking table 
            track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)) # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
            colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
            rownames(track) <- SampleNames
            DadaTables<-list(track)

            #create dada plots
            DadaPlots<-c()
            DadaPlots[[1]]<-plotQualityProfile(fnFs) #example
            DadaPlots[[2]]<-plotQualityProfile(fnRs) #example
            DadaPlots[[3]]<-plotQualityProfile(filtFs) #example
            DadaPlots[[4]]<-plotQualityProfile(filtRs) #example
            DadaPlots[[5]]<-plotErrors(errF, nominalQ=TRUE)
            DadaPlots[[6]]<-plotErrors(errR, nominalQ=TRUE)
                
            #format ESV data 
            ESVtable<-seqtab.nochim
            ESVtable<-t(ESVtable) # make samples columns
            ESVtable<-cbind(paste0("ESV_",1:dim(ESVtable)[1]), rownames_to_column(as.data.frame(ESVtable))) #add esv id column
            names(ESVtable)[1:2]<-c("ESV", "Sequence")
            names(ESVtable)[3:length(names(ESVtable))]<-paste0("Sample_", names(ESVtable)[3:length(names(ESVtable))]) # giving all samples "Sample_" prefix


            return(list("SeqDataTable"=ESVtable, "SecondaryOutputs"=list("DadaPlots"=DadaPlots, "DadaTables"=DadaTables, "SeqLengthDist"=SeqLengthDist)))
        }