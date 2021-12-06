#dada2 wrapped function now optimised for big datasets in accordance with http://benjjneb.github.io/dada2/bigdata.html

#currently assumes all fastqs come from a single sequencing run, and so can have a single learnt error rate applied. it
#doesn't matter if there were other samples included in the sequencing run as 1e8 bases is sufficient to learn error rates, 
# and denoising can take place by sample using this error rate. If multiple sequcing runs need to be combined then this function will need to have a
# multipleruns argument added which will run the steps till denoising on each sequecing run independently and then merge them for chimera removal and production of
#summary outputs, a well as all downstream modules. Further guidance available here: http://benjjneb.github.io/dada2/bigdata.html

        RunDADA2<-function(truncLen=NULL, trimLeft=NULL, maxN=0, maxEE=c(2,2), 
                            truncQ=2, DesiredSequenceLengthRange=NULL, dataname=NULL, multithread, pool,
                            UseCutadapt=FALSE) {
            #stop function if necessary arguments blank                    
            if (is.null(truncLen) | is.null(trimLeft)) stop("You must specify truncLen and trimLeft")

            #get filenames of fastas from right folder depending on if cutadapt run
            if (UseCutadapt==FALSE) {
                # Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
                fnFs<-loadFastq(FASTQ_folder=dataname, pattern="_R1_001.fastq")
                fnRs<-loadFastq(FASTQ_folder=dataname, pattern="_R2_001.fastq")
            } else if (UseCutadapt==TRUE) {
                fnFs<-sort(list.files(file.path(path, "IntermediateOutputs", paste0(dataname,"_CutadaptedSeqs")), pattern="_R1_001.fastq", full.names = TRUE))
                fnRs<-sort(list.files(file.path(path, "IntermediateOutputs", paste0(dataname,"_CutadaptedSeqs")), pattern="_R2_001.fastq", full.names = TRUE))
            }


            # Extract sample names, assuming filenames have format: SAMPLENAME_RN_001.fastq.gz
            SampleNames <- sapply(strsplit(basename(fnFs), "_R1_001.fastq"), `[`, 1)
            names(fnFs)<-SampleNames
            names(fnRs)<-SampleNames

            # Create paths for filtered outputs
            filtFs<-createOutputFilePaths(suffix="_F_filt.fastq.gz", outputDirectoryPrefix="_filteredsequences")
            filtRs<-createOutputFilePaths(suffix="_R_filt.fastq.gz", outputDirectoryPrefix="_filteredsequences")

            #filter and trim
            out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=truncLen,
                                trimLeft=trimLeft, maxN=maxN, maxEE=maxEE, truncQ=truncQ, 
                                rm.phix=TRUE, compress=TRUE, multithread=multithread, verbose=TRUE) 
            errF <- learnErrors(filtFs, multithread=multithread, nbases=1e8, randomize=TRUE)
            errR <- learnErrors(filtRs, multithread=multithread, nbases=1e8, randomize=TRUE)


            # Infer sequence variants
            dds <- vector("list", length(SampleNames))
            names(dds) <- SampleNames
            ddsFs <- dds
            ddsRs <- dds

            for(sam in SampleNames) {
                cat("Processing:", sam, "\n")
                derepFs <- derepFastq(filtFs[[sam]])
                derepRs <- derepFastq(filtRs[[sam]])

                ddsFs[[sam]] <- dada(derepFs, err=errF, multithread=TRUE)
                ddsRs[[sam]] <- dada(derepRs, err=errR, multithread=TRUE)

                #merge pairs
                    dds[[sam]] <- mergePairs(ddsFs[[sam]], derepFs, ddsRs[[sam]], derepRs, verbose=TRUE)

            }

            #make sequence table
            seqtab <- makeSequenceTable(dds)

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
            track <- cbind(out, sapply(ddsFs, getN), sapply(ddsRs, getN), sapply(dds, getN), rowSums(seqtab.nochim)) # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
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