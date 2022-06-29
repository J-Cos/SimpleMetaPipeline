######################################
# Secondary demlutiplexing and unmixing orientations for primary demultiplexed mixed orientation Illumina MiSeq COI data
# Cross-pacific ARMS data
#####################################


# 1) paths and dependencies
    library(tidyverse)
    library(dada2)
    library(ShortRead)

    path<-'/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/2015'
    #contents of this directory should be:
        #Barcodes - this should contain a subdirectory for each run which in turn contains a list of the barcodes related to each run for each adapter
        #RunsSplitByAdapter - this should contain a subdirectory for each run which in turn contain R1 and R2 fastqs for each adapter
    
    #primers
    mlCOIintF<-DNAString("GGWACWGGWTGAACWGTWTAYCCYCC")
    jgHCO2198<-DNAString("TANACYTCNGGRTGNCCRAARAAYCA")



# 2) code
    RunsSplitByAdapter<-list.files(file.path(path, "RunsSplitByAdapter"))
    AllBarcodes<-list.files(file.path(path, "Barcodes"))

    #loop over all runs
    for (r in 1:length(RunsSplitByAdapter)) {

        #get data for a single run
        Run<-RunsSplitByAdapter[r]
        Barcodes<-AllBarcodes[r]

        #get list of fastqs for sequences for which the adapter was determined (i.e. exclude undetermined)
        RunFastqs_withUndetermined<-list.files(file.path(path, "RunsSplitByAdapter", Run), pattern="fastq", full.names = TRUE)
        RunFastqs<-RunFastqs_withUndetermined[-((length(RunFastqs_withUndetermined)-1):length(RunFastqs_withUndetermined))]

        #get list of barcode files
        RunBarcodes<-sort(list.files(file.path(path, "Barcodes", Barcodes), pattern=".txt", full.names = TRUE))

        # if lengths imply there is not one barcode for every pair of fastqs then break out of loop
        if (length(RunFastqs)/2!=length(RunBarcodes)) { break }

        #loop over each adapter in run
        for (a in (1:length(RunBarcodes))) {

            # get data for single adapter
            AdapterFastqs<-RunFastqs[ (a*2-1) : (a*2) ]
            AdapterBarcodesFile<- RunBarcodes[a]

            #make df of barcodes from this adapter
            Barcodes_df<-read.table(AdapterBarcodesFile)
            colnames(Barcodes_df)<-c("Sample", "Barcode")

            # Create paths for prefiltered fastqs
            path_filtered<-str_replace(AdapterFastqs, "RunsSplitByAdapter", "RunsSplitByAdapter_prefiltered")

            #prefilter fastqs
            #The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. Next we are going to “pre-filter” the sequences just to remove those with Ns, but perform no other filtering.
            out <- dada2::filterAndTrim(AdapterFastqs[1], path_filtered[1], AdapterFastqs[2], path_filtered[2], maxN=0, rm.phix=TRUE, compress=TRUE, multithread=TRUE, verbose=TRUE) 

            #create all orients
            FO1<-lapply(paste0(Barcodes_df$Barcode, mlCOIintF), DNAString)
            RO1<-lapply(paste0(Barcodes_df$Barcode, jgHCO2198), DNAString)

            names(FO1)<-paste0(Barcodes_df$Sample, "_FO")
            names(RO1)<-paste0(Barcodes_df$Sample, "_RO")

            BarcodePrimers<-c(FO1, RO1)

            # Create paths for barcodes fastas
            path_barcodePrimerFastas<-str_replace(AdapterBarcodesFile, ".txt", "_withPrimers.fasta")

            #write fastas to this location
            seqinr::write.fasta(sequences= BarcodePrimers, names= names(BarcodePrimers), file.out= file.path(path_barcodePrimerFastas ) , open = "w", nbchar = 60, as.string = TRUE)

            # Create paths for demultiplexed fastqs
            path_demux<-file.path(path, "DemultiplexedFastqs")

            #now run cutadapt
            system2( "cutadapt", args= c( "-e", 0.1, "--no-indels",  "--minimum-length", 1,  # options (10% error rate, no insertion/deletions, and zero length seqs discarded)
                                        "-g", paste0("file:",path_barcodePrimerFastas),  "-G", paste0("file:",path_barcodePrimerFastas),       #barcodes
                                        "-o",   "/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific_2015/{name}_R1_001.fastq", "-p", "/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific_2015/{name}_R2_001.fastq",
                                    path_filtered[[1]] ,  path_filtered[[2]]                       #inputs
                                )
            ) 
        }
    }