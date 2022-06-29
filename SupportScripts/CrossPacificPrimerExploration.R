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
        #Barcodes
        #RunsSplitByAdapter
    
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

            #function to # Create all orientations of the input sequence
            allOrients <- function(primer) {
                require(Biostrings)
                dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
                orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                    RevComp = reverseComplement(dna))
                return(sapply(orients, toString))  # Convert back to character vector
            }

            #create all orients
            #barcode_orients_l<-lapply(Barcodes_df$Barcode, allOrients)


            FO1<-lapply(paste0(Barcodes_df$Barcode, mlCOIintF), DNAString)
            RO1<-lapply(paste0(Barcodes_df$Barcode, jgHCO2198), DNAString)
            RO2<-lapply(FO1, reverseComplement)
            FO2<-lapply(RO1, reverseComplement)

            names(FO1)<-paste0(Barcodes_df$Sample, "_FO")
            names(RO1)<-paste0(Barcodes_df$Sample, "_RO")
            names(RO2)<-paste0(Barcodes_df$Sample, "_RO")
            names(FO2)<-paste0(Barcodes_df$Sample, "_FO")

            BarcodePrimers<-c(FO1, RO1)





            #possibilities_df<- rbind(   expand.grid (barcode_orients_l[[1]], c(allOrients(jgHCO2198), allOrients(mlCOIintF))),
            #                            expand.grid ( c(allOrients(jgHCO2198), allOrients(mlCOIintF)), barcode_orients_l[[1]])
            #    )

            #possibilities<-paste0(possibilities_df[,1], possibilities_df[,2])

            #function to count number of hits of barcode
            Hits <- function(primer, fn) {
                # Counts number of reads in which the primer is found
                nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
                return(sum(nhits > 0))
            }

            #A1_ForwardReads = sapply(possibilities, Hits, fn = path_filtered[[1]])
            #A1_ReverseReads = sapply(possibilities, Hits, fn = path_filtered[[2]])

            #count hits
            #barcodeCounts <-    rbind(  A1_ForwardReads = sapply(barcode_orients_l[[1]], Hits, fn = path_filtered[[1]]), 
            #                            A1_ReverseReads = sapply(barcode_orients_l[[1]], Hits, fn = path_filtered[[2]]),
            #                            A2_ForwardReads = sapply(barcode_orients_l[[2]], Hits, fn = path_filtered[[1]]), 
            #                            A2_ReverseReads = sapply(barcode_orients_l[[2]], Hits, fn = path_filtered[[2]]),
            #                            A3_ForwardReads = sapply(barcode_orients_l[[3]], Hits, fn = path_filtered[[1]]), 
            #                            A3_ReverseReads = sapply(barcode_orients_l[[3]], Hits, fn = path_filtered[[2]]),
            #                            A4_ForwardReads = sapply(barcode_orients_l[[4]], Hits, fn = path_filtered[[1]]), 
            #                            A4_ReverseReads = sapply(barcode_orients_l[[4]], Hits, fn = path_filtered[[2]]),
            #                            A5_ForwardReads = sapply(barcode_orients_l[[5]], Hits, fn = path_filtered[[1]]), 
            #                            A5_ReverseReads = sapply(barcode_orients_l[[5]], Hits, fn = path_filtered[[2]]),
            #                            A6_ForwardReads = sapply(barcode_orients_l[[6]], Hits, fn = path_filtered[[1]]), 
            #                            A6_ReverseReads = sapply(barcode_orients_l[[6]], Hits, fn = path_filtered[[2]]),
            #                            A7_ForwardReads = sapply(barcode_orients_l[[7]], Hits, fn = path_filtered[[1]]), 
            #                            A7_ReverseReads = sapply(barcode_orients_l[[7]], Hits, fn = path_filtered[[2]])
            #                    )

            #count hits
            #barcodePrimerCounts <-    rbind(  ForwardReads = sapply(BarcodePrimers, Hits, fn = path_filtered[[1]]), 
            #                            ReverseReads = sapply(BarcodePrimers, Hits, fn = path_filtered[[2]]),
            #                    )
            #write.table(x=barcodeCounts,
            #            file=file.path("/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific_2015", paste0("BarcodeCount_Adapter_", a, "_", RunsSplitByAdapter[r]))
            #    )

            # Create paths for barcodes fastas
            #path_barcodeFastas<-str_replace(AdapterBarcodesFile, ".txt", ".fasta")

            #write fastas to this location
            #seqinr::write.fasta(sequences= as.list(Barcodes_df$Barcode), names= Barcodes_df$Sample, file.out= file.path(path_barcodeFastas ) , open = "w", nbchar = 60, as.string = TRUE)



            # Create paths for barcodes fastas
            path_barcodePrimerFastas<-str_replace(AdapterBarcodesFile, ".txt", "_withPrimers.fasta")

            #write fastas to this location
            seqinr::write.fasta(sequences= BarcodePrimers, names= names(BarcodePrimers), file.out= file.path(path_barcodePrimerFastas ) , open = "w", nbchar = 60, as.string = TRUE)




            # Create paths for demultiplexed fastqs
            path_demux<-file.path(path, "DemultiplexedFastqs")

            #now run cutadapt
            #system2( "cutadapt", args= c( "-e", 1, "--no-indels",   # options
            #                            "-g", paste0("file:",path_barcodeFastas),  "-G", paste0("file:",path_barcodeFastas),       #barcodes
            #                            "-o",   "/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific_2015/{name}_R1_001.fastq", "-p", "/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific_2015/{name}_R2_001.fastq",
            #                           path_filtered[[1]] ,  path_filtered[[2]]                       #inputs
            #                    )
            #) 


                #alternative cutadapt usign primers
            system2( "cutadapt", args= c( "-e", 0.1, "--no-indels",  "--minimum-length", 1,  # options (10% error rate, no insertion/deletions, and zero length seqs discarded)
                                        "-g", paste0("file:",path_barcodePrimerFastas),  "-G", paste0("file:",path_barcodePrimerFastas),       #barcodes
                                        "-o",   "/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific_2015/{name}_R1_001.fastq", "-p", "/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific_2015/{name}_R2_001.fastq",
                                    path_filtered[[1]] ,  path_filtered[[2]]                       #inputs
                                )
            ) 
        }
    }