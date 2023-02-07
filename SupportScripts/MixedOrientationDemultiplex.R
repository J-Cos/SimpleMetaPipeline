######################################
# Secondary demlutiplexing and unmixing orientations for primary demultiplexed mixed orientation Illumina MiSeq data
# only known from COI data Cross-pacific ARMS data (Emma and Dita)
#####################################


#1) dependencies    
    #installed software (get latest from bioconda e.g.  conda install -c bioconda blast==2.12.0)
        #cutadapt
    #R packages (install from bioconductor or through R devtools as appropriate )
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load( dada2,
                    seqinr,
                    tidyverse,
                    ShortRead,
                    Biostrings
                    )


# 2) code
    primerF<-DNAString(primerF)
    primerR<-DNAString(primerR)

    RunsSplitByAdapter<-list.files(file.path(path, "RunsSplitByAdapter"))
    AllBarcodes<-list.files(file.path(path, "Barcodes"))

    #get number of runs to loop over
        if (!RunSubset) {
            RunsToLoopOver<-1:length(RunsSplitByAdapter)
        } else {
            RunsToLoopOver<-RunSubset
        }

    #loop over runs
    for (r in RunsToLoopOver) {

        #get data for a single run
        Run<-RunsSplitByAdapter[r]
        Barcodes<-AllBarcodes[r]

        #get list of fastqs for sequences for which the adapter was determined (i.e. exclude undetermined)
        RunFastqs_withUndetermined<-list.files(file.path(path, "RunsSplitByAdapter", Run), pattern="fastq", full.names = TRUE)
        RunFastqs<-RunFastqs_withUndetermined[! grepl("Undetermined", RunFastqs_withUndetermined, fixed=TRUE)]

        #get list of barcode files
        RunBarcodes<-sort(list.files(file.path(path, "Barcodes", Barcodes), pattern=".txt", full.names = TRUE))

        # if lengths imply there is not one barcode for every pair of fastqs then break out of loop
        if (length(RunFastqs)/2!=length(RunBarcodes)) { 
            print("Not matching number of barcodes and fastq pairs")
            break 
            }

        #loop over each adapter in run
        for (a in (1:length(RunBarcodes))) {

            # get data for single adapter
            AdapterFastqs<-RunFastqs[ (a*2-1) : (a*2) ]
            AdapterBarcodesFile<- RunBarcodes[a]

            #make df of barcodes from this adapter
            Barcodes_df<-read.table(AdapterBarcodesFile)
            colnames(Barcodes_df)<-c("Sample", "Barcode")

            # Create paths for prefiltered fastqs in ephemeral directory 
            path_filtered<-str_replace(AdapterFastqs, "RunsSplitByAdapter", "RunsSplitByAdapter_prefiltered")
            path_filtered<-str_replace(AdapterFastqs, "home", "ephemeral")


            #prefilter fastqs
            #The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. Next we are going to “pre-filter” the sequences just to remove those with Ns, but perform no other filtering.
            out <- dada2::filterAndTrim(AdapterFastqs[1], path_filtered[1], AdapterFastqs[2], path_filtered[2], maxN=0, rm.phix=TRUE, compress=TRUE, multithread=Ncores, verbose=TRUE) 
            print(paste0("Run ", r, ": Adapter ", a, " filtering complete"))


            FO1<-lapply(paste0(Barcodes_df$Barcode, primerF), DNAString)
            RO1<-lapply(paste0(Barcodes_df$Barcode, primerR), DNAString)
            RO2<-lapply(FO1, reverseComplement)
            FO2<-lapply(RO1, reverseComplement)

            names(FO1)<-paste0(Barcodes_df$Sample, "_FO")
            names(RO1)<-paste0(Barcodes_df$Sample, "_RO")
            names(RO2)<-paste0(Barcodes_df$Sample, "_RO")
            names(FO2)<-paste0(Barcodes_df$Sample, "_FO")

            BarcodePrimers<-c(FO1, RO1)



            # Create paths for barcodes fastas
            path_barcodePrimerFastas<-str_replace(AdapterBarcodesFile, ".txt", "_withPrimers.fasta")

            #write fastas to this location
            seqinr::write.fasta(sequences= BarcodePrimers, names= names(BarcodePrimers), file.out= file.path(path_barcodePrimerFastas ) , open = "w", nbchar = 60, as.string = TRUE)

            # Create paths for demultiplexed fastqs
            path_demux<-file.path(path, "DemultiplexedFastqs")

            #now run cutadapt
            system2( "cutadapt", args= c( "-e", 0, "--no-indels",  "--minimum-length", 1,  # options (0% error rate, no insertion/deletions, and zero length seqs discarded)
                                        "-g", paste0("file:",path_barcodePrimerFastas),  "-G", paste0("file:",path_barcodePrimerFastas),       #barcodes
                                        "-o",   paste0(path,"/../../FASTQs/CrossPacific_COI/",Run, "/{name}_R1_001.fastq.gz"), "-p", paste0(path,"/../../FASTQs/CrossPacific_COI/",Run, "/{name}_R2_001.fastq.gz"),
                                        "-j", Ncores,
                                        path_filtered[[1]] ,  path_filtered[[2]]                       #inputs
                                    )
            ) 
            

            print(paste0("Run ", r, ": Adapter ", a, " cutadapt complete"))
            print(paste0(AdapterBarcodesFile ," matched with ", AdapterFastqs))

        }
    }