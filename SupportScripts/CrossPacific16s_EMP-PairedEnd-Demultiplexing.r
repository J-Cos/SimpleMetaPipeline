######################################
# Demultiplexing EMP paired end reads
#####################################

# 0) HPC set up
    setwd("/rds/general/user/jcw120/home/BioinformaticPipeline_Env") #necessary as it appears different job classes have different WDs.

        #CRAN mirror
            r = getOption("repos")
            r["CRAN"] = "http://cran.us.r-project.org"
            options(repos = r)

# 1a) packages and parameters
    library(tidyverse)
    library(dada2)
    library(ShortRead)
    library(seqinr)

# 1b) functions
    combineFastqSequences<-function(R1, I) {
        numCycles<-length(I)
        seqs<-DNAStringSet(rep("", numCycles))
        quals<-BStringSet(rep("", numCycles))
        ids<-BStringSet(rep("", numCycles))
        for (seq in 1: numCycles) {
            seqs[[seq]]<-c(sread(I)[[seq]], sread(R1)[[seq]])
            quals[[seq]]<-c(quality(I)[[seq]], quality(R1)[[seq]])
            ids[seq]<-id(R1)[seq]
        }
        R1_New<-ShortReadQ( seqs,quals,ids)
        return(R1_New)
    }


# 1c) Directory locations
    RunFolders<-list.dirs(file.path(path, "Runs"), recursive=FALSE, full.names=TRUE)
    AllBarcodes<-list.files(file.path(path, "Barcodes"), full.names=TRUE,  pattern=".csv")

#loop over all runs
    for (Run in 1: length(RunFolders)) {

    # 2 ) combine R1 and I1 fastqs ahead of cutadapt demultiplexing
            FASTQpaths<-list.files( RunFolders[Run], pattern="Undetermined_", recursive=TRUE, full.names=TRUE)    
            I<-ShortRead::readFastq(FASTQpaths[1])
            R1<-ShortRead::readFastq(FASTQpaths[2])
            R1_New<-combineFastqSequences(R1, I)

            #write to ephemeral directory
            CombinedFastqPath<-file.path(RunFolders[Run], "Combined_R1-I1.fastq.gz")
            CombinedFastqPath<-str_replace(CombinedFastqPath, "Runs", "Runs_Combined_R1-I1_Fastqs")
            CombinedFastqPath<-str_replace(CombinedFastqPath, "home", "ephemeral")
            dir.create(CombinedFastqPath, recursive=TRUE)
            writeFastq(object=R1_New, file= CombinedFastqPath )

            print(paste0("Written combined fastq from Run ", Run))


    # 3) Demultiplex

                BarcodesPath<-AllBarcodes[Run]
                R1Path<-CombinedFastqPath
                R2Path<-FASTQpaths[3]

                #barcode prep
                    #make df of barcodes from this adapter
                    Barcodes_df<-read.csv(BarcodesPath)

                #extract barcode lengths to specify that compelte overlap is required in demultiplexing
                if (length(unique(nchar(Barcodes_df$barcode.sequence)))!=1) {
                    print(paste0("ERROR: Multiple Barcode lengths detected. As a result minimum overlap in demultiplexing cannot be set to a single length suitable for all barcode primers. Please note that this is a required to work around as anchoring adapters (using e.g. '-g ^file:pathtofastq' notation) appears impossible through an R system() command "))
                } else {
                    BarcodeLength<-unique(nchar(Barcodes_df$barcode.sequence))
                    print(paste0("All barcodes are of length ", BarcodeLength))
                }

                    #make barcode fasta path
                    BarcodeFastaPath<-str_replace(BarcodesPath, ".csv", ".fasta")
                    #write fastas to this location
                    seqinr::write.fasta(sequences= as.list(Barcodes_df$barcode.sequence), 
                                        names= Barcodes_df$sample, 
                                        file.out= BarcodeFastaPath, 
                                        open = "w", 
                                        nbchar = 200, 
                                        as.string = TRUE)

                #seq prefiltering
                    # Create paths for prefiltered fastqs
                    path_filtered<-str_replace(CombinedFastqPath, "Combined_R1-I1_Fastqs", "Runs_prefiltered")

                    #prefilter fastqs
                    #The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. Next we are going to “pre-filter” the sequences just to remove those with Ns, but perform no other filtering.
                    out <- dada2::filterAndTrim(R1Path, path_filtered[1], R2Path, path_filtered[2], maxN=0, rm.phix=TRUE, compress=TRUE, multithread=TRUE, verbose=TRUE) 

                #demultiplexing
                    # Create paths for demultiplexed fastqs
                    path_demux<-file.path("FASTQs", "CrossPacific_16s")

                    #now run cutadapt
                    system2( "cutadapt", args= c( "-e", 0, "--no-indels",  "--minimum-length", 1,  # options (0% error rate, no insertion/deletions, and zero length seqs discarded)
                                                "--overlap", BarcodeLength, #overlap equal to length of barcodes + primers 
                                                "-g", paste0("^file:",BarcodeFastaPath),       #barcodes
                                                "-o",   file.path(path_demux, "{name}_R1_001.fastq"), "-p", file.path(path_demux, "{name}_R2_001.fastq"),
                                                "-j", Ncores,
                                                path_filtered[[1]] ,  path_filtered[[2]]                       #inputs
                                        )
                    )
                print(paste0("Run ", Run,  " cutadapt complete"))
                print(paste0(RunFolders[Run] ," matched with ", AllBarcodes[Run]))
                 
            }