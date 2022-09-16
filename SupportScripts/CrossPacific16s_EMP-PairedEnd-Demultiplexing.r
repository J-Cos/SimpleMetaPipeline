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


    path <-"../BioinformaticPipeline_Env/FASTQs/EmmaRAWruns/EmmaAll_16s/CrossPacific"
    #contents of this directory should be:
        #Barcodes - this should contain a csv for each run which contains a list of the barcodes and the associated sample name
        #Each Run folder, containing I1, R1 and R2 undetermined fastq,gzs
        #NB - the order of th run folders must match the order of the barcode csvs to ensure they align - 


# 1b) functions
    combineFastqSequences<-function(R1, I) {
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
    RunFolders<-list.dirs(file.path(path), full.names=FALSE, recursive=FALSE)
    RunFolders<-RunFolders[RunFolders!="Barcodes"]
    AllBarcodes<-list.files(file.path(path, "Barcodes"), pattern=".csv")


# 2 ) combine R1 and I1 fastqs ahead of cutadapt demultiplexing
    for (Run in 1: length(RunFolders)) {
        FASTQfiles<-list.files(file.path(path, RunFolders[Run]), pattern="Undetermined_")    
        I<-ShortRead::readFastq(file.path(path, RunFolders[Run], FASTQfiles[1]))
        R1<-ShortRead::readFastq(file.path(path, RunFolders[Run], FASTQfiles[2]))
        R1_New<-combineFastqSequences(R1, I)

        
        CombinedFastqsPath<-file.path(path,"../../..","CrossPacific_16s", "ModifiedMultiplexedFastqs", RunFolders[Run])
        dir.create(CombinedFastqsPath, recursive=TRUE)
        writeFastq(object=R1_New, file= file.path(CombinedFastqsPath, paste0("Combined_R1-I1.fastq.gz")) )

        print(paste0("Written combined fastq from Run ", Run))

    }

# 3) Demultiplex

    #loop over all runs
    for (Run in 1:length(RunFolders)) {
            BarcodesPath<-file.path(path, "Barcodes",AllBarcodes[Run])
            R1Path<-file.path(path,"../../..","CrossPacific_16s", "ModifiedMultiplexedFastqs", RunFolders[Run], "Combined_R1-I1.fastq.gz")
            R2filename<-list.files(file.path(path, RunFolders[Run]), pattern="Undetermined_")[3]
            R2Path<-file.path( path, RunFolders[Run], R2filename)

            #barcode prep
                #make df of barcodes from this adapter
                Barcodes_df<-read.csv(BarcodesPath)

                #make barcode fasta path
                BarcodeFastaPath<-file.path(path, "Barcodes",paste0("Run",Run,".fasta"))
                #write fastas to this location
                seqinr::write.fasta(sequences= as.list(Barcodes_df$barcode.sequence), 
                                    names= Barcodes_df$sample, 
                                    file.out= BarcodeFastaPath, 
                                    open = "w", 
                                    nbchar = 200, 
                                    as.string = TRUE)

            #seq prefiltering
                # Create paths for prefiltered fastqs
                path_filtered<-file.path(path,"../../..","CrossPacific_16s", "ModifiedMultiplexedFastqs", RunFolders[Run], "PrefilteredFastqs", c("Combined_R1-I1.fastq.gz", R2filename))

                #prefilter fastqs
                #The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. Next we are going to “pre-filter” the sequences just to remove those with Ns, but perform no other filtering.
                out <- dada2::filterAndTrim(R1Path, path_filtered[1], R2Path, path_filtered[2], maxN=0, rm.phix=TRUE, compress=TRUE, multithread=TRUE, verbose=TRUE) 

            #demultiplexing
                # Create paths for demultiplexed fastqs
                path_demux<-file.path(path, "../../..", "CrossPacific_16s")

                #now run cutadapt
                system2( "cutadapt", args= c( "-e", 0.1, "--no-indels",  "--minimum-length", 1,  # options (10% error rate, no insertion/deletions, and zero length seqs discarded)
                                            "-g", paste0("file:",BarcodeFastaPath),       #barcodes
                                            "-o",   "/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific_16s/{name}_R1_001.fastq", "-p", "/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific_16s/{name}_R2_001.fastq",
                                        path_filtered[[1]] ,  path_filtered[[2]]                       #inputs
                                    )
            ) 
        }