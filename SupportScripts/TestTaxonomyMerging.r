#test taxonomies succesfully merge - i.e there are not duplicate entries - at higher ranks in particular

#1 set seed for replication
    set.seed(0.1)

#2 dependencies    
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load( tidyverse,
                    seqinr,
                    DECIPHER,
                    Biostrings
                    )

    combinefasta<-"ARMS_classifier_SubsampleTerrestrial_m256_bJAN21_combinedSeqs.fasta"
    path<-"../../BioinformaticPipeline_Env"
    seqs <- readDNAStringSet(file.path(path, "Data", "Classifiers", combinefasta) )


    splits<-strsplit(names(seqs), ";")
    Tax_df<-plyr::ldply(splits, rbind )


#change the number to see different ranks
tail(sort(table(Tax_df[3])), 100)




#descriptive dataframe


    df<-data.frame( Names=c(
                        "Number of Sequences",
                        "Number Domain Assigned",
                        "Number Phyla Assigned",
                        "Number Class Assigned",
                        "Number Order Assigned",
                        "Number Family Assigned",
                        "Number Genus Assigned",
                        "Number Species Assigned"
                ),
                    Answer=NA
    )
   #1 Total number of sequences in reference dataset
   df[1,2] <-   length(seqs)


    #2 How many have a taxonomic assignment at i) Phyla ii) Class iii) Order etc…. species level
    df[2,2] <-  df[1,2] - sum(is.na(Tax_df[2]), na.rm=TRUE)
    df[3,2] <-  df[1,2] - sum(is.na(Tax_df[3]), na.rm=TRUE)
    df[4,2] <-  df[1,2] - sum(is.na(Tax_df[4]), na.rm=TRUE)
    df[5,2] <-  df[1,2] - sum(is.na(Tax_df[5]), na.rm=TRUE)
    df[6,2] <-  df[1,2] - sum(is.na(Tax_df[6]), na.rm=TRUE)
    df[7,2] <-  df[1,2] - sum(is.na(Tax_df[7]), na.rm=TRUE)
    df[8,2] <-  df[1,2] - sum(is.na(Tax_df[8]), na.rm=TRUE)


    # 3. How many sequences available per main phylum? Looking for metazoan phyla in particular 
    # (e.g. Arthropoda, Echinodermata, Mollusca, Bryozoa, Porifera etc…)
    phyla<-sort(table(Tax_df[3])) %>% names
    SeqPerPhylaDf<-data.frame("phyla"=phyla, "Number Of Seqs"=NA )
    for (phylum in phyla) {
            SeqPerPhylaDf[SeqPerPhylaDf$phyla==phylum,2]    <-Tax_df[Tax_df[3]==phylum,] %>% dim %>% '['(1)

    }

    # 4 List of taxonomy assignments at genus and species level alone from the COI reference dataset 
    # (e.g. list of the fasta headers from sequences assigned at those ranks) if possible. I wanted to
    # see if species previously found in Chagos via morphological identification are represented in our reference dataset.
    Tax_df %>% str

#save csvs

    write.csv(x=df, file=file.path(path, "Data", "Classifiers", "ARMS_classifier_SubsampleTerrestrial_combinedSeqs_description.csv"))

    write.csv(SeqPerPhylaDf, file=file.path(path, "Data", "Classifiers", "ARMS_classifier_SubsampleTerrestrial_combinedSeqs_SequencesPerPhyla.csv"))

    write.csv(Tax_df, file=file.path(path, "Data", "Classifiers", "ARMS_classifier_SubsampleTerrestrial_combinedSeqs_AllTaxa.csv"))


