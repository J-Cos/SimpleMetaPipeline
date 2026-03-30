#test taxonomies succesfully merge - i.e there are not duplicate entries - at higher ranks in particular

#1 set seed for replication
    set.seed(0.1)

# parameters
    #fastaname<-"ARMS_classifier_combinedSeqs.fasta"
    #blastdbname<-"MidoriBiocodeBlast"
    fastaname<-"SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta"
    blastdbname<-"SILVA138_BlastDB"

#2 dependencies    
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load( tidyverse,
                    seqinr,
                    DECIPHER,
                    Biostrings
                    )

    seqs <- readDNAStringSet( paste0("/home/j/Dropbox/BioinformaticPipeline_Env/Data/Raw/", fastaname))

#fix taxanames as required
    #Taxanames<-names(seqs)
    Taxanames<-str_split(names(seqs), " ", 2)

     Taxanames<-unlist(lapply(Taxanames, '[', 2)) 


#make taxid
    blastTaxidMap<-data.frame(  SeqId=paste0("Seq", 1:length(seqs)),
                                Taxa=Taxanames)

    
    
    names(seqs)<-blastTaxidMap$SeqId
    writeXStringSet(seqs, paste0("/home/j/Dropbox/BioinformaticPipeline_Env/Data/BlastDBs/", blastdbname, ".fasta"), format="fasta", width=10000)
    write.table(blastTaxidMap, paste0("/home/j/Dropbox/BioinformaticPipeline_Env/Data/BlastDBs/", blastdbname, "_TaxidMap.txt"), quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)



~/miniconda3/bin/makeblastdb -in "/home/j/Dropbox/BioinformaticPipeline_Env/Data/BlastDBs/MidoriBiocodeBlast.fasta" -parse_seqids -blastdb_version 5 -title "MidoriBiocodeBlastDB" -out "MidoriBiocodeBlastDB" -dbtype nucl

~/miniconda3/bin/makeblastdb -in "/home/j/Dropbox/BioinformaticPipeline_Env/Data/BlastDBs/SILVA138_BlastDB.fasta" -parse_seqids -blastdb_version 5 -title "SILVA138_BlastDB" -out "SILVA138_BlastDB" -dbtype nucl




           