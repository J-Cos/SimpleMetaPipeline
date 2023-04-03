#make blastDB for MetaFish
#1 set seed for replication
    set.seed(0.1)

#2 dependencies    
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load( tidyverse,
                    seqinr,
                    DECIPHER,
                    Biostrings
                    )

    seqs <- readDNAStringSet( "../Data/Raw/MetaFish.12s.miya.dada.taxonomy.v253.fasta")

    blastTaxidMap<-data.frame(  SeqId=paste0("Seq", 1:length(seqs)),
                                Taxa=names(seqs))
    
    names(seqs)<-blastTaxidMap$SeqId

    writeXStringSet(seqs, "/home/j/Dropbox/BioinformaticPipeline_Env/Data/BlastDBs/MetafishBlast.fasta", format="fasta", width=10000)
    write.table(blastTaxidMap, "/home/j/Dropbox/BioinformaticPipeline_Env/Data/BlastDBs/MetafishTaxidMap.txt", quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)


#exit R

~/miniconda3/bin/makeblastdb -in "/home/j/Dropbox/BioinformaticPipeline_Env/Data/BlastDBs/MetafishBlast.fasta" -parse_seqids -blastdb_version 5 -title "MetafishBlastDB" -out "MetafishBlastDB" -dbtype nucl





           