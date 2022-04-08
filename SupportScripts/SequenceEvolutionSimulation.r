#script to simulate evolution of some sequences for classifier testing

#parameters
    library(DECIPHER)
    library(tidyverse)

    path<-"../../BioinformaticPipeline_Env/Data/SimulatedData"

    rankSuffixes<-c("-iota", "-ota", "-ia", "-ales", "-aceae", "-us", "-species")
    numranks<-7
    minGroupsPerRank<-2
    maxGroupsPerRank<-6
    mutationRate<-10
    primerLength<-300


    bases<-c("a", "c", "t", "g")

#simulation
    Sequences<-list()
    df<-data.frame(
                ESV=rep(NA, 10^6),
                #Sequence=c(NA),
                Root=rep(NA, 10^6),
                Kingdom=rep(NA, 10^6),
                Phylum=rep(NA, 10^6),
                Class=rep(NA, 10^6),
                Order=rep(NA, 10^6),
                Family=rep(NA, 10^6),
                Genus=rep(NA, 10^6),
                Species=rep(NA, 10^6)
                )

    df$ESV[1]<-1
    #df$Sequence[1]<-paste(sample(bases, 300, replace=TRUE), collapse="")
    df$Root[1]<-"Root"
    df$Kingdom[1]<-NA
    df$Phylum[1]<-NA
    df$Class[1]<-NA
    df$Order[1]<-NA
    df$Family[1]<-NA
    df$Genus[1]<-NA
    df$Species[1]<-NA


    Sequences[[1]]<-sample(bases, primerLength, replace=TRUE)


    nextrow<-2
    for (rank in 1:numranks) {
        #seniorGroups<-unique(df[,rank+1])
        seniorIndices<-which(!is.na(df[,rank+1]))
        for (seniorIndex in seniorIndices) {
            #seniorIndex<-which(df[,rank+1]==seniorGroup  & df[,rank+1]==seniorGroup   )
            numJuniorGroups<-floor(runif(1, minGroupsPerRank,maxGroupsPerRank))

            for (juniorGroup in 1:numJuniorGroups) {
                Sequences[[nextrow]]<-Sequences[[seniorIndex]]
                numMutations<-floor(runif(1, 1, mutationRate))
                Sequences[[nextrow]][sample( primerLength, numMutations, replace=FALSE )]<-sample(bases, numMutations, replace=TRUE)

                df$ESV[nextrow]<-nextrow
                df$Root[nextrow]<-df$Root[seniorIndex]
                df$Kingdom[nextrow]<- df$Kingdom[seniorIndex]
                df$Phylum[nextrow]<- df$Phylum[seniorIndex]
                df$Class[nextrow]<- df$Class[seniorIndex]
                df$Order[nextrow]<- df$Order[seniorIndex]
                df$Family[nextrow]<- df$Family[seniorIndex]
                df$Genus[nextrow]<- df$Genus[seniorIndex]
                df$Species[nextrow]<- df$Species[seniorIndex] 
                if (rank==1) {
                    df[nextrow,rank+2]<-paste0( juniorGroup)
                } else if (rank >1) {
                    df[nextrow,rank+2]<-paste0( df[seniorIndex,rank+1],
                                                juniorGroup)
                }
                nextrow<-nextrow+1
            }
        }
        print(paste0("Rank ", rank, " complete"))
    }

#trim to modern seqs and save
    filled_df<-df[!is.na(df$ESV),]
    modern_Sequences<-Sequences[!is.na(filled_df$Species)]
    modern_df<-filled_df[!is.na(filled_df$Species),]
    # add suffixes
    for (rank in 1:numranks) {
        modern_df[,rank+2]<-paste0(modern_df[,rank+2],rankSuffixes[rank])
    }


    dna <- DNAStringSet ( unlist(lapply(modern_Sequences, paste,collapse="") ))
    names(dna)<-unlist(  unite(modern_df, "combined",2:9, sep = ";") [,2]     )

                    
    writeXStringSet(dna, file.path(path, "AllSimulatedSequences.fasta"), format="fasta", width=10000)
