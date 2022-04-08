#script to assess performaance of different assignment methods on outut of SequenceEvolution Simulation.

set.seed(1)

#parameters
    path<-"../../BioinformaticPipeline_Env/Data/SimulatedData/IntermediateOutputs"

    numranks<-7
    queryfraction<-0.05

    fractions<-c(0.005, 0.01, 0.02, 0.03, 0.04, 0.05)
    BLASTthresholds<-seq(from=85, to=97, by=2)
    IdtaxaThresholds<-seq(from=0, to=100, by =10)

    #either FALSE or the proportion
    #MislabelledSequences<-

#functions
    #duplicated from TrainMidori+Biocode classifier script
        TrainIdtaxaClassifier<-function(seqs, taxid=NULL, maxGroupSize, maxIterations, allowGroupRemoval=FALSE) {
            
            #get some group objects to manipulate
            groups<-names(seqs)
            groupCounts <- table(names(seqs))
            u_groups <- names(groupCounts) # unique groups

            #set maximum number of sequences assigned to each taxa 
            maxGroupSize <- maxGroupSize # max sequences per label (>= 1)
            
            #identify seqs for removal based on this max group size
            remove <- logical(length(seqs))
            for (i in which(groupCounts > maxGroupSize)) {
                index <- which(groups==u_groups[i])
                keep <- sample(length(index),
                                maxGroupSize)
                remove[index[-keep]] <- TRUE   
            }

            #classifier training loop
            probSeqsPrev <- integer() # suspected problem sequences from prior iteration
            for (i in seq_len(maxIterations)) {
                cat("Training iteration: ", i, "\n", sep="")
                
                # train the classifier
                trainingSet <- LearnTaxa(train=seqs[!remove],
                                        taxonomy=names(seqs)[!remove],
                                        verbose=TRUE)
                
                #something prevents this line working on the HPC - it does work locally
                #save(trainingSet, file=file.path(parent.frame()$path, "Data", "Classifiers", paste0(parent.frame()$libraryname, "_IdtaxaClassifier_Iteration:", i, ".Rdata")))
                                        
                # look for problem sequences
                probSeqs <- trainingSet$problemSequences$Index
                if (length(probSeqs)==0) {
                    cat("No problem sequences remaining.\n")
                    break
                } else if (length(probSeqs)==length(probSeqsPrev) && all(probSeqsPrev==probSeqs)) {
                    cat("Iterations converged.\n")
                    break
                }
                if (i==maxIterations)
                    break
                probSeqsPrev <- probSeqs
                
                # remove any problem sequences
                index <- which(!remove)[probSeqs]
                remove[index] <- TRUE # remove all problem sequences
                if (!allowGroupRemoval) {
                    # replace any removed groups
                    missing <- !(u_groups %in% groups[!remove])
                    missing <- u_groups[missing]
                    if (length(missing) > 0) {
                        index <- index[groups[index] %in% missing]
                        remove[index] <- FALSE # don't remove
                    }
                }
            }

            return(list(trainingSet, probSeqs))
        }





QueryIndices<-sample(size=floor(length(dna)*queryfraction), x=1:length(dna), replace=FALSE)
QueryDna<-dna[QueryIndices]



               QueryDnaFasta<-file.path(path, "QuerySimulatedSequences.fasta")
                
                writeXStringSet(QueryDna, QueryDnaFasta, format="fasta", width=10000)





fractions_df<-data.frame( Rank= c(),
                    PercentCorrect= c(),
                    PercentIncorrect= c(),
                    FractionBarcoded= c(),
                    Threshold= c()
                    )

for (barcodedFraction in fractions) {
    RefIndices<-sample(size=floor(length(dna)*barcodedFraction), x=1:length(dna), replace=FALSE)
    RefDna<-dna[RefIndices]


    #BLAST
        BlastRefDna<-RefDna
        RefDnaIdMap<-data.frame(RefSeqTaxa=names(BlastRefDna),
                                RefSeqID=paste0("RefEsv_", 1:length(BlastRefDna) ))
        names(BlastRefDna)<- RefDnaIdMap$RefSeqID
        RefDnaFasta<-file.path(path, paste0( barcodedFraction, "RefSimulatedSequences.fasta"))

        writeXStringSet(BlastRefDna, RefDnaFasta, format="fasta", width=10000)

        system(command= paste0( "~/miniconda3/bin/makeblastdb -in ", RefDnaFasta, " -parse_seqids -blastdb_version 5 -title 'RefSimulatedSequencesBlastDB' -out '", path, "/", barcodedFraction, "RefSimulatedSequencesBlastDB' -dbtype nucl"))


        system(command= paste0( "~/miniconda3/bin/blastn ",
                                "-db ", file.path(path, paste0(barcodedFraction,"RefSimulatedSequencesBlastDB")),
                                " -query ", file.path(QueryDnaFasta),
                                " -out ", file.path(path, paste0("BlastOutput_", barcodedFraction, ".out")),
                                " -outfmt '6 qseqid sseqid pident evalue qcovs' ", 
                                " -evalue 1e-10" #expect value must be higher than this
                                ) 
            )

    
        #read output into r
            blast_output<-read.table( file.path(path, paste0("BlastOutput_", barcodedFraction, ".out")), header=FALSE, sep="\t")

        #adjust names blastoutput
            names(blast_output)<-c("QuerySeqID", "RefSeqID", "Blast_percentIdentical", "Blast_evalue", "Blast_query_coverage")
    
        #combine taxa names and blast_output
            FullBlastOutput<-left_join(blast_output, RefDnaIdMap, by = "RefSeqID")
            FullBlastOutput<-FullBlastOutput[!names(FullBlastOutput)=="SeqID"]

    
        #pick blast top hits
            blastTopHits<-FullBlastOutput[!duplicated(FullBlastOutput$QuerySeqID),]

        #blast thresholding loop
        BLASTthresholds_df<-data.frame( Rank= c(),
                                        PercentCorrect= c(),
                                        PercentIncorrect= c(),
                                        FractionBarcoded= c(),
                                        Threshold= c(),
                                        Assignment=c()
                                        )


        for (threshold in BLASTthresholds) {

            #only take those above blast threshold
            blastTopHitsThresholded<-blastTopHits [blastTopHits$Blast_percentIdentical>threshold , ]


            Taxa_df<-cbind(  plyr::ldply(str_split(blastTopHitsThresholded$RefSeqTaxa, ";"), rbind), 
                    plyr::ldply(str_split(blastTopHitsThresholded$QuerySeqID, ";"), rbind)
            )
            names(Taxa_df)<-c("Root_Assigned", paste0("Rank_", 1:numranks, "_Assigned"), "Root_True",paste0("Rank_", 1:numranks, "_True") )
            Blast_df<-cbind(Taxa_df, blastTopHitsThresholded[,3:5])


            BLASTranks_df<-data.frame(  Rank= paste0(1:numranks, "_Rank"),
                                        PercentCorrect= rep(NA, numranks),
                                        PercentIncorrect= rep(NA, numranks),
                                        FractionBarcoded= rep(barcodedFraction, numranks),
                                        Threshold= rep(threshold, numranks),
                                        Assignment=rep("BLAST", numranks)
                                        )


            for (rank in 1:numranks) {
                BLASTranks_df$PercentCorrect[rank] <- sum ( Blast_df[[paste0("Rank_", rank, "_Assigned")]]==Blast_df[[paste0("Rank_", rank, "_True")]] ) / dim(blastTopHits)[1]
                BLASTranks_df$PercentIncorrect[rank] <-sum ( !Blast_df[[paste0("Rank_", rank, "_Assigned")]]==Blast_df[[paste0("Rank_", rank, "_True")]] ) / dim(blastTopHits)[1]
            }
            
            BLASTthresholds_df<-rbind(BLASTranks_df, BLASTthresholds_df)
        }

    #Idtaxa
        trainingSet<-TrainIdtaxaClassifier(seqs=RefDna, maxGroupSize=Inf, maxIterations=3)



        thresholds_df<-data.frame( Rank= c(),
                            PercentCorrect= c(),
                            PercentIncorrect= c(),
                            FractionBarcoded= c(),
                            Threshold= c(),
                            Assignment=c()
                            )
        for (threshold in IdtaxaThresholds) {
            ids <- IdTaxa(QueryDna,
                                    trainingSet[[1]],
                                    type="extended",
                                    strand="top",
                                    threshold=threshold,
                                    processors=NULL)


            ExtractFromIds<-function(list_item, category){
                    ListItemLength<-length(list_item[[category]])
                    if (ListItemLength < Nranks) {
                    list_item[[category]][(ListItemLength+1):Nranks] <-NA
                    }
                    return(list_item[[category]])
                }

            Nranks<-max( unlist ( lapply ( lapply( ids, '[[', 1), length)))


            TaxonVecList<-lapply(ids, ExtractFromIds, category="taxon")
            TaxonDf<-plyr::ldply(TaxonVecList, rbind)
            names(TaxonDf)<-c("ESV", "Root_Assigned", paste0("Rank_", 1:numranks, "_Assigned"))[1:length(names(TaxonDf))]

            ConfVecList<-lapply(ids, ExtractFromIds, category="confidence")
            ConfDf<-plyr::ldply(ConfVecList, rbind)
            names(ConfDf)<-c("ESV", "Root_Confidence", paste0("Rank_", 1:numranks, "_Confidence"))[1:length(names(ConfDf))]


            TrueIdDf<-plyr::ldply(str_split(names(ids), ";"), rbind)
            TrueIdDf$.id<-names(ids)
            names(TrueIdDf)<-c("Root_True", paste0("Rank_", 1:numranks, "_True"), "ESV")


            IdtaxaDf<-merge(merge(TaxonDf, ConfDf, by="ESV"), TrueIdDf, by="ESV")

                ranks_df<-data.frame( Rank= paste0(1:numranks, "_Rank"),
                            PercentCorrect= rep(NA, numranks),
                            PercentIncorrect= rep(NA, numranks),
                            FractionBarcoded= rep(barcodedFraction, numranks),
                            Threshold= rep(threshold, numranks),
                            Assignment=rep("Idtaxa", numranks)
                            )

                for (ranknum in 1:numranks) {

                    ranks_df$PercentCorrect[ranknum] <- sum( IdtaxaDf[[ paste0("Rank_", ranknum, "_Assigned") ]] == IdtaxaDf[[ paste0("Rank_", ranknum, "_True") ]], na.rm=TRUE) / dim(IdtaxaDf)[1]

                    assignedSeqsIndices<-   !(is.na(IdtaxaDf[[ paste0("Rank_", ranknum, "_Assigned") ]]) | lapply(str_split(IdtaxaDf[[ paste0("Rank_", ranknum, "_Assigned") ]], "_"), '[', 1)  == "unclassified")
                    ranks_df$PercentIncorrect[ranknum] <- sum( IdtaxaDf[[ paste0("Rank_", ranknum, "_Assigned") ]][assignedSeqsIndices] != IdtaxaDf[[ paste0("Rank_", ranknum, "_True") ]][assignedSeqsIndices], na.rm=TRUE) / dim(IdtaxaDf)[1]
            }
            
            thresholds_df<-rbind(ranks_df, thresholds_df)
        }

    fractions_df<- rbind(thresholds_df, BLASTthresholds_df, fractions_df)

}











pdf(file = "/home/j/Desktop/AssignedSimulatedESVs.pdf", width=20, height =8 ) # The height of the plot in inches    

fractions_df$assignmentThreshold<-paste0(fractions_df$Assignment , fractions_df$Threshold)
fractions_df$RankAssignment<-paste0(fractions_df$Assignment , fractions_df$Rank)
fractions_df$ThreholdTightness<-as.factor(fractions_df$Threshold )
levels(fractions_df$ThreholdTightness)<-c(7, 6, 5,4,3,2,1,1,1,7,6,5,1,4,3,2,1,1)
fractions_df$ThreholdTightness<- as.numeric(as.character(fractions_df$ThreholdTightness))




ggplot(fractions_df, aes(x=PercentIncorrect, y=PercentCorrect, color=Rank, shape=Assignment)) +
    geom_path(aes(group=RankAssignment, linetype=Assignment), size=1, alpha=0.7) +
    geom_point(aes(color=Rank, size=ThreholdTightness)) +
    #geom_text(aes(label=as.character(Threshold)),hjust=0,vjust=0, size=2, color="Black") +
    coord_cartesian(xlim = c(0, 1), ylim=c(0, 1)) +
    facet_wrap(~FractionBarcoded) +
    ggtitle("Simulated ESVs assigned using Idtaxa and BLAST for different fractions of lowest rank groups barcoded - larger icons indicate looser thresholds for Idtaxa (0-100 confidence) and BLAST (85-97% similarity)")



dev.off()