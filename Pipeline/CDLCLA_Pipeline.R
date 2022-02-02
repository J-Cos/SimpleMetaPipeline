####################################
# CUTADAPT>DADA2>LULU>CLUSTER(SWARM2/VSEARCH)>LULU>IDTAXA pipeline,  (CDLCLI_Pipeline)

# Functionalised modular pipeline. Each processing step is its own function which outputs the primary result and 
# secondary results of various outputs such as figures and tables detailing the processing step. The primary 
# result of each step is input to an SQL database, the secondary results are saved to a csv or pdf as appropriate.

# DO NOT adjust this pipeline for the purpose of analysing a single dataset. Adjustments here should instead reflect
#optimisation of the general pipeline. 

# to parameterise the pipeline for a specific dataset please create a new ParameterSet script.

#1 set seed for replication
    set.seed(0.1)

#2 dependencies    
    #installed software
        #vsearch
        #swarm v2
    #R packages
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load( dada2,
                    seqinr,
                    DECIPHER,
                    tidyverse,
                    lulu,
                    ggplot2,
                    gridExtra,
                    ShortRead,
                    Biostrings,
                    insect)

    #Functions
        function_files<-list.files(file.path(path, "BioinformaticPipeline", "Pipeline", "Functions"))
        sapply(file.path(path, "BioinformaticPipeline", "Pipeline", "Functions",function_files),source)


#3 Cutadapt
    # 3.0 process inputs
        #none required 
    
    # 3.1 Run Module
    CutadaptOutput<-RunCutadapt(FWD =FWD,  ## CHANGE ME to your forward primer sequence
                REV = REV,  ## CHANGE ME...
                multithread=TRUE,
                dataname=dataname,
                UseCutadapt=UseCutadapt)

        #inspect outputs
            #CutadaptOutput$PrimerCountSamplesList
            #CutadaptOutput$PrimerCountCutadaptedSamplesList

    # 3.2 cache Output
        CacheOutput(CutadaptOutput)

#3 DADA2
    # 3.0 process inputs
        #none required 
    
    # 3.1 Run Module
        DadaOutput<-RunDADA2(   truncLen=truncLen, 
                                trimLeft=trimLeft, 
                                maxN=maxN, 
                                maxEE=maxEE, 
                                truncQ=truncQ,
                                DesiredSequenceLengthRange=DesiredSequenceLengthRange,
                                dataname=dataname,
                                multithread=multithread,
                                pool=pool)

        #inspect outputs
            #DadaOutput$SeqDataTable
            #DadaOutput$SecondaryOutputs$DadaPlots
            #DadaOutput$SecondaryOutputs$DadaTables
            #DadaOutput$SecondaryOutputs$SeqLengthDist

    # 3.2 cache Output
        CacheOutput(DadaOutput)


#4 LULU
    #4.0 process inputs
        ConfirmInputsPresent("DadaOutput")

        OtuTableForLulu<-CreateOtuTableForLulu(Input=DadaOutput$SeqDataTable, clustering = "ESV")

        MatchListForLulu<-CreateMatchlistForLulu(Input=DadaOutput$SeqDataTable ,MatchRate, clustering="ESV", HPC=HPC)
    
    # 4.1 Run Module
        LuluOutput1<-RunLULU(TableToMergeTo=DadaOutput$SeqDataTable, MatchRate=MatchRate, MinRelativeCo=MinRelativeCo, RatioType=RatioType, clustering="ESV") # inputs are Otutable and matchlist created previously

    # 4.2 cache Output
        CacheOutput(LuluOutput1)


#5 Clustering
    #5.0 process inputs
        #in general this will always be checking the previous modules output is present
        # then coverting it to required input formats for next module
        #if DadaOutput not loaded then read it in
            ConfirmInputsPresent("DadaOutput")
            ConfirmInputsPresent("LuluOutput1")

    # 5.1 Run Module in this case you have a choice between vsearch and swarm (algorithms which 
    # efficiently implement complete linkage clustering and single linkage clustering respectively.)
        ClusterOutput<-ClusterOTUs( linkage=linkage, #either 'complete' (vsearch) or 'single' (swarm)
                                    differences=differences,
                                    threads=threads,
                                    TableToMergeTo=LuluOutput1,
                                    clustering="curatedESV",
                                    HPC=HPC,
                                    SimilarityThreshold=SimilarityThreshold)

        #inspect outputs
            #ClusterOutput

    # 5.2 cache Output
        CacheOutput(ClusterOutput)
    
#5 LULU
    #5.0 process inputs
        ConfirmInputsPresent("DadaOutput")
        ConfirmInputsPresent("LuluOutput1")
        ConfirmInputsPresent("ClusterOutput")

        OtuTableForLulu<-CreateOtuTableForLulu(Input=ClusterOutput, clustering="OTU")

        MatchListForLulu<-CreateMatchlistForLulu(Input=ClusterOutput ,MatchRate, clustering="OTU", HPC=HPC)
    
    # 5.1 Run Module
        LuluOutput2<-RunLULU(TableToMergeTo=ClusterOutput, MatchRate=MatchRate, MinRelativeCo=MinRelativeCo, RatioType=RatioType, clustering="OTU") # inputs are Otutable and matchlist created previously

    # 5.2 cache Output
        CacheOutput(LuluOutput2)
    


# 6 IDTAXA
    #6.0 process inputs
        ConfirmInputsPresent("DadaOutput")
        ConfirmInputsPresent("LuluOutput1")
        ConfirmInputsPresent("ClusterOutput")
        ConfirmInputsPresent("LuluOutput2")

        #load reference library

    # 6.1 Run Module
            IdtaxaOutput<-RunIdtaxa(Type=Type, trainingSet=trainingSet, TableToMergeTo=LuluOutput2, SeqsToAssign=SeqsToAssign, threshold=threshold)
            
    # 6.2 cache Output
        CacheOutput(IdtaxaOutput)
    
# 7 Creating final results from intermediate outputs
    #process inputs
        ConfirmInputsPresent("CutadaptOutput")
        ConfirmInputsPresent("DadaOutput")
        ConfirmInputsPresent("LuluOutput1")
        ConfirmInputsPresent("ClusterOutput")
        ConfirmInputsPresent("LuluOutput2")
        ConfirmInputsPresent("IdtaxaOutput")

    # 7.1 save sequences and clusters to results
        #SaveSequenceDataTableToDataBase(Input=IdtaxaOutput)

        #save dada plots            

            glist <- lapply(DadaOutput$SecondaryOutputs$DadaPlots, ggplotGrob)
            ggsave(file.path(path,"Results",paste0(dataname,"_DadaPlots.pdf")), marrangeGrob(glist, nrow = 1, ncol = 1))
            
        #save dada and overall clustering tables
            
            write.csv( CutadaptOutput$PrimerCountSamplesList, file=file.path(path,"Results",paste0(dataname,"_PrimerCounts.csv")) )
            write.csv( CutadaptOutput$PrimerCountCutadaptedSamplesList, file=file.path(path,"Results",paste0(dataname,"_PrimerCountsAfterCutadapt.csv")) )

            write.csv( DadaOutput$SecondaryOutputs$DadaTables, file=file.path(path,"Results",paste0(dataname,"_DadaTable.csv")) )
            
            write.csv( DadaOutput$SecondaryOutputs$SeqLengthDist, file=file.path(path,"Results",paste0(dataname,"_DadaSeqLengthDistribution.csv")) )

            WriteClusteringTable(FinalOutput=LuluOutput2, pipeline="DLSL")

            if ( ! is.null(IdtaxaOutput) ) {
                saveRDS(IdtaxaOutput$SeqDataTable, file=file.path(path, "Results", paste0(dataname,"_SeqDataTable.RDS")))
                write.csv(IdtaxaOutput$SeqDataTable, file=file.path(path, "Results", paste0(dataname,"_SeqDataTable.csv")))
                
                pdf(file = file.path(path, "Results", paste0(dataname, "_TaxaAssignment.pdf")))   # The directory you want to save the file in
                plot(IdtaxaOutput$PlotData)
                dev.off()
            } else {
                saveRDS((LuluOutput2), file=file.path(path, "Results", paste0(dataname,"_SeqDataTable.RDS")))
                write.csv(LuluOutput2, file=file.path(path, "Results", paste0(dataname,"_SeqDataTable.csv")))

            }
