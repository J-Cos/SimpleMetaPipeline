####################################
# DADA2>SWARM2>LULU>IDTAXA pipeline, outputting to SQLite (DSLI_SQL_Pipeline)

# Functionalised modular pipeline. Each processing step is its own function which outputs the primary result and 
# secondary results of various outputs such as figures and tables detailing the processing step. The primary 
# result of each step is input to an SQL database, the secondary results are saved to a csv or pdf as appropriate.

# DO NOT adjust this pipeline for the purpose of analysing a single dataset. Adjustments here should isntead reflect
#optimisation of the general pipeline. 

# to parameterise the pipeline for a specific dataset please create a new ParameterSet script.

#1 set seed for replication
    set.seed(0.1)

#2 dependencies    
    #installed software
        #vsearch
        #swarm v2
    
    #R Packages
        library(dada2)
        library(seqinr)
        library(DECIPHER)
        library(tidyverse)
        library(lulu)
        #library(Biostrings)
        #library(gsubfn)
        library(ggplot2)
        library(gridExtra)

    #Functions
    function_files<-list.files(file.path(path, "BioinformaticPipeline", "Pipeline", "Functions"))
    sapply(file.path(path, "BioinformaticPipeline", "Pipeline", "Functions",function_files),source)

#3 DADA2
    # 3.0 process inputs
        #none required 
    
    # 3.1 Run Module
        DadaOutput<-RunDADA2(   truncLen=truncLen, 
                                trimLeft=trimLeft, 
                                maxN=maxN, 
                                maxEE=maxEE, 
                                truncQ=truncQ, 
                                dataname=dataname,
                                multithread=multithread,
                                pool=pool)

        #inspect outputs
            #DadaOutput$ESVtable
            #DadaOutput$ESVsequences
            #DadaOutput$SecondaryOutputs$DadaPlots
            #DadaOutput$SecondaryOutputs$DadaTables

    # 3.2 cache Output
        CacheOutput(DadaOutput)

#4 SWARM2
    #4.0 process inputs
        #in general this will always be checking the previous modules output is present
        # then coverting it to required input formats for next module
        #if DadaOutput not loaded then read it in
            ConfirmInputsPresent("DadaOutput")
        #convert ESV sequences to fasta with abundances for input into swarmv2    
            CreateFastaWithAbundances(ESVtable=DadaOutput$ESVtable, ESVsequences=DadaOutput$ESVsequences)


    # 4.1 Run Module
            SwarmOutput<-RunSwarm(differences=differences,  # input is FastaWithAbundances created previously
                                    threads =4) 
        #inspect outputs
            #SwarmOutput

    # 4.2 cache Output
        CacheOutput(SwarmOutput)
    
#5 LULU
    #5.0 process inputs
        ConfirmInputsPresent("DadaOutput")
        ConfirmInputsPresent("SwarmOutput")

        OtuTableForLulu<-CreateOtuTableForLulu(Input=SwarmOutput)

        MatchListForLulu<-CreateMatchlistForLulu(Input=SwarmOutput ,MatchlistRate)
    
    # 5.1 Run Module
        LuluOutput<-RunLULU(TableToMergeTo=SwarmOutput) # inputs are Otutable and matchlist created previously
        #inspect outputs

    # 5.2 cache Output
        CacheOutput(LuluOutput)
    


# 6 IDTAXA
    #6.0 process inputs
        ConfirmInputsPresent("DadaOutput")
        ConfirmInputsPresent("SwarmOutput")
        ConfirmInputsPresent("LuluOutput")

        #load reference library
            trainingSet<-GetTrainingSet(Type=Type, RefLibrary= RefLibrary)

    # 6.1 Run Module
            IdtaxaOutput<-RunIdtaxa(trainingSet, TableToMergeTo=LuluOutput, SeqsToAssign=SeqsToAssign, threshold=threshold)
        #inspect outputs


    # 6.2 cache Output
        CacheOutput(IdtaxaOutput)
    
# 7 Creating final results from intermediate outputs
    #process inputs
            ConfirmInputsPresent("DadaOutput")
            ConfirmInputsPresent("SwarmOutput")
            ConfirmInputsPresent("LuluOutput")
            ConfirmInputsPresent("IdtaxaOutput")

    # 7.1 save sequences and clusters to database
        #SaveSequenceDataTableToDataBase(Input=IdtaxaOutput)
        saveRDS(IdtaxaOutput$SeqDataTable, file=file.path(path, "Results", paste0(dataname,"_SeqDataTable.RDS")))


        #make sure measures of each module are output here 

        #save plots            

            glist <- lapply(DadaOutput$SecondaryOutputs$DadaPlots, ggplotGrob)
            ggsave(file.path(path,"Results",paste0(dataname,"_DadaPlots.pdf")), marrangeGrob(glist, nrow = 1, ncol = 1))
            
            pdf(file = file.path(path, "Results", paste0(dataname, "_TaxaAssignment.pdf")))   # The directory you want to save the file in
            plot(IdtaxaOutput$IDTAXAplotdata, IdtaxaOutput$trainingSet)
            dev.off()

        #save tables
            write.csv( DadaOutput$SecondaryOutputs$DadaTables, file=file.path(path,"Results",paste0(dataname,"_DadaTable.csv")) )
            
            ClusteringTable<-data.frame(  
                "Number of ESVs"= dim(IdtaxaOutput$SeqDataTable)[1], 
                "Number of OTUs" = dim(unique(IdtaxaOutput$SeqDataTable[2]))[1], 
                "Number of cOTUs" = dim(unique(IdtaxaOutput$SeqDataTable[7]))[1])

            write.csv( ClusteringTable, file=file.path(path,"Results",paste0(dataname,"_ClusteringTable.csv")) )
