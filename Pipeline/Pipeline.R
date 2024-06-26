####################################
# CUTADAPT>DADA2>LULU>CLUSTER(SWARM2/VSEARCH)>LULU>IDTAXA pipeline,  (CDLCLI_Pipeline)

# Functionalised modular pipeline. Each processing step is its own function which outputs the primary result and 
# secondary results of various outputs such as figures and tables detailing the processing step. The primary 
# result of each step is input to an SQL database, the secondary results are saved to a csv or pdf as appropriate.

# DO NOT adjust this pipeline for the purpose of analysing a single dataset. Adjustments here should instead reflect
# optimisation of the general pipeline. 

# To parameterise the pipeline for a specific dataset please create a new ControlScript.


# set CRAN mirror
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)


# A set seed for replication
set.seed(0.1)

# B dependencies unavailable through conda (see environment.yml for full list) 
#install r packages unavailable through conda
devtools::install_github("tobiasgf/lulu")  
#load packages installed with conda
library(tidyverse)
library(gridExtra)
library(seqinr)
library(ShortRead)
library(Biostrings)
library(dada2)
library(lulu)
library(DECIPHER)
# Functions
function_files<-list.files(file.path(path, "BioinformaticPipeline", "Pipeline", "Functions"))
sapply(file.path(path, "BioinformaticPipeline", "Pipeline", "Functions",function_files),source)


# C identify where the pipeline should start
StartingStep<-IdentifyLastInputPresent(c("CutadaptOutput", "DadaOutput", "LuluOutput1","ClusterOutput", "LuluOutput2", "IdtaxaOutput", "BlastOutput"))

# check if dataPath is set (i.e. fastqs are stored externally), if it is not then default to standard path for pipeline
if (is.null(dataPath)) { dataPath<-path} 

# 1 Cutadapt
if (StartingStep <= 1) {
    # 1.0 process inputs
    # none required 
        
    # 1.1 Run Module
    CutadaptOutput<-RunCutadapt(
        FWD =FWD,  
        REV = REV, 
        multithread=TRUE,
        dataname=dataname,
        UseCutadapt=UseCutadapt)

    # 1.2 cache Output
    CacheOutput(CutadaptOutput)
}

# 2 DADA2
if (StartingStep<=2){
    # 2.0 process inputs
    # none required 
        
    # 2.1 Run Module
    DadaOutput<-RunDADA2(   
        truncLen=truncLen, 
        trimLeft=trimLeft, 
        maxN=maxN, 
        maxEE=maxEE, 
        truncQ=truncQ,
        DesiredSequenceLengthRange=DesiredSequenceLengthRange,
        dataname=dataname,
        multithread=multithread,
        pool=pool,
        MixedOrientation=MixedOrientation,
        NumberOfRuns=NumberOfRuns,
        ReadType=ReadType)

    # 2.2 cache Output
    CacheOutput(DadaOutput)
}

# 3 LULU1
if (StartingStep<=3){
    # 3.0 process inputs
    ConfirmInputsPresent("DadaOutput")
    OtuTableForLulu<-CreateOtuTableForLulu(Input=DadaOutput$SeqDataTable, clustering = "ESV")
    MatchListForLulu<-CreateMatchlistForLulu(Input=DadaOutput$SeqDataTable ,MatchRate1, clustering="ESV", multithread=multithread)

    # 3.1 Run Module
    LuluOutput1<-RunLULU(TableToMergeTo=DadaOutput$SeqDataTable, MatchRate=MatchRate1, MinRelativeCo=MinRelativeCo1, RatioType=RatioType1, clustering="ESV") # inputs are Otutable and matchlist created previously

    # 3.2 cache Output
    CacheOutput(LuluOutput1)
}

# 4 Clustering
if (StartingStep<=4){
    # 4.0 process inputs
    ConfirmInputsPresent("DadaOutput")
    ConfirmInputsPresent("LuluOutput1")

    # 4.1 Run Module in this case you have a choice between vsearch and swarm (algorithms which 
    # efficiently implement complete linkage clustering and single linkage clustering respectively.)
    ClusterOutput<-ClusterOTUs( 
        linkage=linkage, #either 'complete' (vsearch) or 'single' (swarm)
        differences=differences,
        threads=threads,
        TableToMergeTo=LuluOutput1,
        clustering="ESV",
        SimilarityThreshold=SimilarityThreshold,
        multithread=multithread)

    # 4.2 cache Output
    CacheOutput(ClusterOutput)
}

# 5 LULU2
if (StartingStep<=5){
    #5.0 process inputs
    ConfirmInputsPresent("DadaOutput")
    ConfirmInputsPresent("LuluOutput1")
    ConfirmInputsPresent("ClusterOutput")

    OtuTableForLulu<-CreateOtuTableForLulu(Input=ClusterOutput, clustering="OTU")
    MatchListForLulu<-CreateMatchlistForLulu(Input=ClusterOutput ,MatchRate2, clustering="OTU", multithread=multithread)

    # 5.1 Run Module
    LuluOutput2<-RunLULU(TableToMergeTo=ClusterOutput, MatchRate=MatchRate2, MinRelativeCo=MinRelativeCo2, RatioType=RatioType2, clustering="OTU") # inputs are Otutable and matchlist created previously

    # 5.2 cache Output
    CacheOutput(LuluOutput2)
}


# 6 IDTAXA
if (StartingStep<=6){
    # 6.0 process inputs
    ConfirmInputsPresent("DadaOutput")
    ConfirmInputsPresent("LuluOutput1")
    ConfirmInputsPresent("ClusterOutput")
    ConfirmInputsPresent("LuluOutput2")

    # 6.1 Run Module
    IdtaxaOutput<-RunIdtaxa(IDTAXA=IDTAXA, trainingSet=trainingSet, TableToMergeTo=LuluOutput2, SeqsToAssign=SeqsToAssign, threshold=threshold, parallel=parallel, multithread=multithread)
            
    # 6.2 cache Output
    CacheOutput(IdtaxaOutput)
}

#7 Run Blast
if (StartingStep<=7){
    #7.0 process inputs
    ConfirmInputsPresent("CutadaptOutput")
    ConfirmInputsPresent("DadaOutput")
    ConfirmInputsPresent("LuluOutput1")
    ConfirmInputsPresent("ClusterOutput")
    ConfirmInputsPresent("LuluOutput2")
    ConfirmInputsPresent("IdtaxaOutput")

    #7.1 run blast if desired
    if (Blast){
        BlastOutput<-RunBLAST(dbname=dbname, clustering="ESV", TableToMergeTo=IdtaxaOutput$SeqDataTable, assignmentThresholds=assignmentThresholds, Blastdesiredranks= Blastdesiredranks) 
    } else {
        BlastOutput<-NULL
    }

    # 7.2 cache Output
        CacheOutput(BlastOutput)
}

# 8 Creating final results from intermediate outputs

# process inputs
ConfirmInputsPresent("CutadaptOutput")
ConfirmInputsPresent("DadaOutput")
ConfirmInputsPresent("LuluOutput1")
ConfirmInputsPresent("ClusterOutput")
ConfirmInputsPresent("LuluOutput2")
ConfirmInputsPresent("IdtaxaOutput")
ConfirmInputsPresent("BlastOutput")

# 8.1 save sequences and clusters to results

# save dada plots            
glist <- lapply(DadaOutput$SecondaryOutputs$DadaPlots, ggplotGrob)
ggsave(file.path(path,"Results",paste0(dataname,"_DadaPlots.pdf")), marrangeGrob(glist, nrow = 1, ncol = 1))

# save dada and overall clustering tables
write.csv( CutadaptOutput$PrimerCountSamplesList, file=file.path(path,"Results",paste0(dataname,"_PrimerCounts.csv")) )
write.csv( CutadaptOutput$PrimerCountCutadaptedSamplesList, file=file.path(path,"Results",paste0(dataname,"_PrimerCountsAfterCutadapt.csv")) )
write.csv( DadaOutput$SecondaryOutputs$DadaTables, file=file.path(path,"Results",paste0(dataname,"_DadaTable.csv")) )
write.csv( DadaOutput$SecondaryOutputs$SeqLengthDist, file=file.path(path,"Results",paste0(dataname,"_DadaSeqLengthDistribution.csv")) )

WriteClusteringTable(FinalOutput=LuluOutput2, pipeline="DLSL")
       
if ( ! is.null(BlastOutput)) {
    saveRDS(BlastOutput$SeqDataTable, file=file.path(path, "Results", paste0(dataname,"_SeqDataTable.RDS")))
    write.csv(BlastOutput$SeqDataTable, file=file.path(path, "Results", paste0(dataname,"_SeqDataTable.csv")))
    write.csv( BlastOutput$FullBlastOutput, file=file.path(path,"Results",paste0(dataname,"_FullBlastOutput.csv")) )
    pdf(file = file.path(path, "Results", paste0(dataname, "_TaxaAssignment.pdf")))   # The directory you want to save the file in
        plot(IdtaxaOutput$PlotData)
    dev.off()
} else if ( ! is.null(IdtaxaOutput) ) {
    saveRDS(IdtaxaOutput$SeqDataTable, file=file.path(path, "Results", paste0(dataname,"_SeqDataTable.RDS")))
    write.csv(IdtaxaOutput$SeqDataTable, file=file.path(path, "Results", paste0(dataname,"_SeqDataTable.csv")))
    pdf(file = file.path(path, "Results", paste0(dataname, "_TaxaAssignment.pdf")))   # The directory you want to save the file in
        plot(IdtaxaOutput$PlotData)
    dev.off()
} else {
    saveRDS((LuluOutput2), file=file.path(path, "Results", paste0(dataname,"_SeqDataTable.RDS")))
    write.csv(LuluOutput2, file=file.path(path, "Results", paste0(dataname,"_SeqDataTable.csv")))
}
