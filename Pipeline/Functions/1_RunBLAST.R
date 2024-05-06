
# pipeline function to blast query sequences - defaults to ESVs, 
# other clusters may not have fastas auto-output from earlier steps to serve as inputs
# test<-RunBLAST( path="~/Dropbox/BioinformaticPipeline_Env",
#               dbname="BIOCODE_ER")

RunBLAST<-function(
        dbname, 
        clustering="ESV", 
        TableToMergeTo, 
        assignmentThresholds,
        Blastdesiredranks=NULL) {

    # query fasta 
    queryfasta<-paste0(dataname, "_", clustering, "_sequences.fasta")

    # blastn seqs
    system(command = paste0( 
        "blastn ",
        "-db ", 
        file.path(path, "Data", "BlastDBs", dbname, dbname),
        " -query ", 
        file.path(path, "IntermediateOutputs", queryfasta),
        " -out ", 
        file.path(path, "IntermediateOutputs", 
        paste0("BlastOutput_", dbname, ".out")),
        " -outfmt '6 qseqid sseqid pident evalue qcovs' ", 
        " -evalue 1e-10")) # expect value must be higher than this
    
    # read output into r
    blast_output<-read.table( file.path(path,"IntermediateOutputs", paste0("BlastOutput_", dbname, ".out")), header=FALSE, sep="\t")

    # adjust names blastoutput
    names(blast_output)<-c(clustering, "SeqID", "Blast_percentIdentical", "Blast_evalue", "Blast_query_coverage")
    
    # ensure SeqID column has same upper lower format as taxid (blast seems to modify this for certain rows... not sure why)
    blast_output$SeqID<-str_to_title(blast_output$SeqID)

    # map taxonomy df onto blast output
    # load file
    blastTaxidMap<- readLines(file.path(path, "Data", "BlastDBs", paste0(dbname,"_TaxidMap.txt")))
    # seperate id and create base df
    Taxid_components<-str_split(blastTaxidMap, " ", 2)
    Taxid_df<-data.frame(SeqID=unlist(lapply(Taxid_components, '[', 1)))
    # seperate taxa ranks and add into df
    taxa<-unlist(lapply(Taxid_components, '[', 2))
    taxaRanks_vectorList<-str_split(taxa, ";", simplify=TRUE)
    taxaRanks_df<-as.data.frame(taxaRanks_vectorList)
    names(taxaRanks_df)<-Blastdesiredranks
    Taxid_df<-cbind(Taxid_df, taxaRanks_df)

    # combine taxa names and blast_output
    FullBlastOutput<-left_join(blast_output, Taxid_df, by = "SeqID")
    FullBlastOutput<-FullBlastOutput[!names(FullBlastOutput)=="SeqID"]

    # pick blast top hits
    blastTopHits<-FullBlastOutput[!duplicated(FullBlastOutput[clustering]),]


    # use %age match to truncate assignments at appropriate rank
    for (row in 1:dim(blastTopHits)[1]) {
        item<-blastTopHits[row,]
        item[item==""]<-paste0("Unknown_",item[length(item)-sum(item=="")])
        positiveAssignments<-item$Blast_percentIdentical>assignmentThresholds
        item[6:length(item)][ !positiveAssignments]<- paste0("unclassified_",item[5+sum(positiveAssignments)]) # unclassified not capitalised to match standard Idtaxa output
        blastTopHits[row,]<-item
    }

     SeqDataTable<-merge(TableToMergeTo,blastTopHits, by= clustering, all=TRUE)



    return(list("SeqDataTable"=SeqDataTable, "FullBlastOutput"=FullBlastOutput))
}
