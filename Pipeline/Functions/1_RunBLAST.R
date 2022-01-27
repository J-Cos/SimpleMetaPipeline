
# pipeline function to blast query sequences against Biocode (should be generalised long-term)
#test<-RunBLAST( path="~/Dropbox/BioinformaticPipeline_Env",
 #               dbname="BIOCODE_ER",
  #              queryfasta="test.fasta")

    #requires tidyverse


RunBLAST<-function(path, dbname, queryfasta) {

    #blastn seqs
        system(command= paste0( "~/miniconda3/bin/blastn ",
                                "-db ", file.path(path, "Data", "BlastDBs", dbname, dbname),
                                " -query ", file.path(path, "IntermediateOutputs", queryfasta),
                                " -out ", file.path(path, "IntermediateOutputs", paste0("BlastOutput_", dbname, ".out")),
                                " -outfmt '6 qseqid sseqid pident evalue qcovs' ", 
                                " -evalue 1e-10" #expect value must be higher than this
                                ) 
            )

    #read output into r
        blast_output<-read.table( file.path(path,"IntermediateOutputs", paste0("BlastOutput_", dbname, ".out")), header=FALSE, sep="\t")

    #adjust names blastoutput
        names(blast_output)<-c("OTU", "ID", "percent_identical_matches", "evalue", "query_coverage")

    #map taxonomy df onto blast output
        for (i in 1: dim(blast_output)[1]) {
            blast_output$ID[i]<-unlist ( str_split(blast_output$ID[i], "\\|") )[1]
        }

        x<-left_join(blast_output, Tax_df, by = "ID")


        return(x)
    }
