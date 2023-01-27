######################################
# Secondary demlutiplexing and unmixing orientations for primary demultiplexed mixed orientation Illumina MiSeq COI data
# Cross-pacific ARMS data (Emma)
#####################################

    path<-'/rds/general/user/jcw120/home/BioinformaticPipeline_Env/RawFastqs/RansomeAll_COI'
        #contents of this directory should be:
            #Barcodes - this should contain a subdirectory for each run which in turn contains a list of the barcodes related to each run for each adapter
            #RunsSplitByAdapter - this should contain a subdirectory for each run which in turn contain R1 and R2 fastqs for each adapter
    
    Ncores<-30
        #number of cores to use

    primerF="GGWACWGGWTGAACWGTWTAYCCYCC"
        #mlCOIintF

    primerR="TANACYTCNGGRTGNCCRAARAAYCA"
        #jgHCO2198

    RunSubset<-1
        #either FALSE (i.e. all runs) or a vector containing the numbers of the runs desired, e.g. c(1,6,11) or 4:7
    
#run the script
    source(file.path(path, "../../BioinformaticPipeline/SupportScripts/MixedOrientationDemultiplex.R"))