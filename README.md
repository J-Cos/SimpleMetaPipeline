# BioinformaticPipeline
 General Bioinformatic Pipeline, and parameter sets recorded for specific datasets.

# 1) directory structure
 Before you start you need the right directory structure.

 First run this code in R
            setwd() #set this to where you want the pipeline environment to be located
            path <-"../../BioinformaticPipeline_Env"
            dir.create("BioinformaticPipeline_Env")
            dir.create((file.path("BioinformaticPipeline_Env", "FASTQs")))
            dir.create((file.path("BioinformaticPipeline_Env", "Data")))

            dir.create((file.path("BioinformaticPipeline_Env", "IntermediateOutputs")))
            dir.create((file.path("BioinformaticPipeline_Env", "Results")))

This creates a directory structure as follows
                # FASTQs - fill this file with your unmerged multiplexed raw FASTQ files
                # Data - fill this file with your taxonomic reference library trainingsets and metadata
                # IntermediateOutputs - this will be populated by the pipeline as it runs, it will enable the pipeline to be run over multiple sessions as the output from each module is saved here.
                # Results - this is where final results will be saved

Finally copy this downloaded directory into the structure:
                # BioinformaticPipeline - get this from github and then create your own control file from the template - no need to modify anything else


# 2) dependencies    
    installed software - you must install these yourself and ensure they are in usr/bin or usr/local/bin.
        vsearch
        swarm v2
    
    R Packages - will be installed automatically
        library(dada2)
        library(seqinr)
        library(DECIPHER)
        library(tidyverse)
        library(lulu)
        library(ggplot2)
        library(gridExtra)

# 3) now you can start, create a new control script from the template and run the pipeline!
    Remember you will need to find the optimal parameters first - either do this algorithmically or visually for instance in QIIME2.