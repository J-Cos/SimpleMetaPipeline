[![DOI](https://zenodo.org/badge/419331038.svg)](https://zenodo.org/badge/latestdoi/419331038)

# SimplePipeline
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
    dir.create((file.path("BioinformaticPipeline_Env", "Data", "BlastDBs")))
    dir.create((file.path("BioinformaticPipeline_Env", "Data", "Classifiers")))
    dir.create((file.path("BioinformaticPipeline_Env", "Data", "Raw")))

This creates a directory structure as follows:
- FASTQs - fill this directory with directories titled Run1, Run2.. RunN, each filled with the unmerged multiplexed raw FASTQ files from a single sequencing run
- Data - fill this file with your taxonomic reference library trainingsets and metadata (raw fastas in "Raw", IDtaxa classifiers in "Classifiers", and BlastDbs in "BlastDBs")
- IntermediateOutputs - this will be populated by the pipeline as it runs, it will enable the pipeline to be run over multiple sessions as the output from each module is saved here.
- Results - this is where final results will be saved

Finally copy and paste this downloaded directory into the structure:
- BioinformaticPipeline - get this from github and then create your own control file from the template - no need to modify anything else


# 2) dependencies    
installed software - you must install these yourself and ensure they are in usr/bin or usr/local/bin.
- vsearch
- swarm v2
- cutadapt
- BLAST
    
R Packages - will be installed automatically in the pipeline

    library(dada2)
    library(seqinr)
    library(DECIPHER)
    library(tidyverse)
    library(lulu)
    library(ggplot2)
    library(gridExtra)
    library(ShortRead)

# 3) now you can start, create a new control script from the template and run the pipeline!

Remember you will need to find the optimal parameters first - either do this algorithmically or visually for instance in QIIME2.


# Appendix 1
## HPC Instructions
Linux and Mac OS open a Terminal on your desktop and type  ssh -XY username@login.hpc.imperial.ac.uk, substituting in your own college username, and entering in its password when prompted You'll need to be on the College network, or connected to the VPN.

Map drive with "IC" as domain
