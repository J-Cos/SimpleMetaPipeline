[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7740558.svg)](https://doi.org/10.5281/zenodo.7740558)

# SimpleMetaPipeline
 A simple bioinformatic pipeline for metabarcoding and meta-analyses, designed for Linux.

# 1) Download pipeline and create directory structure

 First download the pipeline code from this repository. This contains 4 subdirectories as follows:
- Pipeline: the main directory containing the pipeline itself (Pipeline.R) and a directory (Functions) containing all the functions called by the pipeline. You do not need to modify these scripts.
- ControlScripts: this directory contains a single example control script. First use this example to check your installation is working, and then modify it to create a control script for your own data. This is the only directory where you should modify scripts.
- 2 other directories not required to run the pipeline. These contain supporting functionality that is under develoment.

 Then you need to build the expected directory structure, so the pipeline knows where to find the inputs and can provide outputs in the expected locations. This exact directory structure is required for the pipeline to work.

 First run this code in R

    setwd() #set this to where you want the pipeline environment to be located
    dir.create("BioinformaticsPipeline_Env")
    dir.create((file.path("BioinformaticsPipeline_Env", "FASTQs")))
    dir.create((file.path("BioinformaticsPipeline_Env", "Data")))
    dir.create((file.path("BioinformaticsPipeline_Env", "IntermediateOutputs")))
    dir.create((file.path("BioinformaticsPipeline_Env", "Results")))
    dir.create((file.path("BioinformaticsPipeline_Env", "Data", "BlastDBs")))
    dir.create((file.path("BioinformaticsPipeline_Env", "Data", "Classifiers")))
    dir.create((file.path("BioinformaticsPipeline_Env", "Data", "Raw")))

This creates a directory structure as follows:
- FASTQs - fill this directory with a subdirectory containing all runs you wish to analyse together. Within this subdirectory create further subdirectories for each run titled Run1, Run2.. RunN, each filled with the unmerged multiplexed raw FASTQ files from a single sequencing run
- Data - fill this directory with your taxonomic reference library training sets and metadata (raw training library fastas in "Raw", IDtaxa classifiers in "Classifiers", and BlastDbs in "BlastDBs"). You can download commonly used IDtaxa classifiers here: http://www2.decipher.codes/Downloads.html
- IntermediateOutputs - this will be populated by the pipeline as it runs, it will enable the pipeline to be run over multiple sessions as the output from each step in the pipeline is saved here as it is produced. Note that this means you will need to empty this directory if you wish to start a new pipeline run on data you have previously run through the pipeline.
- Results - this is where final results will be saved. Once you have your results it is recommened you move them out of this directory into a project specific directory. This will keep the pipeline directory structure clean for future pipeline runs.

Finally copy and paste this downloaded directory into the BioinformaticsPipeline_Env directory without renaming it:
- BioinformaticsPipeline

# 2) Install dependencies    
SimpleMetaPipeline aims to make life easy by stiching together the latest bioinformatic tools for metabarcoding data. Of course these tools have to be installed for SimpleMetaPipeline to work. To do this you need to open your terminal navigate to the directory you just created and run one line of code. For example if your directory was in your home directory:
    
    cd BioinformaticsPipeline_Env/BioinformaticsPipeline
    conda env create --prefix ./env --file environment.yml

# 3) Run the example data through the pipeline to confirm your installation is correctly configured.

Download the example data here: https://drive.google.com/drive/folders/1FUALCE8PkWZabMpG4VuMjQqL_6hMmmpZ?usp=sharing

Move the Example directory contained in the downloaded .zip file (along with all contained subdirectories) into the FASTQs directory within the pipeline directory structure you created earlier.

Download a compatible IDtaxa classifier here (this example data is 16s so we need to use a 16s classifier): http://www2.decipher.codes/Classification/TrainingSets/GTDB_r207-mod_April2022.RData

Move the downloaded classifier into the Data/Classifiers directory within the pipeline directory structure you created earlier.

Open terminal and run the following commands:

    conda activate ../BioinformaticsPipeline/env
    Rscript "ControlScripts/ControlScriptExample.r"

Note that (once that environment is activated) this is equivalent to opening ControlScriptExample.R in a code editor (e.g. Visual Studio Code or R Studio) and running each line of the script sequentially.

If installed correctly the pipeline will run and the IntermediateOutputs and Results directories you created earlier will be populated with the outputs of the pipeline.

# 4) Now you can run your own data, create a new control script from the template and run the pipeline from this control script as you did before!

Once your have results take a look at SimpleMetaPackage (https://github.com/J-Cos/SimpleMetaPackage) which provides tools to easily convert pipeline outputs to phyloseq objects, and also enables multi-algorithm agreement tests.

Remember you will need to find the optimal filtering and trimming parameters for your own sequences first. Remember these can differ across multiple runs. See here for an example of how to find these parameters based on vidual inspection of read quality profiles: https://benjjneb.github.io/dada2/tutorial.html

If you want to rerun the pipeline on the same dataset for any reason you will need to delete the intermediate outputs from previous runs on that dataset.
