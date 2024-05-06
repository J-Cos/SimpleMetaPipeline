[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7740558.svg)](https://doi.org/10.5281/zenodo.7740558)

# SimpleMetaPipeline
 A simple bioinformatic pipeline for metabarcoding and meta-analyses, designed for Linux.

# 1) Create directory structure

 First you need to build the expected directory structure so that the pipeline knows where to find the inputs and can provide outputs in the expected locations. This exact directory structure is required for the pipeline to work.

 To do this open a terminal window and run the code below. By default this will create your pipeline directory in home, if you wish to create it somewhere else just change to that directory before running the below code:
 
    mkdir -p BioinformaticPipeline_Env/{FASTQs,Data/{BlastDBs,Classifiers,Raw},IntermediateOutputs,Results}

This creates a directory structure as follows:
- FASTQs - fill this directory with a subdirectory containing all runs you wish to analyse together. Within this subdirectory create further subdirectories for each run titled Run1, Run2... RunN, each filled with the unmerged multiplexed raw FASTQ files from a single sequencing run
- Data - fill this directory with your taxonomic reference library training sets and metadata (raw training library fastas in "Raw", IDtaxa classifiers in "Classifiers", and BlastDbs in "BlastDBs"). You can download commonly used IDtaxa classifiers here: http://www2.decipher.codes/Downloads.html
- IntermediateOutputs - this will be populated by the pipeline as it runs, it will enable the pipeline to be run over multiple sessions as the output from each step in the pipeline is saved here as it is produced. Note that this means you will need to empty this directory if you wish to start a new pipeline run on data that you have previously run through the pipeline.
- Results - this is where final results will be saved. Once you have your results it is recommened that you move them out of this directory into a project specific directory. This will keep the pipeline directory structure clean for future pipeline runs.

# 2) Download the pipeline

Once the directory is set up we need to download the pipeline code, from this repository, into our directory structure. If you have git installed you can run the code below:

    cd BioinformaticPipeline_Env
    git clone https://github.com/J-Cos/SimpleMetaPipeline BioinformaticPipeline

Alternatively you can do this manually by pressing the button marked "< > Code" above and then choosing "Download ZIP" from the dropdown menu. Once donwloaded unzip the file, copy it into the BioinformaticPipeline_Env directory and rename it to:
- BioinformaticPipeline

This repository contains 4 subdirectories as follows:
- Pipeline: the main directory containing the pipeline itself (Pipeline.R) and a directory (Functions) containing all the functions called by the pipeline. You do not need to modify these scripts.
- ControlScripts: this directory contains a single example control script. We will first use this example to check that your installation is working, and you can then modify it to create a control script for your own data. This is the only directory where you should modify scripts.
- 2 other directories not required to run the pipeline. These contain supporting functionality that is under develoment.

# 3) Install dependencies    
SimpleMetaPipeline aims to make life easy by stiching together the latest bioinformatic tools for metabarcoding data. Of course, these tools have to be installed for SimpleMetaPipeline to work. To do this we use conda. 

First you need to install conda. (We recommend miniconda - a streamlined version of conda - if you will only be using conda in the context of SimpleMetaPipeline). Installation instructions are available here:
https://conda.io/projects/conda/en/latest/user-guide/install/index.html

> ⚠️ N.B. though the pipeline is designed for linux if you are installing on macOS note that you will need to install the osx-64 conda installer rather than osx-arm64 as the bioconda channel does not compile for osx-arm64 at this time. See here: https://github.com/conda/conda/issues/11216

Once conda is installed and initialised you need to open a new terminal window, navigate to the pipeline directory you created earlier and run one line of code. This single line tells conda to install everything else you need! For example, if your pipeline directory was created in your home directory:
    
    cd BioinformaticPipeline_Env 
    conda env create --prefix ./BioinformaticPipeline/env --file BioinformaticPipeline/environment.yml 

This creates a conda environment containing all the software dependencies for the pipeline to work. 

# 4) Run the example data through the pipeline to confirm your installation and directory structure are correctly configured.

## 4.1) Running the Example data

Download the example data [here](https://drive.google.com/drive/folders/1dz3LibPnOfTBEOnQucZ_gfvja4y2aS7W?usp=drive_link). Make sure to download the directory called "Example", not just the subdirectories it contains, as the structure of this directory is important.

Unzip the downloaded .zip file and move the contained "Example" directory (along with all contained subdirectories) into the FASTQs directory within the pipeline directory structure you created earlier.

Download a compatible IDtaxa classifier here (this example data is 16s so we need to use a 16s classifier): http://www2.decipher.codes/Classification/TrainingSets/GTDB_r207-mod_April2022.RData

Move the downloaded classifier into the Data/Classifiers directory within the pipeline directory structure you created earlier.

Now let's check everything looks right. Your directory structure should look like this:

```bash
└──BioinformaticPipeline_Env
    ├── BioinformaticPipeline
    │   ├── ControlScripts
    │   ├── env (this was created when you ran "conda env create ..." earlier)
    │   ├── Pipeline
    │   ├── SupportFunctions
    │   ├── SupportScripts
    │   ├── environement.yml
    │   ├── LICENSE
    │   └── README.md
    ├── Data
    │   ├── BlastDBs
    │   ├── Classifiers
    │   │   └── GTDB_r207-mod_April2022.RData
    │   └── Raw
    ├── FASTQs
    │   └── Example
    │       ├── RandomCode_Run1_RandomCode
    │       └── RandomCode_Run2_RandomCode
    ├── IntermediateOutputs
    └── Results
```

If it does then you are ready to run the pipeline. Open a new terminal window and run the following commands:

    cd BioinformaticPipeline_Env
    conda activate ./BioinformaticPipeline/env 
    Rscript "BioinformaticPipeline/ControlScripts/ControlScriptExample.r"
  
Note that (once our conda environment is activated) the final command is equivalent to opening ControlScriptExample.R in a code editor (e.g. Visual Studio Code or R Studio) and running each line of ControlScriptExample.r sequentially.

If you've done everything correctly the pipeline will run and the IntermediateOutputs and Results directories you created earlier will be populated with the outputs of the pipeline. All these files are prefixed with the name of your dataset, in this case "Example".

> [!IMPORTANT]
> Every time you wish to use the pipeline you will first need to run "conda activate ./BioinformaticPipeline/env". This activates the environment in which you have installed all the required software for the pipeline to work. If you try and run the pipeline without first activating your conda environment then it won't work as the software will not be available. When the environment is correctly activated Conda will prepend the path name to the activated environment onto your system command. See the conda documentation for more information on environment management.

## 4.2 Inspecting the Example IntermediateOutputs and Results
IntermediateOutputs files include the output after each algorithm is applied (e.g. Example_DadaOutput.RDS and Example_IdtaxaOutput.RDS)

Results will contain 9 files as follows
- Main results files
    - Example_SeqDataTable.RDS (the main result of the pipeline, use this as the input to SimpleMetaPackage::SeqDataTable2Phyloseq(), see below)
    - Example_SeqDataTable.csv (the same information in csv format, for easy visual inspection)
- Post-run diagnostic files
    - Example_PrimerCounts.csv (if cutadapt is used this contains the number of each primer found, otherwise it will be empty)
    - Example_PrimerCountsAfterCutadapt.csv (if cutadapt is used this will show the numbers of each primer remaining, which should be 0)
    - Example_DadaSeqLengthDistribution.csv (the distribution of ASV lengths as produced by dada2)
    - Example_DadaTable.csv (a tracking table showing the number of reads passing each step in the dada2 process)
    - Example_DadaPlots.pdf (key diagnostic plots produced by dada2)
    - Example_ClusteringTable.csv (a simple table showing the number of ASVs, OTUs or similar clusters produced)
    - Example_TaxaAssignment.pdf (taxonomic assignment pie chart produced by IDTAXA)

# 5) Now you can run your own data, create a new control script from the template and run the pipeline from this control script as you did before!

All parameters in the control script are adjustable to suit your own data. The example control script provided is heavily commented, and includes links to the underlying tools, to guide you in changing these parameters as needed. 

Remember you will need to find the optimal filtering and trimming parameters for your own sequences first. These optimal filtering and trimming parameters will differ across different MiSeq runs. See here for an example of how to find these optimal parameters for each run based on visual inspection of read quality profiles: https://benjjneb.github.io/dada2/tutorial.html.

You will also need to find an appropriate Idtaxa clssifier, BLAST database, or both for your marker gene. A range of pre-trained Idtaxa classifiers are available here: http://www2.decipher.codes/Downloads.html. If required, instructions for training your own Idtaxa classifier / BLAST database are available in the published guidance for each tool, we do not duplicate this here.

Once you have results take a look at SimpleMetaPackage (https://github.com/J-Cos/SimpleMetaPackage), which provides tools to easily convert pipeline outputs to phyloseq objects, and also enables multi-algorithm agreement tests.

If you want to rerun the pipeline on the same dataset for any reason you will need to delete the intermediate outputs (the contents of the "IntermediateOutputs" directory) from previous runs of that dataset. This is because the pipeline checks which intermediate outputs are available and starts running from the furthest point the pipeline previously reached. This feature enables an interrupted run to be start from the last completed step, rather than needing to start from the beginning again.
