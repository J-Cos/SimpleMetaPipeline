######################################
# Testing legacy COI demultiplexxing options
#####################################


# 1) paths and dependencies
    library(tidyverse)
    library(dada2)
    library(ShortRead)
    library(seqinr)

    path<-'/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/LegacyCOI_test'
    #contents of this directory should be:
        #Barcodes - this should contain a subdirectory for each run which in turn contains a list of the barcodes related to each run for each adapter
        #RunsSplitByAdapter - this should contain a subdirectory for each run which in turn contain R1 and R2 fastqs for each adapter
    
    #primers
        mlCOIintF_char<-("GGWACWGGWTGAACWGTWTAYCCYCC")
        jgHCO2198_char<-("TANACYTCNGGRTGNCCRAARAAYCA")
        mlCOIintF<-DNAString(mlCOIintF_char)
        jgHCO2198<-DNAString(jgHCO2198_char)

    #run choice (1 is a test with only 999 reads, 2 is a full adapter)
        run<-2
        subsetsize<-10000 # must be <1000 for run 1

# functions
    MakePrimerBarcodeHeatmap<-function(BarcodePrimers, names, subsetsize=subsetsize) {

        heatmap<-matrix(ncol=length(BarcodePrimers)+1, nrow=length(BarcodePrimers)+1, 0)
        colnames(heatmap)<-rownames(heatmap)<-names

        for (r in 1:length(subsetR1)) {
            x<-y<-c()
            for (bp in 1:length(BarcodePrimers)) {
                if (vcountPattern(BarcodePrimers[[bp]], subsetR1[r], fixed = FALSE) >0 ) {
                    x<-c(x, bp+1)
                }
                if (vcountPattern(BarcodePrimers[[bp]], subsetR2[r], fixed = FALSE) >0 ) {
                    y<-c(y,bp+1)
                }
            }
            if (is.null(y)) {y<-1}
            if (is.null(x)) {x<-1}
            heatmap[x,y]<-heatmap[x,y]+1
        }
        #heatmap[1,1]<-NA
        
        heatmap_proportions<-heatmap/subsetsize*100

        return(heatmap_proportions)
    }
    MakePrimerBarcodeHeatmapPlot<-function(heatmap_proportions, title){
        heatmap_proportions %>% 
                as.data.frame() %>%
                rownames_to_column("BarcodePrimerR1") %>%
                pivot_longer(-c(BarcodePrimerR1), names_to = "BarcodePrimerR2", values_to = "PercentOfReads") %>%
                ggplot(aes(x=BarcodePrimerR2, y=BarcodePrimerR1, fill=PercentOfReads)) + 
                    geom_raster() +
                    ggtitle(title)+
                    geom_text(aes(label=signif(PercentOfReads,2)),size=3)+
                    theme(axis.text.x=element_text(angle=90))
    }
    #function to # Create all orientations of the input sequence
    allOrients <- function(primer) {
        require(Biostrings)
        dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
        orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
            RevComp = reverseComplement(dna))
        return(sapply(orients, toString))  # Convert back to character vector
    }

    #count hits
    Hits_subset <- function(primer, subset) {
        # Counts number of reads in which the primer is found
        nhits <- vcountPattern(primer, subset, fixed = FALSE)
        return(sum(nhits > 0))
    }


# 2) code
    RunsSplitByAdapter<-list.files(file.path(path, "RunsSplitByAdapter"))
    AllBarcodes<-list.files(file.path(path, "Barcodes"))


        #get data for a single run
        Run<-RunsSplitByAdapter[run]
        RunFastqs<-list.files(file.path(path, "RunsSplitByAdapter", Run), pattern="fastq", full.names = TRUE)


        #get subset of reads from fastq
            subsetR1<-sread(readFastq(RunFastqs[[1]]))[1:subsetsize] 
            subsetR2<-sread(readFastq(RunFastqs[[2]]))[1:subsetsize] 


        Barcodes<-AllBarcodes[1]
        RunBarcodes<-sort(list.files(file.path(path, "Barcodes", Barcodes), pattern=".txt", full.names = TRUE))
        AdapterBarcodesFile<-RunBarcodes[run]


        #make df of barcodes from this adapter
        Barcodes_df<-read.table(AdapterBarcodesFile)
        colnames(Barcodes_df)<-c("Sample", "Barcode")
        sampleNum<-dim(Barcodes_df)[1]

        #various unused ways to create versions of barcodes and primers
            #create all orients
            #barcode_orients_l<-lapply(Barcodes_df$Barcode, allOrients)
            
            
           # FO<-lapply(paste0(Barcodes_df$Barcode, mlCOIintF), DNAString)
            #RO<-lapply(paste0(Barcodes_df$Barcode, jgHCO2198), DNAString)
            #FO2<-lapply(paste0(mlCOIintF, Barcodes_df$Barcode), DNAString)
            #RO2<-lapply(paste0(jgHCO2198, Barcodes_df$Barcode), DNAString)
            #FOAll<-lapply(FO1, allOrients) %>% unlist %>% lapply(DNAString)
            #ROAll<-lapply(RO1, allOrients) %>% unlist %>% lapply(DNAString)
            #FO2All<-lapply(FO2, allOrients) %>% unlist %>% lapply(DNAString)
            #RO2All<-lapply(RO2, allOrients) %>% unlist %>% lapply(DNAString)

           # names(FO1)<-paste0(Barcodes_df$Sample, "_FO")
           # names(RO1)<-paste0(Barcodes_df$Sample, "_RO")
           # names(RO2)<-paste0(Barcodes_df$Sample, "_RO")
           # names(FO2)<-paste0(Barcodes_df$Sample, "_FO")

           # BarcodePrimers<-c(FOAll, ROAll, FO2All, RO2All)
           # BarcodePrimers<-c(FO, RO)

           # Sample_for_mlC<- lapply( paste0( Barcodes_df$Barcode , mlCOIintF ) , DNAString)
           # revcomp_mlC<-lapply(Sample_for_mlC, allOrients) %>% unlist %>% lapply(DNAString)
           # Sample_for_JGH<- lapply( paste0( Barcodes_df$Barcode , jgHCO2198 ) , DNAString)
           # revcomp_JGH<-lapply(Sample_for_JGH, allOrients) %>% unlist %>% lapply(DNAString)
           # BarcodePrimers<-c(Sample_for_mlC, revcomp_mlC, Sample_for_JGH, revcomp_JGH)


            #Sample_rev_JGH<- lapply( paste0( lapply(Barcodes_df$Barcode, Biostrings::reverse) , jgHCO2198 ) , DNAString)
            #Sample_rev_mlC <- lapply( paste0( lapply(Barcodes_df$Barcode, Biostrings::reverse) , mlCOIintF ) , DNAString)
            #BarcodePrimers<-c(Sample_for_mlC, Sample_rev_JGH, Sample_for_JGH, Sample_rev_mlC)

        #for barcode primers combo from tripps script
            Sample_for_mlC<- lapply( paste0( Barcodes_df$Barcode , mlCOIintF ) , DNAString)
            Sample_for_JGH<- lapply( paste0( Barcodes_df$Barcode , jgHCO2198 ) , DNAString)
            revcomp_mlC<-lapply(Sample_for_mlC, reverseComplement) %>% unlist %>% lapply(DNAString)
            revcomp_JGH<-lapply(Sample_for_JGH, reverseComplement) %>% unlist %>% lapply(DNAString)
            BarcodePrimers<-c(Sample_for_mlC, Sample_for_JGH, revcomp_mlC, revcomp_JGH)

            matrixnames<-c("None",  paste0("mlCOIint_", rep(1:sampleNum)), paste0("jgHCO2198_", rep(1:sampleNum)), paste0("revcomp_mlCOIint_", rep(1:sampleNum)), paste0("revcomp_jgHCO2198_", rep(1:sampleNum)))
                                
            Heatmap<-MakePrimerBarcodeHeatmap(BarcodePrimers=BarcodePrimers, names=matrixnames, subsetsize=subsetsize)
            p1<-MakePrimerBarcodeHeatmapPlot(Heatmap, title= "Frequencey of all expected primerbarcodess on both reads")

        #for primers only - all orients
            PrimerOrients<- lapply(c( mlCOIintF_char, jgHCO2198_char) , allOrients) %>% unlist() %>%lapply(DNAString)
            matrixnames<-c("None", c("For_mlCOIint", "Comp_mlCOIint", "Rev_mlCOIint", "RevComp_mlCOIint", "For_jgHCO2198", "Comp_jgHCO2198", "Rev_jgHCO2198", "RevComp_jgHCO2198"))
            Heatmap<-MakePrimerBarcodeHeatmap(BarcodePrimers=PrimerOrients, names=matrixnames, subsetsize=subsetsize)
            p2<-MakePrimerBarcodeHeatmapPlot(Heatmap, title= "Frequencey of all orietnts of primers on both reads")

        #for barcodes only
            barcodeOrients<-lapply(Barcodes_df$Barcode, allOrients) %>% unlist %>% lapply(DNAString)
            matrixnames<-c("None", paste0(c("For_", "Comp_", "Rev_", "RevComp_"), unlist(lapply(1:sampleNum, rep,4) )))
            Heatmap<-MakePrimerBarcodeHeatmap(BarcodePrimers=barcodeOrients, names=matrixnames, subsetsize=subsetsize)
            p3<-MakePrimerBarcodeHeatmapPlot(Heatmap, title= "Frequencey of all orients of barcodes on both reads")

        #for all possible combinations but only looking at one sample for speed
            Sample1barcodeOrients<-lapply(Barcodes_df$Barcode[1], allOrients) %>% unlist 
            PrimerOrients<- lapply(c( mlCOIintF_char, jgHCO2198_char) , allOrients) %>% unlist

            Allcombos1sample<-c()
            for (i in 1:length(Sample1barcodeOrients)) {
                Allcombos1sample<-c(Allcombos1sample, 
                                    lapply(paste0(PrimerOrients, Sample1barcodeOrients[i]), DNAString),
                                    lapply(paste0(Sample1barcodeOrients[i], PrimerOrients), DNAString)
                                    )
            }

            matrixnames<-c("None",  paste0(c("For_mlCOIint", "Comp_mlCOIint", "Rev_mlCOIint", "RevComp_mlCOIint", "For_jgHCO2198", "Comp_jgHCO2198", "Rev_jgHCO2198", "RevComp_jgHCO2198"),"_forSample"),
                                    paste0("_forSample", c("For_mlCOIint", "Comp_mlCOIint", "Rev_mlCOIint", "RevComp_mlCOIint", "For_jgHCO2198", "Comp_jgHCO2198", "Rev_jgHCO2198", "RevComp_jgHCO2198")),
                                    paste0(c("For_mlCOIint", "Comp_mlCOIint", "Rev_mlCOIint", "RevComp_mlCOIint", "For_jgHCO2198", "Comp_jgHCO2198", "Rev_jgHCO2198", "RevComp_jgHCO2198"),"_compSample"),
                                    paste0("_compSample", c("For_mlCOIint", "Comp_mlCOIint", "Rev_mlCOIint", "RevComp_mlCOIint", "For_jgHCO2198", "Comp_jgHCO2198", "Rev_jgHCO2198", "RevComp_jgHCO2198")),
                                    paste0("_revSample", c("For_mlCOIint", "Comp_mlCOIint", "Rev_mlCOIint", "RevComp_mlCOIint", "For_jgHCO2198", "Comp_jgHCO2198", "Rev_jgHCO2198", "RevComp_jgHCO2198")),
                                    paste0("_revSample", c("For_mlCOIint", "Comp_mlCOIint", "Rev_mlCOIint", "RevComp_mlCOIint", "For_jgHCO2198", "Comp_jgHCO2198", "Rev_jgHCO2198", "RevComp_jgHCO2198")),
                                    paste0(c("For_mlCOIint", "Comp_mlCOIint", "Rev_mlCOIint", "RevComp_mlCOIint", "For_jgHCO2198", "Comp_jgHCO2198", "Rev_jgHCO2198", "RevComp_jgHCO2198"),"_revcompSample"),
                                    paste0("_revcompSample", c("For_mlCOIint", "Comp_mlCOIint", "Rev_mlCOIint", "RevComp_mlCOIint", "For_jgHCO2198", "Comp_jgHCO2198", "Rev_jgHCO2198", "RevComp_jgHCO2198"))
            )
            
            Heatmap<-MakePrimerBarcodeHeatmap(BarcodePrimers=Allcombos1sample, names=matrixnames, subsetsize=subsetsize)
            Heatmap[1,1]<-NA #as only one sample lots will be none-none so we remove this to make pattern visible

            pALL<-MakePrimerBarcodeHeatmapPlot(Heatmap, title= "Frequencey of all possible primer barcode combos from one sample on both reads")



    pdf(file = file.path(path,paste0("PrimerBarcodePresenceHeatmaps.pdf")), width=20, height =8 ) # The height of the plot in inches        
            egg::ggarrange(p1, p2, p3, widths=c(1,1,1))

            pALL
    dev.off()

#count hits if needed
            #count barcode primers individually by read
            subsetsize<-10
            subset<-sread(readFastq(RunFastqs[[2]]))[1:subsetsize] 
            m<-matrix(ncol=subsetsize, nrow=length(BarcodePrimers))
            for (i in 1:length(BarcodePrimers)) {

            m[i,]<- subset%>%
                    vcountPattern(BarcodePrimers[[i]], ., fixed = FALSE)        
            }


            colSums(m)[colSums(m)>1] %>% length



            # count barcode primers overall

            barcodeCountR1<-sapply(lapply(barcodeOrients, DNAString), Hits_subset, subsetR1)
            barcodeCountR2<-sapply(lapply(barcodeOrients, DNAString), Hits_subset, subsetR2)
            sum(barcodeCountR1)/length(subsetR1)
            sum(barcodeCountR2)/length(subsetR2)

