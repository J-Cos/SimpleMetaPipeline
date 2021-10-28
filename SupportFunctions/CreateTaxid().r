CreateTaxid<-function() {       





        ranks <- readLines("../../PNG_OceanAcidification/Data/microgreen_algaebase.tax")
        taxa <- setNames(c("domain", "phylum", "order", "family", "genus"),
                c("d__", "p__", "o__", "f__", "g__"))
        ranks <- strsplit(ranks, ";", fix=T)
        count <- 1L
        groups <- "Root"
        index <- -1L
        level <- 0L
        rank <- "rootrank"
        pBar <- txtProgressBar(style=3)
        for (i in seq_along(ranks)) {
            for (j in seq_along(ranks[[i]])) {
                rank_level <- taxa[substring(ranks[[i]][j], 1, 3)]
                group <- substring(ranks[[i]][j], 4)
                w <- which(groups==group & rank==rank_level)
                if (length(w) > 0) {
                    parent <- match(substring(ranks[[i]][j - 1], 4),
                        groups)
                if (j==1 || any((parent - 1L)==index[w]))
                    next # already included
                }

                count <- count + 1L
                groups <- c(groups, group)
                if (j==1) {
                    index <- c(index, 0)
                } else {
                    parent <- match(substring(ranks[[i]][j - 1], 4),
                    groups)
                    index <- c(index, parent - 1L)
                }
                level <- c(level, j)
                rank <- c(rank, taxa[j])
            }
            setTxtProgressBar(pBar, i/length(ranks))
        }

        groups <- gsub("^[ ]+", "", groups)
        groups <- gsub("[ ]+$", "", groups)
        taxid <- paste(0:(length(index) - 1L), groups, index, level, rank, sep="*")
        head(taxid, n=10)

        writeLines(taxid, con="asdasfsedfsedfsef")
}