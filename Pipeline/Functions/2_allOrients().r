allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    # The Biostrings works w/ DNAString objects rather than character vectors
    dna <- DNAString(primer) 
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), RevComp = reverseComplement(dna))
    # Convert back to character vector
    return(sapply(orients, toString))  
}