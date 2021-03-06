HBAc <- function(inseq, silent = FALSE){

    if(silent == FALSE) message("species: sc, Chemical-map-based HBA", "\n")
    if(!is(inseq)[1] == "character"){
        if(is(inseq)[1] == "DNAString"){
            if(requireNamespace("Biostrings", quietly = TRUE)){
                inseq <- as.character(inseq)
                if(silent == FALSE){
                message(
                "The class of inseq was changed from DNAString to character")
                }
            }else{
            message("DNAString cannot be changed to a character string")
            out <- NA; names(out) <- "HBA"; return(out)
            }
        }else{
            message("The class of inseq must be DNAString or character")
            out <- NA; names(out) <- "HBA"; return(out)
        }
    }
    if(nchar(inseq) != 147){
        message("Length of inseq: ", nchar(inseq), "bp")
        message("The length of inseq must be 147 bp")
        out <- NA; names(out) <- "HBA"; return(out)
    }

        freqL1 <- nature11142_s2.147.freqL
        tranL1 <- nature11142_s2.147.tranL
        TtranL2 <- nature11142_s2.147.tranL2
        TtranL3 <- nature11142_s2.147.tranL3
        TtranL4 <- nature11142_s2.147.tranL4
        TfreqN4 <- nature11142_s2.147.freqN4SA
        TtranN4 <- nature11142_s2.147.tranN4_SMA
    
    outlist <- .Fortran("HBA_3", 
            inseq = inseq,
            logasc = numeric(length=1), freqL1 = freqL1, 
            tranL1 = tranL1, TtranL2 = TtranL2, TtranL3 = TtranL3, 
            TtranL4 = TtranL4, TfreqN4 = TfreqN4, 
            TtranN4 = TtranN4)[2]

    # N を含んでいるときはすべてを NA にする。
    if(outlist[[1]] == 0){
        # outlist[seq_len(13)] <- as.numeric(NA)    ## 13 を 1 にする！！！！ (20180129)
        outlist[1] <- as.numeric(NA)    
    }
    out <- unlist(outlist)
    names(out) <- "HBA"
    return(out)
}
