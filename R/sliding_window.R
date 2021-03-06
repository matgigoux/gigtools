#' Outputs all peptides consituting longer protein sequence
#'
#' This function generates a list of peptides of desired lengths that compose a longer protein in table format, along with position and length infomration.
#' @param pep.len Peptide length default to 9
#' @export
#' @examples
#' prot.sequence = "ETFNTPAMYVAIQAVLSLYASGRTTGIVMDSGDG"
#' sliding_window(prot.sequence, 8:12)

sliding_window <- function(prot.sequence, pep.len=9){

  if(length(prot.sequence)>1){ ##Only takes first object from prot.sequence and generates warning
    prot.sequence <- prot.sequence[1]
    warning("Sliding window was only applied to first sequence")
  }

  if(xor(typeof(pep.len)!="double",typeof(pep.len)=="integer")){ ##Creates an error if pep.len is not numeric
    stop("pep.len is not numeric")
  }

  if(any(!(pep.len%%1==0))){ ##Creates an error if pep.len is not integer
    stop("pep.len is not integer")
  }

  all.peptides <- NULL

  for(pep.size in pep.len){ ##For every size of required peptide length, start at position 1

    start.pos <- 1

    while(start.pos <= nchar(prot.sequence)-pep.size+1){

      ##Generate peptide
      peptide <- stringr::str_sub(prot.sequence,
                                  start = start.pos,
                                  end = start.pos + pep.size - 1)

      all.peptides <- rbind(all.peptides,
                            c(start.pos, start.pos + pep.size - 1, pep.size, peptide))

      start.pos <- start.pos + 1

    }

  }

  all.peptides <- data.frame(all.peptides)

  colnames(all.peptides) <- c("Start.Pos", "End.Pos", "Pep.Len", "Peptide")

  all.peptides$Start.Pos <- as.numeric(as.character(all.peptides$Start.Pos))
  all.peptides$End.Pos <- as.numeric(as.character(all.peptides$End.Pos))
  all.peptides$Pep.Len <- as.numeric(as.character(all.peptides$Pep.Len))
  all.peptides$Peptide <- as.character(all.peptides$Peptide)

  return(all.peptides)

}
