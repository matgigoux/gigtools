#' Outputs all peptides containing a specific point mutation
#'
#' This function generates a list of peptides of desired lengths that incorporate a specific point mutation, along with position and length information, as well as WT equivalent peptide.
#' @param pep_length Peptide length default to 9
#' @export
#' @examples
#' protein = "ETFNTPAMYVAIQAVLSLYASGRTTGIVMDSGDG"
#' mut = "V15L"
#' peptide.sizes = c(9,10)
#' snv_peptides(prot_sequence = protein, mutation = "V15L", pep_length = peptide.sizes)

snv_peptides <- function(prot_sequence = NULL, mutation = NULL, pep_length = 9){

  if(!is.character(prot_sequence)){stop("Protein sequence improperly formated -- Needs to be character")}
  if(is.null(prot_sequence)){stop("No protein sequence input")}
  if(is.null(mutation)){stop("No mutation input")}
  if(!stringr::str_sub(mutation, start=1, end=1) %in% toupper(letters)){stop("Mutation improperly formated -- WT residue")}
  if(!stringr::str_sub(mutation, start=nchar(mutation), end=nchar(mutation)) %in% toupper(letters)){stop("Mutation improperly formated -- MUTANT residue")}
  if(!is.numeric(as.numeric(stringr::str_sub(mutation, start=2, end=nchar(mutation)-1)))){stop("Mutation improperly formated -- position not numeric")}
  if(!is.numeric(pep_length)){stop("peptide length is not numeric")}

  wt.residue <- stringr::str_sub(mutation, start=1, end=1)
  mutant.residue <- stringr::str_sub(mutation, start=nchar(mutation), end=nchar(mutation))
  position.residue <- as.numeric(stringr::str_sub(mutation, start=2, end=nchar(mutation)-1))

  if(!stringr::str_sub(prot_sequence, start=position.residue, end=position.residue) == wt.residue){stop("WT residue/position in mutation format does not match WT residue in protein sequence")}


  mut_prot_sequence <- paste0(stringr::str_sub(prot_sequence, start = 1, end = (position.residue - 1)),
                              mutant.residue,
                              stringr::str_sub(prot_sequence, start = (position.residue + 1), end = nchar(prot_sequence)))


  all.peptides <- NULL

  for(len in pep_length){

    ##Resets alls values
    temp.protein <- NULL
    adjust <- NULL
    start.pos <- NULL
    end.pos <- NULL
    wt_prot_window <- NULL
    mut_prot_window <- NULL

    if((position.residue - (len-1)) < 1){   ##If the mutation position is near the beginning of the protein, start at position 1
      start.pos <- 1
      adjust <- (position.residue - (len-1)) * -1 + 1
    } else {
      start.pos <- (position.residue - (len-1))
      adjust <- 0
    }


    if((position.residue + (len-1)) > nchar(prot_sequence)){  ##Likewise, if the mutation position is near the end, make the end position the end of the protein
      end.pos <- nchar(prot_sequence)
    } else {
      end.pos <- (position.residue + (len-1))
    }

    wt_prot_window <- stringr::str_sub(prot_sequence, start = start.pos, end = end.pos)  ##Get a smaller peptide that incorporates the mutant residue in all possible sliding window peptide assuming a specific peptide length
    mut_prot_window <- stringr::str_sub(mut_prot_sequence, start = start.pos, end = end.pos)  ##Same as above except this sequence contains the mutation

    temp.protein <- gigtools::sliding_window(mut_prot_window, pep.len = len)  ##Run Sliding window function from gigtools

    temp.protein$Start.Pos <- temp.protein$Start.Pos + position.residue - len  + adjust ##Fix the start and end position
    temp.protein$End.Pos <- temp.protein$End.Pos + position.residue - len + adjust

    colnames(temp.protein)[4] <- "MUT.Peptide"

    temp.protein$WT.Peptide <- gigtools::sliding_window(wt_prot_window, pep.len = len)$Peptide

    all.peptides <- rbind(all.peptides,
                          temp.protein)

  }

  return(all.peptides)

}
