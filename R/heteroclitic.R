#' Generates single amino-acid variant of a peptide.
#'
#' This function mutates a peptide sequence and outputs single amino-acid variants of that peptide.
#' @param variant Default status is to return all possible variants
#' @param case Default is to output upper case unless "lower" is applied
#' @export
#' @examples
#' heteroclitic("LLLALLPPGA")
#'
#' heteroclitic("LLLALLPPGA", variant = "4R", case = "lower")


heteroclitic <- function(peptide, variant = "all", case = "upper"){

  amino.acid <- LETTERS[-match(c("B", "J", "O", "U", "X", "Z"), LETTERS)]

  all.variants <- NULL

  ##If variant = "all" default is applied. Return variant of each amino acid for each position
  if(variant == "all"){

    for(position in 1:nchar(peptide)){

      for(replacement in amino.acid){

        all.variants <- c(all.variants, paste(stringr::str_sub(peptide, start = 1, end = position - 1),
                                              replacement,
                                              stringr::str_sub(peptide, start = position + 1, end = nchar(peptide)),
                                              sep = ""))

      }
    }

  ##Otherwise, run only on selected position
  } else {

    replacement <- stringr::str_sub(variant, start = nchar(variant), end = nchar(variant))

    position <- as.numeric(stringr::str_sub(variant, start = 1, end = nchar(variant) - 1))

    all.variants <- paste(stringr::str_sub(peptide, start = 1, end = position - 1),
                          replacement,
                          stringr::str_sub(peptide, start = position + 1, end = nchar(peptide)),
                          sep = "")

  }


  ##Fix case
  if(case == "lower"){
    all.variants <- tolower(all.variants)
  } else {
    all.variants <- toupper(all.variants)
  }

  ##Remove peptides that are same as original
  all.variants <- all.variants[!all.variants %in% peptide]

  return(all.variants)



}
