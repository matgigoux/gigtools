#' Mutates a nucleotide seuqnce
#'
#' This function mutates a nucleotide sequence based on a specific input, whether it be substitution, duplication, deletion and/or insertion. Splice variants are not tolerated.
#' @param protein.return Default status is to return the nucleotide sequence. Otherwise, return protein sequence.
#' @export
#' @examples
#' mutate_nucl(sequence, "c.497_499delinsAA")
#'
#' mutate_nucl(sequence, "c.487T>A")
#'
#' mutate_nucl(sequence, "c.625_629del")
#'
#' mutate_nucl(sequence, "c.814delG")
#'
#' mutate_nucl(sequence, "c.1042_1043insAAGGCCT")
#'
#' mutate_nucl(sequence, "c.320_324dup")
#'
#' mutate_nucl(sequence, "c.351dup")

mutate_nucl <- function(sequence, mutation, protein.return=F){

  ##There should only be 1 sequence to mutate. If multiple sequences are added, only take the first and send a warning
  if(length(sequence) > 1){
    sequence <- sequence[1]
    warning("Only first sequence was changed.")
  }

  ##Likewise, only 1 mutation can be applied at a time.
  if(length(mutation) > 1){
    mutation <- mutation[1]
    warning("Only first mutation was applied. Consider looping if multiple modifications are required.")
  }

  ##If format of mutation contains splicing nomenclature, "-" or "+", create error
  if(any(grepl("\\+", mutation), grepl("\\-", mutation))){
    stop("This function does not apply to splicing mutations")
  }

  if(!grepl("c.", mutation)){
    warning("Mutation format does not have \"c.\" associated with HGSV nomenclature")
  }

  old.protein <- paste(seqinr::translate(seqinr::s2c(sequence)), collapse = "")



  ##If substitution, do the following
  if(grepl(">", mutation)){

    old.bp <- unlist(stringr::str_split(mutation, ">"))[1]  ##Double-check whether the intendend base pair was actually switched
    old.bp <- rev(unlist(stringr::str_split(old.bp, "")))[1]
    sub.bp <- unlist(stringr::str_split(mutation, ">"))[2]  ##Replace with this nucleotide
    pos.bp <- as.numeric(paste(unlist(stringr::str_split(mutation, ""))[unlist(stringr::str_split(mutation, "")) %in% 0:9], collapse = "")) ##At this position

    ##Double checks whether the wildtype nucleotide is as expected
    if(old.bp != stringr::str_sub(sequence, start = pos.bp, end = pos.bp)){
      warning("Wildtype nucleotide is not as expected")
    }

    ##Create new nucleotide sequence
    new.sequence <- paste(stringr::str_sub(sequence, end = pos.bp - 1),
                          sub.bp,
                          stringr::str_sub(sequence, start = pos.bp + 1), sep = "")

    if(protein.return==T){
      ##Create new protein sequence
      new.protein <- seqinr::translate(seqinr::s2c(new.sequence))
      new.protein <- stringr::str_sub(paste(new.protein, collapse = ""), start = 1, end = which(new.protein=="*")[1] - 1)
      return(new.protein)
    } else {
      return(new.sequence)
    }




  ##check if only an insertion has occured and apply that insertion
  }else if(stringr::str_detect(mutation, "ins") & !stringr::str_detect(mutation, "del")){

    ins.position <- stringr::str_sub(mutation, start = 3, end = stringr::str_locate(mutation, "ins")[1] - 1) ##Isolate the insertion location information
    ins.position <- as.numeric(unlist(stringr::str_split(ins.position, "_"))[1])  ##Isolate the position of the event
    insertion <- stringr::str_sub(mutation, start = stringr::str_locate(mutation, "ins")[2] + 1) ##Isolate the insertion event

    ##Create new nucleotide sequence
    new.sequence <- paste(stringr::str_sub(sequence, end = ins.position),
                          insertion,
                          stringr::str_sub(sequence, start = ins.position + 1),
                          sep = "")

    if(protein.return==T){
      ##Create new protein sequence
      new.protein <- seqinr::translate(seqinr::s2c(new.sequence))
      new.protein <- stringr::str_sub(paste(new.protein, collapse = ""), start = 1, end = which(new.protein=="*")[1] - 1)
      return(new.protein)
    } else {
      return(new.sequence)
    }


  ##check if only a del has occured and apply that del
  } else if(stringr::str_detect(mutation, "del") & !stringr::str_detect(mutation, "ins")){

    del.position <- stringr::str_sub(mutation, start = 3, end = stringr::str_locate(mutation, "del")[1] - 1)  ##isolate the deletion location information
    del.position <- as.numeric(unlist(stringr::str_split(del.position, "_"))[1])  ##Isolate the beginning of the event. The end position is redundant given the length of the deletion

    del.stop <- stringr::str_sub(mutation, start = 3, end = stringr::str_locate(mutation, "del")[1] - 1)  #Isolate the length of the event
    del.stop <- as.numeric(unlist(stringr::str_split(del.stop, "_"))[2])

    if(is.na(del.stop)){del.stop <- del.position}

    ##Create new nucleotide sequence
    new.sequence <- paste(stringr::str_sub(sequence, start = 1, end = del.position - 1),
                          stringr::str_sub(sequence, start = del.stop + 1, end = nchar(sequence)),
                          sep = "")

    if(protein.return==T){
      ##Create new protein sequence
      new.protein <- seqinr::translate(seqinr::s2c(new.sequence))
      new.protein <- stringr::str_sub(paste(new.protein, collapse = ""), start = 1, end = which(new.protein=="*")[1] - 1)
      return(new.protein)
    } else {
      return(new.sequence)
    }


    ##check if a del and an ins has occured and apply that del
  } else if(stringr::str_detect(mutation, "del") & stringr::str_detect(mutation, "ins")){

    del.position <- stringr::str_sub(mutation, start = 3, end = stringr::str_locate(mutation, "del")[1] - 1) ##isolate the deletion location information
    del.position <- as.numeric(unlist(stringr::str_split(del.position, "_"))[1])  ##Isolate the beginning of the event. The end position is redundant given the length of the deletion

    del.stop <- stringr::str_sub(mutation, start = 3, end = stringr::str_locate(mutation, "del")[1] - 1)  #Isolate the length of the event
    del.stop <- as.numeric(unlist(stringr::str_split(del.stop, "_"))[2])

    insertion <- stringr::str_sub(mutation, start = stringr::str_locate(mutation, "ins")[2] + 1) ##Isolate the insertion event

    ##Create new nucleotide sequence
    new.sequence <- paste(stringr::str_sub(sequence, start = 1, end = del.position - 1),
                          insertion,
                          stringr::str_sub(sequence, start = del.stop + 1, end = nchar(sequence)),
                          sep = "")

    if(protein.return==T){
      ##Create new protein sequence
      new.protein <- seqinr::translate(seqinr::s2c(new.sequence))
      new.protein <- stringr::str_sub(paste(new.protein, collapse = ""), start = 1, end = which(new.protein=="*")[1] - 1)
      return(new.protein)
    } else {
      return(new.sequence)
    }


    ##Check if a duplication has occured
  } else if(stringr::str_detect(mutation, "dup")){

    dup.start <- as.numeric(stringr::str_sub(mutation, start = 3, end = stringr::str_locate(mutation, "_")[1] - 1)) ##isolate the deletion location information
    dup.end <- as.numeric(stringr::str_sub(mutation, start = stringr::str_locate(mutation, "_")[2] + 1, end = stringr::str_locate(mutation, "dup")[1] - 1))

    ##If only one nucleotide duplication
    if(is.na(dup.start)){
      dup.start <- as.numeric(stringr::str_sub(mutation, start = 3, end = stringr::str_locate(mutation, "dup")[1] - 1)) ##isolate the deletion location information
      dup.end <- dup.start
    }

    duplication <- stringr::str_sub(sequence, start = dup.start, end = dup.end) ##Isolate the insertion event

    ##Create new nucleotide sequence
    new.sequence <- paste(stringr::str_sub(sequence, end = dup.start - 1),
                          duplication,
                          stringr::str_sub(sequence, start = dup.start),
                          sep = "")

    if(protein.return==T){
      ##Create new protein sequence
      new.protein <- seqinr::translate(seqinr::s2c(new.sequence))
      new.protein <- stringr::str_sub(paste(new.protein, collapse = ""), start = 1, end = which(new.protein=="*")[1] - 1)
      return(new.protein)
    } else {
      return(new.sequence)
    }


  }


}
