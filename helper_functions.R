Fun_check_Primer_InterruptMotif <- function(primer_input_df = primer_df) {
  # check the primer directions and primer number
  # reads with primer information
  
  ## (1) detect concatenated reads
  ## a. if one read contains >1 primers, we see it as too many primers
  primer_input_df$ReadWithTooManyPrimer <- ifelse(
    (primer_input_df[,grep("^FR", colnames(primer_input_df))][,1] > 1 | primer_input_df[,grep("^FR", colnames(primer_input_df))][,2] > 1 | 
       primer_input_df[,grep("^RC", colnames(primer_input_df))][,1] > 1 | primer_input_df[,grep("^RC", colnames(primer_input_df))][,2] > 1), "Yes","No")
  
  ## b. if one read contains both primers from Forward reads and Reverse reads
  ## we define it as conflict primer
  primer_input_df$ReadWithConflictPrimer <- ifelse(
    (primer_input_df[,grep("^FR_", colnames(primer_input_df))][,1] > 0 | primer_input_df[,grep("^FR_", colnames(primer_input_df))][,2] > 0) & 
      (primer_input_df[,grep("^RC_", colnames(primer_input_df))][,1] > 0 | primer_input_df[,grep("^RC_", colnames(primer_input_df))][,2] > 0) ,"Yes","No")
  ## c. annotation for reads with too many primers or conflict primers to keep reads without concatenation
  primer_input_df$Filt.ConcatRead.BasedOn.Primer <- ifelse(
    primer_input_df$ReadWithTooManyPrimer == "No" & primer_input_df$ReadWithConflictPrimer == "No" & 
      rowSums(primer_input_df[,c(grep("^FR_", colnames(primer_input_df)), grep("^RC_", colnames(primer_input_df)))]) <= 4,
    "pass.primer.filter", "not.pass.primer.filter"
  )
  ########### end of detecting concatenated reads ###########
  
  ## (2) decide read direction based on whether it contains Forward primers or Reverse primers
  primer_input_df$ReadWithEitherPrimer <- ifelse(
    (primer_input_df[,grep("^FR_", colnames(primer_input_df))][,1] == 1 | primer_input_df[,grep("^FR_", colnames(primer_input_df))][,2] == 1) |  
      (primer_input_df[,grep("^RC_", colnames(primer_input_df))][,1] == 1 | primer_input_df[,grep("^RC_", colnames(primer_input_df))][,2] == 1),
    "Yes", "No")
  primer_input_df$ReadWithBothPrimer <- ifelse(
    (primer_input_df[,grep("^FR_", colnames(primer_input_df))][,1] == 1 & primer_input_df[,grep("^FR_", colnames(primer_input_df))][,2] == 1) |  
      (primer_input_df[,grep("^RC_", colnames(primer_input_df))][,1] == 1 & primer_input_df[,grep("^RC_", colnames(primer_input_df))][,2] == 1),
    "Yes", "No")
  
  ## (3) check Interrupt motif
  ## only keep InterruptMotif_FR or InterruptMotif_RC == 1
  primer_input_df$INTmotif_check <- ifelse(
    rowSums(primer_input_df[,grep("InterruptMotif_", colnames(primer_input_df))]) == 1, "Yes", "No"
  )
  
  ## (4) use ReadWithEitherPrimer and Interrupt Motif to decide read direction
  ## a. based on interrupt motif
  primer_input_df$FRorRC.based.INTmotif <- ifelse(
    primer_input_df$INTmotif_check == "Yes" & primer_input_df$Filt.ConcatRead.BasedOn.Primer == "pass.primer.filter" & 
      primer_input_df$InterruptMotif_FR == 1, "FR", "TBC"
  )
  primer_input_df[primer_input_df$FRorRC.based.INTmotif == "TBC",]$FRorRC.based.INTmotif <- ifelse(
    primer_input_df[primer_input_df$FRorRC.based.INTmotif == "TBC",]$INTmotif_check == "Yes" & 
      primer_input_df[primer_input_df$FRorRC.based.INTmotif == "TBC",]$Filt.ConcatRead.BasedOn.Primer == "pass.primer.filter" & 
      primer_input_df[primer_input_df$FRorRC.based.INTmotif == "TBC",]$InterruptMotif_RC == 1, "RC", "TBC"
  )
  ## b. based on Either primer and no concatemers
  primer_input_df$FRorRC.based.eitherPrim <- ifelse(
    primer_input_df$ReadWithEitherPrimer == "Yes" & primer_input_df$Filt.ConcatRead.BasedOn.Primer == "pass.primer.filter" & 
      (rowSums(primer_input_df[,grep("^FR_", colnames(primer_input_df))]) <= 2 & rowSums(primer_input_df[,grep("^RC_", colnames(primer_input_df))]) == 0) ,
    "FR", "TBC"
  )
  primer_input_df[primer_input_df$FRorRC.based.eitherPrim == "TBC",]$FRorRC.based.eitherPrim <- ifelse(
    primer_input_df[primer_input_df$FRorRC.based.eitherPrim == "TBC",]$ReadWithEitherPrimer == "Yes" & 
      primer_input_df[primer_input_df$FRorRC.based.eitherPrim == "TBC",]$Filt.ConcatRead.BasedOn.Primer == "pass.primer.filter" & 
      (rowSums(primer_input_df[primer_input_df$FRorRC.based.eitherPrim == "TBC",][,grep("^FR_",colnames(primer_input_df))]) == 0 & 
         rowSums(primer_input_df[primer_input_df$FRorRC.based.eitherPrim == "TBC",][,grep("^RC_",colnames(primer_input_df))]) <= 2) ,
    "RC", "TBC"
  )
  ## c. annotate reads based on ReadWithEitherPrimer and Interrupt Motif
  primer_input_df$FRorRC.based.eitherPrim.INTmotif <- ifelse(
    primer_input_df$FRorRC.based.INTmotif == primer_input_df$FRorRC.based.eitherPrim,
    primer_input_df$FRorRC.based.INTmotif, "TBC"
  )
  
  return(primer_input_df)
}

Fun_Primer_InterruptMotif_Loci <- function(primer_dfchecked, intrr_pattern_iupac_dna, reads_dna, list_primers) {
  # interrupt motif loci
  ## (1) keep clean data
  ## remove reads with problematic primers
  primer_dfchecked_clean <- primer_dfchecked[primer_dfchecked$Filt.ConcatRead.BasedOn.Primer == "pass.primer.filter",]
  ## keep reads with either primer and INT motif information
  primer_dfchecked_clean <- primer_dfchecked_clean[primer_dfchecked_clean$FRorRC.based.eitherPrim.INTmotif != "TBC",]
  
  read_fr_or_rc_vec <- c("FR", "RC") 
  int_pattern_loci_all_df <- NULL
  
  if (nrow(primer_dfchecked_clean) > 1){
    ## (2) separate FR and RC motif
    for (r in seq_along(read_fr_or_rc_vec)){
      ind_read_direction <- read_fr_or_rc_vec[r]
      colname_int_motif <- paste0("InterruptMotif_", ind_read_direction)
      inx_int_motif <- which(colnames(primer_dfchecked_clean) == colname_int_motif)
      
      if (ind_read_direction == "FR"){
        intrr_pattern_iupac <- intrr_pattern_iupac_dna
        primer_5prime <- list_primers$'F'
        primer_3prime <- list_primers$R_RC
      }else{
        intrr_pattern_iupac <- reverseComplement(intrr_pattern_iupac_dna)
        primer_5prime <- list_primers$R
        primer_3prime <- list_primers$F_RC
      }
      read_motif_loci <- primer_dfchecked_clean[primer_dfchecked_clean[,inx_int_motif] == 1,]$read_name
      
      if (length(read_motif_loci) > 0) {
        # check interrupt motif loci
        int_pattern_loci <- vmatchPattern(intrr_pattern_iupac, reads_dna[read_motif_loci], fixed = F)
        int_pattern_loci_df <- as.data.frame(unlist(int_pattern_loci))
        colnames(int_pattern_loci_df)[c(1:3)] <- c("Str.INTmotif", "End.INTmotif", "Width.INTmotif")
        int_pattern_loci_df$InterruptMotif.cat <- colname_int_motif
        
        # check primer loci
        primer_5prime_loci <- vmatchPattern(primer_5prime, reads_dna[read_motif_loci], max.mismatch = 1)
        primer_5prime_loci_df <- as.data.frame(unlist(primer_5prime_loci))
        colnames(primer_5prime_loci_df)[c(1:3)] <- c("Str.5prime", "End.5prime", "Width.5prime")
        
        primer_3prime_loci <- vmatchPattern(primer_3prime, reads_dna[read_motif_loci], max.mismatch = 1)
        primer_3prime_loci_df <- as.data.frame(unlist(primer_3prime_loci))
        colnames(primer_3prime_loci_df)[c(1:3)] <- c("Str.3prime", "End.3prime", "Width.3prime")
        
        primer_prime_loci_df <- merge(primer_5prime_loci_df, primer_3prime_loci_df,
                                      by = "names", all = T)
        
        primer_prime_loci_df$EitherPrimer.cat <- ind_read_direction
        
        int_pattern_loci_df <- merge(int_pattern_loci_df, primer_prime_loci_df,
                                    by = "names", all = T)
      } else {
        int_pattern_loci_df <- NULL
      }
      int_pattern_loci_all_df <- rbind(int_pattern_loci_all_df, int_pattern_loci_df)
    }
    
  } else {
    int_pattern_loci_all_df <- NULL
  }
  return(int_pattern_loci_all_df)
}


###
# user's input primer

Fun_Calculate_Trinucleotide <- function(ind_loci_res, dis.between_5prime.end_cag.tract){
  # #####
  # This function calculate number of trinucleotide between FR primer and interrupting motif CARCARCCRCCR
  # The location of FR primer and the distance using user's input primer
  ####
  ###
  # Input file
  # ind_loci_res: df object from Fun_Primer_InterruptMotif_Loci function output
  ###
  Read_FRorRC_vec <- c("FR","RC")
  Num_Trinucleotides_df <- NULL
  for (r in seq_along(Read_FRorRC_vec)){
    ind_read_dir <- Read_FRorRC_vec[r]
    sub_res_df <- ind_loci_res[ind_loci_res$EitherPrimer.cat==ind_read_dir,]
    if (ind_read_dir=="FR"){
      sub_res_clean_df <- sub_res_df[!is.na(sub_res_df$Str.5prime),]
      # loci for 5' primer should be upstream to INT motif
      sub_res_clean_df <- sub_res_clean_df[sub_res_clean_df$Str.5prime < sub_res_clean_df$Str.INTmotif,]
      sub_res_clean_df$Dist.prim2INTmotif <- sub_res_clean_df$Str.INTmotif - sub_res_clean_df$End.5prime -1
    }else{
      sub_res_clean_df <- sub_res_df[!is.na(sub_res_df$Str.3prime),]
      # for RC read, loci for 3' primer should be downstream to INT motif
      sub_res_clean_df <- sub_res_clean_df[sub_res_clean_df$Str.3prime > sub_res_clean_df$Str.INTmotif,]
      sub_res_clean_df$Dist.prim2INTmotif <- sub_res_clean_df$Str.3prime - sub_res_clean_df$End.INTmotif -1
    }
    # PCR products using FR primer: user's input 5' primer
    sub_res_clean_df$Num.Trinuc <- round((sub_res_clean_df$Dist.prim2INTmotif - dis.between_5prime.end_cag.tract)/3, digits=0)
    
    Num_Trinucleotides_df <- rbind(Num_Trinucleotides_df,sub_res_clean_df)
  }
  return(Num_Trinucleotides_df)
}



