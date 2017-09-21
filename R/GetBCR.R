#' GetBCR
#'
#' Get BCR #'s to subset count data
#' @param alpha alpha code for species of interest
#' @return Either a vector containing the BCR #s of interest or NULL if all BCRs where the species was detected should be used
#' @export

GetBCR	<- function(alpha) {
   bcr_list <- read.csv("inst/spp_list.csv", stringsAsFactors = FALSE)

   spp_list <- dplyr::filter(bcr_list, spp == alpha)
   spp_bcr <- spp_list$bcr

   if(spp_bcr == "NULL"){
     bcr <- NULL
   }else{
     a <- strsplit(spp_bcr, ",")
     a2 <- unlist(a)
     a3 <- as.vector(as.numeric(a2))
     bcr <- a3
   }

   bcr
}


