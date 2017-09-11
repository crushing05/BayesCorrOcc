#' CreateSpp
#'
#' Create folder to store results
#' @export

CreateSpp <- function(alpha){
  spp_test <- dir.exists(paste0("inst/output/", alpha))

  if(!spp_test){
    dir.create(paste0("inst/output/", alpha))


    file.copy(from = "inst/output/.gitignore", to = paste0("inst/output/", alpha, "/.gitignore"))
  }
}

