#'Load example data
#
#'@export
#'
load_example_data <- function(){
    example_data <- load(file.path("inst", "extdata", "SLiPat_data.RData"))
    return(example_data)
 }
