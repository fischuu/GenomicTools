.onAttach <- function(libname, pkgname){
  packageStartupMessage("This is GenomicTools version ", utils::packageVersion("GenomicTools"))
}