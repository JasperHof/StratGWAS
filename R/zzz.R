.onLoad <- function(libname, pkgname) {
  # This is required for RcppParallel to work properly
  # It registers the TBB library path so R can find it at runtime
  RcppParallel::setThreadOptions(numThreads = "auto")
}

.onUnload <- function(libpath) {
  library.dynam.unload("StratGWAS", libpath)
}