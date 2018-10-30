.onLoad <- function(libname, pkgname){
  s <- scala(pkgname)
  scalaLazy(function(s) s + '
            import uk.wellredd.ExtraTrees.ExtraTrees
            import uk.wellredd.Table.Table
            import uk.wellredd.IndividualClaimsModel.{mle, maximumLogLikelihood, meanFutureNumObservations, varianceFutureNumObservations}
            ')
  assign("s",s,envir=parent.env(environment()))
}

.onUnload <- function(libpath){
  close(s)
}

