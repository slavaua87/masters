##' Multilayer optics
##'
##' R/c++ implementation of the dipole emission near a planar multilayer stack
##' 
##' @name integral-package
##' @docType package
##' @useDynLib integral
##' @import Rcpp
##' @title integral
##' @author baptiste Auguie \email{baptiste.auguie@@gmail.com}
##' @references
##' Etchegoin, P. Le Ru, E., Principles of Surface-Enhanced Raman Spectroscopy, Elsevier, Amsterdam (2009).
##' 
##' L. Novotny, E. Hecht, Principles of Nano-optics Cambridge University Press, 2006
##' 
##' H. Raether. Surface Plasmons on Smooth and Rough Surfaces and on Gratings. Springer, 1988.
##' @keywords packagelibrary
##' 
NULL

##' Rcpp module: integral
##' 
##' Exposes C++ function test
##' @name integral
##' @docType data
##' @export
##' @details
##' \itemize{
##'  \item{integrand_gb}{ integrand for gaussian beam excitation at a planar interface}
##' }
##' @examples
##' show( integral )
NULL
