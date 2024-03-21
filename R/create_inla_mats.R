#' @name create_inla_mats
#' @title Create the matrices required for SPDE approximation in CensSpBayes
#'
#' @description
#' Create the mass and stiffness matrices and the index matrices for both the observation locations and the prediction locations
#'
#' @usage create_inla_mats(S,S.pred=NULL,
#'                          offset=NULL,
#'                          cutoff=1e-12,
#'                          max.edge=NULL,...)
#'
#' @param S Matrix of locations where data is observed
#' @param S.pred (optional) Matrix of locations where data is to be predicted. If this is not supplied, the corresponding index matrix will be returned as a NULL object
#' @param offset The automatic extension distance. One or two values, for an inner and an optional outer extension. If negative, interpreted as a factor relative to the approximate data diameter (default=-0.10???). Same as in fmesher::fm_mesh_2d_inla()
#' @param cutoff The minimum allowed distance between points. Point at most as far apart as this are replaced by a single vertex prior to the mesh refinement step. Same as in fmesher::fm_mesh_2d_inla()
#' @param max.edge The largest allowed triangle edge length. One or two values. Same as in fmesher::fm_mesh_2d_inla()
#' @param ... Additional parameters to be passed to the function fmesher::fm_mesh_2d_inla()
#'
#' @import fmesher
#' @importFrom Matrix diag
#' @import INLA
#'
#' @return list of matrices needed for the SPDE approximation. The entries are: c.mat(mass matrix), g1.mat(stiffness matrix), g2.mat(a combination of mass and stiffness matrices), A(index matrix for observation locations), A.pred(index matrix for prediction locations)
#' @export



create_inla_mats <- function(S,S.pred=NULL,
                             offset=NULL,
                             cutoff=1e-12,
                             max.edge=NULL,...){

  Sall <- rbind(S,S.pred)
  sim.mesh <- fmesher::fm_mesh_2d_inla(offset=offset,
                              loc=S, loc.domain = Sall, cutoff=cutoff,
                              max.edge=max.edge,...)
  A <- INLA::inla.spde.make.A(mesh=sim.mesh,loc=S)
  A.pred <- NULL;
  if(!is.null(S.pred)){A.pred <- INLA::inla.spde.make.A(mesh=sim.mesh,loc=S.pred)}
  fem.mesh <- INLA::inla.mesh.fem(sim.mesh,order=2)
  c.mat <- fem.mesh$c0
  g1.mat <- fem.mesh$g1
  g2.mat <- fem.mesh$g2
  inla.mats <- list("c.mat"=c.mat,"g1.mat"=g1.mat,"g2.mat"=g2.mat,"A"=A,"A.pred"=A.pred)
  return(inla.mats)
}
