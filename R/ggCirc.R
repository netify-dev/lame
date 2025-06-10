#' Plot circular network with ggplot2
#' 
#' @param Y (matrix) m by n relational matrix.
#' @param U (matrix) m by 2 matrix of row factors of Y.
#' @param V (matrix) n by 2 matrix of column factors of Y.
#' @param row.names (character vector) names of the row objects. 
#' @param col.names (character vector) names of the columns objects.
#' @param vscale Scaling factor for the V matrix
#' @param prange Range for point size
#' @param lcol Color for the lines
#' @param ltype Line type for the lines
#' @param lsize Line size for the lines
#' @param force Force parameter for ggrepel
#' @param maxIter Maximum iterations for ggrepel
#' @param showActLinks Logical indicating whether to show active links
#' @param geomLabel Logical indicating whether to use geom_label
#' @param geomText Logical indicating whether to use geom_text
#' @param geomPoint Logical indicating whether to use geom_point
#' 
#' @return A ggplot object representing the circular network
#' @author Shahryar Minhas
#' @export
ggCirc = function(
    Y, U=NULL, V=NULL, row.names=rownames(Y), col.names=colnames(Y),
    vscale=.6, prange=c(2,5), lcol='gray85', ltype='dotted', lsize=.5,
    force=2, maxIter = 3e3,
    showActLinks=FALSE, geomLabel=TRUE, geomText=FALSE, geomPoint=TRUE, ...
){
  
  #
  # libs
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggrepel))
  
  if (is.null(U)) {
    a <- rowMeans(Y, na.rm = TRUE)
    b <- colMeans(Y, na.rm = TRUE)
    Y0 <- Y
    Y0[is.na(Y)] <- (outer(a, b, "+"))[is.na(Y)]
    Y0 <- Y0 - mean(Y0)
    if (!all(Y == t(Y), na.rm = TRUE)) {
      sY <- svd(Y0)
      u <- sY$u[, 1:2]
      v <- sY$v[, 1:2]
      mu <- sqrt(apply(u^2, 1, sum))
      mv <- sqrt(apply(v^2, 1, sum))
      u <- diag(1/mu) %*% u
      v <- diag(1/mv) %*% v * vscale
    }
    if (all(Y == t(Y), na.rm = TRUE)) {
      eY <- eigen(Y0)
      bv <- which(abs(eY$val) >= sort(abs(eY$val), decreasing = TRUE)[2])[1:2]
      u <- eY$vec[, bv]
      mu <- sqrt(apply(u^2, 1, sum))
      u <- diag(1/mu) %*% u
      mv <- mu
      v <- u
    }
  }
  if (!is.null(U)) {
    if (is.null(V)) {
      V <- U
      vscale <- 1
    }
    mu <- sqrt(apply(U^2, 1, sum))
    mv <- sqrt(apply(V^2, 1, sum))
    u <- diag(1/mu) %*% U
    v <- diag(1/mv) %*% V * vscale
  }
  
  rsum <- apply(abs(Y), 1, sum, na.rm = TRUE)
  csum <- apply(abs(Y), 2, sum, na.rm = TRUE)
  links <- which(Y != 0, arr.ind = TRUE)
  
  # org df for gg
  uG = data.frame(u*1.2)
  uG$actor = rownames(Y)
  uG$tPch = 0 ; uG$tPch[rsum>0] = (mu[rsum>0])^3
  uG = uG[uG$tPch>0,]
  uG$tPch = uG$tPch
  
  # add v if supplied
  if(!is.null(V)){
    vG = data.frame(v*1.2)
    vG$actor = rownames(Y)
    vG$tPch = 0 ; vG$tPch[csum>0] = (mv[csum>0])^3
    vG = vG[vG$tPch>0,]
    vG$tPch = vG$tPch
    
    uG$eff = 'u' ; vG$eff = 'v'
    uG = rbind(uG, vG)
    ggCirc = ggplot(uG, aes(x=X1, y=X2,color=eff))
  }
  if(is.null(V)){
    ggCirc = ggplot(uG, aes(x=X1, y=X2))
  }
  
  # add segments
  if(showActLinks){
    for(i in 1:nrow(links)){
      ggCirc = ggCirc + geom_segment(
        x=u[links[i,1],1]*1.2, y=u[links[i,1],2]*1.2,
        xend=v[links[i,2],1]*1.2, yend=v[links[i,2],2]*1.2,
        color=lcol, linetype=ltype, size=lsize ) }
  }
  if(geomPoint){ ggCirc = ggCirc + geom_point(aes(...)) }
  if(geomLabel){ ggCirc = ggCirc + geom_label_repel(aes(label=actor, size=tPch, ...),
                                                    force=force, max.iter=maxIter) }
  if(geomText){ ggCirc = ggCirc + geom_text_repel(aes(label=actor, size=tPch, ...),
                                                  force=force, max.iter=maxIter) }
  ggCirc = ggCirc + scale_size(range=prange) +
    theme(
      legend.position='none',
      axis.ticks=element_blank(),
      axis.title=element_blank(),
      axis.text=element_blank(),
      panel.border=element_blank(),
      panel.grid=element_blank()
    )
  return(ggCirc)
}