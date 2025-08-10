#' AB plot of sender/receiver effects
#' @param muEff A vector of sender additive effects
#' @param ylabel Y-axis label for the plot
#' @return A ggplot object
#' @author Shahryar Minhas
#' @export
abPlot=function(muEff, ylabel='Sender Additive Effects'){
  muDf = data.frame(mu=muEff)
  muDf$id = rownames(muDf)
  muDf$id = factor(muDf$id, levels=muDf$id[order(muDf$mu)])
  muDf$ymax = with(muDf, ifelse(mu>=0,mu,0))
  muDf$ymin = with(muDf, ifelse(mu<0,mu,0))
  gg = ggplot(muDf, aes(x=id, y=mu)) +
    geom_point() +
    geom_linerange(aes(ymax=ymax,ymin=ymin)) +
    xlab('') + ylab(ylabel) +
    theme(
      axis.ticks=element_blank(),
      axis.text.x=element_text(angle=45)
    )
  return(gg)
}