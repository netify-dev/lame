#' Plot parameter results from MCMC
#' 
#' @param mcmcData matrix of values from MCMC
#' @return Trace plot and density distribution of 
#' estimated parameters
#' @author Shahryar Minhas
#' @export 
paramPlot <- function(mcmcData){
  
  # libs
  suppressMessages(library(ggplot2))
  suppressMessages(library(plyr))
  suppressMessages(library(reshape2))
  s
  colnames(mcmcData) <- names(data.frame(mcmcData))
  mcmcMelt <- reshape2::melt(mcmcData)
  names(mcmcMelt) <- c('Var1', 'Var2', 'value')
  mcmcMelt$mu <- with(mcmcMelt, ave(value, Var2, FUN=mean))
  mcmcMelt$median <- with(mcmcMelt, ave(value, Var2, FUN=median))
  qts <- c(lo95=0.025,lo90=0.05,hi90=0.95,hi95=0.975)
  for(i in 1:length(qts)){
    mcmcMelt$tmp <- with(mcmcMelt, ave(value, Var2, FUN=function(x){quantile(x,probs=qts[i])}))
    names(mcmcMelt)[ncol(mcmcMelt)] <- names(qts)[i]
  }
  
  ggTrace <- ggplot2::ggplot(mcmcMelt, aes(x=Var1, y=value)) +
    geom_hline(aes(yintercept=mu), color='red') + 
    geom_ribbon(aes(ymin=lo90,ymax=hi90), alpha=.5, fill='grey40') +
    geom_ribbon(aes(ymin=lo95,ymax=hi95), alpha=.3, fill='grey40') +  
    geom_line(lwd=.6) + 
    xlab('Post-Burn Iteration') + 
    ylab('Parameter Value') + 
    facet_wrap(~Var2, ncol=1, scales='free_y') +  
    theme_bw() + 
    theme(
      panel.border=element_blank(),
      axis.ticks=element_blank()
    )
  
  ggMu <- plyr::ddply(mcmcMelt, .(Var2), summarise, mu=mean(value))
  ggDens <- plyr::ddply(mcmcMelt,.(Var2),.fun = function(x){
    tmp <- density(x$value) ; x1 <- tmp$x; y1 <- tmp$y
    q90 <- x1 >= quantile(x$value,0.05) & x1 <= quantile(x$value,0.95)
    q95 <- x1 >= quantile(x$value,0.025) & x1 <= quantile(x$value,0.975)
    data.frame(x=x1,y=y1,q90=q90,q95=q95) })
  
  ggDist <- ggplot2::ggplot(ggDens, aes(x=x)) +
    facet_wrap(~ Var2, scales='free', ncol=1) +
    geom_vline(data=ggMu, aes(xintercept=mu), linetype='solid', size=1, color='red') +
    geom_line(aes(y=y), color='grey40') + ylab('') + 
    xlab('Parameter Value') +
    geom_ribbon(data=subset(ggDens, q90), aes(ymax=y),ymin=0, alpha=.7, fill='grey40') +   
    geom_ribbon(data=subset(ggDens, q95), aes(ymax=y),ymin=0, alpha=.5, fill='grey40') + 
    theme_bw() + 
    theme(
      panel.border=element_blank(),
      axis.ticks=element_blank(),
      axis.text.y=element_blank()
    )
  
  return( gridExtra::grid.arrange(ggTrace, ggDist, ncol=2) )
  
}