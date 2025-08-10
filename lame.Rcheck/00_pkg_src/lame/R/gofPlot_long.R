#' Plot longitudinal goodness of fit results
#' 
#' @param GOF 3D array of GOF statistics from lame model (fit$GOF)
#' @param type string: "actual" or "deviation" for visualization type. "Actual" to show 
#'   observed vs simulated values, or "deviation" to show differences between 
#'   simulated and observed values
#' @param symmetric logical: is the network symmetric?
#' @return ggplot object showing GOF across time periods
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
gofPlot_long <- function(GOF, type = "actual", symmetric = FALSE) {
  
  # Load required libraries
  # (packages are imported via NAMESPACE)
  
  # Validate type parameter
  if (!type %in% c("actual", "deviation")) {
    stop("type parameter must be 'actual' or 'deviation'")
  }
  
  # restructure gof
  gofList = lapply(dimnames(GOF)[[2]], function(t){
    
    # process
    goft = t(GOF[,t,])
    act = goft[1,] # first row is actual
    goft = goft[-1,]
    goft = data.frame(goft)
    goft$time = t
    goft = melt(goft, id='time')
    
    # add in actual values
    stats = unique(goft$variable)
    goft$actual = NA
    for(stat in stats){
      goft$actual[goft$variable==stat] = act[stat] 
    }
    
    return(goft)
  })
  
  # merge all together
  gofDF = do.call('rbind', gofList)
  
  # Remove dyadic dependency if symmetric
  if(symmetric) {
    gofDF = gofDF[gofDF$variable != 'dyad.dep', ]
  }
  
  # Create plot based on type parameter
  if(type == "actual") {
    # viz actual values
    ggplot(gofDF, aes(x=factor(time), y=value)) +
      geom_boxplot() +
      geom_point(
        aes(y=actual), color='grey50', shape='diamond', size=4) +
      labs(
        x = '',
        y = 'Value',
        caption = 'Boxplots show posterior predictive distributions\nDiamond points show observed values'
      ) +
      facet_wrap(~variable, scales='free_y') +
      theme_bw() +
      theme(
        plot.caption = element_text(hjust = 0.5))  # Center the caption
    
  } else if(type == "deviation") {
    # viz deviations from actual
    gofDF$diff = gofDF$value - gofDF$actual
    
    ggplot(gofDF, aes(x=factor(time), y=diff)) +
      geom_boxplot() +
      geom_hline(
        yintercept=0, color='grey50', linetype='dotted') +
      labs(
        y = 'Simulated - Observed',
        x = '',
        caption = 'Boxplots show differences between simulated and observed values'
      ) +
      facet_wrap(~variable, scales='free_y') +
      theme_bw()+
      theme(
        plot.caption = element_text(hjust = 0.5))  # Center the caption
    
  }
}