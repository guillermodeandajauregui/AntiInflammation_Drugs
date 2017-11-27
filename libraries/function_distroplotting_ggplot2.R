#distro plotting ggplots 2 

funcion_cumuldist = function(g, v = V(g)){
  x = 0:max(degree(g, v = v))
  y = degree.distribution(g, v = v, cumulative = TRUE)
  z = data.frame(degree = x, accum.freq = y)
  return(z)
}

cumulative.dist <- function(x, by.breaks = 0.01, maxy = 1){
  breaks = seq(0,maxy, by = by.breaks)
  x.cut = cut(x, breaks, right = FALSE)
  x.freq = table(x.cut)
  x.cumsum = cumsum(x.freq)
  return(x.cumsum)
}

funcion_distroplotting = function(df, 
                                  scale = "",
                                  x_lab = "k",
                                  y_lab = "p(k)",
                                  fit = TRUE,
                                  fit.type = "auto",
                                  fit.color = "red",
                                  titulo = "",
                                  tema = theme_minimal()#,
                                  #postema = ""
                                  )
  {
  p = ggplot(data = df, 
             mapping = aes(degree, accum.freq))
  p = p + geom_point() 
  
  if(scale == "loglog"){
  p = p + scale_y_log10() + scale_x_log10()
  }
  if(scale == "semilog"){
    p = p + scale_y_log10()
  }
  if(fit == TRUE){
  p = p + geom_smooth(se = FALSE, color = "red", method = fit.type)
  }
  p = p + xlab(x_lab) + ylab(y_lab)
  p = p + tema
  p = p + ggtitle(titulo) + theme(plot.title=element_text(hjust=0.5))
  #if(postema!=""){
  #  p = p + theme(postema)
  #}
  return(p)
}
# 
# a = funcion_distroplotting(df, scale = "loglog", fit = FALSE, titulo = "titulo")
# b = funcion_distroplotting(df, scale = "semilog", fit = FALSE)
# 
# a
