cumulative.dist <- function(x, by.breaks = 0.01){
  breaks = seq(0,1, by = by.breaks)
  x.cut = cut(x, breaks, right = FALSE)
  x.freq = table(x.cut)
  x.cumsum = cumsum(x.freq)
  return(x.cumsum)
}

E(tanimoto_g)[which(E(graph = tanimoto_g)$weight==1)]


a = 1 - cumulative.dist(E(tanimoto_g)$weight)/(max(cumulative.dist(E(tanimoto_g)$weight)))
b = seq(0.01,1, by = 0.01)
plot(x = b, y = a)

cut_tanimoto_g = delete.edges(graph = tanimoto_g, 
                              E(tanimoto_g)[weight <= quantile(E(tanimoto_g)$weight, 0.95)]
)

plot(walktrap.community(cut_tanimoto_g), cut_tanimoto_g)
plot(infomap.community(cut_tanimoto_g), cut_tanimoto_g)

funcion_distroplotting(df = funcion_cumuldist(cut_tanimoto_g), 
                       scale = "")

g = cut_tanimoto_g
V(g)
