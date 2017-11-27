#result analysis 

#LOAD NETWORKS
#Correlation, CMap
#fgsea cmap
#hlt cmap
#fgsea lincs
#hlt lincs
#write GMLs (for NX)
#distro plots 
# projecting ADR and PW - cmap
# projecting ADR and PW - lincs

#########################################################################################################
######LOAD NETWORKS#####
#########################################################################################################

#drug drug expression
load(file = "results/drug_drug_networks_correlation_and_footrule_analysis_cmap_lincs_separated.RData")
#drug pathway
load(file = "results/fgsea_nws.RData")
load(file = "results/fgsea_nws_lincs_combined.RData")
#drug ADR 
drug_adr = read.graph(file = "results/drug_adr.gml", format = "gml")
#drug HLT
drug_hlt = read.graph(file = "results/drug_hlt.gml", format = "gml")

##########################
#IMPORTANT################
##########################
#Graphs in this script are pre NX ... which removes disconnected nodes. 

#load correlation matrices

#libraries
library(ggplot2)

#helper function 
fread_nicely = function(x){
  y = fread(x, data.table = FALSE)
  rownames(y) = y[,1]
  y = y[,-1]
  return(y)
}
#
spearman_cmap = fread_nicely("results/cmap_spearman_correlation_matrix.txt")
spearman_lincs = fread_nicely("results/lincs_1_and_2_spearman_correlation_matrix.txt")
tanimoto_matrix = read.csv("results/tanimoto_matrix_complete.txt", sep = "\t")
#########################################################################################################
######Correlation, CMap#####
#########################################################################################################

#more modular than footrule, that's why we choose this one

cmap_correlation_nw_analysis$g
plot(cmap_correlation_nw_analysis$infomap, cmap_correlation_nw_analysis$g)
plot(walktrap.community(cmap_correlation_nw_analysis$g), cmap_correlation_nw_analysis$g)
edge_density(cmap_correlation_nw_analysis$g)
cmap_correlation_nw_analysis$average_path_length
cmap_correlation_nw_analysis$ClusteringCoefficient
cmap_correlation_nw_analysis$GraphDiameter

qplot(x = 0:max(degree(cmap_correlation_nw_analysis$g)), 
      y = degree.distribution(cmap_correlation_nw_analysis$g, 
                              cumulative = TRUE
                              ), 
      log = "xy"
      )
mean(x = degree(cmap_correlation_nw_analysis$g))
1/min(E(cmap_correlation_nw_analysis$g)$weight)

#########################################################################################################
######fgsea cmap#####
#########################################################################################################
nw_fgsea_cmap

table(V(nw_fgsea_cmap)$type)
?bipartite.mapping()
length(E(nw_fgsea_cmap))/(47*304)
degree.distribution(nw_fgsea_cmap, v = V(nw_fgsea_cmap)[type == "PATHWAY"])
degree.distribution(nw_fgsea_cmap, v = V(nw_fgsea_cmap)[type == "DRUG"])
qplot(x = 0:max(degree(nw_fgsea_cmap, v = V(nw_fgsea_cmap)[type == "PATHWAY"])), 
      y = degree.distribution(nw_fgsea_cmap, v = V(nw_fgsea_cmap)[type == "PATHWAY"], 
                              cumulative = TRUE
                              ), 
      log = "y"
      )
mean(degree(nw_fgsea_cmap, v = V(nw_fgsea_cmap)[type == "PATHWAY"]))

mean(degree(nw_fgsea_cmap, v = V(nw_fgsea_cmap)[type == "DRUG"]))
qplot(x = 0:max(degree(nw_fgsea_cmap, v = V(nw_fgsea_cmap)[type == "DRUG"])), 
      y = degree.distribution(nw_fgsea_cmap, v = V(nw_fgsea_cmap)[type == "DRUG"], 
                              cumulative = TRUE
                              ),
      log = "xy"
      )
##remove DRUGS with degree zero
V(nw_fgsea_cmap)[degree(nw_fgsea_cmap)==0]
fgsea_cmap_nonzero = induced_subgraph(graph = nw_fgsea_cmap, vids = V(nw_fgsea_cmap)[degree(nw_fgsea_cmap)>0])
qplot(x = 0:max(degree(fgsea_cmap_nonzero, v = V(fgsea_cmap_nonzero)[type == "DRUG"])), 
      y = degree.distribution(fgsea_cmap_nonzero, v = V(fgsea_cmap_nonzero)[type == "DRUG"]),
      #log = "xy"
)

#########################################################################################################
######hlt-drugs#####
#########################################################################################################
drug_hlt
table(V(drug_hlt)$type)
#########################################################################################################
######hlt cmap#####
#########################################################################################################
#get only those drugs in cmap
drugs_in_cmap = names(V(nw_fgsea_cmap)[type == "DRUG"])
cmap_drug_hlt = induced.subgraph(graph=drug_hlt,
                                         vids=unlist(neighborhood(graph=drug_hlt,
                                                                  order=1,
                                                                  nodes=V(drug_hlt)[name%in%drugs_in_cmap]
                                                                  )
                                                     )
                                 )

V(cmap_drug_hlt)$type = ifelse(V(cmap_drug_hlt)$type == 0, "DRUG", "HLT")
table(V(cmap_drug_hlt)$type)
cmap_drug_hlt
7780/(36*1160)

mean(degree(cmap_drug_hlt, v = V(cmap_drug_hlt)[type == "DRUG"]))
mean(degree(cmap_drug_hlt, v = V(cmap_drug_hlt)[type == "HLT"]))
qplot(x = 0:max(degree(cmap_drug_hlt, v = V(cmap_drug_hlt)[type == "HLT"])), 
      y = degree.distribution(cmap_drug_hlt, v = V(cmap_drug_hlt)[type == "HLT"],
                              cumulative = TRUE
                              ),
      log = "y"
)
qplot(x = 0:max(degree(cmap_drug_hlt, v = V(cmap_drug_hlt)[type == "DRUG"])), 
      y = degree.distribution(cmap_drug_hlt, v = V(cmap_drug_hlt)[type == "DRUG"],
                              cumulative = TRUE
                              ),
      log = "xy"
)

#########################################################################################################
######fgsea lincs#####
#########################################################################################################
nw_fgsea_lincs_combined
table(V(nw_fgsea_lincs_combined)$type)
?bipartite.mapping()
length(E(nw_fgsea_lincs_combined))/(45*302)
degree.distribution(nw_fgsea_lincs_combined, v = V(nw_fgsea_lincs_combined)[type == "PATHWAY"])
degree.distribution(nw_fgsea_lincs_combined, v = V(nw_fgsea_lincs_combined)[type == "DRUG"])
qplot(x = 0:max(degree(nw_fgsea_lincs_combined, v = V(nw_fgsea_lincs_combined)[type == "PATHWAY"])), 
      y = degree.distribution(nw_fgsea_lincs_combined, v = V(nw_fgsea_lincs_combined)[type == "PATHWAY"], 
                              cumulative = TRUE
      ), 
      log = "y"
)
mean(degree(nw_fgsea_lincs_combined, v = V(nw_fgsea_lincs_combined)[type == "PATHWAY"]))

mean(degree(nw_fgsea_lincs_combined, v = V(nw_fgsea_lincs_combined)[type == "DRUG"]))
qplot(x = 0:max(degree(nw_fgsea_lincs_combined, v = V(nw_fgsea_lincs_combined)[type == "DRUG"])), 
      y = degree.distribution(nw_fgsea_lincs_combined, v = V(nw_fgsea_lincs_combined)[type == "DRUG"], 
                              cumulative = TRUE
      ),
      log = "xy"
)
##remove DRUGS with degree zero
V(nw_fgsea_lincs_combined)[degree(nw_fgsea_lincs_combined)==0]
fgsea_lincs_nonzero = induced_subgraph(graph = nw_fgsea_lincs_combined, vids = V(nw_fgsea_lincs_combined)[degree(nw_fgsea_lincs_combined)>0])
qplot(x = 0:max(degree(fgsea_lincs_nonzero, v = V(fgsea_lincs_nonzero)[type == "DRUG"])), 
      y = degree.distribution(fgsea_lincs_nonzero, v = V(fgsea_lincs_nonzero)[type == "DRUG"]),
      #log = "xy"
)
#########################################################################################################
######hlt lincs#####
#########################################################################################################
#get only those drugs in lincs
drugs_in_lincs = names(V(nw_fgsea_lincs_combined)[type == "DRUG"])
lincs_drug_hlt = induced.subgraph(graph=drug_hlt,
                                 vids=unlist(neighborhood(graph=drug_hlt,
                                                          order=1,
                                                          nodes=V(drug_hlt)[name%in%drugs_in_lincs]
                                 )
                                 )
)

V(lincs_drug_hlt)$type = ifelse(V(lincs_drug_hlt)$type == 0, "DRUG", "HLT")
table(V(lincs_drug_hlt)$type)
lincs_drug_hlt
7807/(31*1165)

mean(degree(lincs_drug_hlt, v = V(lincs_drug_hlt)[type == "DRUG"]))
mean(degree(lincs_drug_hlt, v = V(lincs_drug_hlt)[type == "HLT"]))
qplot(x = 0:max(degree(lincs_drug_hlt, v = V(lincs_drug_hlt)[type == "HLT"])), 
      y = degree.distribution(lincs_drug_hlt, v = V(lincs_drug_hlt)[type == "HLT"],
                              cumulative = TRUE
      ),
      log = "y"
)
qplot(x = 0:max(degree(lincs_drug_hlt, v = V(lincs_drug_hlt)[type == "DRUG"])), 
      y = degree.distribution(lincs_drug_hlt, v = V(lincs_drug_hlt)[type == "DRUG"],
                              cumulative = TRUE
      ),
      log = "xy"
)
############################################################################################################
#write GMLs (for NX)
############################################################################################################
write.graph(graph = cmap_correlation_nw_analysis$g, file = "results/gml/nw_corr_cmap.gml", format = "gml")

write.graph(graph = nw_fgsea_cmap, file = "results/gml/nw_fgsea_cmap.gml", format = "gml")
write.graph(graph = cmap_drug_hlt, file = "results/gml/cmap_drug_hlt.gml", format = "gml")

write.graph(graph = nw_fgsea_lincs_combined, file = "results/gml/nw_fgsea_lincs_combined.gml", format = "gml")
write.graph(graph = lincs_drug_hlt, file = "results/gml/lincs_drug_hlt.gml", format = "gml")

write.graph(graph = drug_hlt, file = "results/gml/drug_hlt.gml", format = "gml")
#########################################################################################################
#distro plots 
#########################################################################################################

# drug drug
p = funcion_distroplotting(df = funcion_cumuldist(cmap_correlation_nw_analysis$g), 
                       scale = "", 
                       x_lab = "k", 
                       y_lab = "P(K)", 
                       fit = FALSE, 
                       titulo = "Drug - Drug (perturbation profile similarity)", 
                       tema = theme_gray(base_size = 15),
#                       postema = "plot.title = element_text(size=12, face='bold')"
                       )

p + theme(plot.title = element_text(size = 15))
#fgsea
p = funcion_distroplotting(df = funcion_cumuldist(nw_fgsea_cmap, 
                                              v = V(nw_fgsea_cmap)[type == "DRUG"]), 
                       scale = "loglog", 
                       x_lab = "k", 
                       y_lab = "P(K)", 
                       fit = FALSE, 
      #                 fit.type =  "lm",
                       titulo = "Drug Degree (Drug - Pathway perturbation network)", 
                      tema = theme_gray(base_size = 20)
      )
p + labs(subtitle = "CMap") + theme(plot.title = element_text(size = 15), plot.subtitle = element_text(size = 12, hjust = 0.5))


p = funcion_distroplotting(df = funcion_cumuldist(nw_fgsea_cmap, 
                                              v = V(nw_fgsea_cmap)[type == "PATHWAY"]), 
                       scale = "semilog", 
                       x_lab = "k", 
                       y_lab = "P(K)", 
                       fit = FALSE, 
                       #fit.type =  "lm",
                       titulo = "Pathway Degree (Drug - Pathway perturbation network)", 
                       tema = theme_gray(base_size = 20)
)

p + labs(subtitle = "CMap") + theme(plot.title = element_text(size = 15), plot.subtitle = element_text(size = 12, hjust = 0.5))

p = funcion_distroplotting(df = funcion_cumuldist(nw_fgsea_lincs_combined, 
                                                  v = V(nw_fgsea_lincs_combined)[type == "DRUG"]), 
                           scale = "loglog", 
                           x_lab = "k", 
                           y_lab = "P(K)", 
                           fit = FALSE, 
                           #                 fit.type =  "lm",
                           titulo = "Drug Degree (Drug - Pathway perturbation network)", 
                           tema = theme_gray(base_size = 20)
)

p + labs(subtitle = "LINCS") + theme(plot.title = element_text(size = 15), plot.subtitle = element_text(size = 12, hjust = 0.5))

p = funcion_distroplotting(df = funcion_cumuldist(nw_fgsea_lincs_combined, 
                                                  v = V(nw_fgsea_lincs_combined)[type == "PATHWAY"]), 
                           scale = "semilog", 
                           x_lab = "k", 
                           y_lab = "P(K)", 
                           fit = FALSE, 
                           #fit.type =  "lm",
                           titulo = "Pathway Degree (Drug - Pathway perturbation network)", 
                           tema = theme_gray(base_size = 20)
)

p + labs(subtitle = "LINCS") + theme(plot.title = element_text(size = 15), plot.subtitle = element_text(size = 12, hjust = 0.5))#drug adr

#Drug - ADR
p = funcion_distroplotting(df = funcion_cumuldist(drug_hlt, 
                                                  v = V(drug_hlt)[type == 0]),  #drug
                           scale = "semilog", 
                           x_lab = "k", 
                           y_lab = "P(K)", 
                           fit = FALSE, 
                           #fit.type =  "lm",
                           titulo = "Drug Degree (Drug - ADR network)", 
                           tema = theme_gray(base_size = 20)
)

p + labs(subtitle = "") + theme(plot.title = element_text(size = 15), plot.subtitle = element_text(size = 12, hjust = 0.5))

p = funcion_distroplotting(df = funcion_cumuldist(drug_hlt, 
                                                  v = V(drug_hlt)[type == 1]),  #drug
                           scale = "semilog", 
                           x_lab = "k", 
                           y_lab = "P(K)", 
                           fit = FALSE, 
                           #fit.type =  "lm",
                           titulo = "Adverse Drug Reaction Degree (Drug - ADR network)", 
                           tema = theme_gray(base_size = 20)
)

p + labs(subtitle = "") + theme(plot.title = element_text(size = 15), plot.subtitle = element_text(size = 12, hjust = 0.5))

#using the NX analyzed networks from script 11; keeping these for completeness

p = funcion_distroplotting(df = funcion_cumuldist(cmap_drug_hlt, v = V(cmap_drug_hlt)[type == "DRUG"]), 
                       scale = "semilog", 
                       x_lab = "k", 
                       y_lab = "P(K)", 
                       fit = FALSE, 
                       #fit.type =  "lm",
                       titulo = "Drug Degree (Drug - Adverse Reaction network)", 
                       tema = theme_gray(base_size = 20)
)

p + theme(plot.title = element_text(size = 15))

p = funcion_distroplotting(df = funcion_cumuldist(cmap_drug_hlt, v = V(cmap_drug_hlt)[type == "HLT"]), 
                       scale = "semilog", 
                       x_lab = "k", 
                       y_lab = "P(K)", 
                       fit = FALSE, 
                       #fit.type =  "lm",
                       titulo = "Adverse Reaction Degree (Drug - Adverse Reaction network)", 
                       tema = theme_gray(base_size = 20)
)

p + theme(plot.title = element_text(size = 15))

############################################################################################
# projecting ADR and PW - cmap
############################################################################################

incidence_cmap_pw = get.adjacency(graph = nw_fgsea_cmap, sparse = FALSE)
incidence_cmap_pw =incidence_cmap_pw[rownames(incidence_cmap_pw)%in%names(V(nw_fgsea_cmap)[type == "DRUG"]),
                                     colnames(incidence_cmap_pw)%in%names(V(nw_fgsea_cmap)[type == "PATHWAY"])]

incidence_cmap_adr = get.adjacency(graph = cmap_drug_hlt, sparse = FALSE, attr = "weight")
#incidence_cmap_adr =incidence_cmap_adr[rownames(incidence_cmap_adr)%in%names(V(cmap_drug_hlt)[type == "HLT"]),
#                                     colnames(incidence_cmap_adr)%in%names(V(cmap_drug_hlt)[type == "DRUG"])]
incidence_cmap_adr =incidence_cmap_adr[rownames(incidence_cmap_adr)%in%names(V(cmap_drug_hlt)[type == 1]),
                                       colnames(incidence_cmap_adr)%in%names(V(cmap_drug_hlt)[type == 0])]

secto = intersect(rownames(incidence_cmap_pw), colnames(incidence_cmap_adr))

incidence_cmap_pw = incidence_cmap_pw[rownames(incidence_cmap_pw)%in%secto,]
incidence_cmap_adr =incidence_cmap_adr[,colnames(incidence_cmap_adr)%in%secto]

incidence_cmap_pw = incidence_cmap_pw[order(rownames(incidence_cmap_pw)),
                                      order(colnames(incidence_cmap_pw))]
incidence_cmap_adr = incidence_cmap_adr[order(rownames(incidence_cmap_adr)),
                                      order(colnames(incidence_cmap_adr))]

##actual projection
projection_adr_pw_cmap = incidence_cmap_adr%*%incidence_cmap_pw
projection_adr_pw_cmap[1:5, 1:5]
which(x = projection_adr_pw_cmap == max(projection_adr_pw_cmap), arr.ind = TRUE)
colnames(projection_adr_pw_cmap)[c(116, 118, 137)]

which(x = projection_adr_pw_cmap == max(projection_adr_pw_cmap), arr.ind = TRUE)

a = names(neighbors(graph = nw_fgsea_cmap, v = "REACTOME_CELL_CYCLE"))
b = names(neighbors(graph = nw_fgsea_cmap, v = "REACTOME_CELL_CYCLE_MITOTIC"))
c = names(neighbors(graph = nw_fgsea_cmap, v = "REACTOME_DNA_REPLICATION"))

d = names(neighbors(graph = cmap_drug_hlt, v = "Eye and eyelid infections"))

d%in%c(a,b,c)
d[d%in%c(a,b,c)]

g_adr_pw_cmap = graph_from_incidence_matrix(projection_adr_pw_cmap, directed = FALSE, weighted = TRUE)
table(V(g_adr_pw_cmap)$type)

write.graph(g_adr_pw_cmap, file = "results/gml/projection_adr_pw_cmap.gml", format = "gml")
############################################################################################
# projecting ADR and PW - lincs
############################################################################################

incidence_lincs_pw = get.adjacency(graph = nw_fgsea_lincs_combined, sparse = FALSE)
incidence_lincs_pw =incidence_lincs_pw[rownames(incidence_lincs_pw)%in%names(V(nw_fgsea_lincs_combined)[type == "DRUG"]),
                                     colnames(incidence_lincs_pw)%in%names(V(nw_fgsea_lincs_combined)[type == "PATHWAY"])]

incidence_lincs_adr = get.adjacency(graph = lincs_drug_hlt, sparse = FALSE, attr = "weight")
incidence_lincs_adr =incidence_lincs_adr[rownames(incidence_lincs_adr)%in%names(V(lincs_drug_hlt)[type == "HLT"]),
                                     colnames(incidence_lincs_adr)%in%names(V(lincs_drug_hlt)[type == "DRUG"])]
#incidence_lincs_adr =incidence_lincs_adr[rownames(incidence_lincs_adr)%in%names(V(lincs_drug_hlt)[type == 1]),
#                                       colnames(incidence_lincs_adr)%in%names(V(lincs_drug_hlt)[type == 0])]

secto = intersect(rownames(incidence_lincs_pw), colnames(incidence_lincs_adr))

incidence_lincs_pw = incidence_lincs_pw[rownames(incidence_lincs_pw)%in%secto,]
incidence_lincs_adr =incidence_lincs_adr[,colnames(incidence_lincs_adr)%in%secto]

incidence_lincs_pw = incidence_lincs_pw[order(rownames(incidence_lincs_pw)),
                                      order(colnames(incidence_lincs_pw))]
incidence_lincs_adr = incidence_lincs_adr[order(rownames(incidence_lincs_adr)),
                                        order(colnames(incidence_lincs_adr))]

##actual projection
projection_adr_pw_lincs = incidence_lincs_adr%*%incidence_lincs_pw
projection_adr_pw_lincs[1:5, 1:5]
which(x = projection_adr_pw_lincs == max(projection_adr_pw_lincs), arr.ind = TRUE)
colnames(projection_adr_pw_lincs)[c(116, 118, 137)]

which(x = projection_adr_pw_lincs == max(projection_adr_pw_lincs), arr.ind = TRUE)

a = names(neighbors(graph = nw_fgsea_lincs_combined, v = "REACTOME_CELL_CYCLE"))
b = names(neighbors(graph = nw_fgsea_lincs_combined, v = "REACTOME_CELL_CYCLE_MITOTIC"))
c = names(neighbors(graph = nw_fgsea_lincs_combined, v = "REACTOME_DNA_REPLICATION"))

d = names(neighbors(graph = lincs_drug_hlt, v = "Eye and eyelid infections"))

d%in%c(a,b,c)
d[d%in%c(a,b,c)]

g_adr_pw_lincs = graph_from_incidence_matrix(projection_adr_pw_lincs, directed = FALSE, weighted = TRUE)
table(V(g_adr_pw_lincs)$type)

write.graph(g_adr_pw_lincs, file = "results/gml/projection_adr_pw_lincs.gml", format = "gml")
################################################################################################
#Correlation heatmaps
################################################################################################

library(gplots)

heatmap.2(x = as.matrix(tanimoto_matrix), 
          trace = "none", 
          col = colorpanel(n = 999, low = "white", "lightblue", "blue"),
          cexRow = 0.4,
          cexCol = 0.5,
          breaks = 1000
          )

heatmap.2(x = as.matrix(spearman_cmap), 
          trace = "none", 
          col = "redblue", 
          cexRow = 0.4,
          cexCol = 0.5,
          breaks = 1000
          )

heatmap.2(x = as.matrix(spearman_lincs), 
          trace = "none", 
          col = "redblue", 
          cexRow = 0.4,
          cexCol = 0.5,
          breaks = 1000
          )
