################################################################
#Anti-inflammatory Drug Network                                #
################################################################

################################################################
#05 Drug- Adverse Drug Reaction Network
################################################################

##################################
#For each anti-inflammatory drug
#identify Adverse Drug Reactions
#(ADRs)
#With the following conditions:
#1) FDR =< 0.05 for PRR vs other antiinflammatories
#2) FDR =< 0.05 for PRR All_antiinflammatories vs Non-antiinflammatories
#express as network
##################################

#read normalized FAERS data 
#read anti-inflammatory drug list
#calculate Proportional Reporting Ratios:
#calculate PRR (all anti-inflammatories vs non-antiinflammatories)
#calculate PRR (anti-inflammatories vs anti-inflammatories)
#construct network

##################################
#Libraries
##################################
library(data.table)
library(PhViD)
library(igraph)
##################################
#Inputs
##################################

inputs = scan("input_files.txt", what = "character", quiet = TRUE)

#read normalized FAERS data 
faers = fread(inputs[13])
  #FAERS normalization contains multimapping compounds; remove those
faers_sans4 = faers[V10 != "Type4"]

#read anti-inflammatory drug list
antiinflammatories = fread(inputs[1], header = FALSE)
aif                = antiinflammatories$V2

################################################################
#PRR CALCULATIONS
################################################################

#################
#01 Anti-inflammatories v Anti - inflammatories
#################

#calculate PRR (all anti-inflammatories vs all anti-inflammatories)

dt_aif <- copy(faers_sans4)
dt_aif <-dt_aif[V14%in%aif]

  #make a frequency count of drug - adverse effect 
freq_aif = dt_aif[, .(.N), by=list(V14, V9)]
  #make an object of class PhViD for calculations
phvid_aif = as.PhViD(freq_aif)

  #calculate Proportional Reporting Ratio
res_aif_aif <- PRR(phvid_aif)

#################
#02 All Anti-inflammatories v non Anti-inflammatories
#################
#calculate PRR (all anti-inflammatories vs non-antiinflammatories)

  #convert all anti-inflammatories to a single "ALL_ANTIINFLAMMATORIES"

dt_all <- copy(faers_sans4)
dt_all[V14%in%aif, V14:="ALL_ANTIINFLAMMATORIES"]
  #length(which(dt$V14 == "ALL_ANTIINFLAMMATORIES")) 1451208
  #length(which(faers_sans4$V14%in%aif)) 1451208
  #faers_sans4[which(faers_sans4$V14=="ALL_ANTIINFLAMMATORIES")[1:5]] NA NAs

  #make a frequency count of drug - adverse effect 
freq_all = dt_all[, .(.N), by=list(V14, V9)]
  #make an object of class PhViD for calculations
phvid_all = as.PhViD(freq_all)
  #calculate Proportional Reporting Ratio
res_all <- PRR(phvid_all)

#################
#02-variation All Anti-inflammatories v non Anti-inflammatories
#################
#calculate PRR (all anti-inflammatories vs non-antiinflammatories)

#MAKE NO CHANGES TO faers_sans4

dt_var <- copy(faers_sans4)

#make a frequency count of drug - adverse effect 
freq_var = dt_var[, .(.N), by=list(V14, V9)]
#make an object of class PhViD for calculations
phvid_var = as.PhViD(freq_var)
#calculate Proportional Reporting Ratio
res_var <- PRR(phvid_var)

#################
#02-second variation (dosvar) All Anti-inflammatories v non Anti-inflammatories
#################
#calculate PRR (all anti-inflammatories vs non-antiinflammatories)

#convert all anti-inflammatories to a single "ALL_ANTIINFLAMMATORIES"
#convert all NON anti-inflammatories to a single "NON"

dt_dosvar <- copy(faers_sans4)
dt_dosvar[!(V14%in%aif), V14:="NON"]
dt_dosvar[V14%in%aif, V14:="ALL_ANTIINFLAMMATORIES"]

#length(which(dt$V14 == "ALL_ANTIINFLAMMATORIES")) 1451208
#length(which(faers_sans4$V14%in%aif)) 1451208
#faers_sans4[which(faers_sans4$V14=="ALL_ANTIINFLAMMATORIES")[1:5]] NA NAs

#make a frequency count of drug - adverse effect 
freq_dosvar = dt_dosvar[, .(.N), by=list(V14, V9)]
#make an object of class PhViD for calculations
phvid_dosvar = as.PhViD(freq_dosvar)
#calculate Proportional Reporting Ratio
res_dosvar <- PRR(phvid_dosvar)

################################################################
#Analysis
################################################################

#adverse effects identified in ALL
ADR_all = unique(res_all$SIGNALS$`event effect`) #13645
length(res_all$SIGNALS$`event effect`) #543435
ADR_all_aif = unique(res_all$SIGNALS$`event effect`[res_all$SIGNALS$`drug code` == "ALL_ANTIINFLAMMATORIES"]) #2095
length(res_all$SIGNALS$`event effect`[res_all$SIGNALS$`drug code` == "ALL_ANTIINFLAMMATORIES"])

#adverse effects identified in VAR
ADR_var = unique(res_var$SIGNALS$`event effect`) #13645
length(res_var$SIGNALS$`event effect`)
ADR_var_aif = unique(res_var$SIGNALS$`event effect`[res_var$SIGNALS$`drug code` %in% aif]) #8844
length(res_var$SIGNALS$`event effect`[res_var$SIGNALS$`drug code` %in% aif])

#adverse effects identified in DOSVAR
ADR_dosvar = unique(res_dosvar$SIGNALS$`event effect`) #7820
length(res_dosvar$SIGNALS$`event effect`)
ADR_dosvar_aif = unique(res_dosvar$SIGNALS$`event effect`[res_dosvar$SIGNALS$`drug code` == "ALL_ANTIINFLAMMATORIES"]) #4287
length(res_dosvar$SIGNALS$`event effect`[res_dosvar$SIGNALS$`drug code` == "ALL_ANTIINFLAMMATORIES"])


#adverse effects identified in aif v aif
ADR_aif_aif = unique(res_aif_aif$SIGNALS$`event effect`) #7841
length(res_aif_aif$SIGNALS$`event effect`) # 21623

#intersections of aif_aif and the all approaches

ADR_aifaif_allaif = intersect(ADR_aif_aif, ADR_all_aif)
ADR_aifaif_varaif = intersect(ADR_aif_aif, ADR_var_aif)
ADR_aifaif_2aif = intersect(ADR_aif_aif, ADR_dosvar_aif)

df_a = res_aif_aif$SIGNALS[res_aif_aif$SIGNALS$`event effect`%in%ADR_aifaif_allaif,] #nrow 6445
df_b =res_aif_aif$SIGNALS[res_aif_aif$SIGNALS$`event effect`%in%ADR_aifaif_varaif,] #nrow 20781
df_c =res_aif_aif$SIGNALS[res_aif_aif$SIGNALS$`event effect`%in%ADR_aifaif_2aif,] #nrow 10235

g_a = graph_from_data_frame(df_a, directed = TRUE)
g_b = graph_from_data_frame(df_b, directed = TRUE)
g_c = graph_from_data_frame(df_c, directed = TRUE)
g_0 = graph_from_data_frame(res_aif_aif$SIGNALS, directed = TRUE)

V(g_a)$type = ifelse(V(g_a)$name%in%aif, "DRUG", "ADR")
V(g_b)$type = ifelse(V(g_b)$name%in%aif, "DRUG", "ADR")
V(g_c)$type = ifelse(V(g_c)$name%in%aif, "DRUG", "ADR")
V(g_0)$type = ifelse(V(g_0)$name%in%aif, "DRUG", "ADR")

degree(graph = g_a, v = V(g_a)[type == "DRUG"])
degree(graph = g_b, v = V(g_b)[type == "DRUG"])
degree(graph = g_c, v = V(g_c)[type == "DRUG"])
degree(graph = g_0, v = V(g_0)[type == "DRUG"])

components(g_a)$no
components(g_b)$no
components(g_c)$no
components(g_0)$no

strength(graph = g_a, vids = V(g_a)[type == "DRUG"], weights = E(g_a)$PRR)
strength(graph = g_b, vids = V(g_b)[type == "DRUG"], weights = E(g_b)$PRR)
strength(graph = g_c, vids = V(g_c)[type == "DRUG"], weights = E(g_c)$PRR)
strength(graph = g_0, vids = V(g_0)[type == "DRUG"], weights = E(g_0)$PRR)
strength(graph = g_a, vids = V(g_a)[type == "DRUG"], weights = E(g_a)$PRR) / degree(graph = g_a, v = V(g_a)[type == "DRUG"])
strength(graph = g_b, vids = V(g_b)[type == "DRUG"], weights = E(g_b)$PRR) / degree(graph = g_b, v = V(g_b)[type == "DRUG"])
strength(graph = g_c, vids = V(g_c)[type == "DRUG"], weights = E(g_c)$PRR) / degree(graph = g_c, v = V(g_c)[type == "DRUG"])
strength(graph = g_0, vids = V(g_0)[type == "DRUG"], weights = E(g_0)$PRR) / degree(graph = g_0, v = V(g_0)[type == "DRUG"])

strength(graph = g_0, vids = V(g_0)[type == "ADR"], weights = E(g_0)$PRR) / degree(graph = g_0, v = V(g_0)[type == "ADR"])

plot(g_0, vertex.label = "", vertex.color = ifelse(V(g_0)$type == "DRUG", "green", "yellow"))

plot(g_a, vertex.label = "", vertex.color = ifelse(V(g_0)$type == "DRUG", "green", "yellow"))

bipartite.projection(graph = g_a, types = ifelse(V(g_a)$type=="DRUG", TRUE, FALSE))


# 
# 
# 
# 
# #unique drugs 
# Drugs_Unique = unique(faers_sans4$V14)
# Drugs_Unique = Drugs_Unique[order(Drugs_Unique)]
# 
# #unique adverse events
# Pt_Unique = unique(faers_sans4$V9)
# Pt_Unique = Pt_Unique[order(Pt_Unique)]
# 
# #make a frequency count of drug - adverse effect 
# Drug_Adverse_Frequency = faers_sans4[, .(.N), by=list(V14, V9)]
# DAF_phvid = as.PhViD(Drug_Adverse_Frequency)
# 
# #calculate prr
# res <- PRR(DAF_phvid)
# 
# tail(res$SIGNALS)
# #make network
# 
# 
# aif_prr = res$SIGNALS[res$SIGNALS$`drug code`%in%aif, ]
# head(aif_prr)
# aif_df = aif_prr[,c(1,2,6)]
# 
# 
# #network must be bipartite
# aif_prr_g = graph_from_data_frame(d = aif_df, directed = TRUE)
# drugNodes = V(aif_prr_g)[V(aif_prr_g)$name%in%aif]
# ReactionNodes = V(aif_prr_g)[!(V(aif_prr_g)$name%in%aif)]
# V(aif_prr_g)$type <- ifelse(V(aif_prr_g)%in%drugNodes, "DRUG", "ADR")
# plot(aif_prr_g, 
#      vertex.label = "", 
#      vertex.color = ifelse(V(aif_prr_g)=="DRUG", 
#                            "blue", 
#                            "red")
#      )
# degree(graph = aif_prr_g, v = drugNodes)
# 
# 
# #Drug Nodes -- Adverse Event Nodes ... weighted edges, PRR
# 
# 
# ## version 2, calculate PRR AFTER subsetting 
# 
# Ainf_Adverse_Frequency = Drug_Adverse_Frequency[V14%in%aif]
# Ainf_phvid = as.PhViD(Ainf_Adverse_Frequency)
# Ainf_res <- PRR(Ainf_phvid)
# aif_prr_g_v2 = graph_from_data_frame(d = Ainf_res$SIGNALS[,c(1,2,6)], directed = TRUE)
# drugNodes_v2 = V(aif_prr_g_v2)[V(aif_prr_g_v2)$name%in%aif]
# ReactionNodes_v2 = V(aif_prr_g_v2)[!(V(aif_prr_g_v2)$name%in%aif)]
# 
# max(degree(aif_prr_g_v2, drugNodes_v2))
# min(degree(aif_prr_g_v2, drugNodes_v2))
# max(degree(aif_prr_g_v2, ReactionNodes_v2))
# min(degree(aif_prr_g_v2, ReactionNodes_v2))
# 
# which(degree(aif_prr_g_v2, ReactionNodes_v2)==17)
# which(degree(aif_prr_g, ReactionNodes)==17)
# 
# max(strength(graph = aif_prr_g_v2, vids = drugNodes_v2, weights = E(aif_prr_g_v2)$PRR))
# 
# indice_charro = strength(graph = aif_prr_g_v2, vids = drugNodes_v2, weights = E(aif_prr_g_v2)$PRR)/degree(aif_prr_g_v2, drugNodes_v2)
# indice_charro[order(indice_charro, decreasing = TRUE)]
# 
