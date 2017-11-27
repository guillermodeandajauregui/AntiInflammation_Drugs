source("libraries/libraries.R")
library(ChemmineR)

inputs = scan(file = "input_files.txt", what = "character")

#read structures from the Open SDF file
sdftest_open = read.SDFset(sdfstr = inputs[15])

#function to convert text to "Proper sentence case"
proper=function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))

#list of drugs for cmap
#DRUGS = names(V(graph = nw_fgsea_cmap)[type == "DRUG"])
DRUGS = read.csv(inputs[1], 
                 sep = "\t", 
                 stringsAsFactors = FALSE, 
                 header = FALSE)
DRUGS = DRUGS$V2

#search for drugs, make a SDF list 
lista_sdfs = list()
identificadores = character()
for(i in DRUGS){
  print(i)
  sdf_provisional = grepSDFset(pattern = proper(i), 
                               x = sdftest_open, 
                               field = "datablock", 
                               mode = "subset", 
                               ignore.case = FALSE)
  longines = length(sdf_provisional)
  print(longines)
  lista_sdfs = c(lista_sdfs, sdf_provisional)
  identificadores = c(identificadores, rep(i, longines))
}

#remove duplicates
a = table(identificadores)
b = names(which(table(identificadores)>1))
c = which(identificadores%in%names(which(table(identificadores)>1)))
d = unique(identificadores)
e = c[-c(1,3,12)]
lista_sdfs_prime <- lista_sdfs
lista_sdfs <- lista_sdfs[-e]
identificadores <- identificadores[-e]
rm(a,b,c,d,e)
# #remove item 34, duplicate flurbiprofen 
# lista_sdfs <- lista_sdfs[-34]
# identificadores <- identificadores[-34]

#sdf_ainf_cmap = SDFset(SDFlist = lista_sdfs, ID = identificadores)
sdf_ainf = SDFset(SDFlist = lista_sdfs, ID = identificadores)

#adding the drugs not in open (from SMILES - Wikipedia - PubChem)
drugs_nofound = DRUGS[!(DRUGS%in%identificadores)]
drugs_nofound[order(drugs_nofound)]

drugs_nofound_smiles = c(
                         "ACEMETACIN" = "Clc1ccc(cc1)C(=O)n3c2ccc(OC)cc2c(c3C)CC(=O)OCC(=O)O",
                         "ALMINOPROFEN" = "O=C(O)C(c1ccc(NC/C(=C)C)cc1)C",
                         "ALOXIPRIN" = "[Al+3].O=C(Oc1ccccc1C([O-])=O)C.[OH-].[O-]C(=O)c1ccccc1OC(=O)C",
                         "AMINOPROPIONITRILE" = "NCCC#N",
                         "BENORILATE" = "O=C(C)Oc2ccccc2C(=O)Oc1ccc(NC(C)=O)cc1",
                         "BUFEXAMAC" = "ONC(=O)Cc1ccc(OCCCC)cc1",
                         "BUMADIZONE" = "O=C(O)C(C(=O)N(Nc1ccccc1)c2ccccc2)CCCC",
                         #"CARBASALATE CALCIUM" NON AVAILABLE
                         #"CHOLINE SALICYLATE", NON AVAILABLE
                         "CHONDROITIN SULFATE" ="CC(=O)NC1C(C(C(OC1O)OS(=O)(=O)O)O)OC2C(C(C(C(O2)C(=O)O)O)O)O",
                         "CLOPREDNOL" = "O=C\\1\\C=C/[C@]4(C(=C/1)C(\\Cl)=C/[C@@H]2[C@@H]4[C@@H](O)C[C@@]3([C@@](O)(C(=O)CO)CC[C@@H]23)C)C",
                         "CORTIVAZOL" = "O=C(OCC(=O)[C@@]4(O)[C@H](C)C[C@H]5[C@@H]6/C=C(\\C3=C\\c1c(cnn1c2ccccc2)C[C@@]3([C@H]6[C@@H](O)C[C@]45C)C)C)C",
                         "DERACOXIB" = "O=S(=O)(c3ccc(n1nc(cc1c2ccc(OC)c(F)c2)C(F)F)cc3)N",
                         "DIACEREIN" = "O=C(Oc3cccc2C(=O)c1cc(cc(OC(=O)C)c1C(=O)c23)C(=O)O)C",
                         "DIFENPIRAMIDE" = "O=C(Nc1ncccc1)Cc3ccc(c2ccccc2)cc3",
                         "DIPYROCETYL" = "O=C(Oc1cccc(c1OC(=O)C)C(=O)O)C",
                         "ETHENZAMIDE" = "O=C(c1ccccc1OCC)N",
                         "FENTIAZAC" = "c1ccc(cc1)c2nc(c(s2)CC(=O)O)c3ccc(cc3)Cl", 
                         "FEPRAZONE" = "O=C2N(c1ccccc1)N(C(=O)C2C\\C=C(/C)C)c3ccccc3",
                         "FLUFENAMIC ACID"= "FC(F)(F)c1cc(ccc1)Nc2ccccc2C(=O)O",
                         "FLUMETASONE" = "C[C@@H]1C[C@H]2[C@@H]3C[C@@H](C4=CC(=O)C=C[C@@]4([C@]3([C@H](C[C@@]2([C@]1(C(=O)CO)O)C)O)F)C)F",
                         "FLUNIXIN" = "CC1=C(C=CC=C1NC2=C(C=CC=N2)C(=O)O)C(F)(F)F",
                         "FLUNOXAPROFEN" = "C[C@@H](c1ccc2c(c1)nc(o2)c3ccc(cc3)F)C(=O)O",
                         #"GLUCOSAMINOGLYCAN POLYSULFATE" not available
                         "GRAPIPRANT" = "CCC1=NC2=C(N=C(C=C2N1C3=CC=C(C=C3)CCNC(=O)NS(=O)(=O)C4=CC=C(C=C4)C)C)C",
                         "GUACETISAL" = "O=C(Oc1ccccc1C(=O)Oc2ccccc2OC)C",
                         "IMIDAZOLE SALICYLATE" = "C1=CC=C(C(=C1)C(=O)O)O.C1=CN=CN1",
                         "LONAZOLAC" = "c1ccc(cc1)n2cc(c(n2)c3ccc(cc3)Cl)CC(=O)O",
                         "MAVACOXIB" = "C1=CC(=CC=C1C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F)F",
                         "MEPREDNISONE" = "O=C(CO)[C@@]3(O)[C@@H](C)C[C@H]2[C@@H]4CC\\C1=C\\C(=O)\\C=C/[C@]1(C)[C@H]4C(=O)C[C@@]23C",
                         "MOFEBUTAZONE" = "CCCCC1C(=O)NN(C1=O)C2=CC=CC=C2",
                         #"MORPHOLINE SALICYLATE" not available
                         "NAPROXCINOD" = "C[C@@H](C1=CC2=C(C=C1)C=C(C=C2)OC)C(=O)OCCCCO[N+](=O)[O-]",
                         "NIFLUMIC ACID" = "C1=CC(=CC(=C1)NC2=C(C=CC=N2)C(=O)O)C(F)(F)F",
                         #"ORGOTEIN" not available
                         "OXACEPROL" = "CC(=O)N1C[C@@H](C[C@H]1C(=O)O)O",
                         "OXAMETACIN" = "CC1=C(C2=C(N1C(=O)C3=CC=C(C=C3)Cl)C=CC(=C2)OC)CC(=O)NO",
                         "PENTOSAN POLYSULFATE" = "C1C(C(C(C(O1)OC2COC(C(C2OS(=O)(=O)O)OS(=O)(=O)O)O)OS(=O)(=O)O)OS(=O)(=O)O)O",
                         "PIRPROFEN" = "O=C(O)C(c1cc(Cl)c(cc1)N2C/C=C\\C2)C",
                         "POTASSIUM SALICYLATE" = "[K+].O=C([O-])c1ccccc1O",
                         "PREDNYLIDENE" = "OCC(=O)[C@@]2(O)C(=C)C[C@H]3[C@@H]4CC\\C1=C\\C(=O)\\C=C/[C@]1(C)[C@H]4[C@@H](O)C[C@]23C", 
                         "PROGLUMETACIN" = "Clc1ccc(cc1)C(=O)n3c2ccc(OC)cc2c(c3C)CC(=O)OCCN4CCN(CC4)CCCOC(=O)CCC(C(=O)N(CCC)CCC)NC(=O)c5ccccc5",
                         "PROQUAZONE" = "CC1=CC2=C(C=C1)C(=NC(=O)N2C(C)C)C3=CC=CC=C3",
                         "ROBENACOXIB" = "O=C(O)Cc2c(Nc1c(F)c(F)cc(F)c1F)ccc(c2)CC",
                         "SODIUM SALICYLATE" = "[Na+].O=C([O-])c1ccccc1O",
                         "SUXIBUZONE" = "O=C(O)CCC(=O)OCC2(C(=O)N(c1ccccc1)N(C2=O)c3ccccc3)CCCC",
                         "TENIDAP" = "c1cc(sc1)C(=O)c2c3cc(ccc3n(c2O)C(=O)N)Cl",
                         "TEPOXALIN" = "CN(C(=O)CCc1cc(n(n1)c2ccc(cc2)OC)c3ccc(cc3)Cl)O",
                         "TOLFENAMIC ACID" = "Clc2cccc(Nc1ccccc1C(=O)O)c2C",
                         "VEDAPROFEN" = "O=C(O)C(c2ccc(c1ccccc12)C3CCCCC3)C"
                         )

# cmap_non_found_smiles = c("FLUFENAMIC ACID"= "FC(F)(F)c1cc(ccc1)Nc2ccccc2C(=O)O",
#                           "TOLFENAMIC ACID" = "Clc2cccc(Nc1ccccc1C(=O)O)c2C",
#                           "BUFEXAMAC" = "ONC(=O)Cc1ccc(OCCCC)cc1",
#                           "NIFLUMIC ACID" = "C1=CC(=CC(=C1)NC2=C(C=CC=N2)C(=O)O)C(F)(F)F",
#                           "ACEMETACIN" = "Clc1ccc(cc1)C(=O)n3c2ccc(OC)cc2c(c3C)CC(=O)OCC(=O)O",
#                           "SUXIBUZONE" = "O=C(O)CCC(=O)OCC2(C(=O)N(c1ccccc1)N(C2=O)c3ccccc3)CCCC",
#                           "FLUMETASONE" = "C[C@@H]1C[C@H]2[C@@H]3C[C@@H](C4=CC(=O)C=C[C@@]4([C@]3([C@H](C[C@@]2([C@]1(C(=O)CO)O)C)O)F)C)F",
#                           "FLUNIXIN" = "CC1=C(C=CC=C1NC2=C(C=CC=N2)C(=O)O)C(F)(F)F"
#                           )

#cmap_non_found_sdfset = smiles2sdf(cmap_non_found_smiles)
#cmap_combined_smiles = c(sdf_ainf_cmap, cmap_non_found_sdfset)

non_found_sdfset = smiles2sdf(drugs_nofound_smiles)
combined_smiles = c(sdf_ainf, non_found_sdfset)

# cmap_ap = sdf2ap(cmap_combined_smiles)
AtomPairs = sdf2ap(combined_smiles)

# length(cmap_ap)
# cmp.similarity(a = cmap_ap[1], b = cmap_ap[47], mode = 1)

#make matrix
#cmap_tanimoto_matrix = matrix(nrow = 47, ncol = 47, dimnames = list(cmap_ap@ID, cmap_ap@ID))
tanimoto_matrix = matrix(nrow = length(AtomPairs), 
                         ncol = length(AtomPairs), 
                         dimnames = list(AtomPairs@ID, 
                                         AtomPairs@ID))
#fill with tanimotos
for(i in 1:length(AtomPairs)){
  for(j in 1:length(AtomPairs)){
    tanimoto_matrix[i,j] = cmp.similarity(AtomPairs[i], AtomPairs[j], mode = 1)
  }
}

tanimoto_matrix

tanimoto_g = graph_from_adjacency_matrix(tanimoto_matrix, 
                                         weighted = TRUE, 
                                         mode = "undirected", 
                                         diag = FALSE)

write.table(tanimoto_matrix, file = "results/tanimoto_matrix_complete.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = TRUE)

write.graph(tanimoto_g, file = "results/tanimoto_g_complete.gml", format = "gml")
#plot(cmap_tanimoto_g)
#E(cmap_tanimoto_g)[weight >= 0.50]
quantile(E(cmap_tanimoto_g)$weight, 0.90) #nice cut
quantile(E(cmap_tanimoto_g)$weight, 0.66) #comparable to perturbation similarity network
#quantile(E(cmap_tanimoto_g)$weight, 0.66) #comparable to functional network
#cmap_tanimoto_g = delete.edges(graph = cmap_tanimoto_g, E(cmap_tanimoto_g)[weight <= quantile(E(cmap_tanimoto_g)$weight, 0.66)])

tanimoto_g_cut = funcion_corte_optimo(tanimoto_matrix)$g #0.65
quantile(E(tanimoto_g)$weight, 0.99)
components(tanimoto_g)

#plot(cmap_tanimoto_g)
#components(cmap_tanimoto_g)
cmap_tanimoto_infomap = infomap.community(cmap_tanimoto_g)
cmap_tanimoto_walktrap = walktrap.community(cmap_tanimoto_g)
#cmap_tanimoto_infomap$modularity
#plot(cmap_tanimoto_infomap, cmap_tanimoto_g)
#plot(cmap_tanimoto_walktrap, cmap_tanimoto_g)
# table(degree(cmap_tanimoto_g))
# degree(cmap_tanimoto_g)
# transitivity(cmap_tanimoto_g)
# average.path.length(cmap_tanimoto_g)
# V(cmap_tanimoto_g)
# cmap_tanimoto_g

red_er = erdos.renyi.game(47, p.or.m = 363, type = "gnm")
plot(red_er)
transitivity(red_er)
average.path.length(red_er)
infomap.community(red_er)

save(cmap_tanimoto_matrix, cmap_tanimoto_g, red_er, cmap_tanimoto_infomap, cmap_tanimoto_walktrap, file = "results/structural_similarity.RData")
