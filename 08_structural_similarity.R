library(ChemmineR)
#read structures from the Open SDF file
sdftest_open = read.SDFset(sdfstr = "~/openstructures.sdf")

#function to convert text to "Proper sentence case"
proper=function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))

#list of drugs for cmap
DRUGS = names(V(graph = nw_fgsea_cmap)[type == "DRUG"])

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

#remove item 34, duplicate flurbiprofen 
lista_sdfs <- lista_sdfs[-34]
identificadores <- identificadores[-34]

sdf_ainf_cmap = SDFset(SDFlist = lista_sdfs, ID = identificadores)

#adding the drugs not in open (from SMILES - Wikipedia)

cmap_non_found_smiles = c("FLUFENAMIC ACID"= "FC(F)(F)c1cc(ccc1)Nc2ccccc2C(=O)O",
                          "TOLFENAMIC ACID" = "Clc2cccc(Nc1ccccc1C(=O)O)c2C",
                          "BUFEXAMAC" = "ONC(=O)Cc1ccc(OCCCC)cc1",
                          "NIFLUMIC ACID" = "C1=CC(=CC(=C1)NC2=C(C=CC=N2)C(=O)O)C(F)(F)F",
                          "ACEMETACIN" = "Clc1ccc(cc1)C(=O)n3c2ccc(OC)cc2c(c3C)CC(=O)OCC(=O)O",
                          "SUXIBUZONE" = "O=C(O)CCC(=O)OCC2(C(=O)N(c1ccccc1)N(C2=O)c3ccccc3)CCCC",
                          "FLUMETASONE" = "C[C@@H]1C[C@H]2[C@@H]3C[C@@H](C4=CC(=O)C=C[C@@]4([C@]3([C@H](C[C@@]2([C@]1(C(=O)CO)O)C)O)F)C)F",
                          "FLUNIXIN" = "CC1=C(C=CC=C1NC2=C(C=CC=N2)C(=O)O)C(F)(F)F"
                          )

cmap_combined_smiles = c(sdf_ainf_cmap, cmap_non_found_sdfset)

cmap_ap = sdf2ap(cmap_combined_smiles)
length(cmap_ap)
cmp.similarity(a = cmap_ap[1], b = cmap_ap[47], mode = 1)

#make matrix
cmap_tanimoto_matrix = matrix(nrow = 47, ncol = 47, dimnames = list(cmap_ap@ID, cmap_ap@ID))
#fill with tanimotos
for(i in 1:47){
  for(j in 1:47){
    cmap_tanimoto_matrix[i,j] = cmp.similarity(cmap_ap[i], cmap_ap[j], mode = 1)
  }
}

cmap_tanimoto_matrix

cmap_tanimoto_g = graph_from_adjacency_matrix(cmap_tanimoto_matrix, weighted = TRUE, mode = "undirected", diag = FALSE)
#plot(cmap_tanimoto_g)
#E(cmap_tanimoto_g)[weight >= 0.50]
cmap_tanimoto_g = delete.edges(graph = cmap_tanimoto_g, E(cmap_tanimoto_g)[weight <= 0.25])
#plot(cmap_tanimoto_g)
#components(cmap_tanimoto_g)
cmap_tanimoto_infomap = infomap.community(cmap_tanimoto_g)
cmap_tanimoto_walktrap = walktrap.community(cmap_tanimoto_g)
#cmap_tanimoto_infomap$modularity
plot(cmap_tanimoto_infomap, cmap_tanimoto_g)
plot(cmap_tanimoto_walktrap, cmap_tanimoto_g)
