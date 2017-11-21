library(data.table)
library(PhViD)
setwd(dir = "/extData/HurLab_Projects/AInf_Drugs/results/")

#read the normalized faers data
faers = fread(input = "drug_faers.txt_processed.txt")

#remove Type4 samples 
faers_sans4 = faers[jh_type != "Type4"]

#unique drugs 
Drugs_Unique = unique(faers_sans4$jh_final)
Drugs_Unique = Drugs_Unique[order(Drugs_Unique)]
#unique adverse events
Pt_Unique = unique(faers_sans4$jh_pt_term)
Pt_Unique = Pt_Unique[order(Pt_Unique)]
#make a frequency count of drug - adverse effect 

Drug_Adverse_Frequency = faers_sans4[, .(.N), by=list(jh_final, jh_pt_term)]
DAF_phvid = as.PhViD(Drug_Adverse_Frequency)

#calculate prr

res <- PRR(DAF_phvid)
###


nrow(res$ALLSIGNALS)
#
df = res$SIGNALS
length(unique(df$`event effect`))

max(df$PRR)

#
head(Drug_Adverse_Frequency)

replacements = unique(faers$replacement)
originals = unique(faers$original)
jeyhurs = unique(faers$jh_final)

length(intersect(originals, jeyhurs))


