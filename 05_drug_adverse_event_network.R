library(data.table)
library(PhViD)
inputs = scan(inputs)

#read the normalized faers data
faers = fread(inputs[x])

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


#make network
#Drug Nodes -- Adverse Event Nodes ... weighted edges, PRR