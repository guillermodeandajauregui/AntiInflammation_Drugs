#!/bin/bash
#pipeline for analysis
Rscript scripts/01_drugs_in_resources.R
Rscript scripts/02_extract_expression_profiles.R
Rscript scripts/03_drug_drug_similarity_network.R
Rscript scripts/04_drug_fgsea.R
Rscript scripts/05_drug_adverse_event_network.R
Rscript scripts/06_multigraph.R
Rscript scripts/07_hlt_pt_drug_graph.R
Rscript scripts/08_structural_similarity.R

#analysis of networks with networkx 
# #This section was originally done manually; code untested, added for completeness
# for i in $(ls results/*gml);
#   do
#     a=$(echo $i".nx")
#     tail -n+3 $i > $a
#     sed -i 'N;s/\s\+\[/ \[/g;P;D' $a
#     sed -i 's/name/label/g' $a
#   done
# for i in $(ls results/*.nx);
#   do
#     python libraries/BipartiteAnalysis.py $i $i
#   done
#