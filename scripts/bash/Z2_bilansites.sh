#!/bin/bash

tail -n+2 metapr2_ASVs_selected_abundance_Dinophyceae_2022-03-03.tsv | awk -F"\t" 
'BEGIN{FS=OFS="\t"}{gsub(" 
","_",$23);gsub(" ","_",$10); if($10=="") print $16"%"$17"%noproject%"$23,$1 ; else print 
$16"%"$17"%"$10"%"$23,$1}' | sort | uniq | awk 'NR==1{h=$0;next}{a[$1]=a[$1]"|"$2}END{print h;for (i in 
a){sub("^\\|","",a[i]);print i,a[i]}}' | sed 's/%/\t/g' | sed 's/ /\t/g' | awk -v s=1 'BEGIN{FS=OFS="\t"; 
print "latitude,longitude,project,projectDetails,listSamples,Nbsamples"}{print $1,$2,$3, $4, $5, gsub(/\|/, 
"", $5)+s}' | sed 's/,/\t/g' | awk 'BEGIN{FS=OFS="\t"}{print $1, $2, $3, $4, $6, $5}' > 
Lat_Long_Proj_ProjD_NbS_ListS
