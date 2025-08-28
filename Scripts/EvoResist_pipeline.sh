#!/bin/bash

## This script is used to:
## Find SNP and indel for each drug.

module add r
## snp dir: /work/users/l/i/lingbo1/Mtb/WHO_denovo/source_data/all_ann_convergent_flt_v3.txt

## indels:
#### folder 1: /proj/qliulab/MTB_phy_db/indel_ano/
bash /work/users/l/i/lingbo1/Mtb/WHO_denovo/script/00-indel.sh
#### folder 2: /proj/qliulab/01_Mtb/03_indel/
find /proj/qliulab/01_Mtb/03_indel/ -name "*.indel.ano" -print0 | xargs -0 cat > tmp_second.txt
cat tmp_first.txt tmp_second.txt | cut -f1-5,8 | sort | uniq -c > /work/users/l/i/lingbo1/Mtb/WHO_denovo/source_data/all_indel_100k.txt
sed -i -E 's/^[[:space:]]+//; s/[[:space:]]+/\t/g' "/work/users/l/i/lingbo1/Mtb/WHO_denovo/source_data/all_indel_100k.txt"
rm tmp_second.txt tmp_first.txt


## make input list:
snp_100k_filter_path=/work/users/l/i/lingbo1/Mtb/WHO_denovo/source_data/all_ann_convergent_flt_v3.txt
indel_100k_path=/work/users/l/i/lingbo1/Mtb/WHO_denovo/source_data/all_indel_100k.txt

### RIF
# rpoB +: main region:(759807,763325), upstream 500bp: (759307,763325)
drug_name=RIF
out_dir=/work/users/l/i/lingbo1/Mtb/WHO_denovo/results/${drug_name}/PGC_cutoff3
mkdir ${out_dir}
awk -F' ' '( index($4, "-") == 0 && ($4 >= 759307 && $4 <= 763325))' ${snp_100k_filter_path} > ${out_dir}/denovo_snp.txt
awk -F'\t' '($2 >= 759807 && $2 <= 763325)' ${indel_100k_path} > ${out_dir}/denovo_indel.txt


### INH
## inhA + (1674202,1675011), including promoter:(1673336-1675011)
## katG - (2153889, 2156111), upstream 500bp (2153889, 2156611)
drug_name=INH
out_dir=/work/users/l/i/lingbo1/Mtb/WHO_denovo/results/${drug_name}/PGC_cutoff3
mkdir ${out_dir}
awk -F' ' '( index($4, "-") == 0 && (($4 >= 1673336 && $4 <= 1675011) || ($4 >= 2153889 && $4 <= 2156611)))' ${snp_100k_filter_path} > ${out_dir}/denovo_snp.txt
awk -F'\t' '(($2 >= 1674202 && $2 <= 1675011) || ($2 >= 2153889 && $2 <= 2156111))' ${indel_100k_path} > ${out_dir}/denovo_indel.txt

### EMB
# embB + (4246514, 4249810), upstream 500bp: (4246014, 4249810)
drug_name=EMB
out_dir=/work/users/l/i/lingbo1/Mtb/WHO_denovo/results/${drug_name}/PGC_cutoff3
mkdir ${out_dir}
awk -F' ' '( index($4, "-") == 0 && ($4 >= 4246014 && $4 <= 4249810))' ${snp_100k_filter_path} > ${out_dir}/denovo_snp.txt
awk -F'\t' '($2 >= 4246514 && $2 <= 4249810)' ${indel_100k_path} > ${out_dir}/denovo_indel.txt

### PZA
# pncA - (2288681,2289241), upstream 500bp: (2288681, 2289741)
drug_name=PZA
out_dir=/work/users/l/i/lingbo1/Mtb/WHO_denovo/results/${drug_name}/PGC_cutoff3
mkdir ${out_dir}
awk -F' ' '( index($4, "-") == 0 && ($4 >= 2288681 && $4 <= 2289741))' ${snp_100k_filter_path} > ${out_dir}/denovo_snp.txt
awk -F'\t' '($2 >= 2288681 && $2 <= 2289241)' ${indel_100k_path} > ${out_dir}/denovo_indel.txt

### LFX, MFX
# gyrA + (7302, 9818), gyrB + (5240,7267), whole region(5240-9818), upstream 500bp (4740, 9818)
drug_name=LFX
drug_name=MFX
out_dir=/work/users/l/i/lingbo1/Mtb/WHO_denovo/results/${drug_name}/PGC_cutoff3
mkdir ${out_dir}
awk -F' ' '( index($4, "-") == 0 && ($4 >= 4740 && $4 <= 9818))' ${snp_100k_filter_path} > ${out_dir}/denovo_snp.txt
awk -F'\t' '(($2 >= 5240 && $2 <= 7262) || ($2 >= 7302 && $2 <= 9818))' ${indel_100k_path} > ${out_dir}/denovo_indel.txt

### BDQ
## Rv0678 + (778990,779487), upstream 500bp (778490,779487)
## atpE + (1461045,1461290), upstream 500bp (1460545, 1461290)
## pepQ - (2859300, 2860418), upstream 500bp (2858300, 2860918)
drug_name=BDQ
out_dir=/work/users/l/i/lingbo1/Mtb/WHO_denovo/results/${drug_name}/PGC_cutoff3
mkdir ${out_dir}
awk -F' ' '( index($4, "-") == 0 && (($4 >= 778490 && $4 <= 779487) || ($4 >= 1460545 && $4 <= 1461290) || ($4 >= 2858300 && $4 <= 2860918)))' ${snp_100k_filter_path} > ${out_dir}/denovo_snp.txt
awk -F'\t' '(($2 >= 778990 && $2 <= 779487) || ($2 >= 1461045 && $2 <= 1461290) || ($2 >= 2859300 && $2 <= 2860418))' ${indel_100k_path} > ${out_dir}/denovo_indel.txt


### AMK
## rrs + (1471846, 1473382), upstream 500bp: (1471346, 1473382);
## eis promoter: (2715332, 2715832)
drug_name=AMK
out_dir=/work/users/l/i/lingbo1/Mtb/WHO_denovo/results/${drug_name}/PGC_cutoff3_witheis
mkdir ${out_dir}
awk -F' ' '( index($4, "-") == 0 && (($4 >= 1471346 && $4 <= 1473382) || ($4 >= 2715332 && $4 <= 2715832)))' ${snp_100k_filter_path} > ${out_dir}/denovo_snp.txt
awk -F'\t' '(($2 >= 1471846 && $2 <= 1473382) || ($2 >= 2715332 && $2 <= 2715832))' ${indel_100k_path} > ${out_dir}/denovo_indel.txt

## KAN
## rrs + (1471846, 1473382), upstream 500bp: (1471346, 1473382);
## eis - (2714124,2715332), upstream 500bp: (2714124, 2715832)
drug_name=KAN
out_dir=/work/users/l/i/lingbo1/Mtb/WHO_denovo/results/${drug_name}/PGC_cutoff3_witheis
mkdir ${out_dir}
awk -F' ' '( index($4, "-") == 0 && (($4 >= 1471346 && $4 <= 1473382) || ($4 >= 2714124 && $4 <= 2715832)))' ${snp_100k_filter_path} > ${out_dir}/denovo_snp.txt
awk -F'\t' '(($2 >= 1471846 && $2 <= 1473382) || ($2 >= 2714124 && $2 <= 2715832))' ${indel_100k_path} > ${out_dir}/denovo_indel.txt


## CAP
## rrs + (1471846, 1473382), upstream 500bp: (1471346, 1473382);
## tlyA + (1917940, 1918746), upstream 500bp: (1917440, 1918746)
drug_name=CAP
out_dir=/work/users/l/i/lingbo1/Mtb/WHO_denovo/results/${drug_name}/PGC_cutoff3
mkdir ${out_dir}
awk -F' ' '( index($4, "-") == 0 && (($4 >= 1471346 && $4 <= 1473382) || ($4 >= 1917440 && $4 <= 1918746)))' ${snp_100k_filter_path} > ${out_dir}/denovo_snp.txt
awk -F'\t' '(($2 >= 1471846 && $2 <= 1473382) || ($2 >= 1917940 && $2 <= 1918746))' ${indel_100k_path} > ${out_dir}/denovo_indel.txt

### STM
## rpsL + (781560, 781934), upstream 500bp: (781060, 781934); 
## gid - (4407528, 4408202), upstream 500bp (4407528, 4408703);
## rrs (only 1472359, 1472362)
drug_name=STM
out_dir=/work/users/l/i/lingbo1/Mtb/WHO_denovo/results/${drug_name}/PGC_cutoff3
mkdir ${out_dir}
awk -F' ' '( index($4, "-") == 0 && (($4 >= 781060 && $4 <= 781934) || ($4 >= 4407528 && $4 <= 4408703) || ($4 == 1472359 || $4 == 1472362)))' ${snp_100k_filter_path} > ${out_dir}/denovo_snp.txt
awk -F'\t' '(($2 >= 781560 && $2 <= 781934) || ($2 >= 4407528 && $2 <= 4408202))' ${indel_100k_path} > ${out_dir}/denovo_indel.txt

## ETO
## ethA - (4326004, 4327473), upstream 500bp (4326004, 4327973)
## inhA + (1674202,1675011), including promoter:(1673336-1675011)
drug_name=ETO
out_dir=/work/users/l/i/lingbo1/Mtb/WHO_denovo/results/${drug_name}/PGC_cutoff3
mkdir ${out_dir}
awk -F' ' '( index($4, "-") == 0 && (($4 >= 4326004 && $4 <= 4327973) || ($4 >= 1673336 && $4 <= 1675011)))' ${snp_100k_filter_path} > ${out_dir}/denovo_snp.txt
awk -F'\t' '(($2 >= 4326004 && $2 <= 4327473) || ($2 >= 1674202 && $2 <= 1675011))' ${indel_100k_path} > ${out_dir}/denovo_indel.txt
