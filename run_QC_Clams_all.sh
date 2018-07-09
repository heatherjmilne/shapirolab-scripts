#/bin/bash

# HJM 7/9/18
# Running this pipeline with multiple samples, for the first time.
# The samples are Rachel Turba's clams
# First they will be moved to a new directory so that the other samples beginning with RCH are not included as well.

cp /projects/redser2/raw/180419_M00160_0078_000000000-BRWLV/RCH9* /projects/redser3-notbackedup/projects/hjmilne/Recharge/Rachel_Turba/Clam/Raw_Data
cp /projects/redser2/raw/180419_M00160_0078_000000000-BRWLV/RCH10* /projects/redser3-notbackedup/projects/hjmilne/Recharge/Rachel_Turba/Clam/Raw_Data
cp /projects/redser2/raw/180419_M00160_0078_000000000-BRWLV/RCH11* /projects/redser3-notbackedup/projects/hjmilne/Recharge/Rachel_Turba/Clam/Raw_Data
cp /projects/redser2/raw/180419_M00160_0078_000000000-BRWLV/RCH12* /projects/redser3-notbackedup/projects/hjmilne/Recharge/Rachel_Turba/Clam/Raw_Data
cp /projects/redser2/raw/180419_M00160_0078_000000000-BRWLV/RCH19* /projects/redser3-notbackedup/projects/hjmilne/Recharge/Rachel_Turba/Clam/Raw_Data
cp /projects/redser2/raw/180419_M00160_0078_000000000-BRWLV/RCH22* /projects/redser3-notbackedup/projects/hjmilne/Recharge/Rachel_Turba/Clam/Raw_Data

# Run Pre-Mapping once for each sample
bash /projects/redser3-notbackedup/projects/hjmilne/Scripts/shapirolab-scripts/monica-hjm-v3-premapping.sh -p=RCH -t=clam -dir=/projects/redser3-notbackedup/projects/hjmilne/Recharge/Rachel_Turba/Clam/Raw_Data &> RCH_clam_premapping.log.txt
wait

# Run blastn once for each sample
bash /projects/redser3-notbackedup/projects/hjmilne/Scripts/shapirolab-scripts/monica-hjm-v3-blast.sh -p=RCH4 -t=clam &> RCH_clam_blastn.log.txt

# Run Nuclear Mapping once for each genome
bash /projects/redser3-notbackedup/projects/hjmilne/Scripts/shapirolab-scripts/monica-hjm-v3-mapping.sh -p=RCH -t=clam -ref=/projects/redser2/genomes/Corbicula_fluminea/GCA_001632725.1_ASM163272v1_genomic.fna -refname=CorbFlum_nuclear &> RCH_clam_mapping_CorbFlum.log.txt

bash /projects/redser3-notbackedup/projects/hjmilne/Scripts/shapirolab-scripts/monica-hjm-v3-mapping.sh -p=RCH -t=clam -ref=/projects/redser2/genomes/Bankia_setacea/GCA_001922985.1_ASM192298v1_genomic.fna -refname=BankSeta_nuclear &> RCH_clam_mapping_BankSeta.log.txt

bash /projects/redser3-notbackedup/projects/hjmilne/Scripts/shapirolab-scripts/monica-hjm-v3-mapping.sh -p=RCH -t=clam -ref=/projects/redser2/genomes/Dreissena_polymorpha/GCA_000806325.1_ASM80632v1_genomic.fna -refname=DreiPoly_nuclear &> RCH_clam_mapping_DreiPoly.log.txt

# Run MIA once for each mitogenome
bash  /projects/redser3-notbackedup/projects/hjmilne/Scripts/shapirolab-scripts/monica-hjm-v3-mia.sh -p=RCH -t=clam -mref=/projects/redser2/genomes/mitogenomes/Anodonta_cygnea/Anodonta_cygnea_fullmito_FEMALE_NC_036488.1.fasta -mname=AnoCyg-female-mito &> RCH_clam_AnoCyg_female_mia.log.txt 

bash  /projects/redser3-notbackedup/projects/hjmilne/Scripts/shapirolab-scripts/monica-hjm-v3-mia.sh -p=RCH -t=clam -mref=/projects/redser2/genomes/mitogenomes/Anodonta_cygnea/Anodonta_cygnea_fullmito_MALE_NC_036488.1.fasta -mname=AnoCyg-male-mito &> RCH_clam_AnoCyg_male_mia.log.txt 

bash  /projects/redser3-notbackedup/projects/hjmilne/Scripts/shapirolab-scripts/monica-hjm-v3-mia.sh -p=RCH -t=clam -mref=/projects/redser2/genomes/mitogenomes/Anodonta_cygnea/Anodonta_cygnea_fullmito_NC_036488.1.fasta -mname=AnoCyg-herm-mito &> RCH_clam_AnoCyg_herm_mia.log.txt 

rm -r /projects/redser3-notbackedup/projects/hjmilne/Recharge/Rachel_Turba/Clam/Raw_Data
