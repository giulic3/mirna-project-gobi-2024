# git clone https://github.com/Meffert-Lab/SCRAP.git
# conda activate r4-base
# mamba env create -f SCRAP/SCRAP_environment.yml -n SCRAP
# conda create --clone /home/rmaji/miniconda3/envs/r4-base/envs/SCRAP -n SCRAP
# conda activate SCRAP


#download trascriptome fasta from gencode
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.transcripts.fa.gz


#Building the reference
# bash SCRAP/bin/Reference_Installation.sh \
#     -r SCRAP/ \
#     -m hsa \
#     -g hg38 \
#     -s human


bash SCRAP/bin/SCRAP_isomir.sh \
    -d DATASETS/human/Helwak_GSE50452 \
    -a SCRAP/adapters/CLASH_Human/CLASH_Human_Adapters.txt \
    -p no \
    -f no \
    -r SCRAP/ \
    -m hsa \
    -g hg38

#Peak Calling
bash SCRAP/bin/Peak_Calling.sh \
	-d DATASETS/human/Helwak_GSE50452 \
	-a SCRAP/adapters/CLASH_Human/CLASH_Human_Adapters.txt \
	-c 1 \
	-l 1 \
	-f no \
	-r SCRAP/ \
	-m hsa \
	-g hg38  

#peak annotation
bash SCRAP/bin/Peak_Annotation.sh \
    -p DATASETS/human/Helwak_GSE50452/peaks.bed \
    -r SCRAP/ \
    -s human 