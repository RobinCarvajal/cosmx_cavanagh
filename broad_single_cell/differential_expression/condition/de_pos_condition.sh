#!/bin/bash

# SOFTWARE PATHS
searchlight_dir="/home/robin/Searchlight2"
searchlight_py="${searchlight_dir}/software/Searchlight2.py" # Searchligth script
r_path="/home/robin/miniconda3/envs/DE/bin/Rscript" # R script

# DATA PATHS
analysis_type="broad_single_cell" # type of analysis (e.g. single.cell, niches, celltypes, etc)
main_dir="/mnt/f/cosmx_data/${analysis_type}/differential_expression" # main dir
bg_path="/mnt/f/cosmx_data/de_utils/background/mouse.GRCm39.csv" # background path

de_type="condition" # differential expression variables
position="posterior"
celltypes_dir="${main_dir}/${de_type}/${position}"  # celltypes dir

# celltype list
celltypes_list=(
    #"Astrocytes"
    #"Immune"
	"GABAergic"
	"Glutamatergic"
	"cluster0"
	"Endothelial"
	"cluster1"
	"Neurovascular"
	"Oligodendrocytes"
	"MSNs"
	"Ependymal"   
    )

# Print the array (optional)
echo "${celltypes_list[@]}"

for celltype in "${celltypes_list[@]}"
do
    echo ${celltype}
    out_dir="${celltypes_dir}/${celltype}/SL2_results"
    mkdir -p ${out_dir}
	
    python3 "${searchlight_py}" --r path="${r_path}" --out path="${out_dir}" \
    --bg  file="${bg_path}" \
    --em  file="${celltypes_dir}/${celltype}/em_${celltype}.tsv" \
    --ss  file="${celltypes_dir}/${celltype}/ss_${celltype}.tsv" \
    --de  file="${celltypes_dir}/${celltype}/de_${celltype}.tsv",numerator="pos_A",denominator="pos_C" \
    --ora file="${searchlight_dir}/gene_set_databases/STRING_11.5_mouse.gmt",type=STRING11.5 \
    --ura file="${searchlight_dir}/upstream_regulator_databases/trrust.mouse.tsv",type=TRRUST 

done
