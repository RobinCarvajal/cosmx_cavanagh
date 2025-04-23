#!/bin/bash

# SOFTWARE PATHS
searchlight_py="/home/robin/Searchlight2/software/Searchlight2.py" # Searchligth script
r_path="/home/robin/miniconda3/envs/DE/bin/Rscript" # R script

# DATA PATHS
analysis_type="single.cell" # type of analysis (e.g. single.cell, niches, regions, etc)
main_dir="/mnt/f/cosmx_data/${analysis_type}/differential_expression" # main dir
bg_path="${main_dir}/background/mouse.GRCm39.csv" # background path

de_type="position.condition" # differential expression variables
de_dir="${main_dir}/${de_type}" # folder with de files
ss_dir="${de_dir}/ss/ss_universal.tsv" # sample sheet path 
celltypes_dir="${de_dir}/cell_types" # celltypes dir

# Initialize an empty array
celltypes=()

# Loop through the directories and store their names in the array
for dir in "$celltypes_dir"/*/; do
  if [ -d "$dir" ]; then
    celltypes+=("$(basename "$dir")")
  fi
done

# Print the array (optional)
echo "${celltypes[@]}"

for celltype in "${celltypes[@]}"
do
    echo ${celltype}
    out_dir="${celltypes_dir}/${celltype}/SL2_results"
    mkdir -p ${out_dir}
	
    python3 "${searchlight_py}" --r path="${r_path}" --out path="${out_dir}" \
    --bg  file="${bg_path}" \
    --em  file="${celltypes_dir}/${celltype}/em_${celltype}.tsv" \
    --ss  file="${ss_dir}" \
    --de  file="${celltypes_dir}/${celltype}/de_${celltype}_antC_vs_antA.tsv",numerator="antA",denominator="antC" \
    --de  file="${celltypes_dir}/${celltype}/de_${celltype}_posC_vs_posA.tsv",numerator="posA",denominator="posC" \
    --de  file="${celltypes_dir}/${celltype}/de_${celltype}_antC_vs_posC.tsv",numerator="antC",denominator="posC" \
    --de  file="${celltypes_dir}/${celltype}/de_${celltype}_antA_vs_posA.tsv",numerator="antA",denominator="posA" \
    --ora file="${searchlight_dir}/gene_set_databases/STRING_11.5_mouse.gmt",type=STRING11.5 \
    --ura file="${searchlight_dir}/upstream_regulator_databases/trrust.mouse.tsv",type=TRRUST \
    --mde name=all,numerator="antA"*denominator="antC",numerator="posA"*denominator="posC",numerator="antC"*denominator="posC",numerator="antA"*denominator="posA"

done
