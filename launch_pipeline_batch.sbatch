#!/bin/bash

#SBATCH --job-name=glopp
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --time=72:00:00
#SBATCH --mem=190000000k
#SBATCH --output=mainlogs/%x-%j.out
#SBATCH --error=mainlogs/%x-%j.err

# Setup anaconda3
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/users/astar/gis/ericejs/scratch/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/users/astar/gis/ericejs/scratch/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/users/astar/gis/ericejs/scratch/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/users/astar/gis/ericejs/scratch/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

# Activate environment
conda activate snakemake

# Run 

name="zb_D6331"
total=200

for (( i = 1; i <= $total; i++ ))
do
    snakemake -s single_species_synthetic_gpresets.snk --configfile configs/$name.json --use-conda --conda-prefix /home/users/astar/gis/ericejs/scratch/projects/glopp/conda/ --cores 96 -d res/$name/ --rerun-incomplete --resources ncbi_load=1 --retries 3 --keep-going -p --rerun-triggers mtime --batch all=$i/$total  --scheduler greedy
    rm -rf res/$name/benchmarks/ res/$name/logs/ res/$name/.snakemake/metadata 
    rm -rf res/$name/phasing/glopp/*/inpref/*/longshot/*/long_reads/
done


snakemake -s single_species_synthetic_gpresets.snk --configfile configs/$name.json --use-conda --conda-prefix /home/users/astar/gis/ericejs/scratch/projects/glopp/conda/ --cores 96 -d res/$name/ --rerun-incomplete --resources ncbi_load=1 --retries 3 --keep-going -p --rerun-triggers mtime --scheduler greedy