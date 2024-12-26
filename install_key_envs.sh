#!/bin/bash -l

#Usage: bash install_key_envs.sh CONDA_ENVS_PATH CPUS

#Make sure you have the latest conda version
eval "$(conda shell.bash hook)"
conda activate base
#conda install conda=23.7.3 -y
#python version in base > 3
#avoid conflicts by setting this:


#R_env
eval "$(conda shell.bash hook)"
conda env create -f conda_envs/R_env.yaml -y
conda activate R_env
numberofcpus=$2
Rscript --vanilla support/Install_Rpackages.R $numberofcpus
conda deactivate

#Install snakemake https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
eval "$(conda shell.bash hook)"
conda create -n phylobotl_env -c conda-forge -c bioconda snakemake=7.25.0 Python=3.11.4 -y
conda activate phylobotl_env
conda config --set channel_priority strict
conda deactivate
#EggNOG environment
conda create -n eggnog_mapper_env -c bioconda -c conda-forge eggnog-mapper=2.0.1 -y
conda activate eggnog_mapper_env

CONDA_ENVS_PATH=$1
#export CONDA_ENVS_PATH=$1

#export PATH=$CONDA_ENVS_PATH/eggnog_mapper_env:$CONDA_ENVS_PATH/eggnog_mapper_env/bin:"$PATH"
mkdir -p $CONDA_ENVS_PATH/eggnog_mapper_env/data
#export EGGNOG_DATA_DIR=$CONDA_ENVS_PATH/eggnog_mapper_env/data

download_eggnog_data.py --data_dir $CONDA_ENVS_PATH/eggnog_mapper_env/data -y
conda deactivate


# Create a conda environment for geNomad
eval "$(conda shell.bash hook)"
conda create -n genomad_env -c conda-forge -c bioconda genomad=1.6.1 -y
# Activate the geNomad environment
conda activate genomad_env

mkdir -p $CONDA_ENVS_PATH/genomad_env/data
genomad download-database $CONDA_ENVS_PATH/genomad_env/data
conda deactivate
