#!/bin/bash
#SBATCH --job-name=jupyter
#SBATCH --time=2-00:00:00
#SBATCH --mem=50GB
#SBATCH --output=/net/data.isilon/ag-saez/bq_arifaioglu/home/Projects/CRCDiet/bin/jupyter.log

conda activate crcdiet

cat /etc/hosts
jupyter lab --ip=0.0.0.0 --port=8888
