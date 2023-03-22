import itertools as it
from pathlib import Path
import random
from datetime import date
import os
import csv
today = date.today()
d1 = today.strftime("%d-%m-%Y")


sample_pair_list = ["CD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN", "CD-AOM-DSS-Epi_plus_DN,HFD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN",
                    "CD-AOM-DSS-Immune,LFD-AOM-DSS-Immune", "CD-AOM-DSS-Immune,HFD-AOM-DSS-Immune", "HFD-AOM-DSS-Immune,LFD-AOM-DSS-Immune"]

def generate_generic_job_commands(sample_pair):
    
    
    job_id= "cell_type_prop"
    # out_path = f"./jobs/{job_id}"
    result_path = f"../data/out_data/{job_id}"


    # Path(out_path).mkdir(parents=True, exist_ok=True)
    Path(result_path).mkdir(parents=True, exist_ok=True)
    all_jobs_f = open(f"all_runs.sh", "w")
    all_jobs_f.write("#!/bin/sh\n")
    for i in range(1, 101):

        command_line = f"python utils_test.py -sp {sample_pair} > {result_path}/{sample_pair}_{i}.txt"    
        all_jobs_f.write(f"sbatch --job-name={job_id}_{i} --mem=5g -n 1 --time=7-00:00:00 --output=output_{i} \"{command_line}\"\nsleep 1\n")
    
    all_jobs_f.close()

# python utils_test.py -sp CD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN -o ../data/out_data/ > test.txt
generate_generic_job_commands("CD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN")