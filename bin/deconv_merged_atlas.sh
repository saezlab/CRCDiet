#!/bin/sh


python vis_deconvolute_part2.py -rp ../data/out_data/atlas_cell_type_annot_light_weight.h5ad -ia ../data/out_data/cell2location_atlas/inf_aver.csv -an vis_merged_deconvolution --num_of_cells 15

# sbatch --job-name=deconv_26-01-2023_1 -p gpu --gres=gpu:1 --mem=20 -n 1 --time=7-00:00:00 --output=deconv_atlas_merged.out "deconv_merged_atlas.sh"