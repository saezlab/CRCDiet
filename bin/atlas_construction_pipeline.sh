#!/bin/sh
echo "Merging the samples..."
python sc_merge.py -i ../data/out_data -o ../data/out_data -st atlas -an atlas_merge

echo "Integrating samples by Harmony..."
python sc_integrate.py -i ../data/out_data/atlas_merged.h5ad -o ../data/out_data -st atlas -an atlas_integrate

echo "Clustering..."
python sc_cluster.py -i ../data/out_data/atlas_integrated.h5ad -o ../data/out_data -st atlas  -an atlas_cluster

echo "Calculating DEGs per cluster..."
python sc_cluster_annotate.py -i ../data/out_data/atlas_integrated_clustered.h5ad -o ../data/out_data -an atlas_cluster -st atlas