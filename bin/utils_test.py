from utils import calculate_pval
import argparse
import utils
############################### BOOOORIING STUFF BELOW ############################### 
# Read command line and set args
parser = argparse.ArgumentParser(prog='utils test', description='Run utils')
parser.add_argument('-sp', '--samp_pair', help='Sample pair', required=True)
# parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
# parser.add_argument('-of', '--output_file', help='Output file name', required=False)

args = vars(parser.parse_args())
samp_pair = args['samp_pair']
# output_dir = args['output_dir']
# output_file = args['output_file']

###############################

# p_val_list= []

# for samp_pair in sample_pair_list:
for i in range(100):
    
    c_type_list12, samp_prop_dict12, dict_cell_type_pval12 = calculate_pval(samp_pair, 10000)
    for c_type in dict_cell_type_pval12.keys():
        print(f"{samp_pair}\t{c_type}\t{dict_cell_type_pval12[c_type]}")
        # p_val_list.append([samp_pair, c_type, dict_cell_type_pval12[c_type]])

# python utils_test.py -sp CD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN -o ../data/out_data/ > test.txt