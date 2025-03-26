# IMPORTS
import os
from datetime import datetime
import sys

# DEFINE - file paths
input_path = snakemake.input.work_dir 

# LOAD - params
params_dict = {  
    
    # FOR SNAKEfile
    'mito_percent_thresh'   : snakemake.params.mito_percent_thresh,
    'ribo_percent_thresh'   : snakemake.params.ribo_percent_thresh,
    'doublet_thresh'        : snakemake.params.doublet_thresh,
    'min_genes_per_cell'    : snakemake.params.min_genes_per_cell,
    # 'min_peak_counts'       : snakemake.params.min_peak_counts,
    # 'min_num_cell_by_counts': snakemake.params.min_num_cell_by_counts
    
    # # TEST - for local demo purposes
    # 'mito_percent_thresh'   : 1,
    # 'ribo_percent_thresh'   : 2,
    # 'doublet_thresh'        : 3,
    # 'min_genes_per_cell'    : 4,
    # 'min_peak_counts'       : 5,
    # 'min_num_cell_by_counts': 6
}

# RUN

# CHECK
# LOOP - print all params in dict
for key, value in params_dict.items():
    print(f'{key} = {value}')

# LOGGING

# # Get current date and time
# current_time = datetime.now()
# # Format the date and time
# formatted_time = current_time.strftime("%Y-%m-%d_%H:%M:%S")
# print('formatted_time =', formatted_time)
# TEST - get it from snakemake
formatted_time = snakemake.params.formatted_time

# EXPORT

# DEFINE - file paths
output_path = os.path.join(input_path, 'src', 'output', f'my_params-{formatted_time}.txt')
print('output_path =', output_path)

# Create + write to output file
with open(output_path, 'w+') as file:
    print(formatted_time)
    for key, value in params_dict.items():
        file.write(f'{key} = {value}\n')
        
# Send to Snakemake log file:
with open(snakemake.log[0], "w+") as file:
    sys.stderr = sys.stdout = file
    print(formatted_time)
    for key, value in params_dict.items():
        file.write(f'{key} = {value}\n')