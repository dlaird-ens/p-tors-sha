# To run in the sage shell with the LMFDB API in order to retrieve all the
# modular forms of w/ 1 <= N <= 10000

from lmf import db
from extract_data import *
from collections import defaultdict
import json
import os
import time

def partial_merge(merged, data):
    for key, values in data.items():
        if key not in merged:
            merged[key] = {d['label']: d for d in values}
        else:
            for d in values:
                merged[key][d['label']] = d
    

def merge_local():
    merged = {}

    path = '../nfs_dim_bd_4.json'
    print(path)
    with open(path, 'r') as f:
        data = json.load(f)
    partial_merge(merged, data)
    
    for filename in os.listdir('../local'):
        if filename.endswith('.json'):
            path = os.path.join('../local', filename)
            print(path)
            with open(path, 'r') as f:
                data = json.load(f)
        partial_merge(merged, data)
        
    # Convert back to lists
    merged = {k: list(v.values()) for k, v in merged.items()}
    with open('../nfs_dim_bd_4.json', 'w') as f:
        json.dump(merged, f, indent=2)

        
def main(min_level, max_level):
    os.makedirs("../local/", exist_ok=True)
    table = []
    start = time.time()
    for n in range(min_level, max_level + 1) : 
        group = get_forms(n, 2, 1, deg_bd=4)
        table.append(group)
        if n % 10 == 0: 
            print(f'Done with level {n} : {time.time() - start:.2f}')

        if n % 50 == 0 or n == max_level:
            merged = defaultdict(list)
            for level_group in table:
                for pol, forms in level_group.items():
                    merged[pol].extend(forms)
                    merged_dict = dict(merged)

            filename = f"../local/data_{min_level}_{max_level}.json"
            with open(filename, "w") as f:
                json.dump(merged_dict, f, indent=2)

                
if __name__ == "__main__":
    main(1, 10000)
