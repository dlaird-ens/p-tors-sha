"""
combine_lset_run.sage
=====================
This is a script to combine the hash table computations for each of the
sets L_i. 

When calling this script you should specify where the directory is you have
output the congruence searching into (to allow running on machines with small
~/. storage but some large shared storage area).

Run as for example
> sage combine_lset_run.sage --d "../data/."

"""
load("congruence_search.sage")
import json
from pathlib import Path
import os
from itertools import combinations
import argparse

# If testing reccommend restricting  this list to some [23,29] or something
SMALL_PRIMES = list(primes(5, 4000))

################################################################################
# Take in arguments
parser = argparse.ArgumentParser()
# Target directory
parser.add_argument(
    "--d",
    type=str
)
ARGS, _ = parser.parse_known_args()
if type(ARGS.d) == type(None):
    TARGET_DIR = "../data/"
else:
    if (ARGS.d)[-1] == "/":
        TARGET_DIR = ARGS.d
    else:
        TARGET_DIR = ARGS.d + "/"
################################################################################


def flatten(data):
    """

    """
    data = [cong for licongs in data for cong in licongs]
    merged = []
    for curr in data:
        curr_set = set([x['label'] for x in curr])
        to_merge = []
        for i, m in enumerate(merged):
            if curr_set.intersection(set([x['label'] for x in m])):
                to_merge.append(i)
        
        if len(to_merge) == 0:                              # haven't seen this congruence
            merged.append(curr)
        else:
            union_set = curr_set
            info = curr
            for i in reversed(to_merge):
                union_set = union_set.union(set([x['label'] for x in merged[i]]))
                info += merged[i]
                merged.pop(i)
            ret = []
            for lab in union_set:
                i = [x['label'] for x in info].index(lab)
                ret.append(info[i])
            merged.append(ret)
            
    return merged


def load_congs(p, dim_bd=4, target_dir="../data/"):
    """

    """
    pattern = f"congruences/*/{p}.json"
    files = list(Path(target_dir).glob(pattern))

    data = []
    for filename in files:
        if not "hash" in str(filename):
            with open(filename, "r") as f:
                data.append(json.load(f))

    data = flatten(data)
    data = [[x for x in cong if x['dim'] <= dim_bd] for cong in data]
    data = [cong for cong in data if len(cong) >= 2]

    #now check them for all stored traces
    tmp = []
    for cong in data:
        cc = list_to_verified(cong, p, range(len(PRIMES)))
        for c in cc:
            tmp.append(c)
    data = tmp

    #now sort them by level
    data = [sorted(cong, key=lambda x: x['label']) for cong in data]
    data.sort(key=lambda x: x[0]['label'])
    return data
    

def dump_out_final_congs(p, congs, target_dir="../data/"):
    """

    """
    os.makedirs(target_dir + f"congruences/combined", exist_ok=True)
    cong_file = target_dir + f"congruences/combined/{p}.json"
    with open(cong_file, "w") as f:
        json.dump(congs, f, indent=2)
    cong_file_short = target_dir + f"congruences/{p}.txt"
    with open(cong_file_short, "w") as f:
        short_congs = [[x['label'] for x in c] for c in congs]
        for c in short_congs:
            print(c, file=f)

        

if __name__ == "__main__":
    # Run the code based on user specifications
    for p in SMALL_PRIMES:
        print(f"Doing p={p}")
        congs = load_congs(p, dim_bd=4, target_dir=TARGET_DIR)
        print(f"For p={p} we found {len(congs)} congruences\n\n")
        dump_out_final_congs(p, congs, target_dir=TARGET_DIR)
