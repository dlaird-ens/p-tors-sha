"""
congruence_run.sage
===================
This is a script to compute the hash tables and output duplicated hashes into
the directory specified by the --d flag.

When calling this script you should specify two things
  --d flag : Where the directory is you wish to output the congruence searching into
     (to allow running on machines with small ~/. storage but some large
     shared storage area).
  --l flag : An integer 0 <= l <= 5 specifying the set L_i to be used

Run as for example
> sage congruence_run.sage --d "../data/." --l 0

"""

load("congruence_search.sage")
import json
import os
import argparse

# On a first run you may want to edit this to a small set [23,29] or something
SMALL_PRIMES = list(primes(5, 4000))

################################################################################
# Take in arguments
parser = argparse.ArgumentParser()
# Target directory
parser.add_argument(
    "--d",
    type=str
)
# The index of the set of primes
parser.add_argument(
    "--l",
    type=int
)
ARGS, _ = parser.parse_known_args()
LSET_IND = ARGS.l
if type(ARGS.d) == type(None):
    TARGET_DIR = "../data/"
else:
    if (ARGS.d)[-1] == "/":
        TARGET_DIR = ARGS.d
    else:
        TARGET_DIR = ARGS.d + "/"
################################################################################

def dump_out_congs(p, lset_ind, congs, hs, target_dir="../data/"):
    """Prints the `p`-congruences into a storage file
    
    Prints out `p`-congruences `cong` which we found (for the
    set :math:`\mathcal{L}_{i}` where `i = lset_ind` into a storage file.
    The corresponding hashes for each congruence in `congs` are `hs`.

    Parameters
    ----------
    p : int
        A prime number
    lset_ind : int
        An integer :math:`0 \leq lset_ind \leq 6`
    congs : list
        A list of `p`-congruences
    hs : list
        The corresponding list of hashes for `congs`
    target_dir : str, default="../data/"
        The data directory
    
    Returns
    -------
    None
    
    Note
    ----
    Warning: Congruences are not checked, this is just printing them out.
    """
    os.makedirs(target_dir + f"congruences/{lset_ind}", exist_ok=True)
    os.makedirs(target_dir + f"congruences/{lset_ind}_hash", exist_ok=True)
    cong_file = target_dir + f"congruences/{lset_ind}/{p}.json"
    hash_file = target_dir + f"congruences/{lset_ind}_hash/{p}.json"
    with open(cong_file, "w") as f:
        json.dump(congs, f, indent=2)
    if len(congs) >= 1:
        with open(hash_file, "w") as f:
            mp_z = hs[0][0].parent().modulus().coefficients(sparse=False)
            mp_z = [int(c) for c in mp_z]
            hs = [[hi.polynomial().coefficients(sparse=False) for hi in h] for h in hs]
            hs = [[[int(c) for c in hi] for hi in h] for h in hs]
            json.dump((mp_z, hs), f, indent=2)
    

if __name__ == "__main__":
    print(f"Im the job with the {LSET_IND}-th ell set")

    data = load_forms(filename=TARGET_DIR + "nfs_dim_bd_4.json",
                      dim_bd=4, min_lev=1, max_lev=10000)
    lset = P_SETS[LSET_IND]
    
    os.makedirs(TARGET_DIR + f"sha/{LSET_IND}", exist_ok=True)
    for p in SMALL_PRIMES:
        ht = make_hash_table(data, p, lset)
        congs, hs = get_congs_from_hash_table(ht, p, verify_ish=True)
        print(f"For p={p} we found {len(congs)} congruences\n\n")
        dump_out_congs(p, LSET_IND, congs, hs, target_dir=TARGET_DIR)
