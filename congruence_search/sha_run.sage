"""
sha_run.sage
============
This is a script to compute *potential* examples of III[p]

When calling this script you should specify
  --d flag : Where the directory where you output your congruence searching
     (to allow running on machines with small ~/. storage but some large
     shared storage area).

Run as for example
> sage sha_run.sage --d "../data/."

"""

load("congruence_search.sage")
load("sha_helpers.sage")
import argparse

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

WEIER = load_weier(TARGET_DIR + "label_to_curve.json")

def load_congs(p, dim_bd=4, target_dir="../data/"):
    """Load the precomputed congruences     
    """
    pattern = f"congruences/combined/{p}.json"
    with open(target_dir + pattern, "r") as f:
        data = json.load(f)
    data = [[x for x in cong if x['dim'] <= dim_bd] for cong in data]
    data = [cong for cong in data if len(cong) >= 2]
    return data


if __name__ == "__main__":
    lset = P_SETS[0]
    for p in SMALL_PRIMES:
        congs = load_congs(p, dim_bd=4)
        hs = []
        for cong in congs:
            flag, h = verify_cong_for_primes(cong[0], cong[1], p, lset)
            assert flag
            hs.append(h[0])
        
        vis = get_possible_vis_from_congs(congs, hs, lset, p,
                                      dim_bd=2, weier=WEIER)
        vis.sort(key=sort_sha)
        vis_file = TARGET_DIR + f"sha/{p}.json"
        with open(vis_file, "w") as f:
            json.dump(out_vis(vis, p), f, indent=2)
