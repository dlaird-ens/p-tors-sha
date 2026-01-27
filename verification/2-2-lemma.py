"""
This script verifies the proof of Lemma 2.2
"""

from math import gcd

P_sets = [
    [467, 541, 683, 691, 761, 773, 811, 829, 853, 863, 929, 937, 947, 953, 991],
    [439, 461, 617, 631, 691, 733, 739, 757, 787, 811, 827, 919, 937, 947, 983],
    [439, 503, 587, 647, 683, 701, 709, 727, 739, 757, 797, 829, 839, 853, 929],
    [419, 439, 541, 569, 587, 641, 709, 727, 751, 769, 773, 821, 827, 859, 971],
    [419, 461, 569, 577, 601, 641, 701, 719, 733, 751, 887, 907, 919, 971, 983],
    [419, 431, 577, 587, 599, 617, 733, 773, 823, 859, 877, 883, 887, 941, 983]
]


def is_coprime_to_all(n, primes):
    """
    Return `True` if and only if n is coprime to all primes in the list
    `primes`
    """
    for p in primes:
        if n % p == 0:
            return False
    return True

def is_pair_covered(a, b):
    """
    Return `True` if and only if there exists a set `P` of primes such that
    both `a` and `b` are coprime to `p` for all `p` in `P`.
    """
    for P in P_sets:
        if is_coprime_to_all(a, P) and is_coprime_to_all(b, P):
            return True
    return False


# Check all unordered pairs (a,b) of levels at most 10000
if __name__ == "__main__":
    covered = True
    for a in range(1, 10001):
        for b in range(a, 10001):
            if not is_pair_covered(a, b):
                print(
                    f"PROOF FAILED: For every set P, one of the pair ({a}, {b})"
                    + "is divisible by some p ∈ P."
                )
                covered = False

    if covered:
        print(
            "As required for all pairs (a,b) ∈ {1,...,10000}² there\n"
            + "exists a set P such that a and b are coprime to all p ∈ P."
        )
