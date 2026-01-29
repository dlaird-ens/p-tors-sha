"""
congruence_search.sage
======================

"""

import json
from collections import defaultdict
from itertools import combinations

PRIMES = list(primes(1000))
P_SETS = [
    [467, 541, 683, 691, 761, 773, 811, 829, 853, 863, 929, 937, 947, 953, 991],
    [439, 461, 617, 631, 691, 733, 739, 757, 787, 811, 827, 919, 937, 947, 983],
    [439, 503, 587, 647, 683, 701, 709, 727, 739, 757, 797, 829, 839, 853, 929],
    [419, 439, 541, 569, 587, 641, 709, 727, 751, 769, 773, 821, 827, 859, 971],
    [419, 461, 569, 577, 601, 641, 701, 719, 733, 751, 887, 907, 919, 971, 983],
    [419, 431, 577, 587, 599, 617, 733, 773, 823, 859, 877, 883, 887, 941, 983]
]
P_SETS = [[PRIMES.index(p) for p in pp] for pp in P_SETS]


def load_forms(filename="../data/nfs_dim_bd_4.json",
               dim_bd=4,
               min_lev=1,
               max_lev=10000):
    """
    Loads the forms for the file `filename` with the given dimension bound and
    levels specified.
    """
    with open(filename, "r") as f:
        data = json.load(f)

    to_del = [k for k,v in data.items() if v[0]['dim'] > dim_bd];
    for fms in to_del:
        del data[fms]

    if min_lev != 1 or max_lev != 10000:
        to_del = []
        for k in data.keys():
            data[k] = [f for f in data[k] if f['level'] <= max_lev and f['level'] >=min_lev]
            if len(data[k]) == 0:
                to_del.append(k)
        for fms in to_del:
            del data[fms]
            
    return data


def unpack_traces(form, a):
    """
    Unpack the traces of Frobenius. K = Q(a)
    """
    basis = []
    deg = form['dim']
    if form['hecke_ring_power_basis'] :
        for i in range(deg):
            basis.append(a^i)
    else: 
        for i in range(deg):
            b = 0
            for j in range(deg):
                b += (form['hecke_ring_numerators'][i])[j] * a^j
            b = b / form['hecke_ring_denominators'][i]
            basis.append(b)
            
    coeffs = []
    for i in range(len(form['ap'])):
        c = 0
        for j in range(deg):
            c += form['ap'][i][j] * basis[j]
        coeffs.append(c) 
    return coeffs


def get_trace_and_norm_form(form):
    """
    Given a form (in our format, as a dictionary as loaded by `load_forms`)
    returns the traces and norms of the first few a_p as a list.
    """
    R.<u> = PolynomialRing(QQ)
    K.<a> = NumberField(R(form['field_poly']))
    tt = unpack_traces(form, a)
    return ([t.trace() for t in tt], [t.norm() for t in tt])


def prime_decomposition(K, p):
    """
    Get the decomposition of a prime p in a number field K
    """
    I = []
    J = K.ideal(p)
    for factors in J.factor(): 
            I.append((
                factors[0],                                 # \fp
                factors[1],                                 # e(\fp/p)
                factors[0].residue_class_degree()           # f(\fp/p)
            ))
    return I


def get_reductions(fp, F):
    """
    Gets the embeddings of k in F where k is a resiude field
    """
    K = fp.number_field()
    k = K.residue_field(fp)

    embs = Hom(k, F)

    def reductions(alpha):
        return [emb(k(alpha)) for emb in embs]

    return reductions


def hash_values_sigma_h_L(form, fp, reduction_maps, list_of_l):
    """
    Each of the hash values sigma(h_L(f; fp)) where sigma ranges over
    embeddings
    """
    K = fp.number_field()
    a = K.gens()[0]
    a_l_list = unpack_traces(form, a)
    a_l_list = [a_l_list[i] for i in list_of_l]
    h = [reduction_maps(a_l) for a_l in a_l_list]
    h = [[h_ell[i] for h_ell in h] for i in range(len(h[0]))]
    return h


def make_hash_table(forms_by_field, p, list_of_l):
    """
    Implements the hash table of Cremona--Freitas
     
    """
    R.<u> = PolynomialRing(QQ)
    bad = prod([PRIMES[i] for i in list_of_l])

    # Smallest finite field containing all the relevant ones,
    deg_bd = max([
        len(forms_by_field[poly][0]['field_poly']) - 1
        for poly in forms_by_field
    ])
    bar_Fp.<z> = GF(p^LCM([i for i in range(1, deg_bd + 1)]))
    
    # Creating the hash table 
    hash_table = defaultdict(list)
    for poly in forms_by_field:
        K.<a> = NumberField(R(forms_by_field[poly][0]['field_poly']))
        I = prime_decomposition(K, p)
        for fp, e, f in I:
            if e == 1:                                      #only unramified p
                reds = get_reductions(fp, bar_Fp)
                for form in forms_by_field[poly]:
                    if gcd(form['level'], bad) == 1:
                        hashes = hash_values_sigma_h_L(form, fp, reds, list_of_l)
                        for h in hashes:
                            hash_table[tuple(h)].append(form)

    return dict(hash_table)


def verify_cong_for_primes(form1: dict, form2: dict, p: int, list_of_l: list) -> (bool, list):
    """
    Given a pair of forms `form1` and `form2` a prime number `p` and a list of
    prime numbers `list_of_l`, check whether there exists a pair of prime ideals
    pp and qq above p such that the forms are (pp,qq)-congruent for each l in
    the `list_of_l`
    """
    R.<u> = PolynomialRing(QQ)
    deg_bd = max([form1['dim'], form2['dim']])
    good = [
        l for l in list_of_l
        if GCD(l, form1['level'] * form2['level']) == 1
    ]
    bar_Fp.<z> = GF(p^LCM([i for i in range(1, deg_bd + 1)]))
    K1.<a> = NumberField(R(form1['field_poly']))
    K2.<b> = NumberField(R(form2['field_poly']))
    I1 = prime_decomposition(K1, p)
    I2 = prime_decomposition(K2, p)
    h1 = []
    h2 = []
    for fp, _, _ in I1:
        reds = get_reductions(fp, bar_Fp)
        h1 += hash_values_sigma_h_L(form1, fp, reds, list_of_l)
    for fp, _, _ in I2:
        reds = get_reductions(fp, bar_Fp)
        h2 += hash_values_sigma_h_L(form2, fp, reds, list_of_l)

    match_hashes = [h for h in h1 if h in h2]
    flag = len([h for h in h1 if h in h2]) >= 1
    return flag, match_hashes


def list_to_verified(cong: list, p: int, list_of_l: list):
    """
    Given a list of candidate congruences `cong` for a prime `p` return the 
    pairs which are (pp,qq)-congruence for the test list of primes `list_of_l`
    """
    ret = []
    pairs = list(combinations([x for x in cong], 2))
    for pair in pairs:
        f = pair[0]
        g = pair[1]
        flag, _ = verify_cong_for_primes(f, g, p, list_of_l)
        if flag:
            ret.append([f,g])
    return ret


def get_congs_from_hash_table(hash_table: dict, p: int, verify_ish=False) -> list:
    """
    Given a hash table (of the form returned by `make_hash_table` return the
    congruences within. If the flag `verify_ish` is set to true then the
    congruences are then further checked for each good prime l <= 1000
    """
    # Hashes can sometimes contain twice the same form (for different primes) we
    # don't deal with this case, but it wouldn't be hard to incorporate
    mult_hs = [h for h in hash_table if len(set([x['label'] for x in hash_table[h]])) >= 2]
    # This next loop is doing
    # congs = [list(set(hash_table[h])) for h in mult_hs]
    # but dicts are not hashable (so set breaks)
    congs = []
    for h in mult_hs:
        forms = []
        for form in hash_table[h]:
            if not form['label'] in [f['label'] for f in forms]:
                forms.append(form)
        congs.append(forms)
    print(mult_hs, congs)
    ret = []
    ret_hs = []
    if verify_ish:
        for i in range(len(congs)):
            cong = congs[i]
            these_verified = [
                cc for cc in list_to_verified(
                    cong, p, [i for i in range(len(PRIMES))]
                ) if len(cc) > 1]
            ret += these_verified
            ret_hs += [mult_hs[i] for _ in range(len(these_verified))]
    else:
        ret = congs
        ret_hs = mult_hs
    return ret, ret_hs
