"""
sha_helpers.sage
================

"""

from sage_cluster_pictures.cluster_pictures import Cluster
import subprocess
import signal

################################################################################
# timeout handling for when Magma takes ages to compute component group, we just
# kill it
def timeout_handler(signum, frame):
    subprocess.run(["killall", "magma.exe"])
    signal.alarm(0)
    time.sleep(1)
    raise SystemExit

signal.signal(signal.SIGALRM, timeout_handler)
################################################################################


def load_weier(filename="../data/label_to_curve.json"):
    """

    """
    with open(filename, "r") as f:
        data = json.load(f)
    return data


def mult_primes(N):
    
    facts = list(Integer(N).factor())
    return [p[0] for p in facts if p[1] == 1]


def pdiv(N):
    
    facts = list(Integer(N).factor())
    return [p[0] for p in facts]


def magma_get_component_gp_regmod(C, p):
    
    magma_str = f"C := {C};"
    magma_str += f"p := {p};"
    with open("comp_gp_regmod.m", "r") as f:
        magma_str += f.read()
    return eval(magma.eval(magma_str))


def magma_get_component_gp_kohelstein(fm, ps):
    
    tr, nm = get_trace_and_norm_form(fm)
    magma_str = f"level := {fm['level']};"
    magma_str += f"dim := {fm['dim']};"
    magma_str += f"field_poly := {fm['field_poly']};"
    magma_str += f"ps := {ps};"
    magma_str += f"my_traces := {tr};"
    magma_str += f"my_norms := {nm};"
    with open("comp_gp_kohelstein.m", "r") as f:
        magma_str += f.read()
    return eval(magma.eval(magma_str))


def inverse_divisors_to_string(s):
    """
    Gets integers from the type of group output by genus2reduction
    """
    info = s.split("x")
    info = [i.split("^") for i in info]
    info = [[int(j.strip("()")) for j in i] for i in info]
    ret = []
    for i in info:
        if len(i) == 1:
            ret.append([i[0], 1])
        else:
            ret.append(i)
    return ret


def g2red_component_gp(C, p):
    R.<x> = PolynomialRing(QQ)
    red = genus2reduction(R(C[1]), R(C[0]))
    red = red.local_data
    if p != 2:
        line = (red[p]).split("\n")
        line = line[1]
        line = line.split(",")
        cmp_grp = line[1].strip()
        cmp_grp = inverse_divisors_to_string(cmp_grp)
        return prod([i[0] ** i[1] for i in cmp_grp])
    else:
        return -1


def clusters_component_gp(C, p):
    R.<x> = PolynomialRing(Qp(p, 1000))
    f = R(C[0])
    g = R(C[1])
    C = HyperellipticCurve(f + 1/4*g^2)
    cluster = Cluster.from_curve(C)
    return cluster.tamagawa_number()
    

def good_red_fil(pair, p):
    return GCD(pair[0]['level'] * pair[1]['level'], p) == 1


def rank_fil(pair):
    if pair[0]['analytic_rank'] + 1 < pair[1]['analytic_rank']:
        return True, pair
    elif pair[0]['analytic_rank'] > pair[1]['analytic_rank'] + 1:
        return True, [pair[1], pair[0]]
    else:
        return False, pair


def ell_curve(f):
    N = f['level']
    Es = cremona_curves([N])
    app = [Integer(a[0]) for a in f['ap']][:len(PRIMES)]
    isog_class = []
    for E in Es:
        [E.ap(p) for p in PRIMES][:5]
        if app == [E.ap(p) for p in PRIMES]:
            isog_class.append(E)
    mds = [E.modular_degree() for E in isog_class]
    opt = min(mds)
    i = mds.index(opt)
    return isog_class[i]   


def make_subfield(l):
    """
    Makes a subfield of the parent of the elements of l which is the smallest
    subfield containing all elements of l
    """
    d = [x.minimal_polynomial().degree() for x in l]
    d = LCM(d)
    K = parent(l[0]).subfield(d)
    return K

        
def is_irreducible_from_hash_mod_p(my_hash, pset):
    F = make_subfield(my_hash)
    P.<x> = PolynomialRing(F)
    for i in range(len(my_hash)):
        f = x^2 - F(my_hash[i])*x + pset[i]                 # local L-factor
        if f.is_irreducible():
            return True
    return False    
        

def is_tamagawa_prod_prime_to_p(
        pair: list,
        weier,
        magma_timeout=100) -> (bool, bool):
    """
    First bool is false if we proved (the order of the geometric component
    group) is not coprime to p. If the first bool is true, the second bool is
    true if and only if we proved that the product of tamagawa numbers is
    coprime to p.
    """
    # Reorder because as dimension increases computations get harder
    if pair[0]['dim'] > pair[1]['dim']:
        pair = [pair[1], pair[0]]
        
    for f in pair:
        # start with bigger primes, it's faster
        mult_ps = mult_primes(f['level'])
        mult_ps.sort()
        mult_ps.reverse()
        add_ps = [ell for ell in pdiv(f['level']) if not ell in mult_ps]
        add_ps.sort()
        add_ps.reverse()
        if f['dim'] == 1:
            E = ell_curve(f)
            if GCD(E.tamagawa_product(), p) != 1:
                return False, False
                
        elif len(weier.keys()) != 0 and f['dim'] == 2:
            C = weier[f['label']]
            for ell in mult_ps:
                if ell != 2:
                    cell = clusters_component_gp(C, ell)
                    if GCD(cell, p) != 1:
                        return False, False
                else:
                    try:                                    # First try Liu
                        cell = g2red_component_gp(C, ell)
                        if cell == -1:
                            raise ValueError("genus2reduction failed")
                        if GCD(cell, p) != 1:
                            return False, False
                    except:
                        try:                                # Then try RegularModel
                            signal.alarm(magma_timeout)
                            cell = magma_get_component_gp_regmod(C, ell)
                            signal.alarm(0)
                            if cell == -1:
                                raise ValueError("RegularModel failed")
                            if GCD(cell, p) != 1:
                                return False, False

                        except:
                            try:                            # The try Kohel--Stein
                                signal.alarm(magma_timeout)
                                cell = magma_get_component_gp_kohelstein(f, [ell])[0]
                                signal.alarm(0)
                                if GCD(cell, p) != 1:
                                    return False, False
                            except:
                                print(f"Couldn't compute c_{ell}({f['label']})")
                                return True, False
                                
            if p <= 2*f['dim'] + 1:
                for ell in add_ps:
                    try:                                    # First try Liu
                        cell = g2red_component_gp(C, ell)
                        if cell == -1:
                            raise ValueError("genus2reduction failed")
                        if GCD(cell, p) != 1:
                            return False, False
                    except:
                        try:                                # Then try RegularModel
                            signal.alarm(magma_timeout)
                            cell = magma_get_component_gp_regmod(C, ell)
                            signal.alarm(0)
                            if cell == -1:
                                raise ValueError("RegularModel failed")
                            if GCD(cell, p) != 1:
                                return False, False
                        except:
                            print(f"Couldn't compute c_{ell}({f['label']})")
                            return True, False

        else:                                               # Dimension >= 3
            if (p <= 2*f['dim'] + 1) and (len(add_ps) != 0):
                return True, False
            for ell in mult_ps:
                try:                                            # The try Kohel--Stein
                    signal.alarm(magma_timeout)
                    cell = magma_get_component_gp_kohelstein(f, [ell])[0]
                    signal.alarm(0)
                    if GCD(cell, p) != 1:
                        return False, False
                except:
                    print(f"Couldn't compute c_{ell}({f['label']})")
                    return True, False
            
    return True, True


def get_possible_vis_from_congs(
        congs: list,
        hashes: list,
        pset: list,
        p: int,
        dim_bd=2,
        weier={},
        magma_timeout=100) -> list:
    """
    Filters applied
    """
    # irreducible mod p reps
    irred_congs = []
    for i in range(len(hashes)):
        if is_irreducible_from_hash_mod_p(hashes[i], pset):
            irred_congs.append(congs[i])
    # make pairs
    maybe = []
    for cong in irred_congs:
        pairs = list(combinations(cong, 2))
        maybe += pairs
    # good reduction
    new = []
    for pair in maybe:
        if good_red_fil(pair, p):
            new.append(pair)
    maybe = new
    # rank discrepancy
    new = []
    for pair in maybe:
        flag, new_pair = rank_fil(pair)
        if flag:
           new.append(new_pair)
    maybe = new
    print([[c['dim'] for c in eg] for eg in maybe])
    # dimension bound
    new = []
    for pair in maybe:
        if pair[0]['dim'] <= dim_bd:
            new.append(pair)
    maybe = new
    # associated curve exists
    new = []
    for pair in maybe:
        if len(weier.keys()) != 0 and pair[0]['dim'] == 2:
            if pair[0]['label'] in weier:
                new.append(pair)
        else:
            new.append(pair)
    maybe = new
    # Tamagawa product
    new = []
    for pair in maybe:
        keep, _ = is_tamagawa_prod_prime_to_p(
            pair,
            weier,
            magma_timeout=magma_timeout
        )
        if keep:
            new.append(pair)
    maybe = new
    return maybe


def sort_sha(pair1):
    return pair1[0]['level']


def get_disc(f):
    R.<u> = PolynomialRing(QQ)
    ply = R(f['field_poly'])
    d = ply.disc()
    if d % 4 == 0:
        return d / 4
    else:
        return d

    
def get_cp_prod_genus2_mult(f):
    C = WEIER[f['label']]
    ret = 1
    fail = []
    cells = []
    # use that additive primes have c_ell <= p
    for ell in mult_primes(f['level']):
        if ell != 2:
            cell = clusters_component_gp(C, ell)
            ret *= cell
            cells.append(cell)
        else:
            try:
                signal.alarm(60)
                cell = magma_get_component_gp_regmod(C, ell)
                if cell == -1:
                    raise Exception("Magma failure")
                else:
                    ret *= cell
                    cells.append(cell)
            except:
                try:
                    cell = g2red_component_gp(C, ell)
                    print(f"{f['label']}, {ell}, {cell}")
                    if cell == -1:
                        fail.append(ell)
                    else:
                        ret *= cell
                    cells.append(cell)
                except:
                    print(f"Couldn't compute c_{ell}({f['label']})")
                    fail.append(ell)
                    cells.append(-1)
    return ret, fail, cells


def get_cp_prod_genus2_all(f):
    C = WEIER[f['label']]
    ret = 1
    fail = []
    cells = []
    mps = mult_primes(f['level'])
    for ell in pdiv(f['level']):
        if ell in mps and ell != 2:
            cell = clusters_component_gp(C, ell)
            ret *= cell
            cells.append(cell)
        else:
            try:
                signal.alarm(60)
                cell = magma_get_component_gp_regmod(C, ell)
                if cell == -1:
                    raise Exception("Magma failure")
                else:
                    ret *= cell
                    cells.append(cell)
            except:
                try:
                    cell = g2red_component_gp(C, ell)
                    print(f"{f['label']}, {ell}, {cell}")
                    if cell == -1:
                        fail.append(ell)
                    else:
                        ret *= cell
                    cells.append(cell)
                except:
                    print(f"Couldn't compute c_{ell}({f['label']})")
                    fail.append(ell)
                    cells.append(-1)
    return ret, fail, cells


def out_vis(vis, p):
    out_list = []
    for pair in vis:
        p0 = pair[0]
        p1 = pair[1]
        if p0['dim'] == 1:
            p0['c_ell'] = [int(x) for x in ell_curve(p0).tamagawa_numbers()]
        elif p0['dim'] == 2:
            if p > 2 * p0['dim'] + 1:
                _, _, cells = get_cp_prod_genus2_mult(p0)
            else:
                _, _, cells = get_cp_prod_genus2_all(p0)
            p0['c_ell'] = [int(x) for x in cells]
            
        if p1['dim'] == 1:
            p1['c_ell'] = [int(x) for x in ell_curve(p1).tamagawa_numbers()]
        elif p1['dim'] == 2:
            if p > 2 * p1['dim'] + 1:
                _, _, cells = get_cp_prod_genus2_mult(p1)
            else:
                _, _, cells = get_cp_prod_genus2_all(p1)
            p1['c_ell'] = [int(x) for x in cells]
            
        out_list.append([p0, p1])
    return out_list
