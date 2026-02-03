"""
sha_helpers.sage
================

Functions for searching for examples of III[p] using a precomputed list of
congruences.
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


def load_weier(filename="../data/label_to_curve.json") -> dict:
    """Loads the Weierstrass equations for the data file.

    Returns
    -------
    dict
        Keys are modular form labels and values are Weierstrass equations
    """
    with open(filename, "r") as f:
        data = json.load(f)
    return data


def mult_primes(N) -> list:
    """Returns those prime exactly dividing N
    """
    facts = list(Integer(N).factor())
    return [p[0] for p in facts if p[1] == 1]


def pdiv(N) -> list:
    """Returns the prime divisors of N
    """
    facts = list(Integer(N).factor())
    return [p[0] for p in facts]


def magma_get_component_gp_regmod(C: list, p: int) -> int:
    """Computes the order of the geometric component group using Magma's
    intrinsic `RegularModel`

    Parameters
    ----------
    C : list
        A pair of lists such that in the Magma syntax
        `C := HyperellipticCurve(Polynomial(C[0]), Polynomial(C[1]))`
    p : int
        A prime number

    Returns
    -------
    int
        Either the order of the geometric component group or -1 if failure.
    """
    magma_str = f"C := {C};"
    magma_str += f"p := {p};"
    with open("comp_gp_regmod.m", "r") as f:
        magma_str += f.read()
    return eval(magma.eval(magma_str))


def magma_get_component_gp_kohelstein(fm: dict, ps: list):
    """Computes the order of the geometric component group using Kohel--Stein's
    Magma intrinsics

    Parameters
    ----------
    fm : dict
        Newform
    ps : list
        A prime number

    Returns
    -------
    int
        Either the order of the geometric component group or -1 if failure.
    """
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
    """Gets integers from the type of group output by genus2reduction
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
    """Computes the order of the geometric component group using Liu's
    `genus2reduction`

    Parameters
    ----------
    C : list
        A pair of lists such that in the Magma syntax
        `C := HyperellipticCurve(Polynomial(C[0]), Polynomial(C[1]))`
    p : int
        A prime number

    Returns
    -------
    int
        Either the order of the geometric component group or -1 if failure.
    """
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
    """Computes the order of the geometric component group of a stable model
    of C.

    Parameters
    ----------
    C : list
        A pair of lists such that in the Magma syntax
        `C := HyperellipticCurve(Polynomial(C[0]), Polynomial(C[1]))`
    p : int
        An odd prime number

    Returns
    -------
    int
        The order of the geometric component group or -1 if failure.
    """
    assert p != 2
    R.<x> = PolynomialRing(Qp(p, 1000))
    f = R(C[0])
    g = R(C[1])
    C = HyperellipticCurve(f + 1/4*g^2)
    cluster = Cluster.from_curve(C)
    return cluster.tamagawa_number()
    

def good_red_fil(pair, p):
    """Returns true if and only if both froms have good reduction at p
    """
    return GCD(pair[0]['level'] * pair[1]['level'], p) == 1


def rank_fil(pair):
    """Returns true if there is a rank discrepancy of >= 2
    """
    if pair[0]['analytic_rank'] + 1 < pair[1]['analytic_rank']:
        return True, pair
    elif pair[0]['analytic_rank'] > pair[1]['analytic_rank'] + 1:
        return True, [pair[1], pair[0]]
    else:
        return False, pair


def ell_curve(form: dict):
    """Returns the optimal elliptic curve corresponding to `form`
    """
    N = form['level']
    Es = cremona_curves([N])
    app = [Integer(a[0]) for a in form['ap']][:len(PRIMES)]
    isog_class = []
    for E in Es:
        [E.ap(p) for p in PRIMES][:5]
        if app == [E.ap(p) for p in PRIMES]:
            isog_class.append(E)
    mds = [E.modular_degree() for E in isog_class]
    opt = min(mds)
    i = mds.index(opt)
    return isog_class[i]   


def make_subfield(l : list):
    """Makes a subfield of the parent of the elements of l which is the smallest
    subfield containing all elements of l
    """
    d = [x.minimal_polynomial().degree() for x in l]
    d = LCM(d)
    K = parent(l[0]).subfield(d)
    return K

        
def is_irreducible_from_hash_mod_p(my_hash: list, lset: list) -> bool:
    """Given the hash and the set of ls check whether the hash corresponds to an
    irreducible Galois representation.

    Parameters
    ----------
    my_hash : list
        A list of finite field elements
    lset : list
        A list of prime numbers :math:`\mathcal{L}`. Note that here the
        `lset` is the indices of `\mathcal{L}` in the primes `PRIMES`

    Returns
    -------
    bool
        True only if the hash corresponds to an irreducible Galois representation
        but it may fail to detect an irreducible Galois rep (not enough Frobenii
        checked)
    """
    F = make_subfield(my_hash)
    P.<x> = PolynomialRing(F)
    for i in range(len(my_hash)):
        f = x^2 - F(my_hash[i])*x + lset[i]                 # local L-factor
        if f.is_irreducible():
            return True
    return False    
        

def is_tamagawa_prod_prime_to_p(
        pair: list,
        weier: dict,
        magma_timeout=100
) -> (bool, bool):
    """Returns true if the Tamagawa numbers of pair[0] and pair[1] are
    prime to p
    
    Given a pair of forms (if the dimension is 2 we require to also know a
    Weierstrass equation for a corresponding genus 2 curve) return a pair of
    bools specified as follows:
    - First bool is false if we proved (the order of the geometric component
    group) is not coprime to p.
    - If the first bool is true, the second bool is true if and only if we
    proved that the product of tamagawa numbers is coprime to p.

    Parameters
    ----------
    pair : list
        A pair of modular forms
    weier : dict
        A dictionary specifying Weierstrass equations in dimension 2 from the
        corresponding form label
    magma_timeout : int, default=100
        When using magma timeout the process after this many seconds, then kill
        all magma processes

    Returns
    -------
    (bool, bool)
        Specified above

    Note
    ----
    In dimension 1 we use the sage implementation of Tate's algorithm. In
    dimension 2 we use the algorithms described in Section 5, Filter 4.
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
                try:                                        # The try Kohel--Stein
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
        lset: list,
        p: int,
        dim_bd=2,
        weier={},
        magma_timeout=100
) -> list:
    """Applies the filters in Section 4.1 to give a list of possible III[p]

    Applies the filters from Section 4.1. The output is a list of potential
    elements of III[p]. Note that if we fail to compute a Tamagawa number, we
    still keep the example, so you should not view these examples are proved.

    Parameters
    ----------
    congs : list
        A list of p-congruences
    hashes : list
        The corresponding list of hashes
    lset : list
        A list of prime numbers :math:`\mathcal{L}`. Note that here the
        `lset` is the indices of `\mathcal{L}` in the primes `PRIMES`
    p : int
        A prime number
    dim_bd : int, default=2
        The upper bound on the dimension we search for
    weier : dict, default={}
        A dictionary assocaiting modular form labels to Weierstrass equations
    magma_timeout : int, default=100
        When using magma timeout the process after this many seconds, then kill
        all magma processes

    Returns
    -------
    list
       Possible pairs realising III[p] in the first.
    """
    # irreducible mod p reps
    irred_congs = []
    for i in range(len(hashes)):
        if is_irreducible_from_hash_mod_p(hashes[i], lset):
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


def sort_sha(form: dict) -> int:
    """Gets the level

    Parameters
    ----------
    form : dict
        Modular form

    Returns
    -------
    int
        Level of `form`
    """
    return form[0]['level']


def get_disc(form: dict) -> int:
    """Gets the discriminant of the coefficient field of `form`

    Parameters
    ----------
    form : dict
        Modular form

    Returns
    -------
    int
        The discriminant of K_f
    """
    R.<u> = PolynomialRing(QQ)
    ply = R(form['field_poly'])
    d = ply.disc()
    if d % 4 == 0:
        return d / 4
    else:
        return d

    
def get_cp_prod_genus2_mult(form: dict) -> (int, list, list):
    """Computes the product of the orders of the geometric component groups
    of a Jacobian assocated to `form` at the primes exactly dividing the level

    Parameters
    ----------
    form : dict
        Newform

    Returns
    -------
    int
        The order of the geometric component group
    list
        A list of primes where we failed to compute the component group order
    list
        A list of the computed component group orders
    """
    C = WEIER[form['label']]
    ret = 1
    fail = []
    cells = []
    # use that additive primes have c_ell <= p
    for ell in mult_primes(form['level']):
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
                    print(f"{form['label']}, {ell}, {cell}")
                    if cell == -1:
                        fail.append(ell)
                    else:
                        ret *= cell
                    cells.append(cell)
                except:
                    print(f"Couldn't compute c_{ell}({form['label']})")
                    fail.append(ell)
                    cells.append(-1)
    return ret, fail, cells


def get_cp_prod_genus2_all(form: dict) -> (int, list, list):
    """Computes the product of the orders of the geometric component groups
    of a Jacobian assocated to `form` at the bad primes

    Parameters
    ----------
    form : dict
        Newform

    Returns
    -------
    int
        The order of the geometric component group
    list
        A list of primes where we failed to compute the component group order
    list
        A list of the computed component group orders
    """
    C = WEIER[form['label']]
    ret = 1
    fail = []
    cells = []
    mps = mult_primes(form['level'])
    for ell in pdiv(form['level']):
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
                    print(f"{form['label']}, {ell}, {cell}")
                    if cell == -1:
                        fail.append(ell)
                    else:
                        ret *= cell
                    cells.append(cell)
                except:
                    print(f"Couldn't compute c_{ell}({form['label']})")
                    fail.append(ell)
                    cells.append(-1)
    return ret, fail, cells


def out_vis(vis: list, p: int) -> list:
    """Augments a list of possible pairs of III[p] examples with the keys
    'c_ell' in each of the forms (this is the component group orders at the bad
    primes)
    """
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
