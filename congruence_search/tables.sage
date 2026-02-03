"""
tables.sage
===========
Utilities for producing tables like the ones in the paper.

"""
load("congruence_search.sage")
load("sha_helpers.sage")

SMALL_PRIMES = list(primes(5, 4000))
WEIER = load_weier("../data/label_to_curve.json")

def load_congs(p, dim_bd=4, target_dir="../data/"):
    """

    """
    pattern = f"congruences/combined/{p}.json"
    with open(target_dir + pattern, "r") as f:
        data = json.load(f)
    data = [[x for x in cong if x['dim'] <= dim_bd] for cong in data]
    data = [cong for cong in data if len(cong) >= 2]
    return data


################################################################################
# TABLE B.1
def congs_table(low_p, dim_bd=2, target_dir="../data/"):
    """Outputs in latex format
    """
    ret = ""
    my_p = [p for p in SMALL_PRIMES if p >= low_p]
    data = [load_congs(p, dim_bd=dim_bd) for p in my_p]
    info = []
    for i, congs in enumerate(data):
        for cong in congs:
            info.append((my_p[i], cong[0]['label'], cong[1]['label']))

    n = (len(info) + 1) // 2
    for i in range(n):
        x = info[i]
        ret += f"{x[0]} & \\LMFDBLabelMF{{{x[1]}}} & \\LMFDBLabelMF{{{x[2]}}} & "
        if n + i != len(info):
            y = info[n + i]
            ret += f"{y[0]} & \\LMFDBLabelMF{{{y[1]}}} & \\LMFDBLabelMF{{{y[2]}}} \\\\ \n "
        else:
            ret += f" & & \\\\ \n "

    return ret


def congs_table_txt(low_p, dim_bd=2, target_dir="../data/"):
    """Machine readable format
    """
    my_p = [p for p in SMALL_PRIMES if p >= low_p]
    data = [load_congs(p, dim_bd=dim_bd) for p in my_p]
    info = []
    with open(target_dir + "tmp.txt", "w") as f:
        for i, congs in enumerate(data):
            for cong in congs:
                tr_1, nm_1 = get_trace_and_norm_form(cong[0])
                tr_2, nm_2 = get_trace_and_norm_form(cong[1])
                row = f"{my_p[i]}:"
                row += f"{cong[0]['label']}:{tr_1[:30]}:{nm_1[:30]}:"
                row += f"{cong[1]['label']}:{tr_2[:30]}:{nm_2[:30]}"
                f.write(row + "\n")               


################################################################################
# CONGRUENCES MISMATCHED
# This is not a table we included, but it seems interesting to find congruences
# where the dimensions are different

def congs_mismatched_table(low_p, dim_a, dim_b, target_dir="../data/"):
    """

    """
    ret = ""
    my_p = [p for p in SMALL_PRIMES if p >= low_p]
    data = [load_congs(p, dim_bd=max(dim_a, dim_b)) for p in my_p]
    info = []
    for i, congs in enumerate(data):
        for cong in congs:
            if cong[0]['dim'] == dim_a and cong[1]['dim'] == dim_b:
                info.append((my_p[i], cong[0]['label'], cong[1]['label']))

    n = (len(info) + 1) // 2
    for i in range(n):
        x = info[i]
        ret += f"{x[0]} & \\LMFDBLabelMF{{{x[1]}}} & \\LMFDBLabelMF{{{x[2]}}} & "
        if n + i != len(info):
            y = info[n + i]
            ret += f"{y[0]} & \\LMFDBLabelMF{{{y[1]}}} & \\LMFDBLabelMF{{{y[2]}}} \\\\ \n "
        else:
            ret += f" & & \\\\ \n "

    return ret


################################################################################
# WEIRSTASS TABLE
def make_weier_table(l):
    """Given a list of forms whose associated jacobians are known, output the
    equivalent of Table B.3
    """
    R.<x> = PolynomialRing(QQ)
    ret = ""
    for f in sorted(l):
        C = [R(p) for p in WEIER[f]]
        if C[1] == 0:
            w_eq = f"y^2 = {latex(C[0])}"
        elif C[1].coefficients()[-1] < 0:
            w_eq = f"y^2 - ({latex(-C[1])})y = {latex(C[0])}"
        else:
            w_eq = f"y^2 + ({latex(C[1])})y = {latex(C[0])}"
        ret += f"\\LMFDBLabelMF{{{f}}} & ${w_eq}$ \\\\\n"
    return ret


################################################################################
# SHA TABLE
def print_row(pair, p):
    """
    Prints a row of Table B.2
    """
    ret = f"${p}$ & \\LMFDBLabelMF{{{pair[0]['label']}}} & \\LMFDBLabelMF{{{pair[1]['label']}}} & ${1}$ "
    if pair[0]['dim'] == 2:
        ret += f"& $\\QQ(\\sqrt{{{get_disc(pair[0])}}})$ "
    elif pair[0]['dim'] == 1:
        ret += "& $\\QQ$ "
    else:
        ret += f"& {pair[0]['field_poly']} "
    ret += f"& ${pair[0]['dim'] * pair[0]['analytic_rank']}$ "
    if pair[0]['dim'] == 1:
        ret += f"& ${ell_curve(pair[0]).tamagawa_product()}$ "
    elif pair[0]['dim'] == 2:
        if p > 2 * pair[0]['dim'] + 1:
            cprd, fail = get_cp_prod_genus2_mult(pair[0])
            if len(fail) == 0:
                ret += f"& $\\mathcal{{P}} \\cdot {cprd}$ "
            else:
                ret += "& ??"
        else:
            cprd, fail = get_cp_prod_genus2_all(pair[0])
            if len(fail) == 0:
                ret += f"& ${cprd}$"
            else:
                ret += "& ??"
    else:
        ret += "& ??"
    if pair[1]['dim'] == 2:
        ret += f"& $\\QQ(\\sqrt{{{get_disc(pair[1])}}})$ "
    elif pair[1]['dim'] == 1:
        ret += "& $\\QQ$ "
    else:
        ret += f"& {pair[1]['field_poly']} "
    ret += f"& ${pair[1]['dim'] * pair[1]['analytic_rank']}$ "
    if pair[1]['dim'] == 1:
        ret += f"& ${ell_curve(pair[1]).tamagawa_product()}$ "
    elif pair[1]['dim'] == 2:
        if p > 2 * pair[1]['dim'] + 1:
            cprd, fail = get_cp_prod_genus2_mult(pair[1])
            if len(fail) == 0:
                ret += f"& $\\mathcal{{P}} \\cdot {cprd}$ "
            else:
                ret += "& ??"
        else:
            cprd, fail = get_cp_prod_genus2_all(pair[1])
            if len(fail) == 0:
                ret += f"& ${cprd}$"
            else:
                ret += "& ??"
            
    else:
        ret += "& ??"

    return ret + " \\" + "\\"


def print_txt_row(pair, p):
    """
    Prints a row of Table B.2 in machine readable format
    """
    ret = f"{p}:"
    for fm in pair:
        ret += f"{fm['label']}:"
        if fm['dim'] == 1:
            E = ell_curve(fm)
            ret += f"{list(E.a_invariants())}:"
            ret += "1:"
            ret += f"{fm['analytic_rank']}:"
            ret += f"{E.tamagawa_product()}:"
        else:
            w = WEIER[fm['label']]
            ret += f"{w}:"
            ret += "2:"
            ret += f"{2*fm['analytic_rank']}:"
            if p > 5:
                cprd, fail, _ = get_cp_prod_genus2_mult(fm)
                if len(fail) == 0:
                    ret += f"{cprd}:"
                else:
                    ret += "-1:"
            else:
                cprd, fail, _ = get_cp_prod_genus2_all(fm)
                if len(fail) == 0:
                    ret += f"{cprd}:"
                else:
                    ret += "-1:"
    return ret[:-1]


################################################################################
# Examples of usage

if __name__ == "__main__":
    # Table B.1
    print(congs_table(20, dim_bd=2))
    print("\n\n")

    # Very large p-congruences
    print(congs_table(300, dim_bd=4))
    print("\n\n")
    
    # Dim 2 + Dim 4 congruences with p > 50
    print(congs_mismatched_table(50, 2, 4))
    print("\n\n")

    # # Try to find III[11] and III[13] and then print the row of the table
    # for p in [11,13]:
    #     congs = load_congs(p, dim_bd=2)
    #     lset = P_SETS[0]
    #     hs = []
    #     for cong in congs:
    #         flag, h = verify_cong_for_primes(cong[0], cong[1], p, lset)
    #         assert flag
    #         hs.append(h[0])            
        
    #     vis = get_possible_vis_from_congs(congs, hs, lset, p,
    #                                        dim_bd=2, weier=WEIER)
    #     vis.sort(key=sort_sha)
    #     for pair in vis:
    #         print(print_row(pair, p))
    #         print("\n\n")
