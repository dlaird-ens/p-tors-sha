"""
Author: Sam Frengley and Dylan Laird
Description: A bunch of functions that are useful to import data (modular forms,
elliptic curves ...) from the LMFDB
"""

from collections import defaultdict
from lmf import db

def get_forms(level: int, weight: int, char: int, deg_bd=-1) -> dict:
    """
    Retrives all forms in S_k(Gamma_0(N), chi).

    Retrieves all of the forms of a given weight, level, character with
    coefficients in at most a field of degree deg_bd

    Parameters
    ----------
    level : int
        The level N of the forms.
    weight : int
        The weight k of the forms.
    char : int
        The character order (e.g. 1 for trivial character).
    deg_bd : int, optional
        Max degree of coefficient field. Default -1 means no bound.
    
    Returns
    -------
    dict
        A dictionary whose keys are str of the defining polynomials of the
        coefficient field of the form and whose values are a list of forms
        with N, k, chi and coefficient field degree bounded by deg_bd. 
    """
    N = level
    k = weight
    
    # Importing the required Newforms from the LMFDB
    search_params = {'weight' : k, 'level' : N, 'char_orbit_index' : char}
    output_info = ['label',
                   'level',
                   'field_poly',
                   'ap',
                   'hecke_ring_numerators',
                   'hecke_ring_denominators',
                   'hecke_ring_power_basis'
                   ]
    rawdata = list(db.mf_hecke_nf.search(search_params, output_info))
    if deg_bd == -1:
        data = rawdata
    else:
        data = [rawform for rawform in rawdata
                if len(rawform['field_poly']) - 1 <= deg_bd]

    search_params = {'weight' : k, 'level' : N, 'char_orbit_index' : char}
    output_info = ['label', 'hecke_ring_index', 'analytic_rank', 'dim']
    mydata = list(db.mf_newforms.search(search_params, output_info))
    mylabs = [d['label'] for d in mydata]
    for form in data:
        for info in output_info:
            i = mylabs.index(form['label'])
            form[info] = mydata[i][info]

    #Sort by defining polynomials
    group = defaultdict(list)
    for form in data:
        group[str(form['field_poly'])].append(form)
    return dict(group)


def get_ell_curves(level: int):
    """
    Retrieves all of the elliptic curves over Q with specified level
    """
    N = level
    search_params = {'conductor' : N}
    output_info = ['lmfdb_label', 'lmfdb_iso', 'rank']
    rawdata = list(db.ec_curvedata.search(search_params, output_info))
    sorted_data = defaultdict(list)
    for elem in rawdata: 
        sorted_data[elem['lmfdb_iso']].append(elem)
    return dict(sorted_data)
