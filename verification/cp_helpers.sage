"""
cp_helpers.sage
===============

This module contains functions to help compute component groups of genus 2
Jacobians.

Functions
---------
Selected functions are as follows:
- mult_primes(N)
- pdiv(N)
- magma_get_component_gp_regmod(C, p)
- g2red_component_gp(C, p)
- clusters_component_gp(C, p)

"""
from sage_cluster_pictures.cluster_pictures import Cluster

def mult_primes(N):
    """
    The primes exactly dividing N

    Parameters
    ----------
    N : int or Integer
        An integer.

    Returns
    -------
    list
        The primes l || N
    """
    facts = list(Integer(N).factor())
    return [p[0] for p in facts if p[1] == 1]


def pdiv(N):
    """
    The primes dividing N

    Parameters
    ----------
    N : int or Integer
        An integer.

    Returns
    -------
    list
        The primes l | N
    """
    facts = list(Integer(N).factor())
    return [p[0] for p in facts]


def magma_get_component_gp_regmod(C, p):
    """
    Returns the size of the geometric component group of the Jacobian of C or
    returns -1.

    Parameters
    ----------
    C : list
        A pair consisting of two lists consisting the coefficients of the
        polynomials defining a Weierstrass equation for C
    p : int
        A prime number 

    Returns
    -------
    int 
        The size of the geometric component group of C, or -1 if Magma fails to
        return

    Notes
    -----
    There should be a working Magma installation on your machine.
    """    
    magma_str = f"C := {C};"
    magma_str += f"p := {p};"
    magma_str += "f,g := Explode([Polynomial(c) : c in C]);";
    magma_str += "C := HyperellipticCurve(f,g);"
    magma_str += "C := SimplifiedModel(C);"
    magma_str += "J := Jacobian(C);"
    magma_str += "try print #ComponentGroup(RegularModel(C, p));"
    magma_str += "catch e print -1; end try;"
    return int(magma.eval(magma_str))


def g2red_component_gp(C, p):
    """
    Returns the size of the geometric component group of the Jacobian of C or
    returns -1.

    Parameters
    ----------
    C : list
        A pair consisting of two lists consisting the coefficients of the
        polynomials defining a Weierstrass equation for C
    p : int
        A prime number 

    Returns
    -------
    int 
        The size of the geometric component group of C or -1 genus2reduction
        fails.
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
    """
    Returns the Tamagawa number if C has semistable reduction.

    Parameters
    ----------
    C : list
        A pair consisting of two lists consisting the coefficients of the
        polynomials defining a Weierstrass equation for C
    p : int
        A prime number 

    Returns
    -------
    int 
        The Tamagawa number of C.

    Notes
    -----
    We have not tested this if C does not have semistable reduction, consult
    the documentation of the cluster pictures package.
    """
    R.<x> = PolynomialRing(Qp(p, 1000))
    f = R(C[0])
    g = R(C[1])
    C = HyperellipticCurve(f + 1/4*g^2)
    cluster = Cluster.from_curve(C)
    return cluster.tamagawa_number()
