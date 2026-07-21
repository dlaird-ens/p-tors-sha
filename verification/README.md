This directory contains the computer assistance needed in proofs of the results (except Theorem 1.8). The specific contributions of each file is outlined below.

## File structure 
Files are named so that files associated to, for example, Theorem 1.4 are named `1-4-theorem.*`.

- `1-4-theorem.m` computes the (algebraic) ranks of the abelian varieties recorded in the Theorem.
- `1-4-theorem.sage` computes the Tamagawa numbers of the abelian varieties recorded in the Theorem (note that `Magma` is called inside this script).
- `1-8-theorem.txt` gives the commands which prove the first part of Theorem 1.8 (which should be run in the directory `../congruence_search/`).
- `1-8-theorem.m` gives the commands which prove the congruences in Table B.1 (by checking up to a Sturm bound).
- `1-8-theorem_1-9-remark.m` checks up to a Sturm bound the congruences claimed in Theorem 1.8 (=Table B.1) and Remark 1.9(2).
- `2-2-lemma.py` is a pure python script which verifies Lemma 2.2.
- `2-4-prop.m` checks the congruences up to a Sturm bound, the irreducibility of the Galois representations, computes the geometric endomorphism algebra, and the conductors of the Jacobians.
- `2-4-prop_proof-data.m` stores points on the Elkies--Kumar model of the Hilbert modular surface to verify Lemma 2.5.
- `2-4-prop_star-cases.m` proves the two starred cases of Proposition 2.4 using the methods of [Fis16, Fre24].
- `2-4-prop_star-cases_proof-data.m` stores some data for defining number fields using in `2-4-prop_star-cases.m`.
- `5-section.m` uses the methods of [Fre24] to prove that the claimed twists of the Klein quartic are correct.
- `5-section_proof-data.m` stores some precomputed data for proving the claims in Section 5.

We also provide a script `_2-4-prop_star-cases_torsion-find.m` which allows one to compute the "division polynomials" for the multiplication by $\mathfrak{p}$ on the Jacobian of a genus 2 curve. Since we needed some nontrivial tricks, we include this for the convenience of future researchers. To run this file you need to clone the repo [samfrengley/essential-tricks](https://github.com/SamFrengley/essential-tricks.git) and attach `Attach("/path/to/essential-tricks/tricks.m")`.

## Dependencies
- To run `1-4-theorem.sage` you should have installed [sage-cluster-pictures](https://alexjbest.github.io/cluster-pictures/index.html#how-to-install) (this also implicitly calls `RegularModel` from `Magma` so you need a working `Magma` installation).
- To run `2-3-prop.m` you should have installed and attached [Genus2Conductor](https://github.com/cjdoris/Genus2Conductor) and [genus-2-RM](github.com/SamFrengley/genus-2-RM). Recommend [CHIMP](https://github.com/edgarcosta/CHIMP) for the former, but we have also included it as a submodule. To use `genus-2-RM` you should make sure you run the `setup.py` script in that repository.
- To run `1-4-theorem.m` and `5-section.m`  one should clone the repository [sha-7-examples](https://github.com/SamFrengley/sha-7-examples.git) (we have included it as a submodule).
