This is the GitHub repository containing code related to the paper [*Modular abelian surfaces of small conductor with nontrivial Tate-Shafarevich groups*](https://arxiv.org) by Sam Frengley and Dylan Laird [arXiv:??](https://arxiv.org).

Most of the code in this repository is written in either [SageMath](https://www.sagemath.org) or [Magma](http://magma.maths.usyd.edu.au/magma/). We have tested on Magma version 2.27-8 and SageMath version 9.5 (Python 3.10.12).

## Structure
- The directory `./verifications/` contains code sufficient to prove the claims in the paper (with the exception of Theorem 1.8).
- The directory `./congruence_search/` contains the code to:
  - Prove Theorem 1.8.
  - Search for the examples found in Theorem 1.4.
- The directory `./data/` contains the data from the article. This includes the data in Theorem 1.8 and Tables B.1, B.2, and B.3 in a machine readable format.

## Dependencies
The following capabilities should be installed.

### SageMath / Python
- [sage-cluster-pictures](https://alexjbest.github.io/cluster-pictures/index.html#how-to-install)
- [LMFDB-lite](https://github.com/roed314/lmfdb-lite)

### Magma
- [Genus2Conductor](https://github.com/cjdoris/Genus2Conductor.git) : this is included as a submodule in `./libs/` however we encourage you to install through [CHIMP](github.com/edgarcosta/CHIMP).
- [genus-2-RM](github.com/SamFrengley/genus-2-RM) : this is included as a submodule in `./libs/` and is a simple way of extracting the models for Hilbert modular surfaces given by Elkies--Kumar in [this article](https://arxiv.org/abs/1209.3527).
- [sha-7-examples](github.com/SamFrengley/sha-7-examples) : this is included as a submodule in `./libs/` and is only needed for verifying the claims in Section 5.
