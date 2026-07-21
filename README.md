This is the GitHub repository containing code related to the paper [*Modular abelian surfaces of small conductor with nontrivial Tate-Shafarevich groups*](https://arxiv.org/abs/2602.19813) by Sam Frengley and Dylan Laird [arXiv:2602.19813](https://arxiv.org/abs/2602.19813).

Most of the code in this repository is written in either [SageMath](https://www.sagemath.org) or [Magma](http://magma.maths.usyd.edu.au/magma/). We have tested on Magma version 2.27-8, SageMath version 9.5, and Python 3.10.12.

## Structure
- The directory `./verifications/` contains code sufficient to prove the claims in the paper (with the exception of Theorem 1.8).
- The directory `./congruence_search/` contains the code to:
  - Prove Theorem 1.8.
  - Search for the examples found in Theorem 1.4.
- The directory `./data/` contains the data from the article. This includes the data in Theorem 1.8 and Tables B.1, B.2, and B.3 in a machine readable format.

## Dependencies
The following capabilities should be installed.

### Zenodo data
The simplest way to get the complete dataset needed to *verify* our calculations is to download it from [Zenodo](https://doi.org/10.5281/zenodo.20921990) and extract it into the `./data/` directory (overwriting some of the contents). For example one can run:
```bash
rm -r ./data/ 
curl -fL -o data_for_p-tors-sha.tar.gz 'https://zenodo.org/records/20921990/files/data_for_p-tors-sha.tar.gz?download=1'
tar -xzf data_for_p-tors-sha.tar.gz
rm data_for_p-tors-sha.tar.gz
mv ./data_for_p-tors-sha/ ./data/
```

### SageMath / Python
- [sage-cluster-pictures](https://alexjbest.github.io/cluster-pictures/index.html#how-to-install) : which can be installed via `sage -pip install git+https://github.com/alexjbest/cluster-pictures.git`
- (OPTIONAL) [LMFDB-lite](https://github.com/roed314/lmfdb-lite) : This is not a strict dependency of the project if the larger data files are downloaded from Zenodo as mentioned above, and is used only to populate the file `./data/nfs_dim_bd_4.json` (using `./data/code/run_extraction.py`).

### Magma
- [Genus2Conductor](https://github.com/cjdoris/Genus2Conductor.git) : this is included as a submodule in `./libs/` however we encourage you to install through [CHIMP](https://github.com/edgarcosta/CHIMP).
- [ExactpAdics2](https://github.com/cjdoris/ExactpAdics2.git) : this is a dependency of `Genus2Conductor` and is included as a submodule in `./libs/` again you can install through [CHIMP](https://github.com/edgarcosta/CHIMP).
- [AndrewVSutherland/Magma](https://github.com/AndrewVSutherland/Magma.git) : this is included as a submodule in `./libs/`. It is used in `2-3-prop.m` to convert a Conrey label to a Dirichlet character (for computing modular forms) specifically it depends on [chars.m](https://github.com/AndrewVSutherland/Magma/blob/main/chars.m) (which depends on [utils.m](https://github.com/AndrewVSutherland/Magma/blob/main/utils.m)).
- [genus-2-RM](https://github.com/SamFrengley/genus-2-RM) : this is included as a submodule in `./libs/` and is a simple way of extracting the models for Hilbert modular surfaces given by Elkies--Kumar in [this article](https://arxiv.org/abs/1209.3527). To use `genus-2-RM` you should make sure you run the `setup.py` script in that repository.
- [sha-7-examples](https://github.com/SamFrengley/sha-7-examples) : this is included as a submodule in `./libs/` and is only needed for verifying the claims in Section 5.
