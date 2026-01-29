This directory contains the data for the article. The file `nfs_dim_db_4.json` (which will not be automatically populated) contains the newforms in the LMFBD whose coefficient fields have degree $\leq 4$ (in a format which can be used by the files in `../congruence_search`). This data can also be extracted using `./code/run_extraction.py` (this may take some time and a stable internet connection).

## The (probable) congruences 
The pairs of $(\mathfrak{p}, \mathfrak{q})$-congruent forms claimed in Theorem 1.8 are listed in the files `./congruences/p.txt` where `p` is the prime lying above $\mathfrak{p}$ and $\mathfrak{q}$.

## The tables
The tables from the article are stored in the files labelled `*-*-table.txt` as colon separated values with rows separated by newlines. The individual files have the following formats for each line
- The table `B-1-table.txt` has rows (note that the traces and norms record the first 30 traces and norms of $a_p$ with $p$ *prime*)
  - `p:label_f:tr_f:nm_f:label_g:tr_g:nm_g`
- The table `B-2-table.txt` has rows
  - `p:label_A:eqn_A:dim_A:rk_A:cp_A:label_B:eqn_B:dim_B:rk_B:cp_B`
- The table `B-3-table.txt` has rows
  - `label:disc:eqn`

## Directories for partial output 
The directories`./sha/` and `./local/` are empty but will be populated when the code in `../congruence_search/` is run. 

The directory `./congruences` contains empty subdirectories `./congruences/i/` and `./congruences/combined/` which again will be populated by the code in `../congruence_search/`. In terms of the proof of Theorem 1.8 the directories `./congruences/i/` contain the forms with a hash in common for the set $\mathcal{L}_{i}$.

After the code is run, the directory `./congruences/combined/` will contain the congruences `./congruences/p.txt` except the $a_p$-values are also stored (so the files can be rather large). 
