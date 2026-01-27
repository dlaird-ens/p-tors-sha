This directory contains the code to prove Theorem 1.8. 

## Structure
- `congruence_search.sage` this is the main file containing the functions to perform the seive described in the article. 
- `congruence_run.sage` is the file which is run to perform the sieve proving Theorem 1.8.
- `combine_lset_run.sage` proves Theorem 1.8 (after `congruence_run.sage` is run for all the sets $\mathcal{L}_i$), see below for a sample script.
- `sha_helpers.sage`
  - `comp_gp_kohelstein.m` is a `magma` script to compute the component group of a modular abelian variety using the Kohel--Stein algorithm.
  - `comp_gp_regmod.m` is a `magma` script to compute the component group using `RegularModel`.
- `sha_run.sage`
  - This relies heavily on the contributions of the work [CEH+] in preparation. They give genus 2 curves whose (RM) Jacobians have the corresponding modular forms. This is (in our local copy) recorded in `../data/label_to_curve.json` but this is not our data to share. Please contact us if you wish to run this code. Note that we make no completeness statements about the examples in Theorem 1.4 so this code is only used to *find* the examples, not to prove them.
- `tables.sage` prints the tables from the article into the directory `../data/`

## Proof of Theorem 1.8
This can be done by running the following shell script, which completes the proof. 

*Warning:* This will take a long time and it may be better to do this in the background. For us this needs about 5GB of virtual memory (if you choose to run the commands in parallel rather than in a loop you will need $5n$GB).

```bash
for ii in 0 1 2 3 4 5
do
  # Run the sieve for each of the sets LL_i. This outputs into 
  # the directories `../data/congruences/i/`. If you want for it
  # to output somewhere better for you change the --d flag
  sage congruence_run.sage --d="../data/" --l=$ii
done

# Now combine these runs and output it into 
# `../data/congruences/combined/` and `../data/congruences/p.txt`
sage combine_lset_run.sage --d="../data/"
```
