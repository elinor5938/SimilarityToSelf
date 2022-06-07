# SimilarityToSelf

This script contain two function to get the scores of Needlman-Wunsch algorithm for global alignment for peptides to any target peptiotme.
I used it to compare the human peptidome to my reprence peptides.

The input should be in FASTA format.

send_needle_run_one_pep is faster version (using multiprocessing) of the function send_needle_run .
the fasta peptidome file should be given as a deafult parameter.
