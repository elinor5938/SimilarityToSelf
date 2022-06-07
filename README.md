# SimilarityToSelf

This script contain two function to get the scores of Needlman-Wunsch algorithm for global alignment for peptides to any target peptidome.
I used it to compare the human peptidome to my reference peptides.

The input should be in FASTA format.

send_needle_run_one_pep is faster version (using multiprocessing) of the function send_needle_run .
The fasta peptidome file should be given as a default parameter.

df_creator - return df of the results.


## you should download the EMBOSS from here in your Linux to initiate the program -

ftp://emboss.open-bio.org/pub/EMBOSS/. 
