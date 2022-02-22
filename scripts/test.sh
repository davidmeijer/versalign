#!/usr/bin/env bash
./run_multiple_sequence_alignment.py --fasta ../data/pk_test.fasta --scores ../scoring_matrices/scoring_matrix_similarity.txt > ../_out/msa_similarity.msa
./run_multiple_sequence_alignment.py --fasta ../data/pk_test.fasta --scores ../scoring_matrices/scoring_matrix_match_mismatch.txt > ../_out/msa_match_mismatch.msa
./make_logo_from_msa.py ../_out/msa_match_mismatch.msa --out ../_out/logo_match_mismatch.png
./make_logo_from_msa.py ../_out/msa_similarity.msa --out ../_out/logo_similarity.png

