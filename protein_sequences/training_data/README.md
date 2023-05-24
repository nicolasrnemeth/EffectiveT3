# Collection of training data

### Positive training data

- the sequences of all proteins that have been experimentally or
  sufficiently computationally verified to be secreted by the Type III secretion system of bacteria

### Negative training data

- Downloaded all protein sequences from the ref-seq database for all bacterial species
  that have at least one protein that was experimentally verified to be secreted by the Type III secretion system of bacteria

- Randomly and as much as possible uniformly sampled protein sequences from the proteins in the ref-seq database of each of the species above

- Obtained roughly 9k - 10k protein sequences fulfilling the criteria below that fulfill which where preprocessed and clustered by the cdhit commands below

  - the protein sequence has an N-terminal domain in the range of the first 1-50 amino acids
  - the protein sequence DOES NOT contain a Eukaryotic-like domain

    - -> Assumption 1: proteins containing an N-terminal domain are very likely not secreted
    - -> Assumption 2: Eukaryotic-like domains can be an indication that the protein sequence is secreted by the Type III secretion system

### INTERPROSCAN Version 5.53-87.0

Interproscan in combination with pfam was used to obtain information about the protein domains.
-> this way we figured out which protein sequences contained N-terminal domains and which contained any Eukaryotic-like domains

> ./interproscan.sh -i {path_to_the_fasta_file_with_one_or_more_sequences} -cpu {number_of_cores} -appl pfam -f tsv -dp -t p -T {path_to_temp_folder_where_temporary_data_is_saved} -d {path_to_directory_where_results_are_saved}

(pfam: A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs).)

For a proper description of interproscan and pfam please visit their documentation: https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html
For a proper description of Eukaryotic-like domains and their Domain IDs (PF.....) visit: https://effectors.org/reports/eld_search

---

## CdHit PsiCdHit clustering commands to obtain processed training data with 30% percentage sequence identity cutoff, low redundancy and low evolutionary similarity

the raw (unprocessed) protein sequences can be found in the "unprocessed_sequences" directory

### Clustering T3SE positive dataset (secreted)

> cdhit -i t3se_all_positives.fasta -o t3se_all_positives_90.fasta -c 0.9 -n 5 -g 1 -aS 0.4 -d 0 -p 1 -G 0 -T 0 > t3se_all_positives_90pid.log

> cdhit -i t3se_all_positives_90.fasta -o t3se_all_positives_60pid.fasta -c 0.6 -n 4 -g 1 -aS 0.4 -d 0 -p 1 -G 0 -T 0 > t3se_all_positives_60pid.log

> ./psi-cd-hit_new.pl -i t3se_all_positives_60.fasta -o t3se_all_positives_30.fasta -c 0.3 -P /path/to/blastp/executable/ -prog blastp -aS 0.4 -g 1 -exec local -G 0 -para 4 -blp 1 > t3se_all_positives_30.log

### Clustering T3SE negative dataset (not secreted)

> cdhit -i t3se_all_negatives.fasta -o t3se_all_negatives_90.fasta -c 0.9 -n 5 -g 1 -aS 0.4 -d 0 -p 1 -G 0 -T 0 > t3se_all_negatives_90.log

> cdhit -i t3se_all_negatives_90.fasta -o t3se_all_negatives_60.fasta -c 0.6 -n 4 -g 1 -aS 0.4 -d 0 -p 1 -G 0 -T 0 > t3se_all_negatives_60.log

> ./psi-cd-hit_new.pl -i t3se_all_negatives_60.fasta -o t3se_all_negatives_30.fasta -c 0.3 -P /path/to/blastp/executable/ -prog blastp -aS 0.4 -g 1 -G 0 -exec local -para 4 -blp 1 > t3se_all_negatives_30.log
