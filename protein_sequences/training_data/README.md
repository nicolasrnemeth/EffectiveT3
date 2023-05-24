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
