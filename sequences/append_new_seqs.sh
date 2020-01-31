for FILE in `ls EPI*`; do echo "  " >> gisaid_cov2020_all_sequences.fasta ; cat ${FILE} >> gisaid_cov2020_all_sequences.fasta ; done
