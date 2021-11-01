wget https://raw.githubusercontent.com/klh272/pici-repository/main/BLAST_protein_db.faa -P ./databases/
makeblastdb -in ./databases/BLAST_protein_db.faa -input_type fasta -dbtype prot -out ./databases/PICI_BLAST_DB
