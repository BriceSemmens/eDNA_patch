    # NOTE may need to change the limit of number of files you can open with 'ulimit -n 9999' before running this job

    # Need to move into the directory where the taxdb lives!
    cd /Volumes/CalCOFI_OMICS/Databases/nt_euk
    PATH=$PATH:/usr/local/ncbi-blast-2.16.0+/bin

    BLAST_DB='/Volumes/CalCOFI_OMICS/Databases/nt_euk/nt_euk' #The folder of blast database
    QUERY_FASTA='/Users/nastassiapatin/OneDrive - SCCWRP/Experiments/MiFish optimization/dada2-tests/outputs/seqs_to_annotate.fasta' #The fasta file that you want to blast

    # BLAST PARAMETERS
    PERCENT_IDENTITY="97"
    WORD_SIZE="30"
    EVALUE="1e-30"

    # number of matches recorded in the alignment:
    MAXIMUM_MATCHES="50"
    CULLING="3"
	BLAST_OUTPUT="/Users/nastassiapatin/OneDrive - SCCWRP/Experiments/Intercalibration/MFU/outputs/BLAST_output.txt" #The output file

    # NOTE taxidlist here restricts blast to vertebrates
    blastn -query "${QUERY_FASTA}" -db "${BLAST_DB}" \
    -num_threads 4 \
    -perc_identity "${PERCENT_IDENTITY}" \
    -word_size "${WORD_SIZE}" \
    -evalue "${EVALUE}" \
    -max_target_seqs "${MAXIMUM_MATCHES}" \
    -culling_limit="${CULLING}" \
    -outfmt "6 sscinames scomnames qseqid sseqid pident length mismatch gapopen qcovus qstart qend sstart send evalue bitscore staxids qlen qcovs" \
    -out "${BLAST_OUTPUT}" \
    -taxidlist "/Volumes/CalCOFI_OMICS/Databases/7742.txids"
