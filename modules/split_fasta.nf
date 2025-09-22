process SPLIT_FASTA {
    tag "${params.run_id} - Splitting Fasta File"

    input:
    path target_gene

    output:
    path "split_dir/*.fa", emit: fasta_files

    script:
    """
    mkdir split_dir

    # This awk command reads the multi-fasta file.
    # When it sees a line starting with '>', it creates a new output file in the split_dir.
    # The filename is derived from the fasta header.
    # All subsequent lines are printed to that new file until the next '>' is found.
    awk '/^>/ {
        if (out) close(out);
        // Sanitize the header to make a valid filename
        header = substr(\$0, 2);
        gsub(/[^a-zA-Z0-9_-]/, "_", header);
        out = "split_dir/" header ".fa"
    } { if (out) print \$0 > out }' ${target_gene}
    """
}