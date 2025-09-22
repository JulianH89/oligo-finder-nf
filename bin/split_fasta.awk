# AWK script to split a multi-FASTA file into single-entry FASTA files.
# The output filename for each sequence is derived from its header.

# This block is executed for lines starting with ">" (FASTA headers).
/^>/ {
    # Extract the header, removing the leading ">".
    header = substr($0, 2);

    # Sanitize the header to create a valid filename.
    # Replaces any character that is NOT alphanumeric, underscore, dot, or hyphen with an underscore.
    gsub(/[^a-zA-Z0-9_.-]/, "_", header);

    # If an output file is already open, close it to save changes.
    if (outfile) close(outfile);

    # Define the new output filename based on the sanitized header.
    outfile = header ".fa";
}

# This block is executed for every line in the file.
{
    # If an output file has been defined, append the current line to it.
    if (outfile) print >> outfile;
}