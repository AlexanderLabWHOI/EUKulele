Auxiliary Scripts
=================

A few scripts are included with ``EUKulele`` that can be used outside of the software, or before invoking the program to obtain the required configuration.

Create Protein Table
--------------------

``create_protein_table.py`` is used to generate a taxonomy table and protein map JSON file from a provided taxonomy table tab-delimited file, and a FASTA file containing database protein sequences. The file is invoked using the following arguments:

- ``--infile_peptide``: the FASTA file containing the database peptide sequences
- ``--infile_taxonomy``: the tab-separated file containing the taxonomy of each peptide sequence in the database
- ``--outfile_json``: the protein map JSON file to be written
- ``--output``: the formatted taxonomy table to be written
- ``--delim``: the delimiter that separates tokens in the FASTA headers in the peptide sequence file
- ``--col_source_id``: the column in your tab-separated taxonomy file containing the name of the strain
- ``--taxonomy_col_id``: the column containing the strain taxonomy in the tab-separated taxonomy file
- ``--column``: either a character name for the token in the FASTA headers of the peptide sequence file containing the strain name (matching what is in the taxonomy file), or a numeric value indicating the order of the token
- ``--reformat_tax``: if included, indicates that the taxonomy table should be reformatted instead of left as-is
- ``--euk-prot``: if included, means we are using input from the EukProt database, which has different formatting features.

Download Database
-----------------

Used within ``EUKulele``, ``download_database.sh`` may also be used independently of ``EUKulele`` to download one of the available databases provided with ``EUKulele`` (see Section :ref:`databases`). Invoked via::

    download_database.sh <DATABASE> <REF_FASTA> <REF_TABLE> <REF_FASTA_URL> <REF_TABLE_URL> <REFERENCE_DIR>
    
Where <DATABASE> is the name of the database, <REF_FASTA> is the name of the reference FASTA file to be generated/downloaded, REF_TABLE is the same but for the tab-delimited taxonomy table, <REF_FASTA_URL> and <REF_TABLE_URL> are the URLs to download said databases from (provided in ``reference_url.yaml``), and <REFERENCE_DIR> is where to store the resulting database.
