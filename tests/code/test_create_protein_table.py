import pytest
import sys
from unittest import TestCase

sys.path.insert(1, '..')
sys.path.insert(1, '../src/EUKulele')
import EUKulele
from EUKulele.EUKulele_config import eukulele
import scripts
from scripts.create_protein_table import createProteinTable

import yaml
import os
import subprocess

test_reference = "mmetsp"

def test_create_protein():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")

    string_arguments = ["create_protein_table.py", "--infile_peptide",
                        os.path.join(reference_dir, "reference.pep.fa"),
                        "--infile_taxonomy",
                        os.path.join(reference_dir, "tax-table.txt"),
                        "--outfile_json",
                        os.path.join(reference_dir, "test_table.tsv",
                                     "test-protein-map.json"),
                        "--output", os.path.join(reference_dir,
                                                 "test_table.txt"),
                        "--delim", "/", "--col_source_id",
                        "Source_ID", "--taxonomy_col_id",
                        "taxonomy", "--column",
                        "SOURCE_ID"]

    subprocess.Popen(string_arguments)

    assert os.path.isfile(os.path.join(reference_dir, "test_table.txt"))

    os.system("rm " + os.path.join(reference_dir, "test_table.txt"))
    string_arguments = ["--infile_peptide",
                         os.path.join(reference_dir, "reference.pep.fa"), "--infile_taxonomy",
                         os.path.join(reference_dir, "tax-table.txt"), "--outfile_json",
                         os.path.join(reference_dir, "test-protein-map.json"),
                         "--output", os.path.join(reference_dir, "test_table.txt"), "--delim", "/",
                         "--col_source_id", "Source_ID", "--taxonomy_col_id", "taxonomy",
                         "--column", "SOURCE_ID"]
    createProteinTable(args = string_arguments)


    assert os.path.isfile(os.path.join(reference_dir, "test_table.txt"))

#def test_create_protein_function():
#    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
#    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")
    #os.system("rm " + os.path.join(reference_dir, "test_table.txt"))

#    string_arguments = ["--infile_peptide",
#                         os.path.join(reference_dir, "reference.pep.fa"), "--infile_taxonomy",
#                         os.path.join(reference_dir, "tax-table.txt"), "--outfile_json",
#                         os.path.join(reference_dir, "test-protein-map.json"),
#                         "--output", os.path.join(reference_dir, "test_table.txt"), "--delim", "/",
#                         "--col_source_id", "Source_ID", "--taxonomy_col_id", "taxonomy",
#                         "--column", "SOURCE_ID"]
#    createProteinTable(string_arguments)
#    assert os.path.isfile(os.path.join(reference_dir, "test_table.txt"))
