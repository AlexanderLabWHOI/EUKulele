'''
A second round of metagenomic test cases.
'''

import os
import sys
import yaml
import pytest
from unittest import TestCase

sys.path.insert(1, '..')
sys.path.insert(1, '../src/EUKulele')
import EUKulele
from EUKulele.EUKulele_config import eukulele

test_reference = "mmetsp"

def test_individual():
    '''
    Test running the BUSCO functionality with individually specified
    arguments, rather than a full class/other functional group.
    '''

    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out_E")
    os.system("rm -rf " + output_dir)
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref_MAG")

    string_arguments = " ".join(["all", "--database", "mmetsp",
                                 "--sample_dir", sample_dir,
                                 "--mets_or_mags", "mags",
                                 "--out_dir", output_dir, "-i",
                                 '--organisms', 'Chromera',
                                 '--taxonomy_organisms', 'genus',
                                 "--reference_dir", reference_dir])

    eukulele(string_arguments=string_arguments)
    samplenames = [curr.split(".")[0] for curr in os.listdir(sample_dir)]
    busco_out = os.path.join(output_dir, "busco_assessment",
                             samplenames[0], "individual",
                             "summary_" + samplenames[0] + ".tsv")
    assert os.path.isfile(busco_out)

def test_error_input():
    '''
    Checking that we can catch an expected error in the input formatting.
    '''

    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out_F")
    os.system("rm -rf " + output_dir)
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref_MAGs")

    string_arguments = " ".join(["--database", "mmetsp", "--sample_dir", sample_dir,
                                 "--mets_or_mags", "mmm", "--out_dir", output_dir, "-i",
                                 '--organisms', 'Chromera', '--taxonomy_organisms', 'genus',
                                 "--reference_dir", reference_dir])
    error = 0
    try:
        eukulele(string_arguments=string_arguments)
    except:
        error = 1

    assert error == 1

def test_error_required_input():
    '''
    Tests that we receive an error when a required input is not included.
    '''

    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out_G")
    os.system("rm -rf " + output_dir)
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")

    string_arguments = " ".join(["--database", "mmetsp", "--sample_dir", sample_dir,
                                 "--out_dir", output_dir, "-i",
                                 '--organisms', 'Chromera', '--taxonomy_organisms', 'genus',
                                 "--reference_dir", reference_dir])
    error = 0
    try:
        eukulele(string_arguments=string_arguments)
    except:
        error = 1

    assert error == 1

def test_error_busco_no_orgs():
    '''
    Test that we get an error when we do not specify an organism.
    '''

    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out_H")
    os.system("rm -rf " + output_dir)
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")

    string_arguments = " ".join(["--database", "mmetsp", "--sample_dir", sample_dir,
                                 "--mets_or_mags", "mags", "--out_dir", output_dir, "-i",
                                 "--reference_dir", reference_dir])
    error = 0
    try:
        eukulele(string_arguments=string_arguments)
    except:
        error = 1

    assert error == 1

def test_error_n_extension():
    '''
    Tests that we get an error when we use an improper nucleotide
    extension for the files.
    '''

    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out_I")
    os.system("rm -rf " + output_dir)
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")

    string_arguments = " ".join(["--database", "mmetsp", "--sample_dir", sample_dir,
                                 "--mets_or_mags", "mets", "--out_dir", output_dir, "-i",
                                 "--n_ext", ".hello",
                                 "--reference_dir", reference_dir])
    error = 0
    try:
        eukulele(string_arguments=string_arguments)
    except:
        error = 1

    assert error == 1

def test_error_p_extension():
    '''
    Tests that we get an error when we use an improper protein
    extension for the files.
    '''

    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out_J")
    os.system("rm -rf " + output_dir)
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")

    string_arguments = " ".join(["--database", "mmetsp", "--sample_dir", sample_dir,
                                 "--mets_or_mags", "mags", "--out_dir", output_dir, "-i",
                                 "--p_ext", ".hello",
                                 "--reference_dir", reference_dir])
    error = 0
    try:
        eukulele(string_arguments=string_arguments)
    except:
        error = 1

    assert error == 1

def test_error_busco():
    '''
    Tests that we get an error when we use a BUSCO file that does
    not exist.
    '''

    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out_K")
    os.system("rm -rf " + output_dir)
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")

    string_arguments = " ".join(["--database", "mmetsp", "--sample_dir", sample_dir,
                                 "--mets_or_mags", "mags", "--out_dir", output_dir, "-i",
                                 '--busco_file', os.path.join(base_dir, test_reference,
                                                              "samples_MAGs",
                                                              "busco_file_fake.tsv"),
                                 "--reference_dir", reference_dir])
    error = 0
    try:
        eukulele(string_arguments=string_arguments)
    except:
        error = 1

    assert error == 1

def test_busco_file():
    '''
    Tests the functionality of using a file to specify the organisms
    to explore via BUSCO.
    '''

    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out_K")
    os.system("rm -rf " + output_dir)
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref_MAG")

    string_arguments = " ".join(["--database", "mmetsp", "--sample_dir", sample_dir,
                                 "--mets_or_mags", "mags", "--out_dir",
                                 output_dir, "-i",
                                 '--busco_file', os.path.join(base_dir,
                                                              test_reference,
                                                              "samples_MAGs",
                                                              "test_busco.tsv"),
                                 "--reference_dir", reference_dir])
    error = 0
    eukulele(string_arguments=string_arguments)
    samplenames = [curr.split(".")[0] for curr in os.listdir(sample_dir)]
    busco_out = os.path.join(output_dir, "busco_assessment", samplenames[0], "individual",
                             "summary_" + samplenames[0] + ".tsv")
    out_prefix = samplenames[0]

    assert os.path.isfile(busco_out)

def test_all():
    '''
    Combined test case.
    '''

    base_dir = os.path.join(os.path.dirname(__file__), '..',
                            'aux_data')
    base_config = os.path.join(os.path.dirname(__file__), '..',
                               'aux_data', 'config.yaml')
    base_config_curr = os.path.join(os.path.dirname(__file__), '..',
                                    'aux_data', 'config_O.yaml')
    os.system("cp " + base_config + " " + base_config_curr)
    with open(base_config_curr) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    config["mets_or_mags"] = "mags"
    config["reference"] = os.path.join(base_dir, test_reference,
                                       "sample_ref_MAG")
    config["samples"] = os.path.join(base_dir, "real-world-samples", "MAGs")
    config["subroutine"] = "all"
    config["individual_or_summary"] = "summary"
    config["cutoff"] = os.path.join("tax-cutoffs.yaml")
    config["output"] = os.path.join(base_dir, "test_out_all_K")
    config["database"] = test_reference
    config["organisms"] = ["Chromera"]
    config["taxonomy_organisms"] = ["genus"]
    config["download_reference"] = 0
    config["column"] = "SOURCE_ID"
    config["ref_fasta"] = "reference.pep.fa"
    config["protein_map"] = "prot-map.json"
    config["tax_table"] = "tax-table.txt"

    config_path = os.path.join(base_dir, 'test_configs')
    os.system("mkdir -p " + config_path)
    config_file = os.path.join(config_path, 'curr_config_busco_O.yaml')
    with open(config_file, 'w') as f:
        yaml.dump(config, f)

    eukulele(string_arguments=" ".join(["--config",config_file]))
    samplenames = [curr.split(".")[0] for curr in os.listdir(config["samples"])]
    busco_out = os.path.join(config["output"], "busco_assessment",
                             samplenames[0], "species_combined",
                             "summary_species_" + samplenames[0] + ".tsv")
    out_prefix = samplenames[0]
    mag_file = os.path.join(config["output"], "levels_mags",
                            out_prefix + '.' + "species")
    assert (os.path.isfile(busco_out)) & (os.path.isfile(mag_file))

def test_tester():
    '''
    Tests the setup of the tests.
    '''

    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out")
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")
    os.system("rm -rf " + output_dir)

    string_arguments = " ".join(["setup", "--test", "--database",
                                 "mmetsp", "--sample_dir", sample_dir,
                                 "--mets_or_mags", "mags", "--out_dir",
                                 output_dir, "--ref_fasta",
                                 "reference.pep.fa", "--reference_dir",
                                 reference_dir])

    eukulele(string_arguments=string_arguments)
    assert (not os.path.isdir(output_dir))

#def test_cleanup():
#    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
#    config_path = os.path.join(os.path.dirname(__file__), '..', 'aux_data', 'test_configs')
#    base_configs = [os.path.join(config_path, 'curr_config_alignment.yaml'),\
#                    os.path.join(config_path, 'curr_config_setup.yaml')]

#    successful_test = True
#    for base_config in base_configs:
#        with open(base_config) as f:
#            config = yaml.load(f, Loader=yaml.FullLoader)

#        config["reference"] = os.path.join(base_dir, test_reference)
#       os.system("rm -rf " + os.path.join(config["output"]))

#        successful_test = successful_test & (not os.path.isdir(os.path.join(config["output"])))
#        successful_test = True
#    assert successful_test
