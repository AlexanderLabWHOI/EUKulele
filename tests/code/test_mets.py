'''
Test cases for MET functionality.
'''

import os
import yaml
import sys
import pytest
from unittest import TestCase

sys.path.insert(1, '..')
sys.path.insert(1, '../src/EUKulele')
import EUKulele
import subprocess
from EUKulele.EUKulele_config import eukulele

test_reference = "mmetsp"

def test_setup():
    base_dir = os.path.join(os.path.dirname(__file__), '..',
                            'aux_data')
    base_config = os.path.join(os.path.dirname(__file__), '..',
                               'aux_data', 'config.yaml')
    base_config_curr = os.path.join(os.path.dirname(__file__), '..',
                                    'aux_data', 'config_P.yaml')
    os.system("cp " + base_config + " " + base_config_curr)
    with open(base_config_curr) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    outputdir = os.path.join(base_dir, "test_out_P")
    os.system("rm -rf " + outputdir)

    config["mets_or_mags"] = "mets"
    config["reference"] = os.path.join(base_dir, test_reference,
                                       "sample_ref")
    config["samples"] = os.path.join(base_dir, test_reference,
                                     "samples_METs_small")
    config["subroutine"] = "setup"
    config["output"] = outputdir
    config["database"] = test_reference
    config["download_reference"] = 0
    config["column"] = "SOURCE_ID"
    config["nucleotide_extension"] = ".fasta"
    config["ref_fasta"] = "reference.pep.fa"
    config["protein_map"] = "prot-map.json"
    config["tax_table"] = "tax-table.txt"

    config_path = os.path.join(base_dir, 'test_configs')
    os.system("mkdir -p " + config_path)
    config_file = os.path.join(config_path, 'curr_config_setup_P.yaml')
    with open(config_file, 'w') as f:
        yaml.dump(config, f)

    eukulele(config=config_file)
    assert os.path.isfile(os.path.join(config["reference"], config["tax_table"]))

def test_setup_commandline():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_METs_small")
    output_dir = os.path.join(base_dir, "test_out_Q")
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")
    os.system("rm -rf " + output_dir)
    #subprocess.Popen(["EUKulele", "setup", "--database", "mmetsp", "--sample_dir", sample_dir,
    #                  "--mets_or_mags", "mets", "--out_dir", output_dir,
    #                  "--reference_dir", reference_dir])

    string_arguments = " ".join(["setup", "--database", "mmetsp", "--sample_dir", sample_dir,
                      "--mets_or_mags", "mets", "--out_dir", output_dir, "--ref_fasta",
                      "reference.pep.fa", "--reference_dir", reference_dir])

    eukulele(string_arguments=string_arguments)
    assert os.path.isfile(os.path.join(reference_dir, "tax-table.txt"))

def test_all_commandline():
    """
    Tests that alignment works properly using a string of arguments.
    """

    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_METs_small")
    output_dir = os.path.join(base_dir, "test_out_R")
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")

    #EUKulele alignment --database mmetsp --sample_dir
    # tests/aux_data/mmetsp/samples_METs_small --mets_or_mags mets
    # --out_dir tests/test_out --organisms Chromera --taxonomy_organisms
    # genus --reference_dir tests/aux_data/mmetsp

    string_arguments=" ".join(["alignment", "--database", "mmetsp", "--sample_dir",
                               sample_dir, "--mets_or_mags", "mets", "--out_dir",
                               output_dir, "--organisms", "Chromera", "--ref_fasta",
                               "reference.pep.fa", "--taxonomy_organisms", "genus",
                               "--reference_dir", reference_dir])

    eukulele(string_arguments=string_arguments)
    samplenames = [curr.split(".")[0] for curr in os.listdir(sample_dir)]
    est_out = os.path.join(output_dir, "taxonomy_estimation",
                           samplenames[0] + "-estimated-taxonomy.out")
    assert os.path.isfile(est_out)

def test_all_commandline_busco():
    """
    Tests that BUSCO runs using arguments provided as a string.
    """
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_METs_small")
    output_dir = os.path.join(base_dir, "test_out_S")
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")

    string_arguments = ["setup", "--database", "mmetsp", "--sample_dir", sample_dir,
                      "--mets_or_mags", "mets", "--out_dir", output_dir,
                      "--individual_or_summary","summary",
                      "--ref_fasta", "reference.pep.fa", "--reference_dir", reference_dir]

    eukulele(string_arguments=" ".join(string_arguments))
    string_arguments[0] = "alignment"
    eukulele(string_arguments=" ".join(string_arguments))
    string_arguments[0] = "busco"
    eukulele(string_arguments=" ".join(string_arguments))

    samplenames = [curr.split(".")[0] for curr in os.listdir(sample_dir)]
    busco_out = os.path.join(output_dir, "busco_assessment", samplenames[0],
                             "species_combined", "summary_species_" + samplenames[0] + ".tsv")
    assert os.path.isfile(busco_out)

def test_all_commandline_busco_individual():
    """
    Tests that BUSCO runs using arguments provided as a string.
    """
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_METs_small")
    output_dir = os.path.join(base_dir, "test_out_T")
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")

    string_arguments = ["setup", "--database", "mmetsp", "--sample_dir", sample_dir,
                      "--mets_or_mags", "mets", "--out_dir", output_dir,
                      "--individual_or_summary","individual",'--organisms', 'Chromera',
                      '--taxonomy_organisms', 'genus',"--ref_fasta", "reference.pep.fa",
                      "--reference_dir", reference_dir]

    eukulele(string_arguments=" ".join(string_arguments))
    string_arguments[0] = "alignment"
    eukulele(string_arguments=" ".join(string_arguments))
    string_arguments[0] = "busco"
    eukulele(string_arguments=" ".join(string_arguments))

    samplenames = [curr.split(".")[0] for curr in os.listdir(sample_dir)]
    busco_out = os.path.join(output_dir, "busco_assessment", samplenames[0],
                             "individual", "summary_" + samplenames[0] + ".tsv")
    #busco_out = os.path.join(output_dir, "busco_assessment", samplenames[0],
    #                         "species_combined", "summary_species_" + samplenames[0] + ".tsv")
    assert os.path.isfile(busco_out)

def test_all_force_rerun():
    """
    Tests that BUSCO runs using arguments provided as a string.
    """
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_METs_small")
    output_dir = os.path.join(base_dir, "test_out_U")
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")
    os.system("rm -rf " + output_dir)

    string_arguments = " ".join(["all", "--database", "mmetsp", "--sample_dir", sample_dir,
                                 "--mets_or_mags", "mets", "--out_dir", output_dir,
                                 "--busco_threshold", str(30),
                                 "--individual_or_summary","individual",
                                 '--organisms','Chromera', "--filter_metric",
                                 "bitscore", '--taxonomy_organisms', 'genus',
                                 "--ref_fasta", "reference.pep.fa",
                                 "--reference_dir",
                                 reference_dir])#, "-f"])

    eukulele(string_arguments=string_arguments)

    samplenames = [curr.split(".")[0] for curr in os.listdir(sample_dir)]
    busco_out = os.path.join(output_dir, "busco_assessment", samplenames[0],
                             "individual", "summary_" + samplenames[0] + ".tsv")
    assert os.path.isfile(busco_out)

def test_all_use_counts():
    """
    Tests that BUSCO runs using arguments provided as a string.
    """
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_METs_small")
    salmon_dir = os.path.join(base_dir, test_reference, "samples_METs", "salmon_quant")
    output_dir = os.path.join(base_dir, "test_out_V")
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")
    os.system("rm -rf " + output_dir)

    string_arguments = " ".join(["all", "--database", "mmetsp",
                                 "--sample_dir", sample_dir,
                                 "--mets_or_mags", "mets", "--out_dir",
                                 output_dir, "--busco_threshold", str(30),
                                 "--individual_or_summary","individual",
                                 '--organisms', 'Chromera',
                                 "--filter_metric",
                                 "pid",'--taxonomy_organisms',
                                 'genus', "--ref_fasta",
                                 "reference.pep.fa",
                                 "--reference_dir", reference_dir,
                                 "--use_salmon_counts", "--salmon_dir",
                                 salmon_dir])

    eukulele(string_arguments=string_arguments)

    samplenames = [curr.split(".")[0] for curr in os.listdir(sample_dir)]
    busco_out = os.path.join(output_dir, "busco_assessment", samplenames[0],
                             "individual", "summary_" + samplenames[0] + ".tsv")
    assert os.path.isfile(busco_out)
