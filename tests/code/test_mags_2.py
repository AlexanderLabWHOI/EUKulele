import pytest
import sys
from unittest import TestCase

sys.path.insert(1, '..')
sys.path.insert(1, '../src/EUKulele')
import EUKulele
from EUKulele.EUKulele_config import eukulele
import yaml
import os

test_reference = "mmetsp"

def test_individual():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out")
    os.system("rm -rf " + output_dir)
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")
    
    string_arguments = " ".join(["all", "--database", "mmetsp", "--sample_dir", sample_dir, 
                                 "--mets_or_mags", "mags", "--out_dir", output_dir, "-i",
                                 '--organisms', 'Chromera', '--taxonomy_organisms', 'genus',
                                 "--reference_dir", reference_dir])
    
    eukulele(string_arguments=string_arguments)
    samplenames = [curr.split(".")[0] for curr in os.listdir(sample_dir)]
    busco_out = os.path.join(output_dir, "busco_assessment", samplenames[0], "individual", 
                             "summary_" + samplenames[0] + ".tsv")
    assert os.path.isfile(busco_out)
    
def test_error_input():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out")
    os.system("rm -rf " + output_dir)
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")
    
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
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out")
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
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out")
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
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out")
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
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out")
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
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out")
    os.system("rm -rf " + output_dir)
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")
    
    string_arguments = " ".join(["--database", "mmetsp", "--sample_dir", sample_dir, 
                                 "--mets_or_mags", "mags", "--out_dir", output_dir, "-i",
                                 '--busco_file', os.path.join(base_dir, test_reference, "samples_MAGs", "busco_file_fake.tsv"),
                                 "--reference_dir", reference_dir])
    error = 0
    try:
        eukulele(string_arguments=string_arguments)
    except:
        error = 1
    
    assert error == 1
    
def test_busco_file():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out")
    os.system("rm -rf " + output_dir)
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")
    
    string_arguments = " ".join(["--database", "mmetsp", "--sample_dir", sample_dir, 
                                 "--mets_or_mags", "mags", "--out_dir", output_dir, "-i",
                                 '--busco_file', os.path.join(base_dir, test_reference, "samples_MAGs", "test_busco.tsv"),
                                 "--reference_dir", reference_dir])
    error = 0
    eukulele(string_arguments=string_arguments)
    samplenames = [curr.split(".")[0] for curr in os.listdir(sample_dir)]
    busco_out = os.path.join(output_dir, "busco_assessment", samplenames[0], "individual", 
                             "summary_" + samplenames[0] + ".tsv")
    out_prefix = samplenames[0]
    
    assert (os.path.isfile(busco_out))
    
def test_all():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    base_config = os.path.join(os.path.dirname(__file__), '..', 'aux_data', 'config.yaml')
    with open(base_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
       
    config["mets_or_mags"] = "mags"
    config["reference"] = os.path.join(base_dir, test_reference, "sample_ref")
    config["samples"] = os.path.join(base_dir, "real-world-samples", "MAGs")
    config["subroutine"] = "all"
    config["individual_or_summary"] = "summary"
    config["cutoff"] = os.path.join("tax-cutoffs.yaml")
    config["output"] = os.path.join(base_dir, "test_out_all")
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
    config_file = os.path.join(config_path, 'curr_config_busco.yaml')
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
        
    eukulele(string_arguments=" ".join(["--config",config_file]))
    samplenames = [curr.split(".")[0] for curr in os.listdir(config["samples"])]
    busco_out = os.path.join(config["output"], "busco_assessment", samplenames[0], "species_combined", 
                             "summary_species_" + samplenames[0] + ".tsv")
    out_prefix = samplenames[0]
    mag_file = os.path.join(config["output"], "levels_mags", out_prefix + '.' + "species")
    assert (os.path.isfile(busco_out)) & (os.path.isfile(mag_file)) 
    
def test_tester():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_MAGs")
    output_dir = os.path.join(base_dir, "test_out")
    reference_dir = os.path.join(base_dir, test_reference, "sample_ref")
    os.system("rm -rf " + output_dir)
    
    string_arguments = " ".join(["setup", "--test", "--database", "mmetsp", "--sample_dir", sample_dir, 
                      "--mets_or_mags", "mags", "--out_dir", output_dir, "--ref_fasta", 
                      "reference.pep.fa", "--reference_dir", reference_dir])
    
    eukulele(string_arguments=string_arguments)
    assert (not os.path.isdir(output_dir))
    
def test_cleanup():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    config_path = os.path.join(os.path.dirname(__file__), '..', 'aux_data', 'test_configs')
    base_configs = [os.path.join(config_path, 'curr_config_alignment.yaml'),\
                    os.path.join(config_path, 'curr_config_setup.yaml')]
    
    successful_test = True
    for base_config in base_configs:
        with open(base_config) as f:
            config = yaml.load(f, Loader=yaml.FullLoader)

        config["reference"] = os.path.join(base_dir, test_reference)
        os.system("rm -rf " + os.path.join(config["output"]))

        successful_test = successful_test & (not os.path.isdir(os.path.join(config["output"])))
        successful_test = True      
    assert successful_test
       