import pytest
import sys
from unittest import TestCase

sys.path.insert(1, '..')
sys.path.insert(1, '../src/EUKulele')
import EUKulele
import subprocess
from EUKulele.EUKulele_config import eukulele
import yaml
import os

test_reference = "mmetsp"

def test_setup():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    base_config = os.path.join(os.path.dirname(__file__), '..', 'aux_data', 'config.yaml')
    with open(base_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
       
    outputdir = os.path.join(base_dir, "test_out")
    os.system("rm -rf " + outputdir)
    
    config["mets_or_mags"] = "mets"
    config["reference"] = os.path.join(base_dir, test_reference)
    config["samples"] = os.path.join(base_dir, test_reference, "samples_METs")
    config["subroutine"] = "setup"
    config["output"] = outputdir
    config["database"] = test_reference
    config["download_reference"] = 0
    config["column"] = "SOURCE_ID"
    config["nucleotide_extension"] = ".fasta"
    config["ref_fasta"] = "reference-pep-trunc.pep.faa"
    
    config_path = os.path.join(base_dir, 'test_configs')
    os.system("mkdir -p " + config_path)
    config_file = os.path.join(config_path, 'curr_config_setup.yaml')
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
        
    eukulele(config=config_file)
    assert os.path.isfile(os.path.join(config["reference"], config["tax_table"]))
    
def test_setup_commandline():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_METs")
    output_dir = os.path.join(base_dir, "test_out")
    reference_dir = os.path.join(base_dir, test_reference)
    os.system("rm -rf " + output_dir)
    subprocess.Popen(["EUKulele", "setup", "--database", "mmetsp", "--sample_dir", sample_dir, 
                      "--mets_or_mags", "mets", "--out_dir", output_dir, 
                      "--reference_dir", reference_dir])
    
def test_all_commandline():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_METs")
    output_dir = os.path.join(base_dir, "test_out")
    reference_dir = os.path.join(base_dir, test_reference)
    os.system("rm -rf " + output_dir)
    
    subprocess.Popen(["EUKulele", "all", "--database", "mmetsp", "--sample_dir", sample_dir, 
                      "--mets_or_mags", "mets", "--out_dir", output_dir, "--organisms", "Chromera",
                      "--taxonomy_organisms", "genus", "--reference_dir", reference_dir]).wait()
    
    samplenames = [curr.split(".")[0] for curr in os.listdir(sample_dir)]
    busco_out = os.path.join(output_dir, "busco_assessment", samplenames[0], 
                             "species_combined", "summary_species_" + samplenames[0] + ".tsv")
    assert os.path.isfile(busco_out)
    
def test_all_commandline_busco():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    sample_dir = os.path.join(base_dir, test_reference, "samples_METs")
    output_dir = os.path.join(base_dir, "test_out")
    reference_dir = os.path.join(base_dir, test_reference)
    os.system("rm -rf " + output_dir)
    subprocess.Popen(["EUKulele", "all", "--database", "mmetsp", "--sample_dir", sample_dir, 
                      "--mets_or_mags", "mets", "--out_dir", output_dir, "--individual_or_summary",
                      "summary", "--reference_dir", reference_dir]).wait()
    
    samplenames = [curr.split(".")[0] for curr in os.listdir(sample_dir)]
    busco_out = os.path.join(output_dir, "busco_assessment", samplenames[0], 
                             "species_combined", "summary_species_" + samplenames[0] + ".tsv")
    assert os.path.isfile(busco_out)