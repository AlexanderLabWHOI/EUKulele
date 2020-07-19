import pytest
import sys
from unittest import TestCase

sys.path.insert(1, '..')
import EUKulele
import yaml
import os

test_reference = "mmetsp"

def test_setup():
    base_config = os.path.join(os.getcwd(), '..', 'aux_data', 'config.yaml')
    base_dir = os.path.join('EUKulele', 'tests', 'aux_data')
    base_config = os.path.join('EUKulele', 'tests', 'aux_data', 'config.yaml')
    with open(base_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
       
    config["mets_or_mags"] = "mags"
    config["reference"] = os.path.join(base_dir, test_reference)
    config["samples"] = os.path.join(base_dir, test_reference, "samples_MAGs")
    config["subroutine"] = "setup"
    config["output"] = os.path.join(base_dir, "test_out")
    config["database"] = test_reference
    config["download_reference"] = 0
    config["column"] = "SOURCE_ID"
    config["ref_fasta"] = "reference-pep-trunc.pep.faa"
    config["original_taxonomy"] = "taxonomy-table.txt"
    
    config_path = os.path.join(base_dir, 'test_configs')
    os.system("mkdir -p " + config_path)
    config_file = os.path.join(config_path, 'curr_config_setup.yaml')
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
        
    EUKulele.eukulele(config=config_file)
    assert os.path.isfile(os.path.join(config["reference"],config["tax_table"]))
       
def test_alignment():
    base_config = os.path.join(os.getcwd(), '..', 'aux_data', 'config.yaml')
    base_dir = os.path.join('EUKulele', 'tests', 'aux_data')
    base_config = os.path.join('EUKulele', 'tests', 'aux_data', 'config.yaml')
    with open(base_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
       
    config["mets_or_mags"] = "mags"
    config["reference"] = os.path.join(base_dir, test_reference)
    config["samples"] = os.path.join(base_dir, test_reference, "samples_MAGs")
    config["subroutine"] = "alignment"
    config["create_tax_table"] = 0
    config["cutoff"] = os.path.join("EUKulele","static","tax-cutoffs.yaml")
    config["output"] = os.path.join(base_dir, "test_out")
    config["database"] = test_reference
    config["download_reference"] = 0
    config["column"] = "SOURCE_ID"
    config["ref_fasta"] = "reference-pep-trunc.pep.faa"
    config["original_taxonomy"] = "taxonomy-table.txt"
    
    config_path = os.path.join(base_dir, 'test_configs')
    os.system("mkdir -p " + config_path)
    config_file = os.path.join(config_path, 'curr_config_alignment.yaml')
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
        
    EUKulele.eukulele(config=config_file)
    outprefix = config["output"].split("/")[-1]
    assert os.path.isfile(os.path.join(config["output"],out_prefix + "_all_species_counts.csv"))
    
def test_cleanup():
    base_config = os.path.join(os.getcwd(), '..', 'aux_data', 'config.yaml')
    base_dir = os.path.join('EUKulele', 'tests', 'aux_data')
    base_config = os.path.join('EUKulele', 'tests', 'aux_data', 'config.yaml')
    with open(base_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
        
    config["reference"] = os.path.join(base_dir, test_reference)
    os.system("rm " + os.path.join(config["reference"], "tax-table.txt"))
    os.system("rm " + os.path.join(config["reference"], "protein-map.json"))
              
    successful_test = (not os.path.isfile(os.path.join(config["reference"],"tax-table.txt"))) & \
                      (not os.path.isfile(os.path.join(config["reference"],"protein-map.json")))
              
    assert successful_test
        
    