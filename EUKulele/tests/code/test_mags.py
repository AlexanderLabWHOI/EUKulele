import pytest
import sys
from unittest import TestCase

sys.path.insert(1, '..')
#from EUKuleleconfig import *
#from EUKulele import EUKuleleconfig
import EUKulele
import yaml
import os

def test_setup():
    base_config = os.path.join(os.getcwd(), '..', 'aux_data', 'config.yaml')
    base_dir = os.path.join('EUKulele', 'tests', 'aux_data')
    base_config = os.path.join('EUKulele', 'tests', 'aux_data', 'config.yaml')
    with open(base_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
       
    config["mets_or_mags"] = "mags"
    config["reference"] = os.path.join(base_dir, "mmetsp")
    config["samples"] = os.path.join(base_dir, "mmetsp", "samples_MAGs")
    config["subroutine"] = "setup"
    config["output"] = os.path.join(base_dir, "test_out")
    config["database"] = "mmetsp"
    config["download_reference"] = 0
    config["column"] = "SOURCE_ID"
    config["ref_fasta"] = "reference-pep-trunc.pep.faa"
    config["original_taxonomy"] = "taxonomy-table.txt"
    
    config_path = os.path.join(base_dir, 'test_configs')
    os.system("mkdir -p " + config_path)
    config_file = os.path.join(config_path, 'curr_config.yaml')
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
        
    EUKulele.eukulele(config=config_file)
    assert os.path.isfile(os.path.join(config["reference"],config["tax_table"]))