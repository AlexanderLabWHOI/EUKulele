# Test the database formatting function of EUKulele using a subset database in the test-samples folder.

import pytest
import sys
from unittest import TestCase

sys.path.insert(1, '..')
from EUKuleleconfig import *
#from EUKulele import EUKuleleconfig
import yaml
import os

#pytestmark = pytest.mark.random_order(disabled=True)

def test_generate_tax_table():

    # write base config file with test samples to go into testing directory
    # modify the flags in the config for each of the test cases
    # check whether there is no error returned (try/catch in test case?) as
    # well as whether the returned file is in the correct format
    #os.path.join('EUKulele', 'tests', 'aux_data', 'config.yaml')
    base_config = os.path.join('aux_data', 'config.yaml')
    with open(base_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
        
    print("hello",flush=True)
    print(config)
    config["reference"] = "curr_ref"
    config["subroutine"] = "setup"
    config["database"] = "phylodb"
    config["download_reference"] = 1
    config["ref_fasta"] = "ref-phylodb-trunc.pep.fa"
    
    config_path = os.path.join('EUKulele', 'tests','curr_config.yaml')
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    EUKuleleconfig.eukulele(config=config_path)
    self.assertTrue(os.path.isfile(os.path.join(config["reference"],config["ref_fasta"])))
   

base_config = os.path.join('aux_data', 'config.yaml')
with open(base_config) as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

print("hello",flush=True)
test_generate_tax_table()