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

def test_setup():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    base_config = os.path.join(os.path.dirname(__file__), '..', 'aux_data', 'config.yaml')
    with open(base_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
       
    config["mets_or_mags"] = "mags"
    config["reference"] = os.path.join(base_dir, test_reference, "sample_ref")
    config["samples"] = os.path.join(base_dir, test_reference, "samples_MAGs")
    config["subroutine"] = "setup"
    config["output"] = os.path.join(base_dir, "test_out")
    config["database"] = test_reference
    config["download_reference"] = 0
    config["column"] = "SOURCE_ID"
    config["ref_fasta"] = "reference.pep.fa"
    config["protein_map"] = "prot-map.json"
    config["tax_table"] = "tax-table.txt"
    
    config_path = os.path.join(base_dir, 'test_configs')
    os.system("mkdir -p " + config_path)
    os.system("rm -r " + config["output"])
    config_file = os.path.join(config_path, 'curr_config_setup.yaml')
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
        
    eukulele(config=config_file)
    assert os.path.isfile(os.path.join(config["reference"], config["tax_table"]))
       
def test_alignment():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    base_config = os.path.join(os.path.dirname(__file__), '..', 'aux_data', 'config.yaml')
    with open(base_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
       
    config["mets_or_mags"] = "mags"
    config["reference"] = os.path.join(base_dir, test_reference, "sample_ref")
    config["samples"] = os.path.join(base_dir, test_reference, "samples_MAGs")
    config["subroutine"] = "alignment"
    config["cutoff"] = os.path.join("tax-cutoffs.yaml")
    config["output"] = os.path.join(base_dir, "test_out")
    config["database"] = test_reference
    config["download_reference"] = 0
    config["column"] = "SOURCE_ID"
    config["ref_fasta"] = "reference.pep.fa"
    config["protein_map"] = "prot-map.json"
    config["tax_table"] = "tax-table.txt"
    
    config_path = os.path.join(base_dir, 'test_configs')
    os.system("mkdir -p " + config_path)
    config_file = os.path.join(config_path, 'curr_config_alignment.yaml')
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
        
    eukulele(config=config_file)
    outprefix = config["output"].split("/")[-1]
    assert os.path.isfile(os.path.join(config["output"], "taxonomy_counts", outprefix + "_all_species_counts.csv"))
    
#def test_alignment_blast():
#    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
#    base_config = os.path.join(os.path.dirname(__file__), '..', 'aux_data', 'config.yaml')
#    with open(base_config) as f:
#        config = yaml.load(f, Loader=yaml.FullLoader)
#       
#    config["mets_or_mags"] = "mags"
#    config["reference"] = os.path.join(base_dir, test_reference, "sample_ref_MAG")
#    config["samples"] = os.path.join(base_dir, test_reference, "samples_MAGs")
#    config["subroutine"] = "all"
#    config["alignment_choice"] = "blast"
#    config["cutoff"] = os.path.join("tax-cutoffs.yaml")
#    config["output"] = os.path.join(base_dir, "test_out")
#    config["database"] = test_reference
#    config["download_reference"] = 0
#    config["column"] = "SOURCE_ID"
#    config["ref_fasta"] = "reference.pep.fa"
#    config["protein_map"] = "prot-map.json"
#    config["tax_table"] = "tax-table.txt"
#    
#    config_path = os.path.join(base_dir, 'test_configs')
#    os.system("mkdir -p " + config_path)
#    config_file = os.path.join(config_path, 'curr_config_alignment.yaml')
#    with open(config_file, 'w') as f:
#        yaml.dump(config, f)
#        
#    eukulele(config=config_file)
#    outprefix = config["output"].split("/")[-1]
#    assert os.path.isfile(os.path.join(config["output"], "taxonomy_counts", outprefix + "_all_species_counts.csv"))

def test_setup_blast():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    base_config = os.path.join(os.path.dirname(__file__), '..', 'aux_data', 'config.yaml')
    with open(base_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
       
    config["mets_or_mags"] = "mags"
    config["reference"] = os.path.join(base_dir, test_reference, "sample_ref_MAG")
    config["samples"] = os.path.join(base_dir, test_reference, "samples_MAGs")
    config["subroutine"] = "setup"
    config["alignment_choice"] = "blast"
    config["cutoff"] = os.path.join("tax-cutoffs.yaml")
    config["output"] = os.path.join(base_dir, "test_out")
    config["database"] = test_reference
    config["download_reference"] = 0
    config["column"] = "SOURCE_ID"
    config["ref_fasta"] = "reference.pep.fa"
    config["protein_map"] = "prot-map.json"
    config["tax_table"] = "tax-table.txt"
    
    config_path = os.path.join(base_dir, 'test_configs')
    os.system("mkdir -p " + config_path)
    config_file = os.path.join(config_path, 'curr_config_alignment.yaml')
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
        
    eukulele(config=config_file)
    outprefix = config["output"].split("/")[-1]
    assert os.path.isdir(os.path.join(config["reference"], "blast"))
    
def test_alignment_blast():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    base_config = os.path.join(os.path.dirname(__file__), '..', 'aux_data', 'config.yaml')
    with open(base_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
       
    config["mets_or_mags"] = "mags"
    config["reference"] = os.path.join(base_dir, test_reference, "sample_ref_MAG")
    config["samples"] = os.path.join(base_dir, test_reference, "samples_MAGs")
    config["subroutine"] = "alignment"
    config["alignment_choice"] = "blast"
    config["cutoff"] = os.path.join("tax-cutoffs.yaml")
    config["output"] = os.path.join(base_dir, "test_out")
    config["database"] = test_reference
    config["download_reference"] = 0
    config["column"] = "SOURCE_ID"
    config["ref_fasta"] = "reference.pep.fa"
    config["protein_map"] = "prot-map.json"
    config["tax_table"] = "tax-table.txt"
    
    config_path = os.path.join(base_dir, 'test_configs')
    os.system("mkdir -p " + config_path)
    config_file = os.path.join(config_path, 'curr_config_alignment.yaml')
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
        
    eukulele(config=config_file)
    outprefix = config["output"].split("/")[-1]
    assert os.path.isfile(os.path.join(config["output"], "taxonomy_counts", outprefix + "_all_species_counts.csv"))
    
def test_busco():
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'aux_data')
    base_config = os.path.join(os.path.dirname(__file__), '..', 'aux_data', 'config.yaml')
    with open(base_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
       
    config["mets_or_mags"] = "mags"
    config["reference"] = os.path.join(base_dir, test_reference, "sample_ref")
    config["samples"] = os.path.join(base_dir, test_reference, "samples_MAGs")
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
        
    config["subroutine"] = "alignment"
    eukulele(config=config_file)
        
    config["subroutine"] = "busco"
    eukulele(config=config_file)
    samplenames = [curr.split(".")[0] for curr in os.listdir(config["samples"])]
    busco_out = os.path.join(config["output"], "busco_assessment", samplenames[0], "individual", 
                             "summary_" + samplenames[0] + ".tsv")
    assert os.path.isfile(busco_out)
  