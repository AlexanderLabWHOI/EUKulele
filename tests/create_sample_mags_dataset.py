import os

MMETSP_DIR = "/vortexfs1/omics/alexander/data/mmetsp/reference_dir"
NEW_PATH = "aux_data/mmetsp/reference-pep-trunc.pep.faa"
os.system("reformat.sh in=" + os.path.join(MMETSP_DIR, "reference-pep.fa") + \
          " out=" + NEW_PATH + \
          " samplerate=0.0005 overwrite=true ignorejunk")

os.system("mkdir -p aux_data/mmetsp/samples_MAGs")
for r in range(10):
    CURR_PATH = "aux_data/mmetsp/samples_MAGs/sample_" + str(r) + ".faa"
    os.system("reformat.sh in=" + NEW_PATH + " out=" + CURR_PATH + \
                  " samplerate=0.10 overwrite=true ignorejunk")
