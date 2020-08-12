import os
import sys

mmetsp_dir = "/vortexfs1/omics/alexander/data/mmetsp/reference_dir"
newpath = "aux_data/mmetsp/reference-pep-trunc.pep.faa"
os.system("reformat.sh in=" + os.path.join(mmetsp_dir, "reference-pep.fa") + " out=" + newpath + \
              " samplerate=0.0005 overwrite=true ignorejunk")

os.system("mkdir -p aux_data/mmetsp/samples_MAGs")
for r in range(10):
    currpath = "aux_data/mmetsp/samples_MAGs/sample_" + str(r) + ".faa"
    os.system("reformat.sh in=" + newpath + " out=" + currpath + \
                  " samplerate=0.10 overwrite=true ignorejunk")