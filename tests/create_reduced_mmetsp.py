import os
from random

MMETSP_DIR = "/vortexfs1/omics/alexander/data/mmetsp"
mmetspfiles = os.listdir(MMETSP_DIR)
mmetspfiles = [curr for curr in mmetspfiles if "clean.fasta" in curr]

def generate_metatranscriptome(mmetsp_files, outputdir):
    '''
    Create a combined metatranscriptome using the individual files
    from the MMETSP.
    '''

    samplefiles = [os.path.join(outputdir, "sample_" + str(index_mmetsp) + ".fasta") \
                   for index_mmetsp in range(len(mmetsp_files))]
    for mmetsp_file_ind in range(len(mmetsp_files)):
        mmetsp_file = mmetsp_files[mmetsp_file_ind]
        os.system("reformat.sh in=" + os.path.join(MMETSP_DIR, mmetsp_file) + \
                  " out=" + samplefiles[mmetsp_file_ind] + \
              " samplerate=0.1 overwrite=true ignorejunk")

    metatranscriptome_file = os.path.join(outputdir, "reference-pep.pep.fa")
    os.system("cat " + " ".join(samplefiles) + " > " + metatranscriptome_file)
    [os.system("rm " + curr_file) for curr_file in samplefiles]

outputdir = os.path.join("aux_data","mmetsp","sample_ref")
os.system("mkdir -p " + outputdir)

generate_metatranscriptome(mmetspfiles, outputdir)
