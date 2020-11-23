import os
from random import sample

mmetsp_dir = "/vortexfs1/omics/alexander/data/mmetsp"
mmetspfiles = os.listdir(mmetsp_dir)
mmetspfiles = [curr for curr in mmetspfiles if "clean.fasta" in curr]

def generate_metatranscriptome(mmetsp_files, outputdir, METindex):
    '''
    Create metatranscriptome files from several available sub-files.
    '''

    samplefiles = [os.path.join(outputdir, "sample_" + str(index_mmetsp) \
                                + ".fasta") for index_mmetsp in \
                   range(len(mmetsp_files))]
    for mmetsp_file_ind in range(len(mmetsp_files)):
        mmetsp_file = mmetsp_files[mmetsp_file_ind]
        os.system("reformat.sh in=" + os.path.join(mmetsp_dir, mmetsp_file) \
                  + " out=" + samplefiles[mmetsp_file_ind] + \
              " samplerate=0.0001 overwrite=true ignorejunk")
        
    metatranscriptome_file = os.path.join(outputdir, "met_sample_" + \
                                          str(METindex) + ".fasta")
    os.system("cat " + " ".join(samplefiles) + " > " + metatranscriptome_file)
    [os.system("rm " + curr_file) for curr_file in samplefiles]

number_samples = 10
outputdir = os.path.join("aux_data","mmetsp","samples_METs_small")
os.system("mkdir -p " + outputdir)

for r in range(number_samples):
    mmetsp_files = sample(mmetspfiles, 20)
    generate_metatranscriptome(mmetsp_files, outputdir, r)
