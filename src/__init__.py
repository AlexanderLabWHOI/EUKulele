import os

os.system("export PY_IGNORE_IMPORTMISMATCH=1")
from EUKulele.EUKulele_config import *
from EUKulele.EUKulele_main import *
#rc = os.system("which conda")
#if rc == 0:
#    print("Conda found on system; conda environment will be automatically created and activated.")
#    os.system("conda env create -f EUKulele-env.yaml --force")
#    os.system("conda activate EUKulele")
#else:
#    print("Conda environment not successfully created. Before running the software, you should make sure non-Python dependencies are installed. BLAST and/or DIAMOND (depending on alignment choice), and TransDecoder must be installed.")
