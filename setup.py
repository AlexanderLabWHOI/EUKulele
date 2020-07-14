from setuptools import setup, find_packages

version = open('scripts/VERSION').read().strip()

setup(
    name="EUKulele",
    version=version,
    #distclass=distutils.command.bdist_conda.CondaDistribution,
    conda_buildnum=1,
    url="https://github.com/AlexanderLabWHOI/EUKulele",
    author="Arianna Krinos",
    author_email="akrinos@mit.edu",
    license="MIT",
    include=['code'],
    install_requires=['setuptools','conda',\
        'pandas','numpy','matplotlib','argparse',\
        'yaml','multiprocessing','subprocess','chardet',\
        'shutil','glob','joblib','json','pyyaml'],
    packages=['EUKulele','tests','EUKulele.scripts'],
)