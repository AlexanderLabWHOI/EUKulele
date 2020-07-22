from setuptools import setup, find_packages
import os
import numpy.random
from glob import glob

os.system("export PY_IGNORE_IMPORTMISMATCH=1")
version = open('VERSION').read().strip()

setup(
    name="EUKulele",
    version=version,
    #distclass=distutils.command.bdist_conda.CondaDistribution,
    conda_buildnum=1,
    url="https://github.com/AlexanderLabWHOI/EUKulele",
    author="Arianna Krinos",
    author_email="akrinos@mit.edu",
    #packages=find_packages(exclude=("euk-env",)),
    packages=find_packages('EUKulele'),
    package_dir={'': 'src'},
    #py_modules=[os.path.splitext(os.path.basename(path))[0] for path in glob('src/EUKulele/*.py')],
    #packages=['EUKulele'],
    #package_dir={'EUKulele': 'src'}, #src/EUKulele
    license="MIT",
    include=['code'],
    setup_requires=['pytest-runner'],
    test_suite='tests',
    tests_require=['pytest'],
    scripts=['bin/EUKulele'],
    install_requires=['setuptools','hypothesis',\
        'pandas','numpy','matplotlib','argparse','seaborn',\
        'multiprocess','chardet','biopython',\
        'joblib','ujson','pyyaml'],
    #packages=['EUKulele','EUKulele.tests','EUKulele.scripts'],
)
