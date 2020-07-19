from setuptools import setup, find_packages
import os
import numpy.random#.common
#import numpy.random.bounded_integers
#import numpy.random.entropy

os.system("export PY_IGNORE_IMPORTMISMATCH=1")
version = open('scripts/VERSION').read().strip()

setup(
    name="EUKulele",
    version=version,
    #distclass=distutils.command.bdist_conda.CondaDistribution,
    conda_buildnum=1,
    url="https://github.com/AlexanderLabWHOI/EUKulele",
    author="Arianna Krinos",
    author_email="akrinos@mit.edu",
    packages=find_packages(exclude=("euk-env",)),
    license="MIT",
    include=['code'],
    setup_requires=['pytest-runner'],
    test_suite='tests',#'nose.collector',
    tests_require=['pytest','nose'],
    install_requires=['setuptools','conda','hypothesis',\
        'pandas','numpy','matplotlib','argparse','seaborn',\
        'multiprocess','chardet',\
        'joblib','ujson','pyyaml'],
    #packages=['EUKulele','EUKulele.tests','EUKulele.scripts'],
)