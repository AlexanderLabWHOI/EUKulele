'''
Build and deploy the package.
'''

import os
from setuptools import setup, find_packages

os.system("export PY_IGNORE_IMPORTMISMATCH=1")
os.system("cp VERSION src/EUKulele/static")
version = open('VERSION').read().strip()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name="EUKulele",
    version=version,
    description = ("A package to make the process of taxonomically classifying "
                   "microbial eukaryotes easier."),
    keywords = "eukaryote taxonomy classification",
    conda_buildnum=1,
    url="https://github.com/AlexanderLabWHOI/EUKulele",
    author="Arianna Krinos",
    author_email="akrinos@mit.edu",
    packages=['EUKulele','scripts'],
    package_dir={'EUKulele': 'src/EUKulele'},
    scripts=['bin/EUKulele','scripts/create_protein_table.py',
             'scripts/download_database.sh',\
             'scripts/concatenate_busco.sh','scripts/configure_busco.sh',
             'scripts/run_busco.sh',\
             'scripts/install_dependencies.sh','scripts/after_job.sh',
             'scripts/coordinate_batch.sh'],
    license="MIT",
    include=['code'],
    include_package_data = True,
    setup_requires=['pytest-runner'],
    test_suite='tests',
    tests_require=['pytest'],
    install_requires=required,
    zip_safe=False,
)
