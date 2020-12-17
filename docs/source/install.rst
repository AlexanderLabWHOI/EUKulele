Installation
============

Installing with ``conda``
-------------------------

:code:`EUKulele` may also be downloaded as a :code:`conda` package, which will eventually become the easiest option, as `conda` automatically installs all dependencies for the user. The package can currently be downloaded via::

    conda install -c akrinos -c bioconda -c conda-forge EUKulele
    
Eventually, the :code:`-c akrinos` designation will be replaced, when ``EUKulele`` becomes a ``bioconda`` package, at which point it may be installed from that channel (and this documentation will be updated). 

To ensure you have the most recent version of ``EUKulele``, you can run ``conda update EUKulele``, or you can check the ``EUKulele`` landing page on Anaconda Cloud to check the most recent version number.

If you find that the ``conda`` install is running slowly, you may also, alternatively, follow the steps below to install with ``mamba``.

Installing with ``mamba``
------------------------

:code:`EUKulele` can be installed more rapidly remotely using ``mamba``, a parallelized fast installer that is obtainable via ``conda``. In order to install via ``mamba``, the following steps are necessary. 

- Create a ``conda`` environment
- Install ``mamba``, using::
    
    conda install mamba -c conda-forge

- Install the package as before, but using ``mamba`` as the installer, rather than ``conda``::

    mamba install -c akrinos -c bioconda -c conda-forge EUKulele

This should provide a considerable speedup as compared to ``conda``.

Installing with Pip
-------------------

There are few different ways to install :code:`EUKulele`, but the easiest is::

    pip install EUKulele
    
There may some Python dependencies that are not immediately downloaded when :code:`EUKulele` is pip-installed. If this is the case, the offending dependencies may be installed individually.

The external dependencies can be installed individual, or with conda using::
   
    conda create -n EUKulele
    conda activate EUKulele
    conda install -c bioconda -c conda-forge blast busco=4.0.6 diamond transdecoder
    
As the dependencies external to PyPI are:

- BLAST
- BUSCO
- Diamond
- TransDecoder (if using metatranscriptome samples and the ``--use_transdecoder`` flag)

Cloning the Development Code from GitHub
----------------------------------------

In addition, you can clone :code:`EUKulele` from GitHub (the current development version) using::

    git clone https://github.com/AlexanderLabWHOI/EUKulele
    
And then invoke :code:`EUKulele` either by executing the script :code:`bin/EUKulele` directly, or by installing :code:`EUKulele` as a local package, by calling the following from the :code:`EUKulele` directory::

    python3 -m pip install -e . --user
    python setup.py install  --user
    
Again, in this case, external dependencies may be installed via::
   
    conda create -n EUKulele
    conda activate EUKulele
    conda install -c bioconda -c conda-forge blast busco=4.0.6 diamond transdecoder
    
Or individually for each software.

Invocation
==========

If installed either with pip or ``conda``, ``EUKulele`` can be invoked via::

    EUKulele <arguments>
    
Where the minimal command would be::

    EUKulele --mets_or_mags <choice of data type> --sample_dir <where samples are located>
    
In which case ``EUKulele`` would be run with mostly parameter defaults and using the MMETSP database, by default.

Less common use cases
---------------------

``EUKulele`` may also be run as a module within Python. Include the phrase ``import EUKulele`` in the header of a Python file. Then, you may execute ``EUKulele`` using::

    EUKulele.eukulele(config=config_file)

in the case that you have a configuration file to specify (replace the ``config_file`` variable with this path), or with::

    EUKulele.eukulele(string_arguments=string_of_arguments)

where ``string_of_arguments`` is a string containing the non-default `EUKulele` options you wish to specify.
