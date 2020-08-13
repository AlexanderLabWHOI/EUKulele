Installation and Invocation
===========================

Installing with Pip
-------------------

There are few different ways to install :code:`EUKulele`, but the easiest is::

    pip install EUKulele
    
There may some dependencies that are not immediately downloaded when :code:`EUKulele` is pip-installed. If this is the case, the offending dependencies may be installed individually.

Installing with ``conda``
-------------------------

:code:`EUKulele` may also be downloaded as a :code:`conda` package, which will eventually become the easiest option, as `conda` automatically installs all dependencies for the user. The package can currently be downloaded via::

    conda install -c akrinos -c bioconda -c conda-forge EUKulele
    
Eventually, the :code:`-c akrinos` designation will be replaced, when ``EUKulele`` becomes a ``bioconda`` package, at which point it may be installed from that channel (and this documentation will be updated). 

Cloning the Development Code from GitHub
----------------------------------------

In addition, you can clone :code:`EUKulele` from GitHub using::

    git clone https://github.com/AlexanderLabWHOI/EUKulele
    
And then invoke :code:`EUKulele` either by executing the script :code:`bin/EUKulele` directly, or by installing :code:`EUKulele` as a local package, by calling the following from the :code:`EUKulele` directory::

    python3 -m pip install -e . --user
    python setup.py install  --user
    
Invoking ``EUKulele``
---------------------

If installed either with pip or ``conda``, ``EUKulele`` can be invoked via::

    EUKulele <arguments>
    
Where the minimal command would be::

    EUKulele --mets_or_mags <choice of data type> --sample_dir <where samples are located>
    
In which case ``EUKulele`` would be run with mostly parameter defaults and using the MMETSP database, by default.