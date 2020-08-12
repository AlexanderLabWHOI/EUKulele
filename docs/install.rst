There are few different ways to install :code:`EUKulele`, but the easiest is::

    pip install EUKulele
    
There may some dependencies that are not immediately downloaded when :code:`EUKulele` is pip-installed. If this is the case, the offending dependencies may be installed individually.

:code:`EUKulele` may also be downloaded as a :code:`conda` package, which will eventually become the easiest option, as `conda` automatically installs all dependencies for the user. The package can currently be downloaded via::

    conda install -c akrinos -c bioconda -c conda-forge EUKulele
    
But eventually the :code:`-c akrinos` designation will no longer be required. 

In addition, you can clone :code:`EUKulele` from GitHub using::

    git clone https://github.com/AlexanderLabWHOI/EUKulele
    
And then invoke :code:`EUKulele` either by executing the script :code:`bin/EUKulele` directly, or by installing :code:`EUKulele` as a local package, by calling the following from the :code:`EUKulele` directory::

    python3 -m pip install -e . --user
    python setup.py install  --user