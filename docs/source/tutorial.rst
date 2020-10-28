Trying out ``EUKulele``
=======================

A minimal working example for ``EUKulele`` can be tested using the instructions on this page.

You can obtain sample metagenomic data with simulated MAGs from the following Dropbox link, using the following commands on a Linux system::

    wget https://www.dropbox.com/s/l4kvbpqftdad5ib/sample_eukulele.tar.gz?dl=1
    tar xzf sample_eukulele.tar.gz?dl=1

Then ``cd`` into ``sample_EUKulele``. We have noticed a problem with tarring across systems and through Dropbox. So you may be required to run:: 

    rm samples_MAGs/._sample_*

Now, you should have a clean folder containing only metagenomic samples that you can use for test driving ``EUKulele``!

First, let's create an environment to run ``EUKulele`` in.  Create a ``conda`` environment and activate it, using::

    conda env create -f EUKulele-env.yaml
    conda activate EUKulele

This must also be done inside the directory you installed via Dropbox, which contains a ``conda`` configuration file that will make sure that every dependency of ``EUKulele`` that cannot be installed via pip is available on your system.

You should now download ``EUKulele`` via ``pip`` using::

    python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps EUKulele

If any dependency is not satisfied, you can install it manually using ``pip install <requirement> --user``. If you install ``EUKulele`` via ``conda`` instead of via ``pip``, all of the dependencies are installed for you automatically.

Once everything is setup and all dependencies satisfied, from within the ``sample_EUKulele`` directory, run::
    
    EUKulele --config curr_config.yaml

    
Where ``curr_config.yaml`` is a configuration file in the directory you downloaded which contains all of the flags needed for a basic ``EUKulele`` run. Feel free to open this file if you're curious, and compare it to the full list of parameters available to customize ``EUKulele``. 

Then check the folder ``test_out_23July`` in the current directory.
