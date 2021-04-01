.. include:: links.inc

Installing nimgen
==================


Requirements
^^^^^^^^^^^^

julearn requires the following packages:

* `Python`_ >= 3.8
* `numpy`
* `statsmodels`
* `abagen`

We strongly recommend using virtual environments:

* `venv`_
* `conda env`_

Installing
^^^^^^^^^^

.. _install_development_git:

Local git repository (for developers)
-------------------------------------
First, make sure that you have all the dependencies installed:

.. code-block:: bash

    pip install -U numpy statsmodels abagen

Then, clone `nimgen Github`_ repository in a folder of your choice:

.. code-block:: bash

    git clone https://github.com/juaml/nimgen.git

Finally, install in development mode:

.. code-block:: bash

    cd nimgen
    python setup.py develop

