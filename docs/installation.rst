.. include:: links.inc

Installing nimgen
==================


Requirements
^^^^^^^^^^^^

julearn requires the following packages:

* `Python`_ >= 3.8
* ``numpy``
* ``statsmodels``
* ``abagen``
* ``selenium`` (only for webgestalt)

We strongly recommend using virtual environments:

* `venv`_
* `conda env`_


Using webgestalt
^^^^^^^^^^^^^^^^

To use webgestalt, the respective selenium driver must be installed.

* Chrome: https://sites.google.com/a/chromium.org/chromedriver/downloads
* Firefox: https://github.com/mozilla/geckodriver/releases

Just download the respective file, unzip and place in the users PATH (``/bin`` or ``/usr/bin``)

Note for Macos users: 
Open a Terminal, copy the driver to ``/usr/local/bin`` and then execute. An error message will appear. Press close or cancel. Open system preferences, security and press the button to allow for the execution. Go to the terminal, execute again. Another warning will appear. Allow.

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

