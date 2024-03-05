.. _home-page-install:

**************
How to install
**************

.. autosummary::
   :toctree: generated

Please install nextflow `Nextflow <https://www.nextflow.io/>`_ and either `Singularity <https://sylabs.io/>`_ or `Docker <https://www.docker.com/>`_ before.

For installing Nextflow you need a POSIX compatible system (Linux, OS X, etc). It requires Bash 3.2 (or later) and Java 11 (or later, up to 17). Windows system is supported through WSL. For the installation of Nextflow just run:

.. code-block:: console

  curl -s https://get.nextflow.io | bash

To install the pipeline you need to download the repo:

.. code-block:: console

   git clone --depth 1 --recurse-submodules https://github.com/biocorecrg/MoP3.git


You can use **INSTALL.sh** and the version of Guppy you want to download. 

.. note::
  
  Please consider that the support of VBZ compression of fast5 started with version 3.4.X. 


.. code-block:: console
  
  cd MoP3; bash INSTALL.sh 6.0.1

or for installing the default 3.4.5

.. code-block:: console

  cd MoP3; bash INSTALL.sh

Guppy custom models for RNA basecalling will be downloaded from our repository https://biocore.crg.eu/public/mop3_pub/models.tar and placed automatically within the right path inside the pipeline.

You can install different versions of Guppy but only one will be run during the pipeline execution. For switching among them you need to run INSTALL.sh with the version you prefer.

Testing
============

.. code-block:: console

  cd mop_preprocess

  nextflow run mop_preprocess.nf -params-file params.f5.yaml -with-singularity -bg -profile local > log

.. tip::

  You can replace ```-with-singularity``` with ```-with-docker``` if you want to use the docker engine.



