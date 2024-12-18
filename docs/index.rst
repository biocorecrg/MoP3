.. _home-page-index:

*******************
Welcome to the documentation of Master Of Pores 3
*******************


.. autosummary::
   :toctree: generated

.. image:: ../img/goku3.png
  :width: 600

Master of Pores is a pipeline written in Nextflow DSL2 for the analysis of Nanopore data. It can handle reads from direct RNAseq, cDNAseq, DNAseq etc.

The pipeline is composed by four modules:
   - mop_preprocess: preprocessing
   - mop_mod: detecting chemical modifications. It reads the output directly from mop_preprocess
   - mop_tail: estimating polyA tail size. It reads the output directly from mop_preprocess
   - mop_consensus: it generates a consensus from the predictions from mop_mod. It reads the output directly from mop_mod

.. MoP2 documentation master file, created by
   Luca Cozzuto.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Contents:

.. toctree::
   :maxdepth: 1

   about
   install
   mop_preprocess
   mop_mod
   mop_consensus
   mop_tail
   reporting
   awsbatch
   benchmark
   changelog
   ci
   troubleshooting
