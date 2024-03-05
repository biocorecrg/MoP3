.. _home-page-about:

*******************
About Master of Pores
*******************

.. autosummary::
   :toctree: generated

Master of Pores 3 is a collection of pipelines written in Nextflow DSL2 for the analysis of Nanopore data. It can handle reads from direct RNAseq, cDNAseq, DNAseq etc.

The software is composed by four pipelines:

   - mop_preprocess: preprocessing of input data. Basecalling, demultiplexing, alignment, read counts, and more!
   - mop_mod: detecting chemical modifications. It reads the output directly from mop_preprocess
   - mop_tail: estimating polyA tail size. It reads the output directly from mop_preprocess 
   - mop_consensus: it generates a consensus from the predictions from mop_mod. It reads the output directly from mop_mod

The name is inspired by Metallica's `Master Of Puppets <https://www.youtube.com/watch?v=S7blkui3nQc>`_

.. image:: ../img/goku3.png
  :width: 600  

This is a joint project between `CRG bioinformatics core <https://biocore.crg.eu/>`_ and `Epitranscriptomics and RNA Dynamics research group <https://public-docs.crg.es/enovoa/public/website/index.html>`_.


Reference
======================

If you use this tool, please cite our papers:

`"Nanopore Direct RNA Sequencing Data Processing and Analysis Using MasterOfPores" <https://link.springer.com/protocol/10.1007/978-1-0716-2962-8_13>`__ Cozzuto L, Delgado-Tejedor A, Hermoso Pulido T, Novoa EM, Ponomarenko J. N. Methods Mol Biol. 2023;2624:185-205. doi: 10.1007/978-1-0716-2962-8_13. 

`"MasterOfPores: A Workflow for the Analysis of Oxford Nanopore Direct RNA Sequencing Datasets" <https://doi.org/10.3389/fgene.2020.00211](https://www.frontiersin.org/articles/10.3389/fgene.2020.00211/full>`_ Luca Cozzuto, Huanle Liu, Leszek P. Pryszcz, Toni Hermoso Pulido, Anna Delgado-Tejedor, Julia Ponomarenko, Eva Maria Novoa. Front. Genet., 17 March 2020.





