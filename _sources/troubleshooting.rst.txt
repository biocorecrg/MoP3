.. _home-page-troubleshooting:

*****************
Troubleshooting
*****************

.. autosummary::
   :toctree: generated

Demultiplexing with Guppy
================================
Error:

.. code-block:: console

  Init time: 2515 ms
  
  0%   10   20   30   40   50   60   70   80   90   100%
  |----|----|----|----|----|----|----|----|----|----|
  ***************************************************
  Caller time: 32995 ms, Samples called: 159196780, samples/s: 4.82488e+06
  Finishing up any open output files.
  Basecalling completed successfully.

  Command error:
  rm: cannot remove '*_out/*/*.fastq': No such file or directory

Solution:
Check your barcode kit!. You must indicate the **--barcode_kits** in the tool_opts file at the row **demultiplexing	guppy**. Example:

.. code-block:: console

  demultiplexing guppy  "--flowcell FLO-MIN114 --kit SQK-LSK114 --barcode_kits SQK-NBD114-24"


Other...
================
