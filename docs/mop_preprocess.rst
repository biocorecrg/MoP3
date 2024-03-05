.. _home-page-moprepr:

*******************
MOP_PREPROCESS
*******************

.. autosummary::
   :toctree: generated

This module takes as input the raw fast5 reads - single or multi - and it produces several outputs (basecalled fast5, sequences in fastq format, aligned reads in BAM format etc). The pre-processing module can perform base-calling, demultiplexing (optional), filtering, quality control, mapping to a reference (either a genome or a transcriptome), feature counting, discovery of novel transcripts, and it generates a final report with the performance and results of each of the steps performed. It automatically detects the kind of input fast5 file (single or multi-sequence). In theory, it can also support the new pod5 format but it won't output basecalled fastq useful for the other pipelines. The basecalling can be performed with guppy or dorado and the demultiplexing with either guppy, deeplexicon or seqtagger. Basecalled fastq and Fast5 files can be demultiplexed as well. 

  

Input Parameters
======================

The input parameters are stored in yaml files like the one represented here:

.. literalinclude:: ../mop_preprocess/params.f5.demrna.yaml
   :language: yaml

You can change them by editing this file or using the command line as explained in the next section. 

How to run the pipeline
=============================

Before launching the pipeline, the user can decide which containers to use: either docker or singularity **[-with-docker / -with-singularity]**.

Then, to launch the pipeline, please use the following command by specifying the path of the yaml parameter file:

.. code-block:: console

   nextflow run mop_preprocess.nf -with-singularity -params-file params.yaml > log.txt


You can run the pipeline in the background by adding the nextflow parameter **-bg**:

.. code-block:: console

   nextflow run mop_preprocess.nf -params-file params.yaml -with-singularity -bg > log.txt

You can change the parameters either by changing the yaml config file or by feeding the parameters via command line:

.. code-block:: console

   nextflow run mop_preprocess.nf -with-singularity -params-file params.yaml -bg --output test2 > log.txt


You can specify a different working directory with temporary files:

.. code-block:: console

   nextflow run mop_preprocess.nf -with-singularity -params-file params.yaml -bg -w /path/working_directory > log.txt

You can use different profiles for running the pipeline in different environments. We have one set up for HPC using the SGE scheduler:

.. code-block:: console

   nextflow run mop_preprocess.nf -with-singularity -bg -params-file params.yaml -w /path/working_directory -profile cluster > log.txt

One for HPC using the slurm scheduler

.. code-block:: console

   nextflow run mop_preprocess.nf -with-singularity -bg -params-file params.yaml -w /path/working_directory -profile slurm > log.txt

One for emulating the new M1 processor for Apple:

.. code-block:: console

   nextflow run mop_preprocess.nf -with-singularity -bg -params-file params.yaml -w /path/working_directory -profile m1mac > log.txt


or you can run the pipeline locally:

.. code-block:: console

   nextflow run mop_preprocess.nf -with-singularity -bg -params-file params.yaml -w /path/working_directory -profile local > log.txt


.. note::
 
   * In case of errors you can troubleshoot by seeing the log file (log.txt) for more details. Furthermore, if more information is needed, you can also go to the intermediate directory indicated in the log and check both the `.command.log` and `.command.err` files. 

.. tip::

   Once the error has been solved or if you change a specific parameter, you can resume the execution with the **Netxtlow** parameter **- resume** (only one dash!). If there is an error, the pipeline will resume from the process that had the error and proceed with the rest.  If a parameter is changed, only processes affected by this parameter will be re-run. 


.. code-block:: console
   nextflow run mop_preprocess.nf -with-singularity -params-file params.yaml -bg -resume > log_resumed.txt

   To check whether the pipeline has been resumed properly, please check the log file. If previous correctly executed process are found as   *Cached*, resume worked!

.. code-block:: console

   ...

   [warm up] executor > crg
   [e8/2e64bd] Cached process > baseCalling (RNA081120181_1)
   [b2/21f680] Cached process > QC (RNA081120181_1)
   [c8/3f5d17] Cached process > mapping (RNA081120181_1)
   ...


.. note::
   To resume the execution, temporary files generated previously by the pipeline must be kept. Otherwise, the pipeline will re-start from the beginning. 

Results
====================

Several folders are created by the pipeline within the output directory specified by the **output** parameter:


* **fast5_files**: Contains the basecalled multifast5 files. Each batch contains 4000 sequences. 
* **fastq_files**: Contains one or, in case of demultiplexing, more fastq files.
* **QC_files**: Contains each single QC produced by the pipeline.
* **alignment**: Contains the bam file(s).
* **cram_files**: Contains the cram file(s).
* **counts**: Contains read counts per gene / transcript if counting was performed.
* **assigned**: Contains assignment of each read to a given gene / transcript if counting was performed.
* **report**: Contains the final multiqc report. 
* **assembly**: It contains assembled transcripts.

.. note::
   MOP3 will automatically detect the version of guppy and modify the parameters accordingly. You don't need to add any extra parameter as in MOP2.
   
.. tip::
   You can pass via parameter a custom NAME_tool_opt.tsv file with custom options for each step to customise  to disable the qscore filtering. Some custom files are already available in this package, like drna_tool_splice_m6A_guppy6_opt.tsv
   
   





