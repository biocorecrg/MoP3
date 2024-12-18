# MoP3 - Master of Pores 3
[![Docker Build Status](https://img.shields.io/docker/automated/biocorecrg/nanopore.svg)](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/nanopore/builds)
[![mop2-CI](https://github.com/biocorecrg/MoP3/actions/workflows/build.yml/badge.svg)](https://github.com/biocorecrg/MoP3/actions/workflows/build.yml)[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow version](https://img.shields.io/badge/Nextflow-21.04.1-brightgreen)](https://www.nextflow.io/)
[![Nextflow DSL2](https://img.shields.io/badge/Nextflow-DSL2-brightgreen)](https://www.nextflow.io/)
[![Singularity version](https://img.shields.io/badge/Singularity-v3.2.1-green.svg)](https://www.sylabs.io/)
[![Docker version](https://img.shields.io/badge/Docker-v20.10.8-blue)](https://www.docker.com/)

<br/>



<img align="right" href="https://biocore.crg.eu/" src="https://raw.githubusercontent.com/CRG-CNAG/BioCoreMiscOpen/master/logo/biocore-logo_small.png" />


Master of Pores is a pipeline written in Nextflow DSL2 for the analysis of Nanopore data.
<br/>

It can handle reads from direct RNAseq, cDNAseq, DNAseq etc.

<br/>


![MOP3](https://github.com/biocorecrg/MoP3/blob/master/img/goku3.png?raw=true)

The name is inspired by the Metallica's [Master Of Puppets](https://www.youtube.com/watch?v=S7blkui3nQc)

## Install
Please install nextflow and singularity or docker before.

Then download the repo:

```
git clone --depth 1 --recurse-submodules https://github.com/biocorecrg/MoP3.git
```

You can use INSTALL.sh to download the version 3.4.5 of guppy or you can replace it with the version you prefer. Please consider that the support of VBZ compression of fast5 started with version 3.4.X.

```
cd MoP3; sh INSTALL.sh
```

## Testing
You can replace ```-with-singularity``` with ```-with-docker``` if you want to use the docker engine.

```
cd mop_preprocess
nextflow run mop_preprocess.nf -with-singularity -bg -profile local -params-file params.yaml > log
```

## Upgrading
To upgrade the tool you can type:

```
git pull --recurse-submodules
```

## Documentation
The documentation is available at [https://biocorecrg.github.io/MoP3/](https://biocorecrg.github.io/MoP3/)

## Contact
Please open an [issue](https://github.com/biocorecrg/MOP2/issues) if you encounter any issues/troubles.
However, please go over the previous issues (including closed issues) before opening a new issue, as your same exact question might have been already answered previously. Thank you!


## Reference
If you use this tool, please cite our papers:

["Nanopore Direct RNA Sequencing Data Processing and Analysis Using MasterOfPores"
Cozzuto L, Delgado-Tejedor A, Hermoso Pulido T, Novoa EM, Ponomarenko J. *N. Methods Mol Biol. 2023*;2624:185-205. doi: 10.1007/978-1-0716-2962-8_13.](https://link.springer.com/protocol/10.1007/978-1-0716-2962-8_13)

["MasterOfPores: A Workflow for the Analysis of Oxford Nanopore Direct RNA Sequencing Datasets"
Luca Cozzuto, Huanle Liu, Leszek P. Pryszcz, Toni Hermoso Pulido, Anna Delgado-Tejedor, Julia Ponomarenko, Eva Maria Novoa.
*Front. Genet., 17 March 2020.* https://doi.org/10.3389/fgene.2020.00211](https://www.frontiersin.org/articles/10.3389/fgene.2020.00211/full)
