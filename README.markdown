
# Variant calling pipeline: exome and WGS

## Overview

This is a basic variant-calling and annotation pipeline developed at the 
Life Sciences Computation Centre, University of Melbourne.

It is based around BWA, GATK and ENSEMBL and was originally designed for human (or similar) data. The master branch is configured for WGS data; there is an exome branch configured for variant calling in exome data.

To run the pipeline you will need Rubra: [https://github.com/bjpop/rubra](https://github.com/bjpop/rubra). Rubra uses the python Ruffus library: [http://www.ruffus.org.uk/](http://www.ruffus.org.uk/).

Usage: 
    
    rubra pipeline.py --config <your_config_file> --style {print,run,touchfiles,flowchart}

More command-line options are described in the Rubra documentation.

## Running on VLSCI's clusters (e.g. merri)

On merri we have a version of Rubra installed into Python 2.7.3, which you can load with

    module load python-gcc/2.7.3

To use the flowchart option you will need graphviz, which you can load with

    module load graphviz
