mikan
======

[![Travis](https://travis-ci.org/takayasaito/mikan.svg?branch=master)](https://travis-ci.org/takayasaito/mikan)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/takayasaito/mikan?branch=master&svg=true)](https://ci.appveyor.com/project/takayasaito/mikan)

`mikan` provides an ensemble approach to combine six different miRNA target prediction algorithms.

Algorithms
----------

-   TargetScan
-   PITA
-   RNAhybrid
-   miRanda
-   TargetMiner
-   Two-step SVM

Installation
------------

`mikan` bundles all dependent header libraries and uses [CMake](https://cmake.org/) to compile the source code.

Example on Linux/Unix/OS X
```
mkdir build
cd build
cmake -D CMAKE_INSTALL_PREFIX=/path/to/install -D CMAKE_BUILD_TYPE=Release ../src
make && make install 
```

Programs
--------

The current version of `mikan` provides an ensemble program as well as individual programs for each algorithm. .

| Program              | Description                                                |
|:---------------------|:-----------------------------------------------------------|
| mikan                | mikan ensemble algorithm                                   |
| mkconf               | Show default configuration values to stdout                |
| mk-targetscan        | mikan version of TargetScan implementation                 |
| mk-pita              | mikan version of PITA implementation                       |
| mk-rhahybrid         | mikan version of RNAhybrid implementation                  |
| mk-miranda           | mikan version of miRanda implementation                    |
| mk-targetminer       | mikan version of TargetMiner implementation                |
| mk-twostep-svm       | mikan version of Two-step SVM implementation               |

`mikan` and all `mk-*` programs require four command line arguments. 

| Input                   | Description                                             |
|:------------------------|:--------------------------------------------------------|
| miRNA fasta file name   | micorRNA sequences in Fasta format                      |
| mRNA fasta file name    | mRNA sequences in Fasta foramt                          |
| Output file name 1      | Site level scores are written in this file              |
| Output file name 2      | RNA level scores are written in this file               |

Usage
-----

```
mikan /path/to/mirna.fasta /path/to/mrna.fasta /path/to/output1.txt /path/to/output2.txt
```
