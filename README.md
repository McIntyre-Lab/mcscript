# McIntyre Lab Scripts
This is a collection of scripts that we use for handling various types of genomic data. Additional information can be found at (). 

## Installation
These scripts were created and tested on several Linux environments, mileage may vary on Mac OSX and Windows. The following programs are required:

git

python >=2.7.3

Various python modules may be necessary depending on the script that is being used. Please see the import lines at the top of the script to determine specific dependencies.

To download:

```Shell
git clone https://github.com/McIntyre-Lab/mcscript.git mcscript

git submodule init
git submodule update
```

The `git submodule` statements will install the mclib library which is required by all of the scripts.

## Usage
All of these scripts are designed to be run from the shell or in a Bash script. For example to run `fasta2bed.py` simply type:

```Shell
python fasta2bed.py -f input_file.fa -o output_file.bed
```

Most of the scripts have help documentation that describes the different command line options. 

```Shell
python fasta2bed.py -h
```
