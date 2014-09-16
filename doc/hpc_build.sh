#!/bin/bash - 
#===============================================================================
#         USAGE: ./hpc_build.sh 
# 
#   DESCRIPTION: Build the sphinx documentation in the mcpublic folder.
# 
#===============================================================================

module load python/2.7.6
SOURCE=/scratch/lfs/mcintyre/python.git/doc
TARGET=/bio/mcintyre/mcpublic/mcpython

sphinx-build -E -b html $SOURCE $TARGET
