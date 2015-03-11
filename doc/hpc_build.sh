#!/bin/bash - 
#===============================================================================
#         USAGE: ./hpc_build.sh 
# 
#   DESCRIPTION: Build the sphinx documentation in the mcpublic folder.
# 
#===============================================================================

# Only run this on the hpc
if [[ $HOSTNAME == *ufhpc* ]]; then

    module load python/2.7.6

    SOURCE=$(dirname $0)
    TARGET=/bio/mcintyre/mcpublic/mcpython/mcscript


    # Move to doc folder
    cd $SOURCE

    # Prepare target folder
    if [ ! -e $TARGET ]; then
        mkdir -p $TARGET
    else
        rm -rf $TARGET/*
    fi

    # Make sure build folder is clean
    make clean

    # Build html docs
    make html

    # Move to public
    rsync -av $SOURCE/_build/html/* $TARGET

    find $TARGET -type d -exec chmod 775 {} \;
    find $TARGET -type f -exec chmod 664 {} \;

    # Clean up build folder
    make clean

else
    echo "To run this script you need to be on the hpc."
fi
