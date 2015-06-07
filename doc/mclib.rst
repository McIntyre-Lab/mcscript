McIntyre Python Library
=======================

This is a collection of functions and classes that are generally useful for a
variety of bioinformatic processes.

BAM Files
---------
This is a basic class for dealing with bam files.

.. automodule:: mclib.bam
   :members:

BED Files
---------
This is a basic class for dealing with bed files.

.. automodule:: mclib.bed
   :members:

GFF Files
---------
This is a basic class for dealing with gff files.

.. note::
   Will only work with files from FlyBase.

.. automodule:: mclib.gff
   :members:

VCF Files
---------
This is a basic class for dealing with vcf files.

.. note::
   This library is called ``vcf2`` because it is a wrapper for the ``pyvcf`` wich is called ``vcf``.

.. automodule:: mclib.vcf2
   :members:

Wiggle Plot Helper Functions
----------------------------
This is a set of functions for helping with creating wiggle type plots.

.. automodule:: mclib.wiggle
   :members:

Logging Helper Functions
------------------------
This is a set of functions for helping with consitant logging accross all
McIntyre Lab scripts.

.. automodule:: mclib.logger
   :members:
