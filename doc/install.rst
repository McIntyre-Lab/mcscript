McIntyre Scripts Installation Instructions
==========================================

Install McIntyre Library Repository
-----------------------------------

1. Before installing McIntyre Scripts, you need to install the McIntyre
   Library. Follow these `installation instruction
   <http://bio.rc.ufl.edu/pub/mcintyre/mcpython/mclib/install.html>`_.

Clone Git Repository
--------------------

1. Clone the McIntyre Scripts to the same location as the McIntyre Library repository.::

        $ git clone https://github.com/McIntyre-Lab/mcscript.git $HOME/devel/mcscript

Congrats you now have a local copy of the McIntyre Scripts ``mcscript`` repository.

Setup Path Environment
-------------------------

To make things easier, add this location ``$HOME/devel/mscript`` to your PATH.
This will allow you to run the scripts without having to point to the folder.
You can do this by editing the ``$HOME/.profile``. 

.. note::
    ``.profile`` is loaded at login, so you will need to logout and login for changes to take place.

1. Open ``$HOME/.profile`` with your favorite text editor and add the following ::

    export PATH=$PATH:$HOME/devel/mcscripts

2. Now logout of your computer and log back in. 

3. Test that the McIntyre Library can be found by python. ::
   
   $ bed2fasta.py -h

If you get the help contents then everything is working.

.. note::
    If nothing happens, make sure that permissions are set to executable. Run the following command::

       $ find $HOME/devel/mcscript -type f -iname *.py -exec chmod 775 {} \;
