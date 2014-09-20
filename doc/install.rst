Installation Instructions
=========================

Setup Git
----------

1. First lets set up a hostname for accessing the git repository so things are
   a little easier. Edit the ``.ssh/config`` file in your home directory::

      > vi $HOME/.ssh/config

   Add the following text::
					
    host hpcgit
        hostname git.hpc.ufl.edu
        user mcintyre

   Make sure that you use 'mcintyre' as the user, we all use the same user name.
   Save and exit the file.


2. Before accessing the git repository you need to have your ``id_rsa.pub`` key
   added to it. Contact Alison, Justin, or Alex on how to do this.


3. Now lets setup a local developmental environment. In your home directory,
   create a new directory called ``devel``::

    > mkdir $HOME/devel


3. Lets clone the python.git repository into the ``devel`` folder. Type into
   the command line::
	
    > cd $HOME/devel
    > git clone hpcgit:python.git python.git


Congrats you now have a local copy of the McIntyre lab ``python.git`` repository.

Setup Python Environment
-------------------------

Now lets set up your local environment to be able to use the scripts and
libraries in the python.git folder.

1. First wee need to create some environmental variables. In your home
   directory edit the ``.profile`` file.::

    > vi .profile

   In your text editor, add the following text to the bottom of the file::
					
    export PATH=$PATH:$HOME/devel/python.git
    export PYTHONPATH=$PYTHONPATH:$HOME/devel/python.git
					
.. note::

    You also need to do this step on the HPC, except the file you need to edit
    is ``$HOME/.bash_profile``. Add the following::

        export PATH=$PATH:$HOME/bin:/scratch/lfs/mcintyre/python.git
        export PYTHONPATH=$PYTHONPATH:/scratch/lfs/mcintyre/python.git


2. After you have made these changes you need to logout of your computer and
   log back in. 
   

Congratulations, you now have the python.git repository installed. To test things out, start a python interpreter and type ``import mclib``, if you do not get an error than everything is working.
