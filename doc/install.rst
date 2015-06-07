McIntyre Library Installation Instructions
==========================================

Clone Git Repository
--------------------

1. Decide where you want to clone the McIntyre Library repository. For example
   we will create a development folder in the user's home folder.::

      $ mkdir $HOME/devel

   Then clone the git repository::

        $ git clone https://github.com/McIntyre-Lab/mclib.git $HOME/devel/mclib

Congrats you now have a local copy of the McIntyre library ``mclib`` repository.

Setup Python Environment
-------------------------

In order to use the contents of the McIntyre library, you need to tell python
where these files are. To do this you need to set up a local environmental
variable. On Linux there are several places you can create these variables. I
would suggest adding the variable to a file in the user's home directory called
``.profile``. 

.. note::
    ``.profile`` is loaded at login, so you will need to logout and login for changes to take place.

1. Open ``$HOME/.profile`` with your favorite text editor and add the following ::

    export PYTHONPATH=$PYTHONPATH:$HOME/devel

2. Now logout of your computer and log back in. 

3. Test that the McIntyre Library can be found by python. ::
   
   $ python

   >>> import mclib

If there are no errors then you have successfully installed the library.
