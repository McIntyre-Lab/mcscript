Installation Instructions
=========================


1. In your home directory, create a new directory called 'devel':
		``> mkdir devel``


2. Next, you will need to edit the .ssh config file in your home directory:
          ``> cd /home/user``
		  ``> vi .ssh/config``
	
		In your text editor, add the following text:
					
            host hpcgit
			        hostname git.hpc.ufl.edu
                    user mcintyre

	Make sure that you use 'mcintyre' as the user, we all use the same user name. Save and exit the file.


3. Go back to your home directory and go into the devel folder you created in Step 1 (/home/user/devel). 
	Type into the command line:
	
            ``> git clone hpcgit:python.git python.git``


4. Go back to your home directory. Now you will edit the .profile file. 
        ``> vi .profile``
			In your text editor, add the following text:
					
                export PATH=$PATH:$HOME/devel/python.git
                export PYTHONPATH=$PYTHONPATH:$HOME/devel/python.git
					
	
5. To set this up on the HPC: First login and stay in the home directory. Edit the .bash_profile file.
        ``> vi .bash_profile``
			In your text editor, add the following text:
					
                export PATH=$PATH:$HOME/bin:/scratch/lfs/mcintyre/python.git
                export PYTHONPATH=$PYTHONPATH:/scratch/lfs/mcintyre/python.git
					
		
