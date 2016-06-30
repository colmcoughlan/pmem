The following was taken from the CASA wiki at https://casaguides.nrao.edu/index.php?title=Writing_a_CASA_Task. More information is available in the accompanying pdf.

Building the task

Once you have the XML and Python files, run buildmytasks to compile them into a form that can be imported into CASA. This can be done on the command line (by simply typing "buildmytasks", or within a CASA session:

# In CASA
os.system('buildmytasks')

After buildmytasks has completed, you will see that there are a number of new files, e.g. "newtask.py-e", "newtask_cli.pyc", etc. Also, there is now a file called "mytasks.py", which is the meta-file for all the tasks you might have in this directory. This is what you import into CASA to include these tasks:

# In CASA
execfile('/<path_to_task_directory>/mytasks.py')

Note that if you rebuild a task after having already imported it into CASA, you will need to restart CASA and reimport for the changes to be incorporated.
Running the task

Now that the task has been imported into CASA, it can be run like any other task. For example,

# In CASA
inp newtask
help newtask
myinput = 'Here is some input'
go

That's it -- a new task! 
