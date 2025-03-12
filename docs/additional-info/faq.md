# Frequently asked questions

#### I don't see my files in the terminal. What should I do?

There are a few things you can do to troubleshoot. 

1. Check that you are in the correct directory. You can do this by running `pwd` in the terminal. If you are not in the correct directory, you can change to the correct directory by running `cd /path/to/directory`. 
2. Did your command run correctly? If you are running a command that should produce output, make sure that the command ran correctly. If it did not, you may need to troubleshoot the command. You might see an error message at the end of the output that can help you identify the problem. For example if you see "command not found" you may need to activate the right conda environment.
3. You can try to find the file using the `find` command. For example, if you are looking for a file called `sample1.bam` you can run `find ~/data/ -name sample1.bam`. This will search the entire filesystem for the file. If the file is found, the path to the file will be printed to the terminal.