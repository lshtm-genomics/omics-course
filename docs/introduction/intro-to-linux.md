# Linux

Most bioinformatic software is written for Unix-type systems and is usually operated from the command line. This means that these programs are intended to be used by typing commands and they don't have a graphical user interface (GUI). There are many reasons for this, but two important ones are that command line programs are quicker to develop and also very easy to automate. To get from raw sequences to final results often requires many steps. Imagine processing thousands of isolates: this would potentially involve hundreds of hours and even more mouse clicks if performed with a GUI. Using a command line tool, a workflow can be defined by chaining a series of commands together in a script which can then be rapidly applied to thousands of isolates without much effort. 

Over the duration of this course you will be learning the skills required to independently undertake bioinformatic projects. By the end of this course you should be comfortable to go all the way from raw data to biological insights.

In order to learn how to use bioinformatic programs we must first become comfortable with how to use the terminal. This may be daunting if you have never used a command line environment before, but hopefully by the end of this practical session you should have a basic understanding which will be reinforced over the next few sessions.

Let's start by opening a terminal by clicking on the icon on the sidebar. Your screen should now look like the one below:

![](../img/linux_1.jpg)

The terminal is very similar to the file explorer in other operating systems. It allows you to see files and interact with them. The only difference is that we have to do this by typing a command into the terminal window. To demonstrate this, let's open up Linux's file explorer by clicking on the icon in the side-bar on the left-hand side of the screen. When you open the file explorer it will go to your home folder and show you all the sub-folders and files present. The home folder is simply the root folder which contains all your files and folders.

You can see a number of sub-folders and files present in the file explorer. Let's see if we can find them in the terminal. To do this, we have to type ls into the terminal and hit enter. You should now see a screen similar to the one below. Now you can see the same files and folders in the terminal. 

![](../img/linux_2.jpg)


In the file explorer double click on the data folder. You will see several new folders which hold the data for the proceeding practicals. Let's try find them using the terminal!

When we double clicked on the data folder we effectively moved from the home folder to the data folder. We will have to do the same using the terminal. To do this, type in `cd data` and hit enter. Now try using ls again to produce the same list of folders that you see in the file explorer.

!!! info 
    Folders are often referred to as directories in Linux. We will use directory from now on, but the two words are synonymous. 


You will notice that text in front of the $ sign changed when you moved from home to data. You should see something similar to what is shown below: 

!!! terminal "Terminal output"
    ```
    user@user-VirtualBox:~/data$ 
    ```

The blue text shows the location where the terminal is currently. In this case we are in `~/data`. Another name for this is the current directory. The `~` character is a special character which symbolises your home directory. So, we can interpret this as "We are currently in the data directory which is in the home directory"

!!! question
    Another way of finding the current location of the terminal is by using pwd which stands for print working directory. Try using this command. Does it produce a similar output to what we have seen above? 

We have now successfully used the terminal and the `ls`, `cd` and `pwd` commands. It is as simple as that! Some commands can be used by themselves such as `ls`, however others such as `cd` need additional information called arguments to do something useful. In the example above, we need to tell `cd` that we wanted to go to the `data` directory, so we typed `data` after `cd` separated by a space. While `cd` only took one argument, some programs take many more. You will see examples of this in a bit.

Now that you have mastered using the terminal you can close the file explorer, this is the last time we will use it. We must now face using the terminal! 

## Useful commands/programs

### head

Change to the tb directory and have a look at the files. Hopefully you will be able to see the **tb.fasta** file. The file is just a very large text file which stores the sequence data of the M. tuberculosis reference genome. We can use `head` to take a look at the first few lines of a file. Let's try it: 

```
cd ~/data/tb
ls
head tb.fasta
```

Hopefully you will see something like this: 

!!! terminal "Terminal output"
    ```
    >Chromosome
    TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACC
    CTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAGCAAAGGGCTTG
    GCTCAATCTCGTCCAGCCATTGACCATCGTCGAGGGGTTTGCTCTGTTATCCGTGCCGAGCAGCTTTGTC
    CAAAACGAAATCGAGCGCCATCTGCGGGCCCCGATTACCGACGCTCTCAGCCGCCGACTCGGACATCAGA
    TCCAACTCGGGGTCCGCATCGCTCCGCCGGCGACCGACGAAGCCGACGACACTACCGTGCCGCCTTCCGA
    AAATCCTGCTACCACATCGCCAGACACCACAACCGACAACGACGAGATTGATGACAGCGCTGCGGCACGG
    GGCGATAACCAGCACAGTTGGCCAAGTTACTTCACCGAGCGCCCGCACAATACCGATTCCGCTACCGCTG
    GCGTAACCAGCCTTAACCGTCGCTACACCTTTGATACGTTCGTTATCGGCGCCTCCAACCGGTTCGCGCA
    CGCCGCCGCCTTGGCGATCGCAGAAGCACCCGCCCGCGCTTACAACCCCCTGTTCATCTGGGGCGAGTCC
    ```

We can see that by default `head` prints out the first 10 lines. We can modify this number by providing an additional parameter. 

```
head -5 tb.fasta
```

This command it will print out the first five lines instead. Here we have used an optional parameter (the program will still function properly without it). Optional parameters are usually given by providing two values, in this case `-n` and `5`. The first value specified what parameter we want to provide (e.g. `-n` = number of lines) and the second value is that which the parameter must take (e.g. `5` will tell the command to print the first five lines). 

### less

Another tool that we can use to view text files is the `less` command. Let's try to view the reference file again with this method. 

```
less tb.fasta
```

This will open an interactive viewer which fills up the screen. You can move up and down the file using the up and down keys on the keyboard. When you are finished viewing, hit the `q` key to quit the viewer. 

### cp

If we need to copy a file we can use the `cp` command. We need to give two arguments: the file we want to copy and the filename of the copy. Let's try it: 

```
cp tb.fasta tb_copy
```

!!! info

    Two things to note here:

    1. Spaces are not allowed in file names! Always use an underscore instead of a space.
    2. We have not given our new file an extension. Extensions are useful as they tell us something about how a file is formatted and what it might contain. If we come back a year later, we might have forgotten that `tb_copy` contained sequence data.


### mv

Let's give our file an extension. To do this we have to rename the file. The command for this is `mv` which stands for move. This command allows you to both move and rename files. Try to add an extension to the file by using the following code: 

```
mv tb_copy tb_copy.fasta
```

You can see that this command also required two parameters:

1. The file we want to rename
2. The new file name

### rm

We now have two files with the same content. Let's remove our copy of the reference by using the `rm` command: 

```
rm tb_copy.fasta
```

!!! danger
    Be very careful with this command! Once you remove a file using rm there is no way of getting it back. Always double check you are removing the right file before hitting enter. 

## Pipes

Imagine a cake factory production line. There are many steps that need to be taken to build the final cake (adding different ingredients). Let's say we have 5 people working to make 100 cakes per day and each person specialises in adding a particular ingredient. We could get 100 bowls and add the flour to each one, then add eggs to each one and repeat this with each ingredient until we have the complete mix however we would be losing a lot of time as at any point there would be only one person working while the others wait for their turn. A more efficient way to do this by installing in a conveyor belt and having the employees sequentially add their ingredients to each bowl. This way everyone is working at the same time and the cakes will be made a lot faster. 

![](../img/linux_4.gif)

The same principals apply to bioinformatic analyses. For example, take a look at the steps required to print the 9th and 10th line of the tb reference fasta file. The easiest way to do this is to:

1. first extract the first 10 lines
2. take the last two lines of that extract

There are two jobs to be done and the first job passes its data to the second job. This is a good place to use pipes. We can use `head` to get the first 10 lines as we have done above. We will also introduce `tail` which can be used to print the last lines. Try running this command:

```
head -n 10 tb.fasta | tail -n 2
```

## Troubleshooting

Sometimes an error may occur while trying to run a command. There may be a number of different reasons for this and troubleshooting is an important skill to master. 

There are a number of ways to check if your commands have failed including:

1. **Run time:** Some commands (such as alignment) are expected to run for at least a few minutes. If the command finished instantly and you do not expect it to do so then something has probably gone wrong.
2. **Error messages:** Usually if a command has failed it will print out an error message. Check the last few lines of output from a program and check to see if you see any.

There are some errors which are are commonly made. We have listed a few below. 

### File or folder not present

These types of error may present in different ways depending on the command used. For example, let's try change to a directory that doesn't exist. 

```
cd fake_directory
```

You should see something like this: 

!!! terminal "Terminal output"
    ```
    cd: fake_directory: No such file or directory 
    ```

This tells us that `cd` could not find the directory called `fake_directory`. 

!!! info
    Remember that the terminal is case sensitive and will also not tolerate any spaces in places where they shouldn't be. For example try running this command and see what you get:
    
    ```
    head TB.fasta
    ```

    Can you fix it?

### Missing a required argument

If you run a program and fail to specify an arguement it may give an error. Alternatively, it may display the program usage help and quit. For example, lets run `head` without passing the number of lines to the `-n`` flag: 

```
head -n tb.fasta
```

You should see an error similar to that displayed below: 

!!! terminal "Terminal output"
    ```
    head: invalid number of lines 'tb.fasta' 
    ```

The error is a bit cryptic but from the illegal line count, we can see that whatever parameter we passed to `-n` (which determines how many lines are printed) is not a valid value. It turns out that this value always has to be a number, however since we left out that parameter it parses tb.fasta as the value and that causes the error. 

### Command not found

The different software used for the practicals are installed in "conda environments". These allow to have multiple versions of the same program to be installed. If you run into a `command not found` error, it is always a good idea to check if the correct conda environment for the current practical is activated (the name of the currently active environment is usually printed in parenthesis to the left of your prompt). You can activate an environment with 

```
conda activate environment_name
```

Make sure to replace "environment_name" with your chosen environment. You can check which environments are available with:

```
conda env list
```

!!! question
    Try activate the tb-profiler environment. Did it work?