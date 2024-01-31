<details>
<summary>

## Getting linux terminal on Windows
</summary>

Check if you can use liux commands without installing anything:

1) open the command line or powershell and type some linux command like 'ls' or 'ssh'
2) If it works, great! If not, you may need to proceed with the following steps or do some other kind of trouble-shooting

Installing linux in windows shell:
https://learn.microsoft.com/en-us/windows/wsl/install

If you have a windows10 or newer, open powershell in administrator mode:

1) type powershell in search bar
2) right click powershell
3) click 'open as administrator'

```
wsl --install
```

This will take a few moments. It will not work if you have wsl already. It will prompt you for a username and password. Restart your computer once it is done installing. You will know that it is done when you see the username to set pop up.

Open the powershell, look for your username, and try linux commands now. Hopefully, it worked.

Also, this process will create a linux folder. Explore the file system and find your username in this linux folder. Files you download through scp or sftp will likely end up here by default.  

</details>

<details>
<summary>

## Useful linux commands 
</summary>

Command                                |Description
---------------------------------------|---------------------------------------
pwd                                    |#print working directory
echo hello world                       |#print hello world to the screen 
ls                                     |#list
ls -lath                               |#list all files in long listing format by date modified and with their file size in human readable format
cd                                     |#change directory
cd ..                                  |#go backwards one directory
cd -                                   |#go to the previous directory
cd ~                                   |#go to home directory
cp filename new_filename               |#copy/rename a file
cp -R /path/to/directory /new/path     |#copy a directory recursively 
rm filename                            |#delete filename
rm -r directory                        |#delete a directory with contents
rmdir directory                        |#delete a directory without contents
mkdir directory                        |#create a directory
vi filename                            |#create a file
touch filename                         |#create a file
history                                |#if you want to see all the previous commands you have typed
sbatch script.sh                       |#sumbit a job to the hpc
squeue --me                            |#see the status of a file you have sumbitted to the slurm queue

</details>

<details>
<summary>

## Links for Practice
</summary>
https://sandbox.bio/tutorials/terminal-basics

</details>
