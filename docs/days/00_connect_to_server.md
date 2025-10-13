# Connect to RStudio server

If you are enrolled in the course, you have access to RStudio server. You can find a link to it and your credentials to connect in the Google Docs that has been shared with you.  

During the first part of the course we will be using bash commands and executing different software. For that, we will work on the built-in terminal of RStudio.   

To access the terminal, select the tab **Terminal** on the top left of the panel:

![](../assets/images/select_terminal_tab.png){width=300}

Check your current directory, you should be in your home directory `/home/rstudio` 

```{bash}
pwd
```

In order to activate your conda enviroment and have access to different tools necessary for this course, please run the following command in your terminal:

```{bash}
conda init
source /home/rstudio/.bashrc
conda activate atac_env
```

You should now be able to run tools like samtools, bedtools, macs2, etc. Check with:
```{bash}
bedtools --version
```
  
  
For the second part of the course we will use only R commands. For that, we will move back to RStudio console, by selecting the tab **Console**.