# Setting up NCI Gadi 
This documents talks about various tools and softwares you might need on NCI Gadi to use the AusARG Genome Assembly pipeline 


Table of Contents
=================

- [Conda Environments](#conda-environments)

### Conda Environments

#### Installing Conda in your local folder 
We will use MiniConda3 to locally install conda environments in your local folder 

``` wget URL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh ```

When you run the install script, you should get a prompt asking to confirm the directory. Make sure you provide a local folder, else it will use your home directory as default. Although the home directory install will work fine, it will create problems as you keep using conda. The size of home directory folders on Gadi is 10GB, which is not much considering we will create multiple conda enviroments totalling to hundreds of packages. 

Once you provide the local folder, it should install conda and automatically activate the base environment. 

#### Creating new environments 

``` conda create --name myenv ```

Then, you can activate/deactivate the environment using the following commands 

``` conda activate myenv ```
``` conda deactivate ```
    

#### Running them in a PBS script 
Locate your conda.sh script, and add the following files at the start of your script to activate it

``` source /path/to/miniconda/etc/profile.d/conda.sh ```
``` conda activate trash ```

and that's it! Now you can use conda environments within your PBS scripts. 




