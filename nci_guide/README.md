# Setting up NCI Gadi 
This documents talks about various tools and softwares you might need on NCI Gadi to use the AusARG Genome Assembly pipeline 


Table of Contents
=================

- [Conda Environments](#conda-environments)

### Conda Environments


#### Installing Conda in your local folder 
We will use MiniConda3 to locally install conda environments in your local folder 

``` wget URL ```

When you run the install script, you should get a prompt asking to confirm the directory. Make sure you provide a local folder, else it will use your home directory as default. Although the home directory install will work fine, it will create problems as you keep using conda. The size of home directory folders on Gadi is 10GB, which is not much considering we will create multiple conda enviroments totalling to hundreds of packages. 

Once you provide the local folder, it should install conda and automatically activate the base environment. 

#### Creating new environments 

#### Running them in a PBS script 
1. Locate your conda.sh script 
2. 





