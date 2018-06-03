# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 10:12:03 2018

@author: user
"""

from rpy2.robjects.packages import importr

# Importing the R packages using RPy2
base = importr('base')
utils = importr('utils')
# The following 3 can be commented out once they have ran once
base.source("http://www.bioconductor.org/biocLite.R")
biocInstaller = importr("BiocInstaller")
biocInstaller.biocLite("GEOquery")
# This will throw a RRuntimeError after it has downloaded the GEOquery package.
# Put the GEOquery folder from C:\Program Files\R\<R Version>\library
# into C:\Users\<User name>\Anaconda3\Lib\R\library
# then restart Spyder