# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 13:32:44 2018

@author: Trevor
"""
import pandas as pd
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import r, pandas2ri,RObject
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import DataFrame
import numpy as np

def Grab_Subset(RFrame, IdxA, IdxB, Col):
       DFrame = r.Table(RFrame)
       subset= DFrame[IdxA:IdxB]
       NdArr=pandas2ri.ri2py(subset)
       Matri = np.matrix(NdArr) 
       Cset = Matri[0:(IdxB-IdxA),Col]
       
       return Cset
   
base = importr('base')
utils = importr('utils')

base.source("http://www.bioconductor.org/biocLite.R")
biocInstaller = importr("BiocInstaller")
biocInstaller.biocLite("GEOquery")
GEOquery = r.library("GEOquery")


gds1 = r.getGEO("GDS4280")
gds2 = r.getGEO("GDS3329") 
gds3 = r.getGEO("GDS4182") 
gds4 = r.getGEO("GDS4181")

GDS1_Sub = Grab_Subset(gds1, 2, 26,14823)
GDS2_Sub = Grab_Subset(gds2, 2, 81,14823)
GDS3_Sub = Grab_Subset(gds3, 2, 98,14823)
GDS4_Sub = Grab_Subset(gds4, 2, 82,14823)

#print(GDS1_Sub)
print(GDS1_Sub.shape)
#print(GDS2_Sub) 
print(GDS2_Sub.shape)  
#print(GDS3_Sub)
print(GDS3_Sub.shape) 
#print(GDS4_Sub)
print(GDS4_Sub.shape)

"""
Half1 = np.concatenate(GDS1_Sub,GDS2_Sub)
Half2 = np.concatenate(GDS3_Sub,GDS4_Sub)
INPP4B = np.concatenate(Half1,Half2)
"""

#print(INPP4B)
#print(INPP4B.shape)

"""
gds1df = r.Table(gds1)
subset1= gds1df[2:26]
gds1df2=pandas2ri.ri2py(subset1)
Matri1 = np.matrix(subset1)

Cset1 = Matri1[0:25,14823]
print(Cset1.shape)

gds2df = r.Table(gds2)
subset2= gds2df[2:81]
gds2df2=pandas2ri.ri2py(subset2)
Matri2 = np.matrix(subset2)
Cset2 = Matri2[0:79,14823]
print(Cset2.shape)

gds3df = r.Table(gds3)
subset3= gds3df[2:98]
gds3df3=pandas2ri.ri2py(subset3)
Matri3 = np.matrix(subset3)
Cset3 = Matri3[0:96,14823]
print(Cset3.shape)

gds4df = r.Table(gds4)
subset4= gds4df[2:82]
gds3df4=pandas2ri.ri2py(subset4)
Matri4 = np.matrix(subset4)
Cset4 = Matri4[0:80,14823]
print(Cset4.shape)
"""




