
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import time
from sklearn.cluster import AgglomerativeClustering, KMeans

from rpy2.robjects import r, pandas2ri

#  Manipulate data structures in order to access and save specific elements
def Grab_Subset(RFrame, IdxA, IdxB, Col):
    ''' Put a subset of expression data from GEO into an array

    Use rpy2 to get the Table dataframe from the GDS R Object
    Grab only user specified columns of this Table
    Turn those columns into a matrix in order to get user specified values
    Store these particular values in an array 
    
    Keyword Arguments:
    RFrame -- The GDS R object holding a data frame of expression data
    IdxA -- Starting index of Table data frame samples to retrieve, as in the R Table object, but with Python indexes
    IdxB -- Ending index of Table data frame samples to retrieve, as in the R Table object, but with Python indexes
    Col -- Row index of the gene of interest (INPP4B), as in the R Table object, but with Python indexes
    
    Returns:
    Cset -- Array of specified expression values
    '''
    DFrame = r.Table(RFrame)
    subset= DFrame[IdxA:IdxB]
    NdArr=pandas2ri.ri2py(subset)
    Matri = np.matrix(NdArr) 
    tempM = Matri[0:(IdxB-IdxA),Col]
    Cset = np.array(tempM, ndmin=2)
   
    return Cset

# loading the GEOquery library from R through RPy2
GEOquery = r.library("GEOquery")

# Grabbing specific datasets from GEO
gds1 = r.getGEO("GDS4280")
gds2 = r.getGEO("GDS3329") 
gds3 = r.getGEO("GDS4182") 
gds4 = r.getGEO("GDS4181")


# Changing the data into a matrix from an R DataFrame
GDS1_Sub = Grab_Subset(gds1, 2, 26,14823)
GDS2_Sub = Grab_Subset(gds2, 2, 81,14823)
GDS3_Sub = Grab_Subset(gds3, 2, 98,14823)
GDS4_Sub = Grab_Subset(gds4, 2, 82,14823)

# Concatenating the data into one data set
INPP4B = np.append(GDS1_Sub, GDS2_Sub)
INPP4B = np.append(INPP4B, GDS3_Sub)
INPP4B = np.append(INPP4B, GDS4_Sub)

#Converting to a data frame in order to have specific labels
INPP4Bdf = pd.DataFrame({'INPP4B' : INPP4B, 'Samples' : np.arange(279)})

# Clustering Parameters
numClusters = 2
init = 25

# Ward Clustering
ward = AgglomerativeClustering(n_clusters=numClusters, linkage='ward')

startTime = time.time()
ward.fit(INPP4Bdf)
wardTime = time.time() - startTime
wardLabels = ward.labels_
print("Ward took: %.2fs" % wardTime)

#plot the Ward clusters
plt.scatter(INPP4Bdf['Samples'], INPP4Bdf['INPP4B'], c=wardLabels)
plt.show()

#Kmeans Clustering
startTime = time.time();

#KMeans analysis
kmeans = KMeans(n_clusters=numClusters, n_init=init)
kmeans.fit(INPP4Bdf)
kmeansLabels = kmeans.predict(INPP4Bdf)

kmeansTime = time.time() - startTime
print("K-Means took: %.2fs" % kmeansTime)

#Comparing times
if kmeansTime > wardTime:
    print("K-Means clustering is faster.")
elif kmeansTime < wardTime:
    print("Ward linkage clustering is faster.")  
else:
    print("Algorithms do not differ in time.")

#Display KMeans clusters
plt.scatter(INPP4Bdf['Samples'], INPP4Bdf['INPP4B'], c=kmeansLabels)
plt.show()

#compare if the clusters are the exact same
for i in range(numClusters):
    for j in range(numClusters):
    
        if np.array_equal((wardLabels == j).nonzero(), (kmeansLabels == i).nonzero()):
            print("Found the same clusters")
            print("Ward:", j, " K-Means: ", i)
        else:
            print("Ward Cluster ", j, " is not equal to K-Means Cluster ", i)
      
#total how many points are shared between the clusters
numIn = 0
notIn = 0
for i in range(279):
    if wardLabels[i] != kmeansLabels[i]:
        notIn = notIn + 1
    if wardLabels[i] == kmeansLabels[i]:
        numIn = numIn + 1
print("There are",numIn, "similar data points between the clusters from these algorithms")
print("There are",notIn, "dissimilar data points between the clusters from these algorithms")
        

