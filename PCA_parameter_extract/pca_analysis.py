#!//usr/bin/python3

import matplotlib
matplotlib.rcParams.update({'font.size': 22})

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA, IncrementalPCA

#investigate what the following library does -
#it is used in the ML w/ python book
#from sklearn.preprocessing import StandardScaler

from matplotlib import style
style.use("fivethirtyeight")

import pandas as pd

#------------------------------------------------------------------------
def make_f_intuitive(data,N_part = 3666,eta=0.8,SX=60.0,verbose=False):
    '''
    given naive, unsorted dr values from frame_data file, 
    sort into feature_data, called fI in the Jadrich paper 
    (and also calculated with the analysis code, so this is a checker).

    required parameters:
    
    data: 
    type: numpy array 
    description: containing ALL dr values from all particles w/in one frame

    optional parameters:

    N_part: 
    type: integer
    description: total number of particles from simulation
    default = 3666 

    eta:
    type: float
    description: eta is area fraction (called density in MD simulation), 
    compared w/ rho which is number density,  
    default = 0.7

    SX:
    type: float
    description: sidelength of system.  assumes square simulation box.
    default = 60.0

    verbose: 
    type: boolean
    description: not-so-helpful data dumps
    default = False
    '''
        
    rho=N_part/(SX**2)
    l = 1/np.sqrt(rho)
    
    if verbose == True:
        print("The average separation is: ",l)

    #number of particles came out correctly
    N_nn = int(np.shape(data)[0]/float(N_part))
    
    if verbose == True:
        print("Number of %pagerticles is: ", N_part, 
              "\nPredicted number of neighbors is", N_part-1, 
              "\nActual number of neighbors is: ", N_nn)
    
    #this will split the data into equal size numpy arrays, not pandas dfs
    #grab just the dr value, not the dx or dy values for instance
    data_split = np.array_split(data, N_part)

    #check the shape, rows = #particles, columns = #neighbor measures
    if verbose == True:
        print(np.shape(data_split))

    #As split, each row contains a single feature (but it isn't sorted yet, so we do nothing about this yet)
    
    #set up list for sorting data into.  
    #this is inefficient, but its a start.
    #C would be better - using the qsort() algorithm.  (and exists!)
    delta_r_norm = []

    #SOME DEBUGGING PRINT() STATEMENTS
    if verbose == True:
        print("loop through each row, meaning each feature vector and sort it.  \n There are N_particles = m rows")
        print(np.shape(delta_r_norm))
        print(np.shape(data_split))
        print(np.shape(data_split[i]))

    #loop through each row, meaning each feature vector and sort it.  
    #There are N_particles = m rows

    for i in range(N_part):                                                          
        #we sort the data contained in each row,
        #so a feature is in a row (as of now)
        delta_r_norm.append( np.sort(data_split[i]) ) 

    #The ARRAY SHOULD BE OF LENGTH M=ROWS=3208, N=COLUMNS=3207 
    if verbose == True:
        print(np.shape(delta_r_norm))

    #reassign data to array for efficiency
    delta_r_norm = np.array(delta_r_norm)

    #CALCULATE <r_i>_D in Jadrich
    #calculate the mean of each column in the matrix - 
    #this is the average of the $ith$ nearest neighbor distance 
    r_mean = []
    for i in range(len(delta_r_norm[0])):
        r_mean.append(np.mean(delta_r_norm[:,i]))
    
    if verbose == True:
        #plot to see if r_mean makes sense
        plt.xlabel("$i$") #,fontsize=20)
        plt.ylabel(r"$\langle r_i \rangle$") #,fontsize=20)
        plt.plot(r_mean)  #the value is a function of the counting number 1, 2, 3, 4, ... N_nn
        plt.show()

    #print(len(r_mean))=3207=N_nn

    #if you plotted any single column in the matrix it would look more or less like this
    #take the transpose of the matrix so that the "feature vectors" are column.  
    #This means that the features are in rows and the samples are in columns (opposite of sklearn standard)
    delta_r_norm = delta_r_norm.T

    if verbose == True:
        print("The following matrix should have ascending neighbor distances in columns")
        print("as in Jadrich et al.")
        print(delta_r_norm)
        
        #this is what we want
        #N-1 rows
        #N columns
        print(np.shape(delta_r_norm))
        print(np.shape(delta_r_norm[0]))
        print(np.shape(delta_r_norm[:,0]))
        
    #do the mean subtraction...
    #every column of delta_r_norm is a descending list of neighbor distances of length 3207
    #r_mean is a column vector of length 3207
    for i in range(N_nn):
        #print(delta_r_norm[i], r_mean[i]) # this input is wonderful!
        delta_r_norm[i] = delta_r_norm[i] - r_mean[i]
        #the output on the otherhand...
        
    #normalize the data as in Jadrich
    return delta_r_norm/l

#------------------------------------------------------------------------

if __name__ == "__main__":

    
