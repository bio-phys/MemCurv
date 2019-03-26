#! /usr/bin/env python
#---------------------------------------------------------------------------#
# Function: Fit a spherical surface to a discontinuous lipid patch.         #
# Usage: ./Fit_bicelle_sph_cap.py				            #
# Author: Ramachandra M. Bhaskara (rabhaska@biophys.mpg.de)                 #
# Version: 0.1 (06.06.2018)                                                 #
#---------------------------------------------------------------------------#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
#OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
#MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND 
#NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE 
#DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, 
#WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
#CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#############################################################################
import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.rms
import MDAnalysis.analysis.helanal
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import Bio
#%matplotlib inline
from Bio.Seq import Seq
from Bio import SeqIO
from scipy import optimize
from scipy.spatial import distance
from scipy.spatial.distance import cdist
from scipy import stats
from matplotlib.mlab import PCA
from numpy import (array, dot, arccos, clip)
from numpy.linalg import norm
import csv
from scipy.spatial.distance import cdist
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
import pandas as pd

### define a function to return radius given a set of x, y and z coord
### and center coordinates i.e. xc, yc and zc
def calc_R(xc,yc,zc):
    return np.sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)

### Function to minimize or optimize
### calculate the algebraic distance between the data points and the 
### mean circle centered at point c=(xc, yc, zc) 
def func(c):
    Ri = calc_R(*c)
    return Ri - Ri.mean()

#Load trajectory
u_bicel = mda.Universe('fam_bicelle.gro','fam_bicelle.xtc')
#define output data_matrix
data=np.zeros((len(u_bicel.trajectory),8))
#loop over trajectory
for idx, ts in enumerate(u_bicel.trajectory):
    ## Define midplane from lipid beads
    midplane = u_bicel.atoms.select_atoms('name C3A or name C3B')
    ## define x, y, z for coordinates
    x = midplane.positions[:,0]
    y = midplane.positions[:,1]
    z = midplane.positions[:,2]
    #Define initial center estimate
    center_estimate = np.mean(x), np.mean(y), np.mean(z)
    #Least square optimization usage 
    ## returns final center coordinates and steps by taking in 
    ## the function to optimize along with its input estimate
    final_center, ier = optimize.leastsq(func, center_estimate)
    xc_f, yc_f, zc_f = final_center
    Ri_f       = calc_R(*final_center)
    R_f        = Ri_f.mean()
    residue_f   = sum((Ri_f - R_f)**2)
    po4=u_bicel.atoms.select_atoms('resnum 180-818 and name PO4')
    k=0
    #Assign sign to the bicelle curvature, k is 50% of the DMPC lipids from the upper leaflet
    for j in range(len(po4.positions)):
        dij=np.sqrt((po4.positions[j,0]-xc_f)**2+(po4.positions[j,1]-yc_f)**2+(po4.positions[j,2]-zc_f)**2)
	if (dij >=R_f):
            k=k+1
        
    print k
    if (k <= 319):
        R_f=-1*R_f
#    else:
#        R_f=-1*R_f   
#print residu_2
    data[idx] = [idx, R_f, 1/R_f, xc_f, yc_f, zc_f, residue_f, ier]
    print idx, R_f, 1/R_f, xc_f, yc_f, zc_f, residue_f, ier
#print Ri_2
#print center_2
np.savetxt('bicelle_curv_ts.dat', data, delimiter='\t', fmt='%0.6f')
