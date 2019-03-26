#! /usr/bin/env python
#---------------------------------------------------------------------------#
# Function: Fit a spherical surface to a discontinuous lipid patch.         #
# Usage: ./Get_mem_profile_hcf_abcd.py                                      #
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
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from numpy import exp, arange
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit, least_squares, leastsq, minimize
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from sympy import *
from scipy import linalg
from scipy.misc import derivative
from scipy.integrate import quad
from scipy.integrate import dblquad
from numba import jit
from six.moves import range
from scipy.optimize import check_grad
import time
from scipy.stats.kde import gaussian_kde
from numpy import linspace

###Define a 2D-Fourier series function 
### input variables x-value and y-value
### hcf=[]; matrix with coefficients for height; to be optimized
#define also number of coefficients 
def fs(x, y, acf, bcf, ccf, dcf):
    #Box periodicity
    lx=579.55004883
    ly=281.45999146
    acal = bcal = ccal = dcal = hcal = 0
    ##double summation over the number of wavelength vectors along x and y
    for n in range(len(acf)):
        for m in range(len(acf)):
            acal += acf[n,m]* np.cos(2*np.pi*n*x/lx) * np.cos(2*np.pi*m*y/ly)
    for n in range(len(bcf)):
        for m in range(len(bcf)):
            bcal += bcf[n,m]* np.cos(2*np.pi*n*x/lx) * np.sin(2*np.pi*m*y/ly)
    for n in range(len(ccf)):
        for m in range(len(ccf)):
            ccal += ccf[n,m]* np.sin(2*np.pi*n*x/lx) * np.cos(2*np.pi*m*y/ly)
    for n in range(len(dcf)):
        for m in range(len(dcf)):
            dcal += dcf[n,m]* np.sin(2*np.pi*n*x/lx) * np.sin(2*np.pi*m*y/ly)
    hcal= acal + bcal + ccal + dcal
    return hcal

##Define a cost residual function to minimize the least squares difference
##between the observed and calculated "h"
##hcf is a array of height coefficients to be optimized
##note reshaping input array to compute the function "fs"
def residual(hcf, x, y, h):
    return h - fs(x, y, hcf[0:nc**2].reshape(nc, nc), hcf[nc**2:2*nc**2].reshape(nc,nc), hcf[2*nc**2:3*nc**2].reshape(nc,nc), hcf[3*nc**2:4*nc**2].reshape(nc,nc))

#########################################################################
##Define Suface height function as an optimized fourier series 
##Provides h(x,y) for the periodic analytical surface given any x and y value 
def hf(x,y, a_coef, b_coef, c_coef, d_coef):
    lx=579.55004883
    ly=281.45999146
    acal = bcal= ccal = dcal = hcal = 0
    ##double summation over the number of wavelength vectors along x and y
    for n in range(len(a_coef)):
        for m in range(len(a_coef)):
            acal += a_coef[n,m]* np.cos(2*np.pi*n*x/lx) * np.cos(2*np.pi*m*y/ly)
    for n in range(len(b_coef)):
        for m in range(len(b_coef)):
            bcal += b_coef[n,m]* np.cos(2*np.pi*n*x/lx) * np.sin(2*np.pi*m*y/ly)
    for n in range(len(c_coef)):
        for m in range(len(c_coef)):
            ccal += c_coef[n,m]* np.sin(2*np.pi*n*x/lx) * np.cos(2*np.pi*m*y/ly)
    for n in range(len(d_coef)):
        for m in range(len(d_coef)):
            dcal += d_coef[n,m]* np.sin(2*np.pi*n*x/lx) * np.sin(2*np.pi*m*y/ly)
    hcal= acal + bcal + ccal + dcal
    return hcal
########################################################################
#####Define angle between 2 vectors
def angle_between(p1, p2):
    ang1 = np.arctan2(*p1[::-1])
    ang2 = np.arctan2(*p2[::-1])
    dg= np.rad2deg((ang1 - ang2) % (2 * np.pi))
#     return dg
    return np.deg2rad(dg)

################################################
###Define Area element dA of a monge patch #####
def ae(x, y, a_coef, b_coef, c_coef, d_coef):
    [hf_x, hf_y] = [ 0, 0]
    ###Define partial differential w.r.t x
    for n in range(len(a_coef)):
        for m in range(len(a_coef)):
            hf_x += -((2*np.pi*n/lx)*(a_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))
    for n in range(len(b_coef)):
        for m in range(len(b_coef)):
            hf_x += -((2*np.pi*n/lx)*(b_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))         
    for n in range(len(c_coef)):
        for m in range(len(c_coef)):
            hf_x +=  ((2*np.pi*n/lx)*(c_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))
    for n in range(len(d_coef)):
        for m in range(len(d_coef)):
            hf_x +=  ((2*np.pi*n/lx)*(d_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))
    ###Define partial differential w.r.t y
    for n in range(len(a_coef)):
        for m in range(len(a_coef)):
            hf_y += -((2*np.pi*m/ly)*(a_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))
    for n in range(len(b_coef)):
        for m in range(len(b_coef)):
            hf_y +=  ((2*np.pi*m/ly)*(b_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))         
    for n in range(len(c_coef)):
        for m in range(len(c_coef)):
            hf_y += -((2*np.pi*m/ly)*(c_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))
    for n in range(len(d_coef)):
        for m in range(len(d_coef)):
            hf_y +=  ((2*np.pi*m/ly)*(d_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))
    ###Define area element###
    dA = np.sqrt(1 + hf_x**2 + hf_y**2)
    return dA

######Define shape operator######################
def sho(x, y, a_coef, b_coef, c_coef, d_coef):
    hf = hf_x = hf_y = hf_xx = hf_yy = hf_xy = 0
    lx=579.5500
    ly=281.4600
    
    ###Define partial differential w.r.t x
    for n in range(len(a_coef)):
        for m in range(len(a_coef)):
            hf_x += -((2*np.pi*n/lx)*(a_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))
    for n in range(len(b_coef)):
        for m in range(len(b_coef)):
            hf_x += -((2*np.pi*n/lx)*(b_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))         
    for n in range(len(c_coef)):
        for m in range(len(c_coef)):
            hf_x +=  ((2*np.pi*n/lx)*(c_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))
    for n in range(len(d_coef)):
        for m in range(len(d_coef)):
            hf_x +=  ((2*np.pi*n/lx)*(d_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))
            
    ###Define partial differential w.r.t y
    for n in range(len(a_coef)):
        for m in range(len(a_coef)):
            hf_y += -((2*np.pi*m/ly)*(a_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))
    for n in range(len(b_coef)):
        for m in range(len(b_coef)):
            hf_y +=  ((2*np.pi*m/ly)*(b_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))         
    for n in range(len(c_coef)):
        for m in range(len(c_coef)):
            hf_y += -((2*np.pi*m/ly)*(c_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))
    for n in range(len(d_coef)):
        for m in range(len(d_coef)):
            hf_y +=  ((2*np.pi*m/ly)*(d_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))
            
    ###Define 2nd order partial differential w.r.t x: d^2/dx^2
    for n in range(len(a_coef)):
        for m in range(len(a_coef)):
            hf_xx += -((2*np.pi*n/lx)**2*(a_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))
    for n in range(len(b_coef)):
        for m in range(len(b_coef)):
            hf_xx += -((2*np.pi*n/lx)**2*(b_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))         
    for n in range(len(c_coef)):
        for m in range(len(c_coef)):
            hf_xx += -((2*np.pi*n/lx)**2*(c_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))
    for n in range(len(d_coef)):
        for m in range(len(d_coef)):
            hf_xx += -((2*np.pi*n/lx)**2*(d_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))
            
    ###Define 2nd order partial differential w.r.t y: d^2/dy^2
    for n in range(len(a_coef)):
        for m in range(len(a_coef)):
            hf_yy += -((2*np.pi*m/ly)**2*(a_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))
    for n in range(len(b_coef)):
        for m in range(len(b_coef)):
            hf_yy += -((2*np.pi*m/ly)**2*(b_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))         
    for n in range(len(c_coef)):
        for m in range(len(c_coef)):
            hf_yy += -((2*np.pi*m/ly)**2*(c_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))
    for n in range(len(d_coef)):
        for m in range(len(d_coef)):
            hf_yy += -((2*np.pi*m/ly)**2*(d_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))
            
    ###Define 2nd order artial differential w.r.t x and y d/dy(d/dx(h(x,y)))
    for n in range(len(a_coef)):
        for m in range(len(a_coef)):
            hf_xy +=  ((2*np.pi*n/lx)*(2*np.pi*m/ly)*(a_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))
    for n in range(len(b_coef)):
        for m in range(len(b_coef)):
            hf_xy += -((2*np.pi*n/lx)*(2*np.pi*m/ly)*(b_coef[n,m])*(np.sin(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))         
    for n in range(len(c_coef)):
        for m in range(len(c_coef)):
            hf_xy += -((2*np.pi*n/lx)*(2*np.pi*m/ly)*(c_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.sin(2*np.pi*m*y/ly)))
    for n in range(len(d_coef)):
        for m in range(len(d_coef)):
            hf_xy +=  ((2*np.pi*n/lx)*(2*np.pi*m/ly)*(d_coef[n,m])*(np.cos(2*np.pi*n*x/lx))*(np.cos(2*np.pi*m*y/ly)))

    ###Define shape operator at (x,y)
    S=np.zeros((2,2))
    S[0,0]=((1 + hf_y**2)*hf_xx - (hf_x*hf_y*hf_xy))
    S[0,1]=((1 + hf_x**2)*hf_xy - (hf_x*hf_y*hf_xx))
    S[1,0]=((1 + hf_y**2)*hf_xy - (hf_x*hf_y*hf_yy))
    S[1,1]=((1 + hf_x**2)*hf_yy - (hf_x*hf_y*hf_xy))
    gs=np.sqrt((1 + hf_x**2 + hf_y**2)**3)
    S=S*(1/gs)
    return np.matrix(S)
##############################################
####Defining Mean Curvature field calculations
def mcf(x, y, a_coef, b_coef, c_coef, d_coef):
    S=sho(x,y, a_coef, b_coef, c_coef, d_coef)
    ###Define principal curvatures and principal directions at x,y
    c_vals, c_vecs = linalg.eig(S)
    M=np.zeros((2,2))
    c1=c_vals[0]
    c2=c_vals[1]
    Kg=c1*c2
    H=(c1+c2)/2
    cv1 = c_vecs[:,0]
    cv2 = c_vecs[:,1]
    ux =[1 , 0]
    ###principal directions w.r.t x clockwise
    theta1 = angle_between(ux, cv1)
    theta2 = angle_between(ux, cv2)
    ###define mean curvature tensor at x,y
    M=np.zeros((2,2))
    M[0,0]= (c1*(np.cos(theta1)**2) + c2*(np.sin(theta1)**2))
    M[0,1]= ((c2-c1)*np.cos(theta1)*np.sin(theta1))
    M[1,0]= ((c2-c1)*np.cos(theta1)*np.sin(theta1))
    M[1,1]= (c1*(np.sin(theta1)**2) + c2*(np.cos(theta1)**2))
    return np.matrix(M)
############Define Curvatures at any point########
def curv(x,y, a_coef, b_coef, c_coef, d_coef):
    S=sho(x,y, a_coef, b_coef, c_coef, d_coef)
    ###Define principal curvatures and principal directions at x,y
    c_vals, c_vecs = linalg.eig(S)
    M=np.zeros((2,2))
    c1=c_vals[0]
    c2=c_vals[1]
    cv1 = c_vecs[:,0]
    cv2 = c_vecs[:,1]
    Kg=c1*c2
    H=(c1+c2)/2
    ###principal directions w.r.t x clockwise
    theta1 = angle_between(ux, cv1)
    theta2 = angle_between(ux, cv2)
    return [Kg, H, c1, c2, np.rad2deg(theta1), np.rad2deg(theta2)]
#########################################################################
##########Input trajectory###################
lx=579.55004883
ly=281.45999146
nc = 3
ux=[0, 1]
#px=py=50
u_buckle = mda.Universe('fam_buckle.gro', 'fam_buckle_100ns.xtc')
data=np.zeros((len(u_buckle.trajectory), nc*nc*4))
#Defining the input variables for least square fit from mda frame
for idx, ts in enumerate(u_buckle.trajectory):
    lipids=u_buckle.atoms.select_atoms('name PO4')
    x=lipids.positions[:,0]
    y=lipids.positions[:,1]
    h=lipids.positions[:,2]
#    print u_buckle.trajectory.frame, u_buckle.trajectory.time 
#Define the initial parameters as a 1D array of coefficients
    hcf0=np.zeros(nc*nc*4)
##Perform least squares optimization using scipy.optimize
    popf, pcov, info, mesg, ier = leastsq(residual, hcf0, args=(x,y, h),full_output=True)
#optimized coefficients reshaped into a 2D array from output
#    acf=popf[0:16].reshape(nc,nc)
#    bcf=popf[16:32].reshape(nc,nc)
#    ccf=popf[32:48].reshape(nc,nc)
#    dcf=popf[48:64].reshape(nc,nc)
    print np.absolute(residual(popf, x, y, h)).sum()/len(h)
    data[idx]= popf
    print idx, data[idx]
    hcf0=popf
np.savetxt('popf_1ns_k3_abcd_fam_buckled.dat', data, delimiter='\t', fmt='%0.6f')
