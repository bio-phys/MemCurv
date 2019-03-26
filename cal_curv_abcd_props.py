#! /usr/bin/env python
import numpy as np
import MDAnalysis as mda
import argparse
from numpy import exp, arange
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
from numpy import random
#########################################################################
parser=argparse.ArgumentParser()
parser.add_argument("-lx", help="Float, PBC box vector along x-axis, lx [Angstroms]", type=float)
parser.add_argument("-ly", help="Float, PBC box vector along y-axis, ly [Angstroms]", type=float)
parser.add_argument("-k", help="Integer > 2, Number of nodes used describe the curved surface", type=int)
parser.parse_args()
########################################################################
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
    g=(1 + hf_x**2 + hf_y**2)
    ##Multiplication by (-1) defines the sign of curvatures; Opposite signs correspond to opposite leaflets.
    S=S*(-1/np.sqrt(g**3))
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
    order=np.argsort(c_vals)
    c2, c1=c_vals[order].real
    cv2, cv1 = c_vecs[:, order].transpose()
    M=np.zeros((2,2))
    Kg=np.linalg.det(S)
    H=0.5*np.trace(S)
    ###principal directions w.r.t x clockwise
    theta1 = angle_between(ux, cv1)
    theta2 = angle_between(ux, cv2)
    return [Kg, H, c1, c2, np.rad2deg(theta1), np.rad2deg(theta2), cv1[0], cv1[1]]
##Define principal curvatures
def prin_ax(coord):
    # compute geometric center
    center = np.mean(coord, 0)
    # center with geometric center
    coord = coord - center
    # compute principal axis matrix
    inertia = np.dot(coord.transpose(), coord)
    e_values, e_vectors = np.linalg.eig(inertia)
    order = np.argsort(e_values)
    eval3, eval2, eval1 = e_values[order]
    axis3, axis2, axis1 = e_vectors[:, order].transpose()
    return [axis1, axis2, axis3]

#########################################################################
#########################################################################
##########Input trajectory###################
lx=579.5500
ly=281.4600
nc = 3
ux = [1, 0]
hcf_data=np.genfromtxt("popf_1ns_k3_abcd_fam_buckled.dat")
u_buckle = mda.Universe('fam_buckle.gro', 'fam_buckle_100ns.xtc')
# print len(u_buckle.trajectory)
data=np.zeros((len(u_buckle.trajectory), 26))
hcf0=np.zeros(nc*nc*4)
#Defining the input variables for least square fit from mda frame
for idx, ts in enumerate(u_buckle.trajectory):
    hdx= u_buckle.trajectory.frame
    acf=hcf_data[hdx][0*nc**2:1*nc**2].reshape(nc,nc)
    bcf=hcf_data[hdx][1*nc**2:2*nc**2].reshape(nc,nc)
    ccf=hcf_data[hdx][2*nc**2:3*nc**2].reshape(nc,nc)
    dcf=hcf_data[hdx][3*nc**2:4*nc**2].reshape(nc,nc)
    hcf0=hcf_data[hdx]
#####Defining Protein props#######    
    amp1=u_buckle.atoms.select_atoms('name BB and resid 90:107')
    amp2=u_buckle.atoms.select_atoms('name BB and resid 158:179')
    prot=u_buckle.atoms.select_atoms('name BB')
#####Protein COM position####
    amp1_com = amp1.center_of_mass()
    amp2_com = amp2.center_of_mass()
    prot_com = prot.center_of_mass()
#####Protein principal axis orientation along x and y###
#     print amp1.resids
#     print amp2.resids
#    amp1_v=[ (amp1.positions[0] - amp1.positions[-1])[0], (amp1.positions[0] - amp1.positions[-1])[1]]
#    amp2_v=[ (amp2.positions[0] - amp2.positions[-1])[0], (amp2.positions[0] - amp2.positions[-1])[1]]
#    amp1_rv=[ (amp1.positions[-1] - amp1.positions[0])[0], (amp1.positions[-1] - amp1.positions[0])[1]]
    amp1_v=[prin_ax(amp1.positions)[0][0],prin_ax(amp1.positions)[0][1]]
    amp2_v=[prin_ax(amp2.positions)[0][0],prin_ax(amp2.positions)[0][1]]
#    amp1_rv= -1*amp1_v
    phi1 = angle_between(ux, amp1_v)
#    phir1 =angle_between(ux, amp1_rv)
    phi2 = angle_between(ux, amp2_v)
    gamma = angle_between(amp2_v, amp1_v)
    cuv1= curv(amp1_com[0], amp1_com[1], acf, bcf, ccf, dcf)
    cuv2= curv(amp2_com[0], amp2_com[1], acf, bcf, ccf, dcf)
    cuvp = curv(prot_com[0], prot_com[1], acf, bcf, ccf, dcf)
    k1_v=[cuvp[6], cuvp[7]]
    phic1=angle_between(amp1_v, k1_v)
    phic2=angle_between(amp2_v, k1_v)
#    if np.rad2deg(phi1)>cuvp[4]:
#        phic1=(np.rad2deg(phi1)-cuvp[4])
#    else:
#        phic1=(cuvp[4]-np.rad2deg(phi1)) 
#    if np.rad2deg(phi2)>cuvp[4]:
#        phic2=(np.rad2deg(phi2)-cuvp[4])
#    else:
#        phic2=(cuvp[4]-np.rad2deg(phi2))
#    phic1= (np.rad2deg(phi1)-cuvp[4])
#    phic2= (np.rad2deg(phi2)-cuvp[4])
    data[idx]=[u_buckle.trajectory.time, cuvp[0],cuvp[1],cuvp[2],cuvp[3],cuvp[4],cuvp[5],cuv1[0],cuv1[1],cuv1[2],cuv1[3],cuv1[4],cuv1[5],cuv2[0],cuv2[1],cuv2[2],cuv2[3],cuv2[4],cuv2[5],phi1,phi2,gamma,phic1, phic2, prot_com[0], prot_com[1]]
    print data[idx]
np.savetxt('fam_1ns_abcd_k3_curv_props.dat', data, delimiter='\t', fmt='%0.6f')
