##########################################################################
# program: mdof_initial_nm.py
# author: Tom Irvine
# version: 1.4
# date: May 1, 2012
# description:  multi-degree-of-freedom system subjected to initial
#               condition excitation
#
#               Solution via the Newmark-beta method
#
##########################################################################


from tompy import read_array,GetInteger2,GetInteger3,enter_float

from Newmark import Newmark_free
from generalized_eigen import generalized_eigen

import numpy as np

import matplotlib.pyplot as plt

from ode_plots import ode_plots,ode_plots_sdof

##########################################################################

print("  ")
print(" Select units:      ")
print(" 1=English 2=metric ")

iu=GetInteger2()

print("  ")
if(iu==1):
    print(" units: mass (lbm)")
    print("        damping coefficients (lbf sec/in)")
    print("        stiffness (lbf/in)")
    print("        velocity (in/sec)")
    print("        displacement (in)")
else:
    print(" units: mass (kg)")
    print("        damping coefficients (N sec/m)")
    print("        stiffness (N/m)")
    print("        velocity (m/sec)")
    print("        displacement (m)")

print("  ")

M=read_array("mass matrix")

nlength=len(np.atleast_1d(M)) 


print("  ")

print " Select damping input method "
print "  1=damping coefficient matrix "
print "  2=uniform modal damping "
print "  3=modal damping vector"

idamp=GetInteger3()


if(idamp==1):
    C=read_array("damping coefficient matrix")

if(idamp==2):
    print " Enter uniform damping ratio "
    damp=enter_float()

if(idamp==3):
    dampV=read_array(" Enter damping ratio vector")


K=read_array("stiffness matrix")


if(iu==1):
    M/=386

print(" ")

ndof,FN,MS =generalized_eigen(K,M,1)

MST=MS.T

num=ndof


if(idamp==2 or idamp==3):
    ccc=np.zeros((ndof,ndof),float)

    for i in range(0,ndof):
        if(idamp==2):
            ccc[i,i]=2*damp*(2*np.pi*FN[i])
        else:
            ccc[i,i]=2*dampV[i]*(2*np.pi*FN[i])

    s1=np.dot(MST,M)
    s2=np.dot(ccc,s1)
    s3=np.dot(MS,s2)
    C=np.dot(M,s3)
    print(" ")
    print(" C= ")
    print C
    print(" ")

##########################################################################


max_fn=max(FN)
min_fn=min(FN)

if(min_fn>0):
    period_max=1/min_fn
else:
    period_min=1/max_fn

period_min=1/max_fn


dur=20*period_max
dt=period_min/50
NT=int(round(dur/dt))


t = np.linspace(0, dur, NT)


v0=read_array("initial velocity vector")

d0=read_array("initial displacement vector")

##########################################################################

d,v,a =Newmark_free(M,C,K,v0,d0,dt,NT,ndof)

##########################################################################

if(iu==1):
    a=a/386
else:
    a=a/9.81

fignum=1

if(nlength>1):
    ode_plots(fignum,NT,ndof,t,d,v,a,iu)
else:
    ode_plots_sdof(fignum,NT,ndof,t,d,v,a,iu)

plt.show()
