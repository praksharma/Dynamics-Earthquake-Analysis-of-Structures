##########################################################################
# program: mdof_arbit_force_nm.py
# author: Tom Irvine
# version: 1.1
# date: December 30, 2011
# description:  Calculate the response of an MDOF system to an applied
#               force or forces using the Newmark-beta method.
#
#               The input file must have two columns: time(sec) & force
##########################################################################

from tompy import read_two_columns
from tompy import GetInteger2
from tompy import time_history_plot
from tompy import read_array,enter_int,enter_float

from ode_plots import ode_plots

from Newmark import Newmark_force

from generalized_eigen import generalized_eigen

from numpy import zeros,linspace,interp

import matplotlib.pyplot as plt

########################################################################

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
C=read_array("damping coefficient matrix")
K=read_array("stiffness matrix")


if(iu==1):
    M/=386

print(" ")


ndof,fn,NN =generalized_eigen(K,M,2)

num=ndof

max_fn=max(fn)
min_fn=min(fn)

if(min_fn>0):
    period_max=1/min_fn
else:
    period_min=1/max_fn

period_min=1/max_fn

dtr=period_min/20

########################################################################

v0=read_array("initial velocity vector")

d0=read_array("initial displacement vector")

########################################################################

if(iu==1):
    print "Each force file must have two columns: time(sec) & force(lbf)"
else:
    print "Each force file must two columns: time(sec) & force(N)"


print " "
print " Enter the number of force files "

nff=enter_int()

MAX=100000


tt=zeros((MAX,nff),float)
ff=zeros((MAX,nff),float)
force_dof=zeros(ndof,int)

for i in range (0,ndof):
    force_dof[i]=-999



print "  "
print " Note: the first dof is 1 "


for i in range (0,nff):
    ii=i+1
    print " "
    print " Enter force file %d " % ii
    a,b,num=read_two_columns()
    L=len(a)
    if(L>MAX):
        L=MAX
    tt[0:L,i]=a[0:L]
    ff[0:L,i]=b[0:L]

    print " "
    print " Enter the number of dofs at which this force is applied"

    nfa=enter_int()


    for j in range (0,nfa):
        print " "

        if(j==0 and nfa==1):
            print " Enter the dof number for this force "

        if(j==0 and nfa>1):
            print " Enter the first dof number for this force "

        if(j>0 and nfa>1):
            print " Enter the next dof number for this force "

        nn=enter_int()
        nn=nn-1
        force_dof[nn]=i


#*******************************************************************************

print " "
print " Enter the duration (sec) "
dur=enter_float()

srr=20*max_fn

print " "
print " Enter the sample rate (recommend > %8.4g) " % srr
sr=enter_float()

dt=1/sr

nt=int(dur/dt)

t = linspace(0., dur, nt)

#*******************************************************************************
#
#  interpolate force
#

print " begin interpolation"

FFI=zeros((nt,nff),float)

for i in range (0,nff):
    tstart=tt[0,i]
    tt[:,i]=tt[:,i]-tstart
#    print max(t),max(tt[:,i])

    last=MAX
    for j in range (1,MAX):
        if(tt[j,i]<=tt[(j-1,i)]):
            last=j
            break

#    print last


    tint=tt[0:last,i]
    fint=ff[0:last,i]

    FFI[:,i] = interp(t,tint,fint)

print " end interpolation"

nrows=FFI.shape[0]
ncolumns=FFI.shape[1]

print nrows,ncolumns

#*******************************************************************************
#

x=zeros((nt,ndof),float)
y=zeros((nt,ndof),float)
a=zeros((nt,ndof),float)


x[i,:]=d0
y[i,:]=v0



##########################################################################

d,v,a =Newmark_force(M,C,K,FFI,force_dof,v0,d0,dt,nt,ndof)

##########################################################################

if(iu==1):
    a=a/386
else:
    a=a/9.81

fignum=1
ode_plots(fignum,nt,ndof,t,d,v,a,iu)

plt.show()
