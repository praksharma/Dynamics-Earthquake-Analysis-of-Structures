##########################################################################
# program: sdof_initial_nm.py
# author: Tom Irvine
# version: 1.2
# date: December 13, 2011
# description:  single-degree-of-freedom system subjected to initial 
#               condition excitation via Newmark-beta method
##########################################################################


from tompy import enter_damping,enter_fn
from tompy import enter_float,time_history_plot

from Newmark import Newmark_free

import numpy as np

import matplotlib.pyplot as plt

##########################################################################

iu,fn,omegan,period=enter_fn()

damp,Q=enter_damping()


if(iu==1):
    print(" Enter initial velocity (in/sec)")
else:
    print(" Enter initial velocity (m/sec)")

V0=enter_float()


if(iu==1):
    print(" Enter initial displacement (in)")
else:
    print(" Enter initial displacement (m)")

D0=enter_float()


##########################################################################

omegan2=omegan**2

dur=20*period
dt=period/50
NT=int(round(dur/dt))


cdm = 2*damp*omegan

print ' '
print ' omegan=%8.4g rad/sec fn=%8.4g Hz' %(omegan,fn)

print ' cdm=%8.4g (rad/sec) ' %(cdm)

t = np.linspace(0, dur, NT)

##########################################################################
            
omegad=omegan*np.sqrt(1-damp**2)
domegan=damp*omegan

M=1
C=2*damp*omegan
K=omegan**2

##########################################################################
ndof=1


d,v,a =Newmark_free(M,C,K,V0,D0,dt,NT,ndof)
 
##########################################################################

dd=np.zeros(NT,float)
vv=np.zeros(NT,float)
aa=np.zeros(NT,float)


for i in range (0,NT):
    dd[i]=d[0,i]
    vv[i]=v[0,i]    
    aa[i]=a[0,i]
    
if(iu==1):
    aa=aa/386
else:
    aa=aa/9.81

if(fn <=10):
    fna=round(fn,2)

if(fn>10 and fn <=100):
    fna=round(fn,1)

if(fn>100):
    fna=round(fn)


dtitle_string='Displacement Response fn='+str(fna)+' Hz  Q='+str(Q)
vtitle_string='Velocity Response fn='+str(fna)+' Hz  Q='+str(Q)
atitle_string='Acceleration Response fn='+str(fna)+' Hz  Q='+str(Q)


for i in range(1,200):
    if(Q==float(i)):
        dtitle_string='Displacement Response fn='+str(fna)+' Hz  Q='+str(i)
        vtitle_string='Velocity Response fn='+str(fna)+' Hz  Q='+str(i)
        atitle_string='Acceleration Response fn='+str(fna)+' Hz  Q='+str(i)
        break;  


if(iu==1):
    time_history_plot(t,dd,1,'Time(sec)','Disp(in)',dtitle_string,'disp')
    time_history_plot(t,vv,2,'Time(sec)','Vel(in/sec)',vtitle_string,'vel')
else:
    time_history_plot(t,dd,1,'Time(sec)','Disp(n)',dtitle_string,'disp')
    time_history_plot(t,vv,2,'Time(sec)','Vel(m/sec)',vtitle_string,'vel')    
  
time_history_plot(t,aa,3,'Time(sec)','Accel(G)',atitle_string,'accel')  

print " "
print " view plots"

plt.show()