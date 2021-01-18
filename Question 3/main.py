import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.optimize import fsolve

# Input Parameters
a=0.35               # Dimension of square column
area=np.square(a)    # Area of C/S of column in m^2
I=np.power(a,4)/12   # Moment of inertia in m^4
ndof=3               # Number of degree of freedom in the system
h=3                  # Interfloor space in m
E=25e9               # Young's modulus of concrete in N/m^2

# Mass of each floor
m1=400000
m2=300000
m3=200000

# Loading the Response Spectra from Question 2
PeriodArray=np.load('DispResponseSpectra.npy')[:,0]
PeakDispArray=np.load('DispResponseSpectra.npy')[:,1]
PeakAcclnArray=np.load('AccResponseSpectra.npy')[:,1]

# Lumped mass matrix
M=np.diag([m1,m2,m3])   # Mass matrix
print('Mass matrix (kg):\n'+str(M)+'\n')

# Lateral stiffness
k=12*E*I/np.power(h,3)  # Stiffness for one column per floor
k=6*k                   # Stiffness per floor
K=np.matrix([[2*k,-k,0],[-k,2*k,-k],[0,-k,k]])  #Stiffness
print('Stiffness matrix (N/m):\n'+str(K)+'\n')
print('Moment of inetria (m^4): '+str(I)+'\n')

# MODULE 1: using eigenvalue solution-Exact solution------------------
print('**************************EIGENVALUE SOLUTION*********************\n' )
OmegaSquare,EigVector=np.linalg.eig(K*np.linalg.inv(M))
OmegaExact=np.sqrt(OmegaSquare)          # Natural frequency
print('Omega (1/s): '+str(OmegaExact))

V1=EigVector[:,0]/EigVector[0,0]    # Scale the modal shape
V2=EigVector[:,1]/EigVector[0,1]
V3=EigVector[:,2]/EigVector[0,2]
# np.insert transpose the matrix, will use transpose again
V1plot=np.transpose(np.insert(V1,0,0)) #inserting zero at the start of array for plotting
V2plot=np.transpose(np.insert(V2,0,0))
V3plot=np.transpose(np.insert(V3,0,0))

xArray=np.arange(np.shape(V1plot)[0])

# Mode plots
fig,ax=plt.subplots(1, 3,sharey=True)

ax[0].grid(color='k', linestyle='-', linewidth=1)
ax[0].plot(V3plot,xArray)
ax[0].set_yticklabels([])
ax[0].set_xlabel('1st mode')

ax[1].grid(color='k', linestyle='-', linewidth=1)
ax[1].plot(V2plot,xArray)
ax[1].set_xlabel('2nd mode')

ax[2].grid(color='k', linestyle='-', linewidth=1)
ax[2].plot(V1plot,xArray)
ax[2].set_xlabel('3rd mode')

plt.savefig('Modes_exact.jpg',dpi=200)
OmegaExact=np.flip(OmegaExact)
T=2*np.pi/OmegaExact                 # time period 
 
# Displacement calculation
mBar2=np.zeros((ndof,1))
mBar2[0,0]=(np.matmul(np.matmul(np.transpose(V1),M),V1))    # Modal mass  
mBar2[1,0]=(np.matmul(np.matmul(np.transpose(V2),M),V2))    # Modal mass  
mBar2[2,0]=(np.matmul(np.matmul(np.transpose(V3),M),V3))
    
Sa2=np.zeros(ndof) 
Su2=np.zeros(ndof) 
T=2*np.pi/OmegaExact                 # time period 

# Using linear interpolation to find displacement corresponding to time period (T)
for j in range(0,int(np.shape(Sa2)[0])):   # for each Sa2
    for i in range(0,int(np.shape(PeriodArray)[0])):  # searching over period
         if PeriodArray[i]>T[j]:    # Value after T i.e. T(i+2)
            T3=PeriodArray[i]
            T1=PeriodArray[i-1]
            Y3=PeakDispArray[i]
            Y1=PeakDispArray[i-1]
            A3=PeakAcclnArray[i]
            A1=PeakAcclnArray[i-1]
            Su2[j]=Y1+(Y3-Y1)/(T3-T1)*(T[j]-T1)     # Peak displacement corresponding to T in the response spectra      
            Sa2[j]=A1+(A3-A1)/(T3-T1)*(T[j]-T1)
            break

# Load participation factor
d=np.matrix([[1],[1],[1]])              # Earthquake input direction vector
l2=np.matmul(np.matmul(np.transpose(V1),M),d)/mBar2
print('Load participation factor: \n'+str(l2)+'\n')


# Maximum Displacement vectors
uMax2_1=l2[0,0]*Su2[0]*V1     
uMax2_2=l2[1,0]*Su2[1]*V2
uMax2_3=l2[2,0]*Su2[2]*V3


# Total maximum displacement using SRSS
uMaxExact=np.zeros(ndof)
uMaxExact[0]=np.square(uMax2_1[0,0])+ np.square(uMax2_2[0,0])+np.square(uMax2_3[0,0])
uMaxExact[1]=np.square(uMax2_1[1,0])+ np.square(uMax2_2[1,0])+np.square(uMax2_3[1,0])
uMaxExact[2]=np.square(uMax2_1[2,0])+ np.square(uMax2_2[2,0])+np.square(uMax2_3[2,0])
print('Maximum Displacment (m): '+str(uMaxExact)+'\n')


# Maximum floor force vector
F=np.zeros(ndof)
F2_1=float(l2[0,0]*Sa2[0])*np.matmul(M,V1)
F2_2=float(l2[1,0]*Sa2[1])*np.matmul(M,V2)
F2_3=float(l2[2,0]*Sa2[2])*np.matmul(M,V3)
# Using SRSS
F[0]=np.square(F2_1[0,0])+ np.square(F2_2[0,0])+np.square(F2_3[0,0])
F[1]=np.square(F2_1[1,0])+ np.square(F2_2[1,0])+np.square(F2_3[1,0])
F[2]=np.square(F2_1[2,0])+ np.square(F2_2[2,0])+np.square(F2_3[2,0])

print('Shear forces (N): '+str(F)+'\n')

# Base shear
VbExact=np.sum(F)
print('Base shear force (N): '+str(VbExact)+'\n')

# Overturning moment
z=np.arange(h,h*ndof+1,h)        # Height of floors   
MbExact=np.sum(z[:]*F[:])
print('Overturning moment (N-m): '+str(MbExact)+'\n')
                          








# MODULE 2: Using linearly increasing mode---------------------------
print('*********************LINEARLY INCREASING MODE*********************\n' )
V1=np.matrix([[1],[2],[3]])     # Linearly increasing mode

print('V1:\n'+str(V1)+'\n' )
mBar=float(np.matmul(np.matmul(np.transpose(V1),M),V1))    # Modal mass
kBar=float(np.matmul(np.matmul(np.transpose(V1),K),V1))    # Modal stiffness

Omega=np.sqrt(kBar/mBar)        # Omega approx
T=2*np.pi/Omega                 # time period 
                           
# Using linear interpolation to find displacement corresponding to time period (T)
for i in range(0,int(np.shape(PeriodArray)[0])):
     if PeriodArray[i]>T:    # Value after T i.e. T(i+2)
        T3=PeriodArray[i]
        T1=PeriodArray[i-1]
        Y3=PeakDispArray[i]
        Y1=PeakDispArray[i-1]
        A3=PeakAcclnArray[i]
        A1=PeakAcclnArray[i-1]
        break
Su=Y1+(Y3-Y1)/(T3-T1)*(T-T1)     # Peak displacement corresponding to T in the response spectra      
Sa=A1+(A3-A1)/(T3-T1)*(T-T1)
        
# Load participation factor
d=np.matrix([[1],[1],[1]])       # Earthquake input direction vector
l1=float(np.matmul(np.matmul(np.transpose(V1),M),d)/mBar)
print('Load participation factor: \n'+str(l1)+'\n')
# Maximum Displacement
uMax=l1*Su*V1                    # Umax for each floor
print('Maximum Displacment(m): '+str(uMax)+'\n')    
# Base shear
totalMass=m1+m2+m3               # Total mass of the structure
Vb=totalMass*Sa                  # Base Shear force
print('Base shear force (N): '+str(Vb)+'\n')

# Floor shear force

z=np.arange(h,h*ndof+1,h)        # Height of floors
F1=np.zeros(ndof)
m=np.array([m1,m2,m3])           # Array is mass
denominator=np.dot(m,z)
 
for i in range(0,int(np.shape(F)[0])):
   F1[i]=Vb*m[i]*z[i]/denominator
print('Shear forces (N): '+str(F1)+'\n')   
# Overturning moment
  
Mb1=np.sum(z[:]*F1[:])
print('Overturning moment (N-m): '+str(Mb1)+'\n')   










# MODULE 3: Using two ritz vector------------------------
print('********************TWO RITZ VECTOR WITH SRSS********************\n' )
r1= np.matrix([[1],[2],[3]])     # Linearly increasing mode
r2=np.matrix([[1],[4],[9]])      # Quadratically increasing mode 
print('R1:\n'+str(r1)+'\n')
print('R2:\n'+str(r2)+'\n')
 
R=np.append(r1,r2,1)
M_Hat=np.matmul(np.matmul(np.transpose(R),M),R)
K_Hat=np.matmul(np.matmul(np.transpose(R),K),R)

OmegaSquare2,EigVector2=np.linalg.eig(K_Hat*np.linalg.inv(M_Hat))
Omega2=np.sqrt(OmegaSquare2)          # Natural frequency
x1=EigVector2[:,0]/EigVector2[0,0]    # Scale the modal shape
x2=EigVector2[:,1]/EigVector2[0,1]

V1=np.matmul(R,x1)
V2=np.matmul(R,x2)

mBar2=np.zeros((2,1))
mBar2[0,0]=(np.matmul(np.matmul(np.transpose(V1),M),V1))    # Modal mass  
mBar2[1,0]=(np.matmul(np.matmul(np.transpose(V2),M),V2))    # Modal mass  
    
Sa2=np.zeros(2) 
Su2=np.zeros(2) 
T=2*np.pi/Omega2                 # time period 

# Using linear interpolation to find displacement corresponding to time period (T)
for j in range(0,int(np.shape(Sa2)[0])):
    for i in range(0,int(np.shape(PeriodArray)[0])):
         if PeriodArray[i]>T[j]:    # Value after T i.e. T(i+2)
            T3=PeriodArray[i]
            T1=PeriodArray[i-1]
            Y3=PeakDispArray[i]
            Y1=PeakDispArray[i-1]
            A3=PeakAcclnArray[i]
            A1=PeakAcclnArray[i-1]
            Su2[j]=Y1+(Y3-Y1)/(T3-T1)*(T[j]-T1)     # Peak displacement corresponding to T in the response spectra      
            Sa2[j]=A1+(A3-A1)/(T3-T1)*(T[j]-T1)
            break

# Load participation factor
d=np.matrix([[1],[1],[1]])              # Earthquake input direction vector
l2=np.matmul(np.matmul(np.transpose(V1),M),d)/mBar2
print('Load participation factor: \n'+str(l2)+'\n')

# Maximum Displacement vectors
uMax2_1=l2[0,0]*Su2[0]*V1     
uMax2_2=l2[1,0]*Su2[1]*V2


# Total maximum displacement using SRSS

uMax2=np.zeros(ndof)
uMax2[0]=np.square(uMax2_1[0,0])+ np.square(uMax2_2[0,0])
uMax2[1]=np.square(uMax2_1[1,0])+ np.square(uMax2_2[1,0])
uMax2[2]=np.square(uMax2_1[2,0])+ np.square(uMax2_2[2,0])
print('Maximum Displacment (m): '+str(uMax2)+'\n')

# Maximum floor force vector
F2=np.zeros(ndof)
F2_1=float(l2[0,0]*Sa2[0])*np.matmul(M,V1)
F2_2=float(l2[1,0]*Sa2[1])*np.matmul(M,V2)
# Using SRSS
F2[0]=np.square(F2_1[0,0])+ np.square(F2_2[0,0])
F2[1]=np.square(F2_1[1,0])+ np.square(F2_2[1,0])
F2[2]=np.square(F2_1[2,0])+ np.square(F2_2[2,0])
print('Shear forces (N): '+str(F2)+'\n')

# Base shear
Vb_2=np.sum(F)
print('Base shear force (N): '+str(Vb_2)+'\n')

# Overturning moment
z=np.arange(h,h*ndof+1,h)        # Height of floors   
Mb2=np.sum(z[:]*F2[:])
print('Overturning moment (N-m): '+str(Mb2)+'\n')   