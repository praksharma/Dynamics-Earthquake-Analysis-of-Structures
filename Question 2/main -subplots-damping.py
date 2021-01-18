import numpy as np
import matplotlib.pyplot as plt
import scipy as sy
import scipy.fftpack as syfp

# Importing the data from file
time = np.loadtxt( 'elcentro.dat' )[:,0]
ag=np.loadtxt( 'elcentro.dat' )[:,1]    # Force=m*ag=ag (mass=1)
plt.close("all")
# plt.figure(0)
# plt.plot(time,ag)
# plt.xlabel('Time')
# plt.ylabel('Ground acceleration')
# plt.savefig('Given_data.jpg',dpi=200)


# Initial conditions
uInitial=0
udInitial=0

# Input parameters
#omega=1         # Angular frequency
m=1             # Mass
gamma=1/2       # for Newmark method
beta=1/4 
chiArray=np.array([0.01,0.05,0.1]) 
#chi=0.05
T=5
timeStep=time[2]-time[1]

#Initializing subplots
fig,ax=plt.subplots(int(np.shape(chiArray)[0]),1)


for j in range(0,len(chiArray)):
    chi=chiArray[j]
    # Calculating parameters
    omega=2*np.pi/T
    c=2*chi*(omega)     # Damping when mass is 1
    k=np.power(omega, 2)
    
    
    # Using Newmark's method
    
    u0=np.zeros(np.shape(time))
    ud0=np.zeros(np.shape(time))
    udd=np.zeros(np.shape(time))
    ud=np.zeros(np.shape(time))
    u=np.zeros(np.shape(time))
    

    for i in range(np.shape(time)[0]-1):
        C1=u[i]+timeStep*ud[i]+0.5*np.power(timeStep,2)*(1-2*beta)*udd[i]
        C2=ud[i]+timeStep*(1-gamma)*udd[i]
        
        udd[i+1]=(ag[i+1]-c*C2-k*C1)/(1+gamma*c*timeStep+beta*k*np.power(timeStep,2))
        ud[i+1]=ud[i]+timeStep*((1-gamma)*udd[i]+gamma*udd[i+1])
        u[i+1]=u[i]+timeStep*ud[i]+0.5*np.power(timeStep,2)*((1-2*beta)*udd[i]+2*beta*udd[i+1])
    ax[j].plot(time,u,linestyle='solid') 
    ax[j].title.set_text(r'$\xi =$'+str(chi))

# error=ag-udd
# plt.figure(1)
# plt.plot(time,u)
# plt.xlabel('Time')
# plt.ylabel('Displacement (u)')  

# plt.savefig('FIG1_u.jpg',dpi=200)

# plt.figure(2)
# plt.plot(time,ud)
# plt.xlabel('Time')
# plt.ylabel('Velocity ($\\dot{u}$)')
# plt.savefig('FIG2_ud.jpg',dpi=200)

# plt.figure(3)
# plt.plot(time,udd)
# plt.xlabel('Time')
# plt.ylabel('Acceleration ($\\ddot{u}$)')
# plt.savefig('FIG3_udd.jpg',dpi=200)


# Tarray=np.arange(0.1,10+0.1,0.1) # Array of time period
# u0=np.zeros(np.shape(time))
# ud0=np.zeros(np.shape(time))
# udd=np.zeros(np.shape(time))
# ud=np.zeros(np.shape(time))
# u=np.zeros(np.shape(time))

# peakDisp=np.zeros(np.shape(Tarray)[0])
# peakAcceleration=np.zeros(np.shape(Tarray)[0])

# # Here np.shape(Tarray)[0]-1 is not needed because j does not need 
# # to access the j+1 index in the loop. so np.shape(Tarray)[0] will work
# # Also see the size of variable
# for j in range(0,np.shape(Tarray)[0]): # Loop for different timeperiod
#     T=Tarray[j]
#     # Calculating parameters
#     omega=2*np.pi/T
#     c=2*chi*(omega)     # Damping when mass is 1
#     k=np.power(omega, 2)

#     for i in range(np.shape(time)[0]-1): # newmark solver
#         C1=u[i]+timeStep*ud[i]+0.5*np.power(timeStep,2)*(1-2*beta)*udd[i]
#         C2=ud[i]+timeStep*(1-gamma)*udd[i]
        
#         udd[i+1]=(ag[i+1]-c*C2-k*C1)/(1+gamma*c*timeStep+beta*k*np.power(timeStep,2))
#         ud[i+1]=ud[i]+timeStep*((1-gamma)*udd[i]+gamma*udd[i+1])
#         u[i+1]=u[i]+timeStep*ud[i]+0.5*np.power(timeStep,2)*((1-2*beta)*udd[i]+2*beta*udd[i+1])
    

#     peakDisp[j]=np.max(u)
#     peakAcceleration[j]=np.max(udd)
    
# plt.figure(4)
# plt.plot(Tarray,peakDisp)
# plt.plot(Tarray,peakAcceleration)
# plt.xlabel('Time period (T)')
# plt.ylabel('Peak values')
# plt.legend(['Peak Displacement','Peak Acceleration'])
# plt.savefig('FIG4_responce_spectra.jpg',dpi=200)
    
# # Saving response spectra for question 3
# Tarray=np.reshape(Tarray,(len(Tarray),1))
# peakDisp=np.reshape(peakDisp,(len(peakDisp),1))
# ResponseSpectraDisp=np.concatenate((Tarray,peakDisp),axis=1)
# np.save('DispResponseSpectra',ResponseSpectraDisp)

# peakAcceleration=np.reshape(peakAcceleration,(len(peakAcceleration),1))
# ResponseSpectraAcc=np.concatenate((Tarray,peakAcceleration),axis=1)
# np.save('AccResponseSpectra',ResponseSpectraAcc)

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1)

# left  = 0.125  # the left side of the subplots of the figure
# right = 0.9    # the right side of the subplots of the figure
# bottom = 0.1   # the bottom of the subplots of the figure
# top = 0.9      # the top of the subplots of the figure
# wspace = 0.2   # the amount of width reserved for blank space between subplots
# hspace = 0.2   # the amount of height reserved for white space between subplots

# x and y common labels https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.text.html
fig.text(0.5, 0.02, 'Time', ha='center')
#fig.text(0.04, 0.5, 'Displacement (u)', va='center', rotation='vertical')
ax[1].set_ylabel('Displacement (u)')



# plt.xlabel('Time Steps')
# plt.ylabel('Maximum kinetic energy ')
plt.savefig('FIG1.jpg', dpi=200)
plt.figure(2)
fft1=np.fft.fft(ag)
plt.plot(fft1)