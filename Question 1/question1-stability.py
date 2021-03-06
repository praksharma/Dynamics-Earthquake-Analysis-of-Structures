import numpy as np
import matplotlib.pyplot as plt
import scipy as sy
import scipy.fftpack as syfp


plt.close("all")
# Initial conditions
uInitial=0      # Inital displacement
udInitial=1    # Intial velocity
 
# Input parameters
omega=5         # Angular frequency
m=1             # Mass
t=np.array([0.1,0.5,1,5])      # time step array
dT=np.zeros(int(np.shape(t)[0]))
tArray=np.zeros(int(np.shape(t)[0]))
gamma=1/2       # for Newmark method
beta=1/12
# maxKEArray=np.zeros(int(np.shape(t)[0]),drange.size)     # For storing the max KE for each timestep
fig,ax=plt.subplots(int(np.shape(t)[0]),1)

for j in range(0,len(t)) :
    #beta=betaArray[k]
    timeStep=t[j]           #time step
    # Exact solution
    T=2*np.pi/omega # Cyclic frequency
    #t=np.linspace(0,200*T,2)  # Time interval
    multiple=100
    drange=np.arange(0,multiple*T,timeStep) 
    umax=pow(pow(uInitial,2)+pow(udInitial/omega,2),0.5)
    uExact=umax*np.sin(omega*drange)    # Exact solution
    
    # plt.figure(1)
    # plt.grid(True)                                  # Enable grids
    # plt.rc('grid', linestyle="--", color='black')   # Darker grids
    #plt.plot(drange,uExact)
    
    # Newmark's solution
    
    # initializing vectors
    u=np.zeros(drange.size)
    ud=np.zeros(drange.size)
    udd=np.zeros(drange.size)
    KE=np.zeros(drange.size)
    
    # Applying initial conditions
    u[0]=uInitial
    ud[0]=udInitial
    KE[0]=0.5*m*np.power(ud[0],2)
    
    
    # python for loop does not support the use of float step-size
    #drange=np.arange(0,100*T,timeStep)  # limits for the loop
    
    # t.size-1 is used otherwise will exceed bounds limit for (i+1)
    for i in range(0,drange.size-1):
        C1=u[i]+timeStep*ud[i]+0.5*np.power(timeStep,2)*(1-2*beta)*udd[i]
        C2=ud[i]+timeStep*(1-gamma)*udd[i]
        
        udd[i+1]=-np.power(omega,2)*C1/(1+beta*np.power(omega,2)*np.power(timeStep,2))
        ud[i+1]=ud[i]+timeStep*((1-gamma)*udd[i]+gamma*udd[i+1])
        u[i+1]=u[i]+timeStep*ud[i]+0.5*np.power(timeStep,2)*((1-2*beta)*udd[i]+2*beta*udd[i+1])
        KE[i+1]=0.5*m*np.power(ud[i+1],2)
        
    
    ax[j].plot(drange,KE,linestyle='solid') 
    ax[j].title.set_text(r'$\Delta t=$'+str(timeStep))
    #print(len(KE))
    # maxKE[j]=np.amax(KE)
    
    # #plt.plot(drange,u)
    # # plt.xlabel('time (seconds)')
    # # plt.ylabel('Displacement (u)')
    # # plt.legend(['Exact solution','Newmark\'s solution'])
    # # plt.figure(2)
    # # #plt.plot(drange,KE)
    # # plt.xlabel('time (seconds)')
    # # plt.ylabel('Kinetic energy')
    
    # # Calculating the period error
    
    # freqdist = np.fft.fftfreq(int(len(u)))
    # #freqdist=(1/(timeStep*int(np.shape(drange)[0])))*np.arange(int(np.shape(drange)[0]))
    
    # # plt.figure(3)
    # FFTNewmark=np.fft.fft(u)
    # FFTNewmark=FFTNewmark[:int(np.shape(FFTNewmark)[0]/2)]
    # #plt.plot(abs(FFTNewmark))
    # maxFreqNewmarkIndex=np.argmax(abs(FFTNewmark))
    # maxFreqNewmark=freqdist[maxFreqNewmarkIndex]
    # TNewmark=abs(1/maxFreqNewmark)
    
    
    # FFTExact=np.fft.fft(uExact)
    # FFTExact=FFTExact[:int(np.shape(FFTExact)[0]/2)]
    # #plt.plot(abs(FFTExact))
    # # plt.legend(['Newmark\'s solution','Exact solution'])
    
    # maxFreqExactIndex=np.argmax(abs(FFTExact))
    # maxFreqExact=freqdist[maxFreqExactIndex]
    # TExact=abs(1/maxFreqExact)
    
    # cycles=drange[-1]/TNewmark
    # TExactEnd=abs(TExact)*cycles
    # #TExactEnd=abs(TExact)*cycles
    # deltaT=(drange[-1]-TExactEnd)
    # dT[j]=deltaT/T
    # #dT[j]=deltaT/TExact
    # tArray[j]=timeStep/T
    # #tArray[j]=timeStep/TExact
    # #print('t=',timeStep)
    # #print(deltaT/TExact)
    # rate=maxFreqExactIndex/maxFreqNewmarkIndex
    # TNewmark=T/rate
    # cycles=drange[-1]/TNewmark
    # TExactEnd=abs(T)*cycles
    # #deltaT=(drange[-1]-TExactEnd)
    # deltaT=(T-TNewmark)
    # dT[j]=deltaT/T

# plt.figure(4)    
# plt.plot(tArray,dT)
#plt.legend(['$\beta$=0','$\beta$=1/12','$\beta$=1/6','$\beta$=1/4','$\beta$=1/3'])

# plt.grid(True)                                  # Enable grids
# plt.rc('grid', linestyle="--", color='black')
# plt.legend(['$\\beta$=0','$\\beta$=1/12','$\\beta$=1/6','$\\beta$=1/4','$\\beta$=1/3'])
# plt.xlabel(r'$\frac{\Delta t}{T}$')
# plt.ylabel('$\\frac{\Delta T}{T}$')
#plt.savefig('Period error.jpg', dpi=200)

# plt.figure(5)
# plt.plot(drange,KE,linestyle='solid')

# Preventing overlaps between title
#fig.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.09)

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1)

# left  = 0.125  # the left side of the subplots of the figure
# right = 0.9    # the right side of the subplots of the figure
# bottom = 0.1   # the bottom of the subplots of the figure
# top = 0.9      # the top of the subplots of the figure
# wspace = 0.2   # the amount of width reserved for blank space between subplots
# hspace = 0.2   # the amount of height reserved for white space between subplots

# x and y common labels https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.text.html
fig.text(0.5, 0.02, 'Time integration interval', ha='center')
fig.text(0.04, 0.5, 'Kinetic Energy', va='center', rotation='vertical')




# plt.xlabel('Time Steps')
# plt.ylabel('Maximum kinetic energy ')
plt.savefig('Stability_'+str(multiple)+'T_Beta='+str(beta)+'.jpg', dpi=200)