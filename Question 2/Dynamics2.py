import numpy as np
import math as ma
import matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import time
# find peak value of displacement, veloctiy, total acceleration for different periods.  
start = time.process_time()

filename = 'elcentro.dat'
indata = np.loadtxt(filename, usecols=(0,1))
time = (indata[:,0])
excit_F = indata[:,1]
# earth_data = np.sin(time)
earth_data = np.multiply((indata[:,1]), 9.8)

# define paremeter

def trange(start, final, dt):
	trange = np.arange(start, final, dt)
	# print(len(trange), 'with', dt,'and', final)
	return trange

def compute_amp(u_c, u_cd, u_cdd, dt, finaltime, c, k, beta, gamma):
	timerange = np.arange(0, finaltime, dt)
	amp = [u_c]
	velocity = [u_cd]
	acceleration = [u_cdd]
	for i in range(0,len(timerange)-1):
		u0_np1 = u_c + dt*u_cd + 0.5*(dt**2)*(1-2*beta)*u_cdd
		u0_dnp1 = u_cd + dt*(1 - gamma)*u_cdd

		u_ddnp1	=(excit_F[i] - c * u0_dnp1 - k*u0_np1)/(1 + gamma * c*dt + beta * k * dt**2)

		u_np1 = u_c + dt*u_cd + 0.5*(dt**2)*((1-2 * beta) * u_cdd + 2*beta*u_ddnp1)
		u_dnp1 = u_cd + dt*((1 - gamma)* u_cdd + gamma*u_ddnp1)

		amp.append(u_np1)
		velocity.append(u_dnp1)
		acceleration.append(u_ddnp1)
		u_c, u_cd, u_cdd = u_np1,u_dnp1, u_ddnp1
	
	return amp,	velocity, acceleration
def abc(T):
	# for k in r:
	frequency = (2*np.pi)/T
	beta = 1/4
	gamma = 1/2
	dt = 0.02
	m = 1
	omega = frequency
	chi = 50
	k = np.square(omega)
	c = chi * omega
	# dk = k + (gamma/(beta*dt)) * c + (1/ (beta* np.square(dt)))* m
	u_0 = 0
	ud_0 = 0
	udd_0 = 0
	
	displacement, veloctiy, acceleration = compute_amp(u_0, ud_0, udd_0, dt, 31.18, c, k, beta, gamma)
	displacement = np.array(displacement)
	veloctiy = np.array(veloctiy)
	acceleration = np.array(acceleration)
    
	dis = np.amax(displacement)
	vel = np.amax(veloctiy)
	acc =np.amax(acceleration)
	return dis, vel, acc

# r = np.arange(10, 0.1 ,-0.1)
max_dis= []
max_vel = []
max_acc = []
T = np.arange(1,10,0.1)
for DT in T:
	dis = abc(T)[0]
	vel = abc(T)[1]
	acc = abc(T)[2]
	max_dis.append(dis)
	max_vel.append(vel)
	max_acc.append(acc)


plt.figure(1)
plt.plot(max_dis)

plt.xlabel('time')
plt.ylabel('Displacement')

plt.figure(2)
plt.plot(max_vel)
plt.xlabel('time')
plt.ylabel('velocity')

plt.figure(3)
plt.plot(max_acc)
plt.xlabel('time')
plt.ylabel('acceleration')




plt.figure(4)
plt.plot(earth_data)
plt.show()




