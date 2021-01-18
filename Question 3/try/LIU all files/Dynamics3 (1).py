#dynamics 3
import numpy as np
import math as ma
import matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.fftpack import fft,ifft
import time

m1 = 400
m2 = 300
m3 = 200
diemension_column = 0.35
E = 25*1000000000 # GN/M^2
height = 3
h_cubic = np.power(height, 3)
I = np.power(diemension_column, 4)/12
k = 12*(E*I)*6/h_cubic
K_matrix = np.zeros((3,3))
K_matrix[0,0] = k*2
K_matrix[1,1] = k*2
K_matrix[2,2] = k
K_matrix[0,1] = -k
K_matrix[1,0] = -k
K_matrix[0,2] = -k
K_matrix[2,0] = -k

# mass matrix
M_matrix = np.zeros((3,3))
M_matrix[0,0] = 4
M_matrix[1,1] = 3
M_matrix[2,2] = 2
M_matrix = M_matrix*100000
# Ritz Vector
r = [[1], [4], [9]]
r = np.array(r)
#
M_hat = np.dot(np.dot(r.transpose(), M_matrix), r)
K_hat = np.dot(np.dot(r.transpose(), K_matrix), r)
omega = np.sqrt(K_hat[0,0]/M_hat[0,0])
T = 2*ma.pi/omega
# print(omega)
print(T)
acc_spectrum = np.load('acc_spectrum.npy')
dis_spectrum = np.load('dis_spectrum.npy')
period = np.load('period.npy')

## interplotation
def find_peak_dis():
	x1 = period[41]
	y1 = dis_spectrum[41]
	x3 = period[42]
	y3 = dis_spectrum[42]
	x2 = T
	y2 = y1 + (x2-x1) * (y3 - y1)/(x3 - x1)
	return y2
def find_peak_acc():
	x1 = period[41]
	y1 = acc_spectrum[41]
	x3 = period[42]
	y3 = acc_spectrum[42]
	x2 = T
	y2 = y1 + (x2-x1) * (y3 - y1)/(x3 - x1)
	return y2
peak_dis = find_peak_dis()
peak_acc = find_peak_acc()

V = M_hat*peak_acc
sum_M = (M_matrix[0,0]*3 + M_matrix[1,1]*6 + M_matrix[2,2]*9)
F1 = V * (M_matrix[0,0]*3)/sum_M
F2 = V * (M_matrix[1,1]*6)/sum_M
F3 = V * (M_matrix[2,2]*9)/sum_M
F = np.hstack((F1,F2,F3))
print(peak_dis)

## question 2
r1= [[1], [4], [9]]
r1= np.array(r1)
M_hat_2 = np.dot(np.dot(r1.transpose(), M_matrix), r1)
K_hat_2 = np.dot(np.dot(r1.transpose(), K_matrix), r1)
