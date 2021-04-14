import numpy as np 
import scipy 
import math 
import matplotlib.pyplot as plt 
from matplotlib.widgets import Slider 
from scipy.interpolate import interp1d 
from scipy.optimize import curve_fit, least_squares 

def funcxn(r, xn): 
    return 4 * r * xn * (1 - xn)

r = 0.958 
x0 = 0.54 
n = 100 
N = np.arange(0, n, 1) 
X = np.zeros([n]) 
X[0] = x0 
for i in range(1, n): 
    X[i] = funcxn(r, X[i - 1]) 

fig = plt.figure(figsize=(12, 6)) 
plt.plot(N, X) 

plt.title('Отображение Xn+1 = 4*r*Xn*(1-Xn)', fontsize=18) 
plt.ylabel('xn', fontsize=18) 
plt.xlabel('n', fontsize=18) 
ax = plt.gca() 
ax.tick_params(labelsize=14) 
plt.show() 

x01 = 0.76 
N = np.arange(0, n, 1) 
X1 = np.zeros([n]) 
X1[0] = x01 
for i in range(1, n): 
    X1[i] = funcxn(r, X1[i - 1]) 

def funcdefxn(xn, xn1): 
    return abs(xn - xn1) * 100 / xn 

f = np.zeros([n]) 
f[0] = abs(X[0] - X[0]) / X[0] 
for i in range(0, n): 
    f[i] = funcdefxn(X[i], X1[i]) 

fig = plt.figure(figsize=(12, 6)) 
plt.plot(N, f) 
   
plt.title('Стандартное отклонение Xn при изменении начального значения', fontsize=18) 
plt.ylabel('|xn - x1n|/xn * 100', fontsize=18) 
plt.xlabel('n', fontsize=18) 
ax = plt.gca() 
ax.tick_params(labelsize=14) 
plt.show() 
   
import random 
   
xn = [0.001] 
a = np.linspace(0.70, 0.98, 400000) 
for i in range(1, a.shape[0]): 
    xn.append(a[i] * np.sin(np.pi * xn[i - 1]) + 0.001 * random.random()) 

fig = plt.figure(figsize=(12, 6)) 
plt.scatter(a, xn, s=0.1) 
plt.title('Бифуркационная диаграмма', fontsize=18) 
plt.ylabel('xn', fontsize=18) 
plt.xlabel('r', fontsize=18) 
ax = plt.gca() 
ax.tick_params(labelsize=14) 
plt.show()
