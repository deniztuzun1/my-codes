import numpy as np
from math import*

tau1 = float(input("tau1 value:"))
tau3 = float(input("tau3 value:"))
r2x = float(input("r2x value:"))
r2y = float(input("r2y value:"))
r2z = float(input("r2z value:"))
r2_dotx = float(input("r2_dotx value:"))
r2_doty = float(input("r2_doty value:"))
r2_dotz = float(input("r2_dotz value:"))
flag = (input("function, 3rd-order series, or 4th order series:"))
r2_mag = (r2x**2 + r2y**2 + r2z**2)**0.5
r2_dot_mag = (r2_dotx**2+r2_doty**2+r2_dotz**2)**0.5
r2 = np.array([r2x, r2y, r2z])
r2_dot = np.array([r2_dotx, r2_doty, r2_dotz])
mu = 1
u = mu / r2_mag**3
z = (np.dot(r2, r2_dot)) / r2_mag**2
q = ((np.dot(r2_dot,r2_dot)) / r2_mag**2) - u 
a = ((2/r2_mag) - np.dot(r2_dot_mag,r2_dot_mag))**-1
h= np.cross(r2, r2_dot) 
h_mag= (h[0]**2+ h[1]**2+h[2]**2)**0.5
e = (1-(h_mag**2/a))**0.5
n = 1/ a**1.5

def delta_E(r2, r2_dot, tau):
    sign_match = ((np.dot(r2,r2_dot))/(n*a**2))*cos((n*tau)-((np.dot(r2,r2_dot))/(n*a**2))) + ((1-(r2_mag/a))*sin((n*tau)-((np.dot(r2,r2_dot))/(n*a**2))))
    if e <= 0.1:
        x0 = n*tau
    elif sign_match<0:
        x0 = n*tau-0.85*e-((np.dot(r2,r2_dot))/(n*a**2))
    elif sign_match>0:
        x0 = n*tau-0.85*e+((np.dot(r2,r2_dot))/(n*a**2))
    xn= x0
    x0= 0
    while abs(xn-x0) > 1*10**-12:
        x0 = xn
        f_x = x0 - ((1-(r2_mag/a))*sin(x0)) + (((np.dot(r2,r2_dot))/(n*a**2))*(1-cos(x0))) - (n*tau)
        f_dx = 1 - ((1-(r2_mag/a))*cos(x0)) + ((np.dot(r2,r2_dot))/(n*a**2)*sin(x0))
        xn = x0 - (f_x/ f_dx)
    return xn 

def f(tau, r2, r2_dot , flag):
    z = (np.dot(r2, r2_dot)) / r2_mag**2
    if flag == '3rd-order series':
        f = 1 - (0.5* u * tau**2) + (0.5* u * z * tau**3) 
    elif flag == '4th order series': 
        f = 1 - (0.5* u * tau**2) + (0.5* u * z * tau**3) + (1/24* ((3*u*q)-(15*u*z**2)+(u**2)) * tau**4)
    elif flag == 'function':
        if tau == tau1:
            DE1 = delta_E(r2, r2_dot, tau)
            f = 1 - (a/r2_mag)*(1- cos(DE1))
        elif tau == tau3:
            DE3 = delta_E(r2, r2_dot, tau)
            f = 1 - (a/r2_mag)*(1- cos(DE3))
    return f

def g(tau, r2, r2_dot , flag):
    z = (np.dot(r2, r2_dot)) / r2_mag**2
    if flag == '3rd-order series':
        g = tau - ((1/6)*u*tau**3)
    elif flag == '4th order series': 
        g = tau - ((1/6)*u*tau**3) + (0.25* u*z*tau**4)
    elif flag == 'function':
        if tau == tau1:
            DE1 = delta_E(r2, r2_dot, tau)
            g = tau + (1/n)*(sin(DE1)-DE1)
        elif tau == tau3:
            DE3 = delta_E(r2, r2_dot, tau)
            g = g = tau + (1/n)*(sin(DE3)-DE3)
    return g

f1 = f(tau1,r2,r2_dot,flag)
f3 = f(tau3,r2,r2_dot,flag)
g1 = g(tau1,r2,r2_dot,flag)
g3 = g(tau3,r2,r2_dot,flag)

print(f'f1 = {f1}')
print(f'f3 = {f3}')
print(f'g1 = {g1}')
print(f'g3 = {g3}')