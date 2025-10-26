import numpy as np
from math import*

ob = radians(23.43829194) 
r_vec = np.array([0.61080832,0.24293337, -0.05423979])
r_dot_vec = np.array([-0.90753098, -0.22522702, 0.14404739])
t = 2459399.66766
r_mag= (r_vec[0]**2 + r_vec[1]**2 + r_vec[2]**2)**0.5
r_dot_mag= (r_dot_vec[0]**2+ r_dot_vec[1]**2+ r_dot_vec[2]**2)**0.5

a = ((2/r_mag) - np.dot(r_dot_mag,r_dot_mag))**-1
print(f'a = {a}')

h= np.cross(r_vec, r_dot_vec) 
h_mag= (h[0]**2+ h[1]**2+h[2]**2)**0.5
e = (1-(h_mag**2/a))**0.5
print(f'e = {e}')

matrix_ec = np.array([[1,0,0],[0,cos(ob),sin(ob)], [0, -sin(ob), cos(ob)]])   
h_ec = np.matmul(matrix_ec, h)  #converted h to ecliptic
r_ec = np.matmul(matrix_ec, r_vec)  #converted r to ecliptic
r_dot_ec = np.matmul(matrix_ec, r_dot_vec) #converted r dot to ecliptic
r_ec_mag= (r_ec[0]**2+r_ec[1]**2+r_ec[2]**2)**0.5
r_dot_ec_mag= (r_dot_ec[0]**2+r_dot_ec[1]**2+r_dot_ec[2]**2)**0.5
i= acos((h_ec[2])/h_mag)
print(f'i = {degrees(i)}')

sin_o= h_ec[0]
cos_o= -h_ec[1]
Ω= atan2(sin_o,cos_o)
print(f'Ω = {degrees(Ω)}')

sinf_w= r_ec[2] / (r_ec_mag*sin(i))
cosf_w= (1/cos(Ω))*(r_ec[0]/r_ec_mag+(cos(i)*sinf_w*sin(Ω)))
f_w= atan2(sinf_w,cosf_w)

sinf= (((np.dot(r_ec, r_dot_ec))/(e*r_ec_mag))*((a*(1-e**2))**0.5))
cosf = (1/e)*(((a*(1-e**2))/r_ec_mag)-1)
f= atan2(sinf,cosf)
w= f_w-f
print(f'w = {degrees(w)}')

E = acos((1/e)*(1-(r_mag/a)))
if f < 0:
    E1= E*-1
if degrees(E) < 360:
    E2 = degrees(E1)+360

M = radians(E2)- e*sin(radians(E2))
print(f'M = {degrees(M)}')
print(f'E = {E2}')

#collab: ari and theodor