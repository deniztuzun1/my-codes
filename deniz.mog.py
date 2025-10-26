import numpy as np
from math import*

file= open("2004LU3test.txt")
year = []
month = []
day = []
time = []
date = []
ra = []
dec = []
Rx = []
Ry = []
Rz = []
for i in file:
    year.append(int(i.split(' ')[0]))
    month.append(int(i.split(' ')[1]))
    day.append(int(i.split(' ')[2]))
    time.append(i.split(' ')[3])
    ra.append(i.split(' ')[4])
    dec.append(i.split(' ')[5])
    Rx.append(i.split(' ')[6])
    Ry.append(i.split(' ')[7])
    Rz.append(i.split(' ')[8][:-1])  #made lists with each value in the text file

y1, y2 , y3 = year[0] , year[1], year[2]
m1, m2 , m3 = month[0], month[1], month[2]
d1, d2, d3 = day[0], day[1], day[2]
ra1, ra2, ra3 = ra[0], ra[1], ra[2]
dec1, dec2, dec3 = dec[0], dec[1], dec[2]
Rx1, Rx2, Rx3 = Rx[0], Rx[1], Rx[2]
Ry1, Ry2, Ry3 = Ry[0], Ry[1], Ry[2]
Rz1, Rz2, Rz3 = Rz[0], Rz[1], Rz[2]
R1 = np.array([float(Rx1), float(Ry1), float(Rz1)])
R2 = np.array([float(Rx2), float(Ry2), float(Rz2)])  #made R vectors
R3 = np.array([float(Rx3), float(Ry3), float(Rz3)])
R2_mag = (R2[0]**2+R2[1]**2+R2[2])**0.5
k = 0.0172020989484
mu = 1
ob = radians(23.43829194) 
C = 173.142168257

time_dec=[]
for i in time:
    hour = float(i.split(':')[0])
    min = float(i.split(':')[1])*10/6
    sec = float(i.split(':')[2])/36
    time_dec.append(hour+min+sec)  #calculated time in decimals
time1= time_dec[0]
time2= time_dec[1]
time3= time_dec[2]

def julian(year, month, day, time):
    J= (367*year)-int(7/4*(year+int((month+9)/12)))+int(275/9*month)+day+1721013.5  #converted to julian date
    JD=(J+(time/24))
    return JD 

t1 = julian(y1,m1,d1,time1)
t2 = julian(y2,m2,d2,time2)
t3 = julian(y3,m3,d3,time3)

t1_orig = t1
t2_orig = t2 #made original time values to use after inside the while loop for the light travel time correction
t3_orig = t3

tau3 = k*(t3-t2)
tau1 = k*(t1-t2)
tau = tau3 - tau1 

def convert(x):
    x_hour= float(x.split(':')[0])
    if int(x_hour)<0:
        x_min2hour= float(x.split(':')[1])/-60  #if the dec is negative divide to -60 and -3600 to have min and sec same sign with hour
        x_sec2hour= float(x.split(':')[2])/-3600
        x_dec= x_hour + x_min2hour + x_sec2hour
    else:
        x_min2hour= float(x.split(':')[1])/60
        x_sec2hour= float(x.split(':')[2])/3600
        x_dec= x_hour + x_min2hour + x_sec2hour
    return x_dec

ra1dec= radians(convert(ra1)*15) #multiplied with 15 to convert hours to degrees
ra2dec= radians(convert(ra2)*15)
ra3dec= radians(convert(ra3)*15)
dec1dec= radians(convert(dec1))
dec2dec= radians(convert(dec2))
dec3dec= radians(convert(dec3))

def p_hat(ra, dec):
    p_hat = np.array([cos(ra)*cos(dec), sin(ra)*cos(dec), sin(dec)]) 
    return p_hat

p_hat1 = p_hat(ra1dec,dec1dec) #rho hat calculation using the function
p_hat2 = p_hat(ra2dec,dec2dec)
p_hat3 = p_hat(ra3dec,dec3dec)

matrix_ec = np.array([[1,0,0],[0,cos(ob),sin(ob)], [0, -sin(ob), cos(ob)]]) #conversion matrix: to convert to ecliptic
p_hat1_ec = np.matmul(matrix_ec, p_hat1)
p_hat2_ec = np.matmul(matrix_ec, p_hat2)
p_hat3_ec = np.matmul(matrix_ec, p_hat3)
R1_ec = np.matmul(matrix_ec, R1)
R2_ec = np.matmul(matrix_ec, R2)
R3_ec = np.matmul(matrix_ec, R3)
R2_ec_mag = (R2_ec[0]**2+R2_ec[1]**2+R2_ec[2]**2)**0.5

D0 = np.dot(p_hat1_ec,(np.cross(p_hat2_ec,p_hat3_ec)))
D11 = np.dot((np.cross(R1_ec, p_hat2_ec)), p_hat3_ec)
D12 = np.dot((np.cross(R2_ec, p_hat2_ec)), p_hat3_ec)
D13 = np.dot((np.cross(R3_ec, p_hat2_ec)), p_hat3_ec)
D21 = np.dot((np.cross(p_hat1_ec, R1_ec)), p_hat3_ec)
D22 = np.dot((np.cross(p_hat1_ec, R2_ec)), p_hat3_ec)
D23 = np.dot((np.cross(p_hat1_ec, R3_ec)), p_hat3_ec)
D31 = np.dot(p_hat1_ec, (np.cross(p_hat2_ec,R1_ec)))
D32 = np.dot(p_hat1_ec, (np.cross(p_hat2_ec,R2_ec)))
D33 = np.dot(p_hat1_ec, (np.cross(p_hat2_ec,R3_ec)))

A1 = tau3/tau
B1 = (A1/6)*((tau**2)-(tau3**2))
A3 = (-tau1) / tau
B3 = (A3/6)*((tau**2)-(tau1**2))

A= ((A1*D21) - D22 + (A3*D23))/ (-D0)
B = ((B1*D21)+(B3*D23))/ (-D0)
E = -2*(np.dot(p_hat2_ec, R2_ec))
F = R2_ec_mag**2

a1 = -((A**2) + (A*E) + F)
b = -mu*((2*A*B) + (B*E))
c = -(mu**2)*(B**2)

roots = np.polynomial.polynomial.polyroots([c, 0, 0, b, 0, 0, a1, 0, 1])

correct_roots = []
for i in roots:
    if i>0 and np.isreal(i):
        correct_roots.append(np.real(i))  #finding the real and positive roots
print(f'roots = {correct_roots}')

choice= int(input("which root do you want to use (1, 2 or 3):"))
r2_mag= correct_roots[choice-1] # choice-1 to find the index of the choosen root 

u = mu / r2_mag**3

f1 = 1 - (0.5*u*tau1**2)
f3 = 1 - (0.5*u*tau3**2)   #truncated f and g functions
g1 = tau1 - (1/6*u*tau1**3)
g3 = tau3 -(1/6*u*tau3**3)

c1 = g3 / (f1*g3 - g1*f3)
c2 = -1
c3 = -g1 / (f1*g3 - g1*f3)
d1 = -f3 / (f1*g3 - f3*g1)
d3 = f1 / (f1*g3 - f3*g1)

p1_mag = (c1*D11 + c2*D12 + c3*D13) / (c1*D0)
p2_mag = (c1*D21 + c2*D22 + c3*D23) / (c2*D0)
p3_mag = (c1*D31 + c2*D32 + c3*D33) / (c3*D0)

p1 = p1_mag*p_hat1_ec
p2 = p2_mag*p_hat2_ec #rho vector calculation
p3 = p3_mag*p_hat3_ec

r1_vec = p1 - R1_ec
r2_vec = p2 - R2_ec  #calculated r vectors using fundamental triangle 
r3_vec = p3 - R3_ec

r2_dot_vec = d1*r1_vec + d3*r3_vec
r2_dot_mag = (r2_dot_vec[0]**2+r2_dot_vec[1]**2+r2_dot_vec[2]**2)**0.5

r_new = 0
r_old = r2_mag
counter = 0
flag = (input("function, 3rd-order series, or 4th order series (3, 4 or f):"))

while abs(r_old-r_new)> 1*10**-12:

    t1 = t1_orig -p1_mag/C
    t2 = t2_orig -p2_mag/C  #light travel time correction 
    t3 = t3_orig -p3_mag/C

    tau3 = k*(t3-t2)
    tau1 = k*(t1-t2)
    tau = tau3 - tau1

    z = (np.dot(r2_vec, r2_dot_vec)) / r2_mag**2
    q = ((np.dot(r2_dot_vec,r2_dot_vec)) / r2_mag**2) - u 
    a = ((2/r2_mag) - np.dot(r2_dot_vec,r2_dot_vec))**-1
    h= np.cross(r2_vec, r2_dot_vec) 
    h_mag= (h[0]**2+ h[1]**2+h[2]**2)**0.5
    e = (1-(h_mag**2/a))**0.5
    n = 1/ a**1.5
    u = mu/ r2_mag**3

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
        if flag == '3':
           f = 1 - (0.5* u * tau**2) + (0.5* u * z * tau**3) 
        elif flag == '4': 
           f = 1 - (0.5* u * tau**2) + (0.5* u * z * tau**3) + (1/24* ((3*u*q)-(15*u*z**2)+(u**2)) * tau**4)
        elif flag == 'f':
            if tau == tau1:
              DE1 = delta_E(r2, r2_dot, tau)
              f = 1 - (a/r2_mag)*(1- cos(DE1))
            elif tau == tau3:
              DE3 = delta_E(r2, r2_dot, tau)
              f = 1 - (a/r2_mag)*(1- cos(DE3))
        return f

    def g(tau, r2, r2_dot , flag):
        z = (np.dot(r2, r2_dot)) / r2_mag**2
        if flag == '3':
           g = tau - ((1/6)*u*tau**3)
        elif flag == '4': 
           g = tau - ((1/6)*u*tau**3) + (0.25* u*z*tau**4)
        elif flag == 'f':
            if tau == tau1:
              DE1 = delta_E(r2, r2_dot, tau)
              g = tau + (1/n)*(sin(DE1)-DE1)
            elif tau == tau3:
              DE3 = delta_E(r2, r2_dot, tau)
              g = tau + (1/n)*(sin(DE3)-DE3)
        return g

    f1 = f(tau1,r2_vec,r2_dot_vec,flag)
    f3 = f(tau3,r2_vec,r2_dot_vec,flag)
    g1 = g(tau1,r2_vec,r2_dot_vec,flag)
    g3 = g(tau3,r2_vec,r2_dot_vec,flag)
    
    c1 = g3 / (f1*g3 - g1*f3)
    c2 = -1
    c3 = -g1 / (f1*g3 - g1*f3)
    d1 = -f3 / (f1*g3 - f3*g1)
    d3 = f1 / (f1*g3 - f3*g1)

    p1_mag = (c1*D11 + c2*D12 + c3*D13) / (c1*D0)
    p2_mag = (c1*D21 + c2*D22 + c3*D23) / (c2*D0)
    p3_mag = (c1*D31 + c2*D32 + c3*D33) / (c3*D0)

    p1 = p1_mag*p_hat1_ec
    p2 = p2_mag*p_hat2_ec
    p3 = p3_mag*p_hat3_ec

    r1_vec = p1 - R1_ec
    r2_vec = p2 - R2_ec
    r3_vec = p3 - R3_ec

    r2_dot_vec = d1*r1_vec + d3*r3_vec
    r2_dot_mag = (r2_dot_vec[0]**2+r2_dot_vec[1]**2+r2_dot_vec[2]**2)**0.5
    r2_mag = (r2_vec[0]**2+r2_vec[1]**2+r2_vec[2]**2)**0.5
    
    r_old = r_new  
    r_new = r2_mag 

    counter += 1  #iteration counter

print(f'iterations: {counter}'), 
print(f'r2= {r2_vec} = {r2_mag} au')
print(f'r2dot = {r2_dot_vec*k} (au/d) = {r2_dot_mag*k} au/day')
print(f'rho2 = {p2_mag} au')

#*************baby od*******************************
a = ((2/r2_mag) - np.dot(r2_dot_vec,r2_dot_vec))**-1
print(f'a = {a} au')

h= np.cross(r2_vec, r2_dot_vec) 
h_mag= (h[0]**2+ h[1]**2+h[2]**2)**0.5
e = (1-(h_mag**2/a))**0.5
print(f'e = {e}')

i= acos((h[2])/h_mag)
print(f'i = {degrees(i)} degrees')

sin_o= h[0]
cos_o= -h[1]
Ω= atan2(sin_o,cos_o)
print(f'Ω = {degrees(Ω)} degrees')

sinf_w= r2_vec[2] / (r2_mag*sin(i))
cosf_w= (1/cos(Ω))*(r2_vec[0]/r2_mag+(cos(i)*sinf_w*sin(Ω)))
f_w= atan2(sinf_w,cosf_w)

sinf= (((np.dot(r2_vec, r2_dot_vec))/(e*r2_mag))*((a*(1-e**2))**0.5))
cosf = (1/e)*(((a*(1-e**2))/r2_mag)-1)
f= atan2(sinf,cosf)
w= f_w-f
print(f'w = {degrees(w)} degrees')

E = acos((1/e)*(1-(r2_mag/a)))
if f < 0:
    E1= E*-1
if degrees(E) < 360:
    E2 = degrees(E1)+360

M = radians(E2)- e*sin(radians(E2))
print(f'E = {E2}')
print(f'M = {degrees(M)} degrees at JD = 2459037.74868')

Jt0= (367*2023)-int(7/4*(2023+int((7+9)/12)))+int(275/9*7)+23+1721013.5
t0=(Jt0+(6/24))

n = ((k**2) / (a**3))**0.5
M0 = M + n*(t0-t2)
if degrees(M0) > 360:
   M0 = degrees(M0) - 360
print(f'M = {M0} degrees at 2023 July 23 6:00 UTC')