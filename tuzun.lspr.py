from math import *
import numpy as np

x_ast= input("asteroid's position x axis:")
y_ast= input("asteroid's position y axis:")

file= open("LSPRtestinput1.txt")
x_values = []
y_values = []
ra_values= []
dec_values = []
for i in file:
    x_values.append(i.split(' ')[0])
    y_values.append(i.split(' ')[1])
    ra_values.append(i.split(' ')[2])
    dec_values.append(i.split(' ')[3][:-1])

ra_dec_values = []
for ra in ra_values:
    ra_hour= float(ra.split(':')[0])
    ra_min2hour= float(ra.split(':')[1])/60
    ra_sec2hour= float(ra.split(':')[2])/3600  
    ra_dec= ra_hour + ra_min2hour + ra_sec2hour
    ra_dec_values.append(ra_dec)
# print(ra_dec_values)

dec_dec_values = []
for dec in dec_values:
    dec_deg= float(dec.split(':')[0])
    dec_min2deg= float(dec.split(':')[1])/60
    dec_sec2deg= float(dec.split(':')[2])/3600  
    dec_dec= dec_deg + dec_min2deg + dec_sec2deg
    dec_dec_values.append(dec_dec)
# print(dec_dec_values)
# print(dec_values)
sum_x=0
for x in x_values:
    sum_x += float(x)
# print(sum_x)

sum_y=0
for y in y_values:
    sum_y += float(y)

sum_ra=0
for ra in ra_dec_values:
    sum_ra += float(ra)
# print(sum_ra)

sum_dec=0
for dec in dec_dec_values:
    sum_dec += float(dec)
# print(sum_dec)
# print(dec_dec_values)

sum_x_x=0
for x in x_values:
    sum_x_x += float(x)*float(x)
# print(sum_x_x)

sum_y_y=0
for y in y_values:
    sum_y_y += float(y)*float(y)

sum_x_y=0
for num in range(12):
    sum_x_y += float(x_values[num])*float(y_values[num])    #theodor helped me with multiplying two list   
# print(sum_x_y)

sum_ra_x=0
for num in range(12):
    sum_ra_x += float(ra_dec_values[num])*float(x_values[num])

sum_ra_y=0
for num in range(12):
    sum_ra_y += float(ra_dec_values[num])*float(y_values[num])

sum_dec_x=0
for num in range(12):
    sum_dec_x += float(dec_dec_values[num])*float(x_values[num])

sum_dec_y=0
for num in range(12):
    sum_dec_y += float(dec_dec_values[num])*float(y_values[num])
# print(sum_dec)
# print(sum_dec_x)
# print(sum_dec_y)

matrix_A= np.array([[12, sum_x, sum_y],[sum_x, sum_x_x, sum_x_y], [sum_y, sum_x_y, sum_y_y]])
matrix_ra= np.array([sum_ra, sum_ra_x, sum_ra_y])
matrix_dec= np.array([sum_dec, sum_dec_x, sum_dec_y])
# print(matrix_dec)

matrix_A_inv=np.linalg.inv(matrix_A)

matrix_1= np.matmul(matrix_A_inv, matrix_ra)
matrix_2= np.matmul(matrix_A_inv, matrix_dec)
# print(matrix_ra)
# print(matrix_dec)

print('***************')
print('plate constants')
print('***************')
b1= matrix_1[0]*15
print(f'{b1} deg') # theodor showed me how to add deg at the end
b2= matrix_2[0]
print(f'{b2} deg')
a11= matrix_1[1]*15
print((f'{a11} deg/pix'))
a12= matrix_1[2]*15
print((f'{a12} deg/pix'))
a21= matrix_2[1]
print((f'{a21} deg/pix'))   #got different value
a22= matrix_2[2]    
print((f'{a22} deg/pix'))

print('************')
print('uncertainity')
print('************')


sum_1=0
for num in range(12):
    sum_1 += (float(ra_dec_values[num])*15-float(b1)-(float(a11)*float(x_values[num]))-(float(a12)*float(y_values[num])))**2  #multiplied re_dec_values with 15 to convert hour to deg
# print(sum_1)

uncertainty_ra= round((((sum_1/9)**0.5)*3600) , 2)
print((f'{uncertainty_ra} arcsec'))

sum_2=0
for num in range(12):
    sum_2 += (float(dec_dec_values[num])-float(b2)-(float(a21)*float(x_values[num]))-(float(a22)*float(y_values[num])))**2
# print(sum_2)

uncertainty_dec= round((((sum_2/9)**0.5)*3600) , 2)  
print((f'{uncertainty_dec} arcsec'))



print('***************************************')
print(f'astrometry for (x,y) = ({x_ast}, {y_ast})')
print('***************************************')


matrix_b= np.array([b1,b2])
matrix_a= np.array([[a11, a12],[a21, a22]])
matrix_x_y = np.array([float(x_ast), float(y_ast)])

m_axm_x_y= np.matmul(matrix_a, matrix_x_y)
# print(m_axm_x_y)

matrix_ra_dec= matrix_b + m_axm_x_y
# print(matrix_ra_dec)

ra_hours= matrix_ra_dec[0]/15
ra_HH= int(ra_hours)
# print(ra_HH) 
ra_min= int((ra_hours-ra_HH)*60)
# print(ra_min)
ra_sec= round((((ra_hours-ra_HH)*60-ra_min)*60), 3)
# print(ra_sec)
print(f'RA = {ra_HH}:{ra_min}:{ra_sec}')

dec_degs= matrix_ra_dec[1]
dec_deg= int(dec_degs)
# print(dec_deg) 
dec_min= int((dec_degs-dec_deg)*60)
# print(dec_min)
dec_sec= round((((dec_degs-dec_deg)*60-dec_min)*60), 2)
# print(dec_sec)
print(f'Dec = +{dec_deg}:{dec_min}:{dec_sec}')

