from astropy.io import fits
import numpy as np
from math import*

image = fits.open("aptest.fit")
image = np.array(image[0].data)

x_1 = int(input("x value:")) -1 #- 1 because fits file starts from 1 but np array starts from 0
y_1 = int(input("y value:")) -1
r_ap = int(input("radius of the apeture:"))
r_in = int(input("radius of the inner circle:"))
r_out = int(input("radius of the outer circle:"))

newimage= image[y_1 - r_out: y_1 + r_out + 1,x_1 - r_out: x_1 + r_out + 1 ]

x_2= int(len(newimage[0])/2)
y_2= int(len(newimage[1])/2)

pixcount= 0
sum= 0

for i in range(len(newimage)):
    for j in range(len(newimage)):
        r= ((i-x_2)**2+(j-y_2)**2)**0.5
        if r > r_in and r < r_out:
            sum += newimage[i,j]
            pixcount +=1

avg_bkgd= sum/ pixcount
# print(avg_bkgd)

newimage2 = newimage[y_2- r_ap: y_2 + r_ap + 1, x_2 - r_ap: x_2 + r_ap + 1 ]  
# print(newimage2)

image2copy = newimage2.copy()

x_3= int(len(image2copy[0])/2) 
y_3= int(len(image2copy[1])/2)
# print(image2copy)

for i in range(len(image2copy)):
    for j in range(len(image2copy)):
        image2copy[i,j]= newimage2[j,i]
        r= ((i-x_3)**2+(j-y_3)**2)**0.5
        if r >= 5:
            image2copy[i,j]= 0
        elif image2copy[i,j]- avg_bkgd > 0:
            image2copy[i,j]= image2copy[i,j]- avg_bkgd
        else:
            image2copy[i,j]= 0

# print(image2copy)

sumweight= 0
weight_col= 0
for i in range(len(image2copy)):
    sumcol = 0
    for j in range(len(image2copy)):
        sumcol += image2copy[i,j]
        sumweight += sumcol
        weight_col += sumcol*i

x_= weight_col/sumweight

# print(x_)

x= x_+ (x_1 +1 - r_ap)

print(x)


sumweight= 0
weight_col= 0
for j in range(len(image2copy)):
    sumcol = 0
    for i in range(len(image2copy)):
        sumcol += image2copy[i,j]
        sumweight += sumcol
        weight_col += sumcol*j
y_= weight_col/sumweight

y= y_+ (y_1 +1 - r_ap)

print(y)