from astropy.io import fits
import numpy as np
from math import*

image = fits.open("obs1img1.fits")
image = np.array(image[0].data)

dark = fits.open("obs1masterdark.fits")
dark = np.array(dark[0].data)

darkpix = 0
sumdark = 0
for i in range(len(dark)):
    for j in range(len(dark)):
            sumdark += dark[i,j]
            darkpix +=1


x_1 = int(input("x value:")) -1 #- 1 because fits file starts from 1 but np array starts from 0
y_1 = int(input("y value:")) -1
r_ap = int(input("radius of the apeture:"))
r_in = int(input("radius of the inner circle:"))
r_out = int(input("radius of the outer circle:"))
gain = 0.8
d_e = ((sumdark / darkpix) / 60) *gain
read= 11
p_2 = read**2 + ((gain**2)/12)
newimage= image[y_1 - r_out: y_1 + r_out + 1,x_1 - r_out: x_1 + r_out + 1 ]
print(d_e)
x_2= int(len(newimage[0])/2)
y_2= int(len(newimage[1])/2)

pixcount_an= 0
sumADU_an= 0

for i in range(len(newimage)):
    for j in range(len(newimage)):
        r= ((i-x_2)**2+(j-y_2)**2)**0.5
        if r > r_in and r < r_out:
            sumADU_an += newimage[i,j]
            pixcount_an +=1

avg_bkgd= sumADU_an/ pixcount_an
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
        if r >= r_ap:
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

x= round((x_+ (x_1 +1 - r_ap)),3)

# print(x)


sumweight= 0
weight_col= 0
for j in range(len(image2copy)):
    sumcol = 0
    for i in range(len(image2copy)):
        sumcol += image2copy[i,j]
        sumweight += sumcol
        weight_col += sumcol*j
y_= weight_col/sumweight

y= round((y_+ (y_1 +1 - r_ap)),3)

# print(y)

sumcol3x = 0
sum_x= 0
for i in range(len(image2copy)):
    sumcol2x = 0
    for j in range(len(image2copy)):
        sumcol2x += image2copy[i,j] 
        sumcol3x += image2copy[i,j] 
    sum_x += sumcol2x*((i-x_)**2)
uncertainty_x = round(((sum_x/(sumcol3x*(sumcol3x-1)))**0.5),3)
# print(uncertainty_x)
# print(image2copy)

sumcol3y = 0
sum_y= 0
for j in range(len(image2copy)):
    sumcol2y = 0
    for i in range(len(image2copy)):
        sumcol2y += image2copy[i,j] 
        sumcol3y += image2copy[i,j] 
    sum_y += sumcol2y*((j-y_)**2)
uncertainty_y = round(((sum_y/(sumcol3y*(sumcol3y-1)))**0.5),4) #prisha helped me debug x and y uncertainty
# print(uncertainty_y)


pixcount_ap= 0
sumADU_ap= 0
for i in range(len(image2copy)):
    for j in range(len(image2copy)):
        image2copy[i,j]= newimage2[j,i]
        r= ((i-x_)**2+(j-y_)**2)**0.5
        if r <= r_ap:
            pixcount_ap += 1
            sumADU_ap += image2copy[i,j]
# print(pixcount_ap, sumADU_ap)

signal= round(sumADU_ap - pixcount_ap*avg_bkgd)
m_inst = round((-2.5*log10(signal)),9)
# print(signal,m_inst)

SNR= round(((signal**0.5)/ (1+(pixcount_ap*(1+ pixcount_ap/ pixcount_an)*((avg_bkgd+d_e+p_2)/signal)))**0.5),1)
# print(SNR)

uncertainty_m_inst= round((1.0857 / SNR),3)
uncertainty_signal = round((signal / SNR))
# print(uncertainty_m_inst, uncertainty_signal)

print((f'x= {x} ± {uncertainty_x}; y= {y} ± {uncertainty_y}; aperture contains {pixcount_ap} pixels'))
print(f'S={signal} ± {uncertainty_signal} (SNR={SNR})')
print(f'm_inst= {m_inst} ± {uncertainty_m_inst}')

#collab: theodor
