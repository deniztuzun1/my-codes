import numpy as np
from math import*
import mog4 
import matplotlib.pyplot as plt
import matplotlib as mpl

file= open("d.txt")

ra = []
dec = []
for i in file:
        ra.append(i.split(' ')[4])
        dec.append(i.split(' ')[5])

ra1, ra2, ra3 = ra[0], ra[1], ra[2]
dec1, dec2, dec3 = dec[0], dec[1], dec[2]

def convert(x):
    x_hour= float(x.split(':')[0])
    if int(x_hour)<0:
        x_min2hour= float(x.split(':')[1])/-60 
        x_sec2hour= float(x.split(':')[2])/-3600
        x_dec= x_hour + x_min2hour + x_sec2hour
    else:
        x_min2hour= float(x.split(':')[1])/60
        x_sec2hour= float(x.split(':')[2])/3600
        x_dec= x_hour + x_min2hour + x_sec2hour
    return x_dec

ra1dec= radians(convert(ra1)*15) 
ra2dec= radians(convert(ra2)*15)
ra3dec= radians(convert(ra3)*15)
dec1dec= radians(convert(dec1))
dec2dec= radians(convert(dec2))
dec3dec= radians(convert(dec3))

number = 1000000
ra1list = np.random.normal(ra1dec, 0.07*0.00000485, number)
ra2list = np.random.normal(ra2dec, 0.04*0.00000485, number )
ra3list = np.random.normal(ra3dec, 0.08*0.00000485, number )
dec1list = np.random.normal(dec1dec, 0.15*0.00000485, number ) 
dec2list = np.random.normal(dec2dec, 0.06*0.000000485, number )
dec3list = np.random.normal(dec3dec, 0.018*0.000000485, number )

counter= 0
a_values= []
e_values = []
i_values = []
o_values = []
w_values = []
M_values = []
M0_values = []

for i in range(0, len(ra1list)) :
    a,e,i,Ω,w,M,M0 = mog4.mog(ra1list[i], ra2list[i], ra3list[i],dec1list[i], dec2list[i], dec3list[i])
    if not a==None:
        a_values.append(a)
        e_values.append(e)
        i_values.append(i)
        o_values.append(Ω)
        w_values.append(w)
        M_values.append(M)
        M0_values.append(M0)
        counter += 1
    print(counter)

plt.title("Histogram of $a$", fontsize=16, fontname="Times New Roman")
plt.xlabel("Semi-major Axis (AU)", fontsize=14, fontname="Times New Roman")
plt.ylabel("Frequency", fontsize=14, fontname="Times New Roman")

n_bins = 50

n, bins, patches = plt.hist(a_values, bins=n_bins, density=True)
frequency_fractions = n/n.max()
normalized_color = mpl.colors.Normalize(frequency_fractions.min(), frequency_fractions.max())
for thisfrac, thispatch in zip(frequency_fractions, patches):
    color = plt.cm.Spectral(normalized_color(thisfrac))
    thispatch.set_facecolor(color)
        
plt.show()

plt.title("Histogram of $e$", fontsize=16, fontname="Times New Roman")
plt.xlabel("Eccentricity", fontsize=14, fontname="Times New Roman")
plt.ylabel("Frequency", fontsize=14, fontname="Times New Roman")

n_bins = 50

n, bins, patches = plt.hist(e_values, bins=n_bins, density=True)
frequency_fractions = n/n.max()
normalized_color = mpl.colors.Normalize(frequency_fractions.min(), frequency_fractions.max())
for thisfrac, thispatch in zip(frequency_fractions, patches):
    color = plt.cm.Spectral(normalized_color(thisfrac))
    thispatch.set_facecolor(color)
        
plt.show()

plt.title("Histogram of $i$", fontsize=16, fontname="Times New Roman")
plt.xlabel("Inclination (radians)", fontsize=14, fontname="Times New Roman")
plt.ylabel("Frequency", fontsize=14, fontname="Times New Roman")

n_bins = 50

n, bins, patches = plt.hist(i_values, bins=n_bins, density=True)
frequency_fractions = n/n.max()
normalized_color = mpl.colors.Normalize(frequency_fractions.min(), frequency_fractions.max())
for thisfrac, thispatch in zip(frequency_fractions, patches):
    color = plt.cm.Spectral(normalized_color(thisfrac))
    thispatch.set_facecolor(color)
        
plt.show()

plt.title("Histogram of $Ω$", fontsize=16, fontname="Times New Roman")
plt.xlabel("Longitude of the ascending node (radians)", fontsize=14, fontname="Times New Roman")
plt.ylabel("Frequency", fontsize=14, fontname="Times New Roman")

n_bins = 50

n, bins, patches = plt.hist(o_values, bins=n_bins, density=True)
frequency_fractions = n/n.max()
normalized_color = mpl.colors.Normalize(frequency_fractions.min(), frequency_fractions.max())
for thisfrac, thispatch in zip(frequency_fractions, patches):
    color = plt.cm.Spectral(normalized_color(thisfrac))
    thispatch.set_facecolor(color)
        
plt.show()

plt.title("Histogram of $w$", fontsize=16, fontname="Times New Roman")
plt.xlabel("Longitude of the perihelion (radians)", fontsize=14, fontname="Times New Roman")
plt.ylabel("Frequency", fontsize=14, fontname="Times New Roman")

n_bins = 50

n, bins, patches = plt.hist(w_values, bins=n_bins, density=True)
frequency_fractions = n/n.max()
normalized_color = mpl.colors.Normalize(frequency_fractions.min(), frequency_fractions.max())
for thisfrac, thispatch in zip(frequency_fractions, patches):
    color = plt.cm.Spectral(normalized_color(thisfrac))
    thispatch.set_facecolor(color)
        
plt.show()

plt.title("Histogram of $M$", fontsize=16, fontname="Times New Roman")
plt.xlabel("M (radians)", fontsize=14, fontname="Times New Roman")
plt.ylabel("Frequency", fontsize=14, fontname="Times New Roman")

n_bins = 50

n, bins, patches = plt.hist(M_values, bins=n_bins, density=True)
frequency_fractions = n/n.max()
normalized_color = mpl.colors.Normalize(frequency_fractions.min(), frequency_fractions.max())
for thisfrac, thispatch in zip(frequency_fractions, patches):
    color = plt.cm.Spectral(normalized_color(thisfrac))
    thispatch.set_facecolor(color)
        
plt.show()

plt.title("Histogram of $M0$", fontsize=16, fontname="Times New Roman")
plt.xlabel("M0 (radians)", fontsize=14, fontname="Times New Roman")
plt.ylabel("Frequency", fontsize=14, fontname="Times New Roman")

n_bins = 50

n, bins, patches = plt.hist(M0_values, bins=n_bins, density=True)
frequency_fractions = n/n.max()
normalized_color = mpl.colors.Normalize(frequency_fractions.min(), frequency_fractions.max())
for thisfrac, thispatch in zip(frequency_fractions, patches):
    color = plt.cm.Spectral(normalized_color(thisfrac))
    thispatch.set_facecolor(color)
        
plt.show()

mean_a =np.mean(a_values)
mean_e =np.mean(e_values)
mean_i =np.mean(i_values)
mean_o =np.mean(o_values)
mean_w =np.mean(w_values)
mean_M =np.mean(M_values)
mean_M0 =np.mean(M0_values)

print("mean", mean_a, mean_e, mean_i, mean_o, mean_w, mean_M, mean_M0)

uncertainty_a =np.std(a_values)
uncertainty_e =np.std(e_values)
uncertainty_i =np.std(i_values)
uncertainty_o =np.std(o_values)
uncertainty_w =np.std(w_values)
uncertainty_M =np.std(M_values)
uncertainty_M0 =np.std(M0_values)

print( "uncertainty",uncertainty_a, uncertainty_e,uncertainty_i,uncertainty_o,uncertainty_w,uncertainty_M,uncertainty_M0)

