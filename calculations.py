a = 3.844*(10**8)
m1 = 7.34767309*(10**22) #kg
m2 = 5.97219*(10**24) #kg
M = m1 + m2 #kg
re = 6371*(10**3) #meters
rm = 1737.4*(10**3) #meters
Im = (2/5)*m1*(rm**2)
Ie = (2/5)*m2*(re**2)
T = 86164.091 #sec
Tdot = 1.6 /100/365.25/24/3600/1000  #converted to sec per sec (isabella helped with the conversion)
 
adot1 = ((T*(M**2)*(-1*T**-2)*Tdot*(Ie+Im))/((-2*a)*((m1*(m2**2))+(m2*(m1**2)))))
adot2 = ((T*(M**2)*(a**2)*(-1*T**-2)*Tdot)/(2*a*(M**2)))
adot = adot1 - adot2 # meters per sec
print(adot*3600*100*365.25*24) #isabella helped me with the conversion

