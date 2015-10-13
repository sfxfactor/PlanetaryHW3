import numpy as np
import matplotlib.pyplot as plt
import diskModel as dM

h = 6.626e-27 #erg s
c = 2.99792458e10 #cm s^-1
k = 1.380658e-16 #erg K^-1
Jy = 1e-23 # erg s^-1 cm^-2 Hz^-1
AU = 1.496e13 #cm
Rsun = 6.9599e10 #cm

######## 1 ######
Fnu10AU=dM.calcF(8590,1.842*Rsun,10.*AU)
l = 10**(np.linspace(-1,3,100)-4.)
nu=c/l

plt.plot(np.log10(nu),np.log10(Fnu10AU(nu)),label="10 AU")
plt.xlabel(r"$\log\nu$ [Hz]")
plt.ylabel(r"$\log F_\nu$ [Jy]")

Fnu130AU=dM.calcF(8590,1.842*Rsun,130.*AU)
l = 10**(np.linspace(-1,3,100)-4.)
nu=c/l

plt.plot(np.log10(nu),np.log10(Fnu130AU(nu)),label="130 AU")
plt.legend(loc=2)
plt.savefig("FomFnu.pdf")

######## 2 ########
print "Pin for 0.1 um grain at 10 AU: ", dM.Pin(Fnu10AU,0.1)
print "Pin for 1 um grain at 10 AU: ", dM.Pin(Fnu10AU,1.)
print "Pin for 10 um grain at 10 AU: ", dM.Pin(Fnu10AU,10.)
print "Pin for 1 mm grain at 10 AU: ", dM.Pin(Fnu10AU,1000.)

print "Pin for 0.1 um grain at 130 AU: ", dM.Pin(Fnu130AU,0.1)
print "Pin for 1 um grain at 130 AU: ", dM.Pin(Fnu130AU,1.)
print "Pin for 10 um grain at 130 AU: ", dM.Pin(Fnu130AU,10.)
print "Pin for 1 mm grain at 130 AU: ", dM.Pin(Fnu130AU,1000.)
