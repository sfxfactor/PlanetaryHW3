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
plt.clf()

######## 2 ########
rg = np.array([0.1,1.,10.,1000.])
Pin10AU = np.array([dM.Pin(Fnu10AU,i) for i in rg])
Pin130AU = np.array([dM.Pin(Fnu130AU,i) for i in rg])

print "Pin for 0.1 um grain at 10 AU: ", Pin10AU[0]
print "Pin for 1 um grain at 10 AU: ", Pin10AU[1]
print "Pin for 10 um grain at 10 AU: ", Pin10AU[2]
print "Pin for 1 mm grain at 10 AU: ", Pin10AU[3]

print "Pin for 0.1 um grain at 130 AU: ", Pin130AU[0]
print "Pin for 1 um grain at 130 AU: ", Pin130AU[1]
print "Pin for 10 um grain at 130 AU: ", Pin130AU[2]
print "Pin for 1 mm grain at 130 AU: ", Pin130AU[3]

####### 3 ########
Teq10AU = np.array([dM.Teq(Pin10AU[i],rg[i]) for i in range(4)])
Teq130AU = np.array([dM.Teq(Pin130AU[i],rg[i]) for i in range(4)])

Fnu10AUp1 = 
Fnu10AU1 = 
Fnu10AU10 = 
Fnu10AU1mm = 

plt.plot(np.log10(nu),np.log10(Fnu10AUp1(nu)),label=r"$a=0.1\mu\mathrm{m}")
plt.plot(np.log10(nu),np.log10(Fnu10AU1(nu)),label=r"$a=1\mu\mathrm{m}")
plt.plot(np.log10(nu),np.log10(Fnu10AU10(nu)),label=r"$a=10\mu\mathrm{m}")
plt.plot(np.log10(nu),np.log10(Fnu10AU1mm(nu)),label=r"$a=1\mathrm{mm}")


Fnu130AUp1 = 
Fnu130AU1 = 
Fnu130AU10 = 
Fnu130AU1mm = 

plt.plot(np.log10(nu),np.log10(Fnu130AUp1(nu)),label=r"$a=0.1\mu\mathrm{m}")
plt.plot(np.log10(nu),np.log10(Fnu130AU1(nu)),label=r"$a=1\mu\mathrm{m}")
plt.plot(np.log10(nu),np.log10(Fnu130AU10(nu)),label=r"$a=10\mu\mathrm{m}")
plt.plot(np.log10(nu),np.log10(Fnu130AU1mm(nu)),label=r"$a=1\mathrm{mm}")
