import numpy as np
import matplotlib.pyplot as plt
import diskModel as dM

h = 6.626e-27 #erg s
c = 2.99792458e10 #cm s^-1
k = 1.380658e-16 #erg K^-1
sigsb = 5.67051e-5 #erg cm^-2 K^-4 s^-1
G = 6.67259e-8 #cm^3 gram^-1 s^-2
Jy = 1e-23 # erg s^-1 cm^-2 Hz^-1
AU = 1.496e13 #cm
Rsun = 6.9599e10 #cm
Msun = 1.989e33 #g

######## 1 ######
Tfom = 8590
Rfom = 1.842*Rsun
Fnu10AU=dM.calcF(Tfom,Rfom,10.*AU)
l = 10**np.linspace(-5,1,1000)
nu=c/l

plt.plot(np.log10(nu),np.log10(Fnu10AU(nu)),label="10 AU")
plt.xlabel(r"$\log\nu$ [Hz]")
plt.ylabel(r"$\log F_\nu$ [Jy]")

Fnu130AU=dM.calcF(Tfom,Rfom,130.*AU)

plt.plot(np.log10(nu),np.log10(Fnu130AU(nu)),label="130 AU")
plt.legend(loc=2)
plt.savefig("FomFnu.pdf")
plt.clf()

######## 2 ########
rg = np.array([0.1,1.,10.,1000.])
Pin10AU = np.array([dM.Pin(Fnu10AU,i) for i in rg])
PinBB10AU = np.array([(np.pi*(Rfom**2)*sigsb*(Tfom**4)*((i*1e-4)/(10.*AU))**2) for i in rg])
Pin130AU = np.array([dM.Pin(Fnu130AU,i) for i in rg])
PinBB130AU = np.array([(np.pi*(Rfom**2)*sigsb*(Tfom**4)*((i*1e-4)/(130.*AU))**2) for i in rg])

print "Pin for 0.1 um grain at 10 AU: ", Pin10AU[0]
print "Pin BB for 0.1 um grain at 10 AU: ", PinBB10AU[0]
print "Pin for 1 um grain at 10 AU: ", Pin10AU[1]
print "Pin BB for 1 um grain at 10 AU: ", PinBB10AU[1]
print "Pin for 10 um grain at 10 AU: ", Pin10AU[2]
print "Pin BB for 10 um grain at 10 AU: ", PinBB10AU[2]
print "Pin for 1 mm grain at 10 AU: ", Pin10AU[3]
print "Pin BB for 1 mm grain at 10 AU: ", PinBB10AU[3]

print "Pin for 0.1 um grain at 130 AU: ", Pin130AU[0]
print "Pin BB for 0.1 um grain at 130 AU: ", PinBB130AU[0]
print "Pin for 1 um grain at 130 AU: ", Pin130AU[1]
print "Pin BB for 1 um grain at 130 AU: ", PinBB130AU[1]
print "Pin for 10 um grain at 130 AU: ", Pin130AU[2]
print "Pin BB for 10 um grain at 130 AU: ", PinBB130AU[2]
print "Pin for 1 mm grain at 130 AU: ", Pin130AU[3]
print "Pin BB for 1 mm grain at 130 AU: ", PinBB130AU[3]

####### 3 ########
Teq10AU = np.array([dM.Teq(Pin10AU[i],rg[i]) for i in range(4)])
TeqBB10AU = (280*(Rfom/Rsun)**0.5*(10**-0.5)*(Tfom/5800))
print Teq10AU
print TeqBB10AU
Teq130AU = np.array([dM.Teq(Pin130AU[i],rg[i]) for i in range(4)])
TeqBB130AU = (280*(Rfom/Rsun)**0.5*(130**-0.5)*(Tfom/5800))
print Teq130AU
print TeqBB130AU

D=7.7*3.0857e18

Fnu10AUp1 = dM.calcFQ(Teq10AU[0],rg[0],D)
Fnu10AU1 = dM.calcFQ(Teq10AU[1],rg[1],D)
Fnu10AU10 = dM.calcFQ(Teq10AU[2],rg[2],D)
Fnu10AU1mm = dM.calcFQ(Teq10AU[3],rg[3],D)

plt.plot(np.log10(nu),np.log10(Fnu10AUp1(nu)),label=r"$a=0.1\mu\mathrm{m}$")
plt.plot(np.log10(nu),np.log10(Fnu10AU1(nu)),label=r"$a=1\mu\mathrm{m}$")
plt.plot(np.log10(nu),np.log10(Fnu10AU10(nu)),label=r"$a=10\mu\mathrm{m}$")
plt.plot(np.log10(nu),np.log10(Fnu10AU1mm(nu)),label=r"$a=1\mathrm{mm}$")
plt.xlabel(r"$\log\nu$ [Hz]")
plt.ylabel(r"$\log F_\nu$ [Jy]")
plt.ylim(-50,-15)
plt.xlim(9,15)
plt.legend(loc=8)
plt.savefig('Fnu10Au.pdf')
plt.clf()


Fnu130AUp1 = dM.calcFQ(Teq130AU[0],rg[0],D)
Fnu130AU1 = dM.calcFQ(Teq130AU[1],rg[1],D)
Fnu130AU10 = dM.calcFQ(Teq130AU[2],rg[2],D)
Fnu130AU1mm = dM.calcFQ(Teq130AU[3],rg[3],D)

plt.plot(np.log10(nu),np.log10(Fnu130AUp1(nu)),label=r"$a=0.1\mu\mathrm{m}$")
plt.plot(np.log10(nu),np.log10(Fnu130AU1(nu)),label=r"$a=1\mu\mathrm{m}$")
plt.plot(np.log10(nu),np.log10(Fnu130AU10(nu)),label=r"$a=10\mu\mathrm{m}$")
plt.plot(np.log10(nu),np.log10(Fnu130AU1mm(nu)),label=r"$a=1\mathrm{mm}$")
plt.xlabel(r"$\log\nu$ [Hz]")
plt.ylabel(r"$\log F_\nu$ [Jy]")
plt.ylim(-50,-15)
plt.xlim(9,15)
plt.legend(loc=8)
plt.savefig('Fnu130Au.pdf')


########## Part 4 #######
FpeakFom = np.array([6e2,1e4])*(1e-3)
Fpeak10AU = np.array([np.max(Fnu10AUp1(nu)),np.max(Fnu10AU1(nu)),np.max(Fnu10AU10(nu)),np.max(Fnu10AU1mm(nu))])
Fpeak130AU = np.array([np.max(Fnu130AUp1(nu)),np.max(Fnu130AU1(nu)),np.max(Fnu130AU10(nu)),np.max(Fnu130AU1mm(nu))])

N10AU = FpeakFom[0]/Fpeak10AU
N130AU = FpeakFom[1]/Fpeak130AU

print N10AU, " grains are needed at 10AU for 0.1, 1, 10, 1000 um grains."
print N130AU, " grains are needed at 130AU for 0.1, 1, 10, 1000 um grains."

M10AU = N10AU * 2. * (4./3.)*np.pi*(rg*1e-4)**3
M130AU = N130AU * 2. * (4./3.)*np.pi*(rg*1e-4)**3

print M10AU, " grams of dust at 10AU for 0.1, 1, 10, 1000 um grains."
print M130AU, " grams of dust at 130AU for 0.1, 1, 10, 1000 um grains."


######### Part 5 ##########
Mfom = 1.92*Msun
Frad10AU = Pin10AU/c
Frad130AU = Pin130AU/c

B10AU = Frad10AU*((10.*AU)**2)/(G*Mfom*2.*(4./3.)*np.pi*(rg*1e-4)**3)
B130AU = Frad130AU*((130.*AU)**2)/(G*Mfom*2.*(4./3.)*np.pi*(rg*1e-4)**3)
print B10AU
print B130AU

Fpr10AU = Pin10AU*np.sqrt(G*Mfom/(10*AU))/c**2
Fpr130AU = Pin130AU*np.sqrt(G*Mfom/(130*AU))/c**2

T10AU = (4e2/(Mfom/Msun))*(10.**2)/B10AU
T130AU = (4e2/(Mfom/Msun))*(130.**2)/B130AU
print T10AU
print T130AU

#Bpr10AU = Fpr10AU*((10.*AU)**2)/(G*Mfom*2.*(4./3.)*np.pi*(rg*1e-4)**3)
#Bpr130AU = Fpr130AU*((130.*AU)**2)/(G*Mfom*2.*(4./3.)*np.pi*(rg*1e-4)**3)
#print Bpr10AU
#print Bpr130AU
