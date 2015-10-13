import numpy as np
from astropy.io import ascii

h = 6.626e-27 #erg s
c = 2.99792458e10 #cm s^-1
k = 1.380658e-16 #erg K^-1
Jy = 1e-23 # erg s^-1 cm^-2 Hz^-1
AU = 1.496e13 #cm
Rsun = 6.9599e10 #cm

def calcF(T,Rs,Dp):
    def Fnu(nu):
        Inu = ((2.*h*nu**3)/(c**2))*(1./(np.exp((h*nu)/(k*T))-1.))
        Omg = np.pi*(Rs/Dp)**2
        return Inu*Omg/Jy
    return Fnu

def Pin(Fnu,rg):
    if rg == 0.1:
        Qdat=ascii.read("suvSil_21",header_start=2436,data_start=2437,data_end=2677)
    elif rg == 1.:
        Qdat=ascii.read("suvSil_21",header_start=3651,data_start=3656,data_end=3892)
    elif rg == 10.:
        Qdat=ascii.read("suvSil_21",header_start=4866,data_start=4867,data_end=5107)
    elif rg == 1000.:
        Qdat=ascii.read("suvSil_21",header_start=4866,data_start=4867,data_end=5107)
        Qdat['Q_abs']=1.
    else:
        raise Error("Grain radius must be 0.1, 1, 10, or 1000 um")

    Pin = (np.pi*(rg*1e-4)**2)*np.trapz(Qdat['Q_abs']*Fnu(c/Qdat['w(micron)'])*Jy,x=c/Qdat['w(micron)'])

    return Pin

def Teq(Pin,rg):
    #Pin = Pout = (4.*np.pi*rg**2)*Q*sigma*T**4
    return 0.

