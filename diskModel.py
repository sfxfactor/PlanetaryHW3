import numpy as np
from astropy.io import ascii
import scipy.optimize as opt
import scipy.interpolate as intp

h = 6.626e-27 #erg s
c = 2.99792458e10 #cm s^-1
k = 1.380658e-16 #erg K^-1
Jy = 1e-23 # erg s^-1 cm^-2 Hz^-1
AU = 1.496e13 #cm
Rsun = 6.9599e10 #cm

def Bnu(nu,T):
    return ((2.*h*nu**3)/(c**2))*(1./(np.exp((h*nu)/(k*T))-1))


def calcF(T,Rs,Dp):
    def Fnu(nu):
        Omg = np.pi*(Rs/Dp)**2
        return Bnu(nu,T)*Omg/Jy
    return Fnu

def getQdat(rg):
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
    return Qdat

def Pin(Fnu,rg):
    Qdat=getQdat(rg)

    Pin = (np.pi*(rg*1e-4)**2)*np.trapz(Qdat['Q_abs']*Fnu(c/(Qdat['w(micron)']*1e-6))*Jy,x=c/(Qdat['w(micron)']*1e-6))

    return Pin

def Teq(Pin,rg):
    Qdat = getQdat(rg)
    #Pin = Pout = (4.*np.pi*rg**2)*Q*sigma*T**4
    def f(T):
        return (4.*np.pi*rg**2)*np.trapz(Qdat['Q_abs']*Bnu(c/(Qdat['w(micron)']*1e-6),T)*np.pi,x=c/(Qdat['w(micron)']*1e-6))-Pin

    Teq = opt.root(f,300.).x[0]
    
    return Teq

def calcFQ(T,R,D):
    Qdat=getQdat(R)
    Qabsnu = intp.UnivariateSpline(c/(Qdat['w(micron)']*1e-6),Qdat['Q_abs'])

    def Fnu(nu):
        Omg = np.pi*(R/D)**2
        return Bnu(nu,T)*Omg*Qabsnu(nu)/Jy
    return Fnu
