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
    '''
    Returns the specific intensity of a black body given by the Planck function.

    :param nu:
    the frequency in Hz at which the specific intensity is calculated

    :param T:
    the temperature in K at which the specific intensity is calculated
    '''
    return ((2.*h*nu**3)/(c**2))*(1./(np.exp((h*nu)/(k*T))-1))


def calcF(T,Rs,Dp):
    '''
    Returns F_nu, the flux density of a black-body in Jy of a given temperature and size and at a certain distance. F_nu is a function of nu.

    :param T:
    the temperature in K at which the flux density is calculated

    :param Rs:
    the radius in cm of the black-body which is emitting

    :param Dp:
    the distance in cm at which to calculate the flux density.
    '''
    def Fnu(nu):
        Omg = np.pi*(Rs/Dp)**2
        return Bnu(nu,T)*Omg/Jy
    return Fnu

#test for sun, should be ~1.24e6 Jy
#print calcF(5800,7e10,1.49598e13)(10e9)

def getQdat(rg):
    '''
    Reads in the Q data of a given size grain from the "suvSil_21" file.

    :param rg:
    the radius in um of the grain
    '''
    if rg == 0.1:
        Qdat=ascii.read("suvSil_21",header_start=2436,data_start=2437,data_end=2678)
    elif rg == 1.:
        Qdat=ascii.read("suvSil_21",header_start=3651,data_start=3652,data_end=3893)
    elif rg == 10.:
        Qdat=ascii.read("suvSil_21",header_start=4866,data_start=4867,data_end=5108)
    elif rg == 1000.:
        Qdat=ascii.read("suvSil_21",header_start=4866,data_start=4867,data_end=5108)
        Qdat['Q_abs']=1.
    else:
        raise Error("Grain radius must be 0.1, 1, 10, or 1000 um")
    return Qdat

def exQabsnu(Qdat,Rg):
    '''
    Returns a Qabs as a function of nu which is interpolated from Qdat and extrapolated to longer wavelengths using a 1/l^2 power law.

    :param Qdat:
    The Q data to extrapolate

    :param Rg:
    the radius in um of the grain, this is used to ignore 1mm perfect absorbers
    '''
    def Qabsnu(nu):
        if Rg == 1000.:
            exQabsnu = np.ones(np.size(nu))
        else: 
            exQabs = intp.interp1d(c/(Qdat['w(micron)']*1e-4),Qdat['Q_abs'],bounds_error=0,fill_value=-1.)
            exQabsnu = exQabs(nu)
            # check to make sure only use power law for small nu
            for i in np.where((exQabsnu==-1) & (nu<(c/0.01)))[0]:
                exQabsnu[i] = exQabs(c/1e-1) * ((nu[i]*1e-1)/c)**2
            if np.size(np.where(exQabsnu==-1)[0]>0):
                print "no Q data available for wavelengths shorter than 1 nm"
        return exQabsnu
    return Qabsnu

def Pin(Fnu,rg):
    '''
    Returns the power absorbed by a grain

    :param Fnu:
    the flux density in Jy of the radiation field the grain is in

    :param rg:
    radius in um of the grain
    '''
    Qdat=getQdat(rg)
    Qabs = exQabsnu(Qdat,rg)
    l = np.append(np.logspace(4,3,20),Qdat['w(micron)'])*1e-4
    nu = c/l

    Pin = (np.pi*(rg*1e-4)**2)*np.trapz(Qabs(nu)*Fnu(nu)*Jy,x=nu)

    return Pin

def Teq(Pin,rg):
    '''
    Returns the equilibrium temperature of a grain by balancing the power absorbed and emitted.

    :param Pin:
    the power absorbed by the grain

    :param rg:
    the radius in um of the grain
    '''
    Qdat = getQdat(rg)
    Qabs = exQabsnu(Qdat,rg)
    l = np.append(np.logspace(4,3,20),Qdat['w(micron)'])*1e-4
    nu = c/l
    #Pin = Pout = (4.*np.pi*rg**2)*Q*sigma*T**4
    #define a function to find the root of
    def f(T):
        return Pin-(4.*np.pi*(rg*1e-4)**2)*np.trapz(Qabs(nu)*Bnu(nu,T)*np.pi,x=nu)

    Teq = opt.root(f,100.).x[0]
    
    return Teq

def calcFQ(T,rg,D):
    '''
    Returns Fnu, the spectrum of a grain, in Jy as a function of frequency

    :param T:
    the temperature in K of the grain

    :param rg:
    the radius in um of the grain

    :param D:
    the distance to the observer
    '''
    Qdat=getQdat(rg)
    #spline does not do a good job of interpolation- use power law for low frequencies
    Qabsnu=exQabsnu(Qdat,rg)
    def Fnu(nu):
        Omg = np.pi*(rg/D)**2
        return Bnu(nu,T)*Omg*Qabsnu(nu)/Jy
    return Fnu
