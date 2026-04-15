from scipy.constants import hbar
from scipy.constants import e as eV

materialclass = "dielectric"


def f_lorentz(w, fi, wi, gi):
    '''
    w is the energy of the photon in eV
    Lorentz model for dielectric matter. 
    see: Zangwill page 635 <- !
    see: Zangwill pages 624-626 + 631
    see: Casimir physics pages 98-99
    '''
    #return fi*( wi**2 /(wi**2 - w**2 - 1j*gi*w ) )
    return fi*( wi**2 /(wi**2 + w**2 + gi*w ) )


def debye(w, epsD, epsInf, tau):
    '''
    only reference I found (+refs. inside...)
    ver: https://en.wikipedia.org/wiki/Dielectric#Debye_relaxation
    '''
    #return (epsD-epsInf) /(1 - 1j*w*tau)
    return (epsD-epsInf) /(1 + w*tau)

class eps_f_lorentz_n1_debye:
    def __init__(self, f1, w1, g1, epsD, epsInf, tau):
        self.f1 = f1
        self.w1 = w1
        self.g1 = g1
        self.epsD = epsD
        self.epsInf = epsInf
        self.tau = tau

    def __call__(self, w):
        return ( 1 + f_lorentz(w, self.f1, self.w1, self.g1) 
                   + debye(w, self.epsD, self.epsInf, self.tau) )



epsilon = eps_f_lorentz_n1_debye(0.77, 2.79545e+16, 2.05101e+16, 78.3, 1.8430, 7.76690e-12)
#print(epsilon(0.))
