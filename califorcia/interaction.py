import numpy as np
from scipy.integrate import quad
from math import sqrt, exp, inf, log1p, pi
from scipy.constants import c

def k_integrand_energy(k, k0, d, epsm, rL, rR):
    """
    Integrand of the radial part of in-plane wave vector for the Casimir energy.

    Parameters
    ----------
    k : float
        in-plane wave vector
    k0 : float
        vacuum wave number
    d : float
        separation between the two plates
    epsm : float
        dielectric function of the medium evaluated at the vacuum wave number
    rL, rR : float
        reflection coefficient of the left and right plate, respectively

    Returns
    -------
    (float, float)
        result for TE and TM polarization
    """
    kappa = sqrt(epsm * k0 ** 2 + k ** 2)
    rTM_L, rTE_L = rL(k0, k)
    rTM_R, rTE_R = rR(k0, k)
    res_TE = k / 2 / pi * log1p(- rTE_L * rTE_R * exp(-2 * kappa * d))
    res_TM = k / 2 / pi * log1p(- rTM_L * rTM_R * exp(-2 * kappa * d))
    return res_TE, res_TM

def k0_func_energy(k0, d, epsm_func, rL, rR, epsrel=1.e-8, epsabs=0.0):
    """
    Casimir free energy contribution at a given wave number k0.

    Parameters
    ----------
    k0 : float
        vacuum wave number
    d : float
        separation between the two plates
    epsm_func : function
        dielectric function of the medium
    rL, rR : float
        reflection coefficient of the left and right plate, respectively

    Returns
    -------
    numpy array of two floats
        result for TE and TM polarization
    """
    f_TE = lambda t: k_integrand_energy(t/d, k0, d, epsm_func(k0 * c), rL, rR)[0]/d
    res_TE = quad(f_TE, 0, inf, epsrel=epsrel, epsabs=epsabs)[0]
    f_TM = lambda t: k_integrand_energy(t/d, k0, d, epsm_func(k0 * c), rL, rR)[1]/d
    res_TM = quad(f_TM, 0, inf, epsrel=epsrel, epsabs=epsabs)[0]
    return np.array([res_TE, res_TM])

def k_integrand_pressure(k, k0, d, epsm, rL, rR):
    """
    Integrand of the radial part of in-plane wave vector for the Casimir pressure.

    Parameters
    ----------
    k : float
        in-plane wave vector
    k0 : float
        vacuum wave number
    d : float
        separation between the two plates
    epsm : float
        dielectric function of the medium evaluated at the vacuum wave number
    rL, rR : float
        reflection coefficient of the left and right plate, respectively

    Returns
    -------
    (float, float)
        result for TE and TM polarization
    """
    kappa = sqrt(epsm * k0 ** 2 + k ** 2)
    rTM_L, rTE_L = rL(k0, k)
    rTM_R, rTE_R = rR(k0, k)
    res_TE = -2 * k * kappa / 2 / pi * rTE_L * rTE_R * exp(-2 * kappa * d) / (1 - rTE_L * rTE_R * exp(-2 * kappa * d))
    res_TM = -2 * k * kappa / 2 / pi * rTM_L * rTM_R * exp(-2 * kappa * d) / (1 - rTM_L * rTM_R * exp(-2 * kappa * d))
    return res_TE, res_TM

def k0_func_pressure(k0, d, epsm_func, rL, rR, epsrel=1.e-8, epsabs=0.0):
    """
    Casimir pressure contribution at a given wave number k0.

    Parameters
    ----------
    k0 : float
        vacuum wave number
    d : float
        separation between the two plates
    epsm_func : function
        dielectric function of the medium
    rL, rR : float
        reflection coefficient of the left and right plate, respectively

    Returns
    -------
    numpy array of two floats
        result for TE and TM polarization
    """
    f_TE = lambda t: k_integrand_pressure(t / d, k0, d, epsm_func(k0 * c), rL, rR)[0] / d
    res_TE = quad(f_TE, 0, inf, epsrel=epsrel, epsabs=epsabs)[0]
    f_TM = lambda t: k_integrand_pressure(t / d, k0, d, epsm_func(k0 * c), rL, rR)[1] / d
    res_TM = quad(f_TM, 0, inf, epsrel=epsrel, epsabs=epsabs)[0]
    return np.array([res_TE, res_TM])

def k_integrand_pressuregradient(k, k0, d, epsm, rL, rR):
    """
    Integrand of the radial part of in-plane wave vector for the Casimir pressure gradient.

    Parameters
    ----------
    k : float
        in-plane wave vector
    k0 : float
        vacuum wave number
    d : float
        separation between the two plates
    epsm : float
        dielectric function of the medium evaluated at the vacuum wave number
    rL, rR : float
        reflection coefficient of the left and right plate, respectively
    Returns
    -------
    (float, float)
        result for TE and TM polarization
    """
    kappa = sqrt(epsm * k0 ** 2 + k ** 2)
    rTM_L, rTE_L = rL(k0, k)
    rTM_R, rTE_R = rR(k0, k)
    res_TE = 4 * k * kappa ** 2 / 2 / pi * rTE_L * rTE_R * exp(-2 * kappa * d) / (1 - rTE_L * rTE_R * exp(-2 * kappa * d))**2
    res_TM = 4 * k * kappa ** 2 / 2 / pi * rTM_L * rTM_R * exp(-2 * kappa * d) / (1 - rTM_L * rTM_R * exp(-2 * kappa * d))**2
    return res_TE, res_TM

def k0_func_pressuregradient(k0, d, epsm_func, rL, rR, epsrel=1.e-8, epsabs=0.0):
    """
    Casimir pressure gradient contribution at a given wave number k0.

    Parameters
    ----------
    k0 : float
        vacuum wave number
    d : float
        separation between the two plates
    epsm_func : function
        dielectric function of the medium
    rL, rR : float
        reflection coefficient of the left and right plate, respectively

    Returns
    -------
    numpy array of two floats
        result for TE and TM polarization
    """
    f_TE = lambda t: k_integrand_pressuregradient(t / d, k0, d, epsm_func(k0 * c), rL, rR)[0] / d
    res_TE = quad(f_TE, 0, inf, epsrel=epsrel, epsabs=epsabs)[0]
    f_TM = lambda t: k_integrand_pressuregradient(t / d, k0, d, epsm_func(k0 * c), rL, rR)[1] / d
    res_TM = quad(f_TM, 0, inf, epsrel=epsrel, epsabs=epsabs)[0]
    return np.array([res_TE, res_TM])
