from math import sqrt, exp
from scipy.constants import c


def _combine_reflection_coefficients(r_left, r_right, kappa_layer, thickness):
    phase = exp(-2 * kappa_layer * thickness)
    rTM_left, rTE_left = r_left
    rTM_right, rTE_right = r_right
    rTM = (rTM_left + rTM_right * phase) / (1 + rTM_left * rTM_right * phase)
    rTE = (rTE_left + rTE_right * phase) / (1 + rTE_left * rTE_right * phase)
    return rTM, rTE

def def_reflection_coeff(medium, materials, thicknesses):
    """
    Define reflection coefficients of the plane by specifying medium and materials of the plane and thickness of the
    coating layers.

    Parameters
    ----------
    medium : object
    materials : list of objects
    thicknesses : list of floats

    Returns
    -------
    function(float, float)->(float, float)
        reflection coefficients taking vacuum wavenumber k0 and in-plane wave number k and returning the values of the
        reflection coefficients for TM and TE polarization, respectively

    """
    Nlayers = len(materials)
    interface_coefficients = [def_fresnel_coefficients(medium, materials[0])]
    interface_coefficients.extend(
        def_fresnel_coefficients(materials[idx], materials[idx + 1])
        for idx in range(Nlayers - 1)
    )

    if Nlayers == 1:
        return interface_coefficients[0]

    def reflection_coeff(k0, k):
        effective_reflection = interface_coefficients[-1](k0, k)

        for idx in range(Nlayers - 2, 0, -1):
            effective_reflection = _combine_reflection_coefficients(
                interface_coefficients[idx](k0, k),
                effective_reflection,
                kappa(materials[idx], k0, k),
                thicknesses[idx],
            )

        return _combine_reflection_coefficients(
            interface_coefficients[0](k0, k),
            effective_reflection,
            kappa(materials[0], k0, k),
            thicknesses[0],
        )

    return reflection_coeff

def kappa(mat, k0, k):
    if k0 == 0.:
        if mat.materialclass == "dielectric":
            return k
        elif mat.materialclass == "drude":
            return k
        elif mat.materialclass == "plasma":
            Kp = mat.wp/c
            return sqrt(Kp**2 + k**2)
    else: # k0 > 0
        return sqrt(mat.epsilon(k0*c)*k0**2 + k**2)


def def_fresnel_coefficients(mat1, mat2):
    """
    Defines Fresnel reflection coefficients for a halfspace.

    Parameters
    ----------
    mat1 : object
        material representing the medium from which the planewave is incident
    mat2 : object
        material of the halfspace

    Returns
    -------
    function(float, float)->(float, float)
        reflection coefficients taking vacuum wavenumber k0 and in-plane wave number k and returning the values of the
        reflection coefficients for TM and TE polarization, respectively
    """
    def fresnel_coefficients(k0, k):
        rTM, rTE = 0., 0.
        if mat1.materialclass == "pec":
            return -1., 1.
        if mat2.materialclass == "pec":
            return 1., -1.
        if k0 == 0.:
            if mat1.materialclass == "dielectric":
                if mat2.materialclass == "dielectric":
                    eps1 = mat1.epsilon(0.)
                    eps2 = mat2.epsilon(0.)
                    rTM = (eps2 - eps1)/(eps2 + eps1)
                    rTE = 0.
                elif mat2.materialclass == "drude":
                    rTM = 1.
                    rTE = 0.
                elif mat2.materialclass == "plasma":
                    Kp = mat2.wp/c
                    rTM = 1.
                    rTE = (k - sqrt(Kp**2 + k**2))/(k + sqrt(Kp**2 + k**2))
            elif mat1.materialclass == "drude":
                if mat2.materialclass == "dielectric":
                    rTM = -1.
                    rTE = 0.
                elif mat2.materialclass == "drude":
                    wp1 = mat1.wp
                    gamma1 = mat1.gamma
                    wp2 = mat2.wp
                    gamma2 = mat2.gamma
                    rTM = (gamma1*wp2**2 - gamma2*wp1**2)/(gamma1*wp2**2 + gamma2*wp1**2)
                    rTE = 0.
            elif mat1.materialclass == "plasma":
                if mat2.materialclass == "dielectric":
                    Kp = mat1.wp/c
                    rTM = -1.
                    rTE = -(k - sqrt(Kp**2 + k**2))/(k + sqrt(Kp**2 + k**2))
                elif mat2.materialclass == "plasma":
                    Kp1 = mat1.wp/c
                    Kp2 = mat2.wp/c
                    q1 = sqrt(Kp1**2 + k**2)
                    q2 = sqrt(Kp2**2 + k**2)
                    rTM = (Kp2**2*q1 - Kp1**2*q2)/(Kp2**2*q1 + Kp1**2*q2)
                    rTE = (q1 - q2)/(q1 + q2)
            
        else: # k0 > 0
            eps1 = mat1.epsilon(k0*c)
            eps2 = mat2.epsilon(k0*c)
            q1 = sqrt(eps1*k0**2 + k**2)
            q2 = sqrt(eps2*k0**2 + k**2)
            rTM = (eps2*q1 - eps1*q2)/(eps1*q2 + eps2*q1)
            rTE = (q1 - q2)/(q2 + q1)

        return rTM, rTE
    return fresnel_coefficients
