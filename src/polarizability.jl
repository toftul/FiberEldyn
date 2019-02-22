# by Ivan Toftul
# 2018/07

using SpecialFunctions
include("./consts.jl")
include("./MieTheory.jl")

export 
    polarizabilityDipole, 
    polarizabilityDipoleRadCorr, 
    polarizabilityMie,
    polarizabilityAtom
    

function polarizabilityDipole(a, epsilon_p, epsilon_m)
    #=
    Calculates the polarizability of a sphere in a dipole approximation (SI units)
    
    Parameters
    ----------
        a : float
            sphere radius;
        epsilon_p, epsilon_m : complex
            epsilon of particle and media
    
    Returns
    -------
        alpha : complex
            polarizability of a sphere in a dipole approx
    =#
    return(4pi * epsilon0 * epsilon_m * a^3 * (epsilon_p - epsilon_m) / (epsilon_p + 2epsilon_m))
end
    
    
function polarizabilityDipoleRadCorr(k0, a, epsilon_p, epsilon_m)
    #=
    Calculates the exact polarizability of a sphere (SI units)
    with radiation correnctions
    
    For reference see Radiative correction in approximate treatments of electromagnetic scattering by point and body scatteres
    by Eric C. Le Ru and others
    
    Parameters
    ----------
        a : float
            sphere radius;
        epsilon_p, epsilon_m : complex
            epsilon of particle and media
    
    Returns
    -------
        alpha : complex
            polarizability of a sphere in a dipole approx 
            with radiation corrections
    =#
    alpha0 = polarizabilityDipole(a, epsilon_p, epsilon_m)
    k1 = k0 * sqrt(epsilon_m)
    return(alpha0 / (1 - 1im * alpha0 * k0^2/epsilon0 * k1/6pi))
end


function polarizabilityMie(k0, a, epsilon_p, epsilon_m)
    #=
    Calculates the exact polarizability of a sphere (SI units)
    
    Parameters
    ----------
        k0 : float
            wavelength in vacuum,
            k0 = omega/c;
        a : float
            sphere radius;
        epsilon_p, epsilon_m : complex
            epsilon of particle and media
    
    Returns
    -------
        alpha : complex
            polarizability of a sphere
    =#
    x = sqrt(epsilon_m) * k0 * a
    m = sqrt(epsilon_p / epsilon_m + 0im)
    return(4pi * epsilon0 * epsilon_m * 1.5im * Mie_an(1, m, x) * (a / x)^3)
end


# atomic polarizability
function polarizabilityAtom(d, Delta, Gamma)
    #=
    from 
    Force of light on a two-level atom near an ultrathin optical fiber
    Fam le Kien and others
    eq. (44)
    
    Parameters
    ----------
        - `d` : complex vector column of matrices elements of the dipole moment operator;
        - `Delta` : omega_field - omega_atom is the detuning frequency;
        - `Gamma` : rate of spontaneous emission
    
    Returns
    -------
        tensor polarizability of a atom
    =#
    # conj(d * d') - is the outer product d^* \otimes d
    return -conj(d * d') / (hbar * (Delta + 0.5im * Gamma))
end
