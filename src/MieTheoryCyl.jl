# by Ivan Toftul
# 2018/10

using SpecialFunctions
using LinearAlgebra

export 
    MieCyl_an,
    MieCyl_bn,
    MieCylEs

function MieCyl_an(n, m, x)
    #= 
        scattering coefficents in cylindrical symmetry
        from Bohren eq. (8.32)

        Parameters
        ----------
            n - order
            m - relative refractive index
            x - size parameter (=ka)

        Returns
        -------
            an - complex number, dimentionless
    =# 
    mx = m*x

    Jn_x   = besselj(n, x)
    Jn_mx  = besselj(n, mx)
    Jnp_x  = 0.5( besselj(n-1,  x) -  besselj(n+1,  x))
    Jnp_mx = 0.5( besselj(n-1, mx) -  besselj(n+1, mx))
    H1np   = 0.5(hankelh1(n-1,  x) - hankelh1(n+1,  x))
    H1n    = hankelh1(n, x)

    return (m * Jnp_x * Jn_mx - Jn_x * Jnp_mx) / (m * Jn_mx * H1np - Jnp_mx * H1n)
end


function MieCyl_bn(n, m, x)
    #= 
        scattering coefficents in cylindrical symmetry
        from Bohren eq. (8.30)

        Parameters
        ----------
            n - order
            m - relative refractive index
            x - size parameter (=ka)

        Returns
        -------
            bn - complex number, dimentionless
    =# 
    mx = m*x

    Jn_x   = besselj(n, x)
    Jn_mx  = besselj(n, mx)
    Jnp_x  = 0.5( besselj(n-1,  x) -  besselj(n+1,  x))
    Jnp_mx = 0.5( besselj(n-1, mx) -  besselj(n+1, mx))
    H1np   = 0.5(hankelh1(n-1,  x) - hankelh1(n+1,  x))
    H1n    = hankelh1(n, x)

    return (Jn_mx * Jnp_x - m * Jnp_mx * Jn_x) / (Jn_mx * H1np - m * Jnp_mx * H1n)
end


function MieCylEs(r, phi, k0, Rf, epsf, epsm, E0, case; nmax=5000, rtol=1e-6)    
    #=
    Scattered el-field of a lin.pol. plane wave scattered on an infinite cylinder 
       in cylindrical coodinates. Plane wave is coming towards -x direction.
    
    See theory from Craig F Bohren & Donald R Huffman
    Absorbtion and scattering of light by Small Particles
    Ch 8 par 4
       
    Case I: a Transverse Magnetic (TM) mode. 
        The magnetic field of the incident wave is 
        perpendicular to the cylinder axis.
    Case II: a Transverse Electric (TE) mode.
        The electric field is perpendicular 
        to the cylinder axis.

    Geometry (Fig. 8.2 and Fig. 8.3 in Bohren)
        * field incidence direction: -x 
        * cylinder axis: z
        * polar angle: phi -- angle with x-axis
    
    Only for normal incident light (zeta = pi/2) 
    which simplifies to 
    h = - k cos(zeta) = 0, l = k
    a_nI = 0, b_nII = 0

    Parameters
    ----------
        r, phi, z : float
            reciever coordinates;
        k0 : float
            k0-vector value, 2pi/\lambda_0 = \omega/c;
        R : float
            cylynder radius;
        epsf : float 
            fiber epsilon
        epsm : float 
            medium epsion
        E0 : float
            the el-field amplitude of incident wave;
        nmax : int
            max indexe in the sum;
        case : string
            TM (case I) or TE (case II)
    
    Returns
    -------
        Es : vector of scattered field
             in cylindrical coordinates
    
    =#
    # h = 0
    m = sqrt(epsf/epsm)
    k = k0 * sqrt(epsm)
    
    rho = r * k
    x = k * Rf
    mx = m * x
    
    # n = 0
    Zn = hankelh1(0, rho)
    Jn_x = besselj(0, x)
    Jn_mx = besselj(0, mx)
    Jnp_x = .5(besselj(-1, x) - besselj(1, x))
    Jnp_mx = .5(besselj(-1, mx) - besselj(1, mx))
    H1np = .5(hankelh1(-1, x) - hankelh1(1, x))
    H1n = hankelh1(0, x)

    EsPLUS = 0
    
    if case == "TM"
        # Case I
        b_nI = (Jn_mx * Jnp_x - m * Jnp_mx * Jn_x) / (Jn_mx * H1np - m * Jnp_mx * H1n)
        EsPLUS = b_nI * Zn * [0, 0, 1]
    elseif case == "TE"
        # Case II
        Znp = .5(hankelh1(-1, rho) - hankelh1(1, rho))
        a_nII = (m * Jnp_x * Jn_mx - Jn_x * Jnp_mx) / (m * Jn_mx * H1np - Jnp_mx * H1n)
        EsPLUS = a_nII * [0, -Znp, 0]
    end

    Es = EsPLUS
    
    timesBeforeBreak = 2
    j = 0
    # terms with |n| > 0
    for n in range(1, stop=nmax)
        # n > 0
        Zn = hankelh1(n, rho)
        Jn_x = besselj(n, x)
        Jn_mx = besselj(n, mx)
        Jnp_x = .5(besselj(n-1, x) - besselj(n+1, x))
        Jnp_mx = .5(besselj(n-1, mx) - besselj(n+1, mx))
        H1np = .5(hankelh1(n-1, x) - hankelh1(n+1, x))
        H1n = hankelh1(n, x)

        # n < 0
        Zn_ = Zn * (-1)^n

        if case == "TM"
            # Case I
            b_nI = (Jn_mx * Jnp_x - m * Jnp_mx * Jn_x) / (Jn_mx * H1np - m * Jnp_mx * H1n)
            EsPLUS = (-1.0im)^n * exp(im * n * phi) * b_nI * Zn * [0, 0, 1]

            b_nI_ = b_nI
            EsPLUS_ = (-1.0im)^(-n) * exp(-im * n * phi) * b_nI_ * Zn_ * [0, 0, 1]
            
            EsPLUS += EsPLUS_
        elseif case == "TE"
            # Case II
            Znp = .5(hankelh1(n-1, rho) - hankelh1(n+1, rho))
            a_nII = (m * Jnp_x * Jn_mx - Jn_x * Jnp_mx) / (m * Jn_mx * H1np - Jnp_mx * H1n)
            # factor i E0 added at the end
            EsPLUS = (-1.0im)^n * exp(im * n * phi) * a_nII * [im * n * Zn / rho, -Znp, 0]

            Znp_ = Znp * (-1)^n
            a_nII_ = a_nII
            # factor i E0 added at the end
            EsPLUS_ = (-1.0im)^(-n) * exp(-im * n * phi) * a_nII_ * [-im * n * Zn_ / rho, -Znp_, 0]
            
            EsPLUS += EsPLUS_
        end
        if norm(EsPLUS)/norm(Es) < rtol
            j += 1
        end
        if j == timesBeforeBreak
            break
        end
        if n == 4999
            println("Desired accuracy can not be reached!")
        end
        Es += EsPLUS
    end

    if case == "TM"
        return(-E0 * Es)
    elseif case == "TE"
        # Bohren & Huffman book has a typo in Es for Case II (= TE)
        # one (-1) is missing
        return(-im*E0 * Es)
    else
        error("Invalid case!")
        return [0, 0, 0]
    end
end
    