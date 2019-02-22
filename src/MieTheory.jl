# by Ivan Toftul
# 2018 / 08

using SpecialFunctions
using HCubature  # works faster than Cubature
using StaticArrays  # for faster integration

include("./consts.jl")
include("./supportFunctions.jl")


export 
    Mie_an,
    Mie_bn,
    MieCrossSectionSCn,
    MieCrossSectionEXTn,
    MieCrossSectionABSn,
    MieCrossSectionSC,
    MieCrossSectionEXT,
    MieCrossSectionABS,
    generalMieEHrFromPol,
    generalMieAlm,
    generalMieBlm,
    generalMieBuildAlm,
    generalMieBuildBlm,
    generalMieFxFy,
    generalMieFz,
    generalMieTx,
    generalMieTy,
    generalMieTz,
    vectorHarmonicsMnm,
    vectorHarmonicsNnm


# electric Mie coefficent
function Mie_an(n, m, x, mu=1)
    #=
    For detatails see Bohren p. 100
    
    Arguments:
        - `m = N1/N = n_p / n_m` : relative refractive index;
        - `x = n_m k0 a` : size parameter;
        - `n` : 2^n multipole order;
        - `mu = mu_p / mu_m` : relative magnetic permittivity
    =#
    mx = m * x
    jnmx = sphericaljn(n, mx)
    jnx = sphericaljn(n, x)
    h1nx = sphericalh1(n, x)
    xjnx_p = jnx + x * sphericaljnp(n, x)
    mxjnmx_p = jnmx + mx * sphericaljnp(n, mx)
    xh1nx_p = h1nx + x * sphericalh1p(n, x)
    
    return (m^2 * jnmx * xjnx_p - mu * jnx * mxjnmx_p) / (m^2 * jnmx * xh1nx_p - mu * h1nx * mxjnmx_p)
end

# magnetic Mie coefficent
function Mie_bn(n, m, x, mu=1)
    #=
    For detatails see Bohren p. 100
    
    Arguments:
        - `m = N1 / N = n_p / n_m` : relative refractive index;
        - `x = n_m k0 a` : size parameter;
        - `n` : 2^n multipole order;
        - `mu = mu1 / mu =  mu_p / mu_m` : relative magnetic permittivity
    =#
    mx = m * x
    jnmx = sphericaljn(n, mx)
    jnx = sphericaljn(n, x)
    h1nx = sphericalh1(n, x)
    xjnx_p = jnx + x * sphericaljnp(n, x)
    mxjnmx_p = jnmx + mx * sphericaljnp(n, mx)
    xh1nx_p = h1nx + x * sphericalh1p(n, x)

    return (mu * jnmx * xjnx_p - jnx * mxjnmx_p) / (mu * jnmx * xh1nx_p - h1nx * mxjnmx_p)
end


function MieCrossSectionSCn(n, x, k0, epsm, epsp; mu_m=1.0, mu_p=1.0, part="both")
    #=
        single term in the sum of scattering cross section 
        based on Bohren & Huffman book (1984)
        eq. (4.61)

        part = "both", "electric", "magnetic"
    =#
    if part == "electric"
        k = sqrt(epsm * mu_m) * k0
        m = sqrt(epsp/epsm)
        mu = sqrt(mu_p / mu_m)
        an = Mie_an(n, m, x, mu)
        bn = 0
        return 2pi/k^2 * (2n+1) * (abs2(an) + abs2(bn))
    elseif part == "magnetic"
        k = sqrt(epsm * mu_m) * k0
        m = sqrt(epsp/epsm)
        mu = sqrt(mu_p / mu_m)
        an = 0
        bn = Mie_bn(n, m, x, mu)
        return 2pi/k^2 * (2n+1) * (abs2(an) + abs2(bn))
    elseif part == "both"
        k = sqrt(epsm * mu_m) * k0
        m = sqrt(epsp/epsm)
        mu = sqrt(mu_p / mu_m)
        an = Mie_an(n, m, x, mu)
        bn = Mie_bn(n, m, x, mu)
        return 2pi/k^2 * (2n+1) * (abs2(an) + abs2(bn))
    else
        error("invalid part")
    end


end

function MieCrossSectionEXTn(n, x, k0, epsm, epsp; mu_m=1.0, mu_p=1.0, part="both")
    #=
        single term in the sum of scattering cross section 
        based on Bohren & Huffman book (1984)
        eq. (4.62)

        part = "both", "electric", "magnetic"
    =#
    if part == "electric"
        k = sqrt(epsm * mu_m) * k0
        m = sqrt(epsp/epsm)
        mu = sqrt(mu_p / mu_m)
        an = Mie_an(n, m, x, mu)
        bn = 0
        return 2pi/k^2 * (2n+1) * real(an + bn)
    elseif part == "magnetic"
        k = sqrt(epsm * mu_m) * k0
        m = sqrt(epsp/epsm)
        mu = sqrt(mu_p / mu_m)
        an = 0
        bn = Mie_bn(n, m, x, mu)
        return 2pi/k^2 * (2n+1) * real(an + bn)
    elseif part == "both"
        k = sqrt(epsm * mu_m) * k0
        m = sqrt(epsp/epsm)
        mu = sqrt(mu_p / mu_m)
        an = Mie_an(n, m, x, mu)
        bn = Mie_bn(n, m, x, mu)
        return 2pi/k^2 * (2n+1) * real(an + bn)
    else
        error("invalid part")
    end
end

function MieCrossSectionABSn(n, x, k0, epsm, epsp; mu_m=1.0, mu_p=1.0, part="both")
    #=
        single term in the sum of scattering cross section 
        based on Bohren & Huffman book (1984)
        \sigma^{\text{ext}} = \sigma^{\text{abs}} + \sigma^{\text{sc}}
    =#
    return MieCrossSectionEXTn(n, x, k0, epsm, epsp; mu_m=mu_m, mu_p=mu_p, part=part) - MieCrossSectionSCn(n, x, k0, epsm, epsp; mu_m=mu_m, mu_p=mu_p, part=part)
end

function MieCrossSectionSC(x, k0, epsm, epsp; mu_m=1.0, mu_p=1.0, rtol=1e-6, part="both")
    #=
        single term in the sum of scattering cross section 
        based on Bohren & Huffman book (1984)
        eq. (4.61)
    =#
    sum = 0.0
    nsum = 0.0
    err = 1.0
    n = 1
    while ((err > rtol) && (n < 10000))
        nsum = MieCrossSectionSCn(n, x, k0, epsm, epsp; mu_m=mu_m, mu_p=mu_p)
        sum += nsum

        err = abs(nsum / sum)
        n += 1
    end
    return sum
end

function MieCrossSectionEXT(x, k0, epsm, epsp; mu_m=1.0, mu_p=1.0, rtol=1e-6, part="both")
    #=
        single term in the sum of scattering cross section 
        based on Bohren & Huffman book (1984)
        eq. (4.62)
    =#
    sum = 0.0
    nsum = 0.0
    err = 1.0
    n = 1
    while ((err > rtol) && (n < 10000))
        nsum = MieCrossSectionEXTn(n, x, k0, epsm, epsp; mu_m=mu_m, mu_p=mu_p, part=part)
        sum += nsum

        err = abs(nsum / sum)
        n += 1
    end
    return sum
end


function MieCrossSectionABS(x, k0, epsm, epsp; mu_m=1.0, mu_p=1.0, rtol=1e-6, part="both")
    #=
        single term in the sum of scattering cross section 
        based on Bohren & Huffman book (1984)
        \sigma^{\text{ext}} = \sigma^{\text{abs}} + \sigma^{\text{sc}}
    =#
    sum = 0.0
    nsum = 0.0
    err = 1.0
    n = 1
    while ((err > rtol) && (n < 10000))
        nsum = MieCrossSectionABSn(n, x, k0, epsm, epsp; mu_m=mu_m, mu_p=mu_p, part=part)
        sum += nsum

        err = abs(nsum / sum)
        n += 1
    end
    return sum
end


function generalMieEHrFromPol(r, theta, varphi, rG0, Apol)
    #=
        projection of the field in polar to spherical with the origin rG0

        based on Generalized Lorenz-Mie theories by Gerard Gouesbet
        Fig. 3.1

        Arguments
            - r, theta, varphi -- spherical coordinates for the surface integration
            - rG0 = [u0, v0, w0] -- position of the particle 
                                    in Cartesian coordinates
            - Apol(rho, phi, z) = [Arho, Aphi, Az] -- el./mag. field funciton
    =#    
    rG = [rG0[1] + r * sin(theta) * cos(varphi),
          rG0[2] + r * sin(theta) * sin(varphi),
          rG0[3] + r * cos(theta)]
    
    rho = sqrt(rG[1]^2 + rG[2]^2)
    phi = atan(rG[2], rG[1])
    z = rG[3]

    A = Apol(rho, phi, z)
    
    return(A[1] * sin(theta) * cos(varphi - phi) + 
           A[2] * sin(theta) * sin(varphi - phi) + 
           A[3] * cos(theta))
end


# beam shape coefficients
function generalMieAlm(Er, l, m, Rsurface, km; rtol=1e-6, maxevals=Int(1e4))
    #=
        Based on PHYSICAL REVIEW A 88, 063845 (2013)
        eq. (B3)
    
        Arguments:
            Er(r, theta, phi) -- function of the mode el. field
                                 in particle's spherical coordinate system
            l, m -- l = 1, Inf; m = -l, ..., +l
            Rsurface -- radius of the integration surface
            km -- wave vector in the media around the particle
    =#
    
    # q[1] -- theta
    # q[2] -- phi
    #qmin = SVector(0., 0.)
    #qmax = SVector(pi, 2pi)
    #q::StaticArrays.SArray{Tuple{2},Float64,1,2}
    qmin = [0., 0.]
    qmax = [pi, 2pi]

    integrand(q) = Er(Rsurface, q[1], q[2]) * conj(SphericalHarmonicY(l, m, q[1], q[2])) * sin(q[1])
    
    sintegral = hcubature(integrand, qmin, qmax; rtol=rtol, maxevals=maxevals)
    return Rsurface^2 / (l*(l+1) * riccatiPsi(l, km * Rsurface)) * sintegral[1]
    #return sintegral[1]
end


function generalMieBlm(Hr, l, m, Rsurface, km; rtol=1e-6, maxevals=Int(1e4))
    #=
        Based on PHYSICAL REVIEW A 88, 063845 (2013)
        eq. (B3)
    
        Arguments:
            Er(r, theta, phi) -- function of the mode el. field
                                 in particle's spherical coordinate system
            l, m -- l = 1, Inf; m = -l, ..., +l
            Rsurface -- radius of the integration surface
            km -- wave vector in the media around the particle
    =#
    
    # q[1] -- theta
    # q[2] -- phi
    #qmin = SVector(0., 0.)
    #qmax = SVector(pi, 2pi)
    #q::StaticArrays.SArray{Tuple{2},Float64,1,2}
    qmin = [0., 0.]
    qmax = [pi, 2pi]
    integrand(q) = Hr(Rsurface, q[1], q[2]) * conj(SphericalHarmonicY(l, m, q[1], q[2])) * sin(q[1])
    
    sintegral = hcubature(integrand, qmin, qmax; rtol=rtol, maxevals=maxevals) 
    # here the vacuum impedance is used
    # vacuumImpedance = Z0 = 1 / (c eps0) ~ 377 Ohm
    return Rsurface^2 * vacuumImpedance / (l*(l+1) * riccatiPsi(l, km * Rsurface)) * sintegral[1]
end


function generalMieBuildAlm(Er, lmax, Rsurface, km; rtol=1e-6, maxevals=Int(1e4))
    #=
        construct the matrix of coefficents Alm
    =#
    Alm = zeros(Complex{Float64}, lmax, 2lmax + 1)
    for l = 1:lmax, m = -l:l
        Alm[l, lmax+1 + m] = generalMieAlm(Er, l, m, Rsurface, km; rtol=rtol, maxevals=maxevals)
    end
    return Alm
end

function generalMieBuildBlm(Hr, lmax, Rsurface, km; rtol=1e-6, maxevals=Int(1e4))
    #=
        construct the matrix of coefficents Blm
    =#
    Blm = zeros(Complex{Float64}, lmax, 2lmax + 1)
    for l = 1:lmax, m = -l:l
        Blm[l, lmax+1 + m] = generalMieBlm(Hr, l, m, Rsurface, km; rtol=rtol, maxevals=maxevals)
    end
    return Blm
end


function generalMieFxFy(Rp, k0, epsp, epsm, AlmMATRIX, BlmMATRIX)
    #=
        based on PHYSICAL REVIEW A 88, 063845 (2013)
        eq. (14)
    
        Arguments:
            - Er, Er -- functions f(r, theta, phi) in particle's spherical coordinate system
            - Rp -- particle radius
            - k0 -- vacuum wave vector 
            - epsp, epsm -- electric permittivities (p - particle, m - media)
            - Alm, Blm -- field shape calculated matrix of coefficents
                            shapeof(A) = [lmax, 2lmax + 1]
    =#
    
    FxFy = 0.0im

    RTOL = 1e-6

    if length(AlmMATRIX[:, 1]) != length(BlmMATRIX[:, 1])
        println("Alm and Alm has different dimensions.")
        error()
        return 0.0im
    else
        lmax = length(AlmMATRIX[:, 1])
    end
    relativeIndex = sqrt(epsp / epsm)
    # size parameter
    km = k0 * sqrt(epsm)
    x = km * Rp
    
    function fAlm(l, m)
        if l > lmax
            return 0.0im
        elseif abs(m) > l
            return 0.0im
        else
            return AlmMATRIX[l, lmax+1 + m]
        end
    end

    function fBlm(l, m)
        if l > lmax
            return 0.0im
        elseif abs(m) > l
            return 0.0im
        else
            return BlmMATRIX[l, lmax+1 + m]
        end
    end

    err = 1.0
    for l = 1:lmax
        FxFy_l = 0.0im
        for m = -l:l
            # Mie coefficents
            al  = Mie_an(  l, relativeIndex, x)
            al1 = Mie_an(l+1, relativeIndex, x)
            bl  = Mie_bn(  l, relativeIndex, x)
            bl1 = Mie_bn(l+1, relativeIndex, x)
            
            # generalized Mie coefficents
            Alm    = fAlm(  l, m)
            Al1m1  = fAlm(l+1, m+1)
            Al1m_1 = fAlm(l+1, m-1)
            Alm1   = fAlm(  l, m+1)
            
            Blm    = fBlm(  l, m)
            Bl1m1  = fBlm(l+1, m+1)
            Bl1m_1 = fBlm(l+1, m-1)
            Blm1   = fBlm(  l, m+1)
            
            alm    = - al  * Alm
            al1m1  = - al1 * Al1m1
            al1m_1 = - al1 * Al1m_1
            alm1   = - al  * Alm1
            
            blm    = - bl  * Blm
            bl1m1  = - bl1 * Bl1m1
            bl1m_1 = - bl1 * Bl1m_1
            blm1   = - bl  * Blm1
            
            
            term1 = epsm * (2alm * conj(al1m1) + alm * conj(Al1m1) + Alm * conj(al1m1))
            term1 += 2blm * conj(bl1m1) + blm * conj(Bl1m1) + Blm * conj(bl1m1)
            
            term2 = epsm * (2al1m_1 * conj(alm) + al1m_1 * conj(Alm) + Al1m_1 * conj(alm))
            term2 += 2bl1m_1 * conj(blm) + bl1m_1 * conj(Blm) + Bl1m_1 * conj(blm)
            
            term3 = sqrt(epsm) * sqrt((l+m+1) * (l - m)) * (2alm * conj(blm1) + alm * conj(Blm1) + Alm * conj(blm1) -
                                                            2blm * conj(alm1) - blm * conj(Alm1) - Blm * conj(alm1))
            
            FxFy_l += l * (l + 2) * (sqrt((l+m+2)*(l+m+1) / ((2l+1)*(2l+3))) * term1 + sqrt((l-m+2)*(l-m+1) / ((2l+1)*(2l+3))) * term2) + term3
        end
        FxFy += FxFy_l
        
        err = abs(FxFy_l) / abs(FxFy)

        if err < RTOL
            break
        end
    end
    #@printf("RELERR = %.1e for x = %.2f \n", err, k0*sqrt(epsm)*Rp)
        
    return 0.25im * epsilon0 * km^2 * FxFy 
end


function generalMieFz(Rp, k0, epsp, epsm, AlmMATRIX, BlmMATRIX)
    #=
        based on PHYSICAL REVIEW A 88, 063845 (2013)
        eq. (15)
    
        Arguments:
            - Rp -- particle radius
            - k0 -- vacuum wave vector 
            - epsp, epsm -- electric permittivities (p - particle, m - media)
            - AlmMATRIX, BlmMATRIX -- coefficents of field decompositions
    =#
    
    RTOL = 1e-6

    if length(AlmMATRIX[:, 1]) != length(BlmMATRIX[:, 1])
        error("Alm and Alm has different dimensions.")
        return 0.0
    else
        lmax = length(AlmMATRIX[:, 1])
    end

    function fAlm(l, m)
        if l > lmax
            return 0.0im
        elseif abs(m) > l
            return 0.0im
        else
            return AlmMATRIX[l, lmax+1 + m]
        end
    end

    function fBlm(l, m)
        if l > lmax
            return 0.0im
        elseif abs(m) > l
            return 0.0im
        else
            return BlmMATRIX[l, lmax+1 + m]
        end
    end

    Fz = 0.0

    relativeIndex = sqrt(epsp / epsm)
    # size parameter
    km = k0 * sqrt(epsm)
    x = km * Rp
    
    err = 1.0
    for l = 1:lmax
        Fz_l = 0.0
        for m = -l:l
            # Mie coefficents
            al  = Mie_an(  l, relativeIndex, x)
            al1 = Mie_an(l+1, relativeIndex, x)
            bl  = Mie_bn(  l, relativeIndex, x)
            bl1 = Mie_bn(l+1, relativeIndex, x)
            
            # generalized Mie coefficents
            Alm  = fAlm(  l, m)
            Al1m = fAlm(l+1, m)
            
            Blm  = fBlm(  l, m)
            Bl1m = fBlm(l+1, m)

            
            alm  = - al  * Alm
            al1m = - al1 * Al1m
            
            blm  = - bl  * Blm
            bl1m = - bl1 * Bl1m
            
            term1 = l * (l + 2) * sqrt((l-m+1) * (l+m+1) / ((2l+1)*(2l+3))) * (
                        epsm * (2al1m * conj(alm) + al1m * conj(Alm) + Al1m * conj(alm)) 
                        + 2bl1m * conj(blm) + bl1m * conj(Blm) + Bl1m * conj(blm))
            term2 = m * sqrt(epsm) * (2alm * conj(blm) + alm * conj(Blm) + Alm * conj(blm))
            
            Fz_l += imag(term1 + term2)
        end
        Fz += Fz_l
        
        err = abs(Fz_l) / abs(Fz)
        if err < RTOL
            break
        end
    end
    #@printf("RELERR = %.1e for x = %.2f \n", err, k0*sqrt(epsm)*Rp)
        
    return -0.5 * epsilon0 * km^2 * Fz
end


function generalMieTx(Rp, k0, epsp, epsm, AlmMATRIX, BlmMATRIX)
    #=
        based on eq. (10) from
        Barton, J. P., D. R. Alexander, and S. A. Schaub. 
        "Theoretical determination of net radiation force 
        and torque for a spherical particle illuminated by a 
        focused laser beam." Journal of Applied Physics 66.10 (1989)
        
        but converted to the SI units
    
        Arguments:
            - Rp -- particle radius
            - k0 -- vacuum wave vector 
            - epsp, epsm -- electric permittivities (p - particle, m - media)
            - AlmMATRIX, BlmMATRIX -- coefficents of field decompositions
    =#

    if length(AlmMATRIX[:, 1]) != length(BlmMATRIX[:, 1])
        println("Alm and Alm has different dimensions.")
        error()
        return 0.0
    else
        lmax = length(AlmMATRIX[:, 1])
    end

    function fAlm(l, m)
        if l > lmax
            return 0.0im
        elseif abs(m) > l
            return 0.0im
        else
            return AlmMATRIX[l, lmax+1 + m]
        end
    end

    function fBlm(l, m)
        if l > lmax
            return 0.0im
        elseif abs(m) > l
            return 0.0im
        else
            return BlmMATRIX[l, lmax+1 + m]
        end
    end

    Tx = 0.0

    RTOL = 1e-6

    relativeIndex = sqrt(epsp / epsm)
    # size parameter
    km = k0 * sqrt(epsm)
    x = km * Rp
    
    err = 1.0
    for l = 1:lmax
        Tx_l = 0.0
        for m = -l:l
            # Mie coefficents
            al = Mie_an(l, relativeIndex, x)
            bl = Mie_bn(l, relativeIndex, x)
            
            # generalized Mie coefficents
            Alm  = fAlm(l, m)
            Alm1 = fAlm(l, m+1)
            
            Blm  = fBlm(l, m)
            Blm1 = fBlm(l, m+1)

            
            alm  = - al * Alm
            alm1 = - al * Alm1
            
            blm  = - bl * Blm
            blm1 = - bl * Blm1
            
            tmp = l*(l+1) * sqrt((l-m) * (l+m+1)) * (
                epsm * alm * conj(alm1) + blm * conj(blm1) + 
                0.5 * (epsm * alm * conj(Alm1) + epsm * alm1 * conj(Alm) +
                      blm * conj(Blm1) + blm1 * conj(Blm))
            )
            
            Tx_l += real(tmp)
        end
        Tx += Tx_l
        
        err = abs(Tx_l) / abs(Tx)
        if err < RTOL
            break
        end
    end
    #@printf("RELERR = %.1e for x = %.2f \n", err, k0*sqrt(epsm)*Rp)
        
    return - epsilon0/2 * k0 * sqrt(epsm) * Tx
end


function generalMieTy(Rp, k0, epsp, epsm, AlmMATRIX, BlmMATRIX)
    #=
        based on eq. (11) from
        Barton, J. P., D. R. Alexander, and S. A. Schaub. 
        "Theoretical determination of net radiation force 
        and torque for a spherical particle illuminated by a 
        focused laser beam." Journal of Applied Physics 66.10 (1989)
        
        but converted to the SI units
    
        Arguments:
            - Rp -- particle radius
            - k0 -- vacuum wave vector 
            - epsp, epsm -- electric permittivities (p - particle, m - media)
            - AlmMATRIX, BlmMATRIX -- coefficents of field decompositions
    =#
    RTOL = 1e-6

    if length(AlmMATRIX[:, 1]) != length(BlmMATRIX[:, 1])
        error("Alm and Alm has different dimensions.")
        return 0.0
    else
        lmax = length(AlmMATRIX[:, 1])
    end

    function fAlm(l, m)
        if l > lmax
            return 0.0im
        elseif abs(m) > l
            return 0.0im
        else
            return AlmMATRIX[l, lmax+1 + m]
        end
    end

    function fBlm(l, m)
        if l > lmax
            return 0.0im
        elseif abs(m) > l
            return 0.0im
        else
            return BlmMATRIX[l, lmax+1 + m]
        end
    end

    Ty = 0.0

    relativeIndex = sqrt(epsp / epsm)
    # size parameter
    km = k0 * sqrt(epsm)
    x = km * Rp
    
    err = 1.0
    for l = 1:lmax
        Ty_l = 0.0
        for m = -l:l
            # Mie coefficents
            al = Mie_an(l, relativeIndex, x)
            bl = Mie_bn(l, relativeIndex, x)
            
            # generalized Mie coefficents
            Alm  = fAlm(l, m)
            Alm1 = fAlm(l, m+1)
            
            Blm  = fBlm(l, m)
            Blm1 = fBlm(l, m+1)

            
            alm  = - al * Alm
            alm1 = - al * Alm1
            
            blm  = - bl * Blm
            blm1 = - bl * Blm1
            
            tmp = l*(l+1)*sqrt((l-m) * (l+m+1)) * (
                epsm * alm * conj(alm1) + blm * conj(blm1) + 
                0.5 * (epsm * alm * conj(Alm1) - epsm * alm1 * conj(Alm) + 
                      blm * conj(Blm1) - blm1 * conj(Blm))
            )
            
            Ty_l += imag(tmp)
        end
        Ty += Ty_l
        
        err = abs(Ty_l) / abs(Ty)
        #@printf("RELERR = %.1e for x = %.2f \n", err, k0*sqrt(epsm)*Rp)
        if err < RTOL
            break
        end
    end
    #@printf("RELERR = %.1e for x = %.2f \n\n", err, k0*sqrt(epsm)*Rp)
        
    return - epsilon0/2 * k0 * sqrt(epsm) * Ty
end


function generalMieTz(Rp, k0, epsp, epsm, AlmMATRIX, BlmMATRIX)
    #=
        based on eq. (12) from
        Barton, J. P., D. R. Alexander, and S. A. Schaub. 
        "Theoretical determination of net radiation force 
        and torque for a spherical particle illuminated by a 
        focused laser beam." Journal of Applied Physics 66.10 (1989)
        
        but converted to the SI units
    
        Arguments:
            - Rp -- particle radius
            - k0 -- vacuum wave vector 
            - epsp, epsm -- electric permittivities (p - particle, m - media)
            - AlmMATRIX, BlmMATRIX -- coefficents of field decompositions
    =#

    if length(AlmMATRIX[:, 1]) != length(BlmMATRIX[:, 1])
        println("Alm and Alm has different dimensions.")
        error()
        return 0.0
    else
        lmax = length(AlmMATRIX[:, 1])
    end

    function fAlm(l, m)
        if l > lmax
            return 0.0im
        elseif abs(m) > l
            return 0.0im
        else
            return AlmMATRIX[l, lmax+1 + m]
        end
    end

    function fBlm(l, m)
        if l > lmax
            return 0.0im
        elseif abs(m) > l
            return 0.0im
        else
            return BlmMATRIX[l, lmax+1 + m]
        end
    end

    Tz = 0.0

    RTOL = 1e-6

    relativeIndex = sqrt(epsp / epsm)
    # size parameter
    km = k0 * sqrt(epsm)
    x = km * Rp
    
    err = 1.0
    for l = 1:lmax
        Tz_l = 0.0
        for m = -l:l
            # Mie coefficents
            al = Mie_an(l, relativeIndex, x)
            bl = Mie_bn(l, relativeIndex, x)
            
            # generalized Mie coefficents
            Alm  = fAlm(l, m)
            Blm  = fBlm(l, m)

            alm  = - al * Alm
            blm  = - bl * Blm
            
            tmp = (l+1) * m * (epsm * abs2(alm) + abs2(blm) + real(epsm * alm * conj(Alm) + blm * conj(Blm)))

            Tz_l += tmp
        end
        Tz += Tz_l
        
        err = abs(Tz_l) / abs(Tz)
        if err < RTOL 
            break
        end
    end
    #@printf("RELERR = %.1e for x = %.2f \n", err, k0*sqrt(epsm)*Rp)
        
    return - epsilon0/2 * k0 * sqrt(epsm) * Tz
end


function vectorHarmonicsMnm(n, m, rho, theta, phi, superscript=1)
    #=
        complex vector Harmonics in spherical coordinates
        from Bohren & Huffman, eq. (4.17), (4.18)

        sin m \phi -> - i exp(i m \phi)
        cos m \phi -> exp(i m \phi)

        M^even = Re M
        M^odd  = Im M

        Input
        -----
            rho = k r, k = k0 n_media 
            theta = [0, pi] -- azimuthal angle
            phi = [0, 2pi] -- polar angle
            n = 0, 1, ...
            m = 0, ..., n ?????
            superscript = 1, 2, 3, 4
    =#
    Mt = 0
    Mp = 0

    zn = 0
    if superscript == 1
        # zn = jn -- spherical Bessel 1st kind 
        zn = sphericaljn(n, rho)
    elseif superscript == 2
        # zn = yn -- sherical Bessel 2nd kind
        zn = sphericalyn(n, rho)
    elseif superscript == 3
        # zn = h(1) -- spherical Hankel 1st
        zn = sphericalh1(n, rho)
    elseif superscript == 4
        # zn = h(2) -- spherical Hankel 2nd
        zn = sphericalh2(n, rho)
    else
        error("Superscript invalid!")
    end

    Mt = im * m / sin(theta) * LegendrePnm(n, m, cos(theta))

    # lazy to write a proper derivative
    dtheta = 1e-12
    dPnm = (LegendrePnm(n, m, cos(theta + dtheta)) - LegendrePnm(n, m, cos(theta - dtheta))) / 2dtheta

    Mp = dPnm

    return [0, Mt, Mp] * exp(im * m * phi) * zn
end


function vectorHarmonicsNnm(n, m, rho, theta, phi, superscript=1)
    #=
        complex vector Harmonics in spherical coordinates
        from Bohren & Huffman, eq. (4.17), (4.18)

        sin m \phi -> - i exp(i m \phi)
        cos m \phi -> exp(i m \phi)

        M^even = Re M
        M^odd  = Im M

        Input
        -----
            rho = k r, k = k0 n_media 
            theta = [0, pi] -- azimuthal angle
            phi = [0, 2pi] -- polar angle
            n = 0, 1, ...
            m = 0, ..., n ?????
            superscript = 1, 2, 3, 4
    =#
    Nr = 0
    Nt = 0
    Np = 0

    zn = 0
    znp = 0
    if superscript == 1
        # zn = jn -- spherical Bessel 1st kind 
        zn = sphericaljn(n, rho)
        znp = sphericaljnp(n, rho)
    elseif superscript == 2
        # zn = yn -- sherical Bessel 2nd kind
        zn = sphericalyn(n, rho)
        znp = sphericalynp(n, rho)
    elseif superscript == 3
        # zn = h(1) -- spherical Hankel 1st
        zn = sphericalh1(n, rho)
        znp = sphericalh1p(n, rho)
    elseif superscript == 4
        # zn = h(2) -- spherical Hankel 2nd
        zn = sphericalh2(n, rho)
        znp = sphericalh2p(n, rho)
    else
        error("Superscript invalid!")
    end

    Pnm = LegendrePnm(n, m, cos(theta))
    dtheta = 1e-12
    dPnm = (LegendrePnm(n, m, cos(theta + dtheta)) - LegendrePnm(n, m, cos(theta - dtheta))) / 2dtheta

    Nr = zn/rho * n * (n+1) * Pnm
    Nt = dPnm * (zn/rho + znp)
    Np = im * m * Pnm / sin(theta) * (zn/rho + znp)

    return [Nr, Nt, Np] * exp(im * m * phi)
end