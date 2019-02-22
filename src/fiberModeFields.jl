# by Ivan Toftul
#
# created          2018/07
# last edit        2018/12

using SpecialFunctions
using LinearAlgebra
include("./consts.jl")

export
    Mode,
    coefHybridElectric,
    coefHybridMagnetic,
    coefTEElectric,
    coefTMMagnetic,
    stokes2Jones,
    P2A,
    eHybrid,
    eHybridPrime,
    eTE,
    eTEPrime,
    eTM,
    eTMPrime,
    eFiberProfile,
    ePrimeFiberProfile,
    eFiber,
    ePrimeFiber,
    intensity,
    eFiberCustomPol,
    ePrimeFiberCustomPol,
    hHybrid,
    hTE,
    hTM,
    hFiber,
    hFiberCustomPol,
    setUpPrm4Fast,
    eHybridFast,
    eTEFast,
    eTMFast,
    eFiberFast,
    eFiberCustomPolFast,
    hHybridFast,
    hTEFast,
    hTMFast,
    hFiberFast,
    hFiberCustomPolFast

# ## fucntions below are based on the Fam's paper Phys Rev A 96, 2017

# custom type wich contains all the informaiton about the mode
struct Mode
    modeType  # HE, EH, TE, TM
    l  # azimuthal number
    m  # radial mode order
    f  # propagation direction
    p  # phase circulation direction (+1 -- CCW (R), -1 -- CW (L))
    beta  # propagation constant
end


# structures for the faster versions of the field profiles
struct coefHybridElectric
    # hybrid modes coefficents
    Ahless_r1
    Ahless_r2
    Ahless_p1
    Ahless_p2
    Ahless_z1

    Ahmore_r1
    Ahmore_r2
    Ahmore_p1
    Ahmore_p2
    Ahmore_z1
end

struct coefHybridMagnetic
    # hybrid modes coefficents
    Bhless_r1
    Bhless_r2
    Bhless_p1
    Bhless_p2
    Bhless_z1

    Bhmore_r1
    Bhmore_r2
    Bhmore_p1
    Bhmore_p2
    Bhmore_z1
end

struct coefTEElectric
    # TE modes coefficents
    ATEless_p1

    ATEmore_p1
end

struct coefTEMagnetic
    # TE modes coefficents
    BTEless_r1
    BTEless_z1

    BTEmore_r1
    BTEmore_z1
end

struct coefTMElectric
    # TM modes coefficents
    ATMless_r1
    ATMless_z1

    ATMmore_r1
    ATMmore_z1
end

struct coefTMMagnetic
    # TM modes coefficents
    BTMless_p1
    BTMmore_p1
end


# polarization control
function stokes2Jones(S)
    # see details here
    # https://physics.stackexchange.com/questions/238957/converting-stokes-parameters-to-jones-vector

    # p = the degree of polarization from the Stokes vector so you can do
    # something with the unpolarized part separately if you like
    p = sqrt(S[2]^2 + S[3]^2 + S[4]^2) / S[1]

    # Normalize the Stokes parameters (first one will be 1, of course)
    Q = S[2] / (S[1] * p)
    U = S[3] / (S[1] * p)
    V = S[4] / (S[1] * p)

    # And construct the 2 Jones components
    A = sqrt((1+Q) / 2)
    B = 1.0
    if A != 0.0
       B = (U - im * V) / 2A
    end

    Jv = sqrt(S[1] * p) * [A, B]

    return Jv
end


function P2A(modeType, a, n1, n2, k0, l, beta, P)
    #=
        finds the const A which is respectful for the field amplitude
        based on PHYSICAL REVIEW A 96, 023835 (2017) Appendix C: Power

        # Arguments
        - 'modeType' : HE, EH, TE or TM.
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `l::Int`: Azymuthal mode number.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.
    =#
    h = sqrt((n1*k0)^2 - beta^2)
    q = sqrt(beta^2 - (n2*k0)^2)
    ha = h * a
    qa = q * a

    if ((modeType == "HE") || (modeType == "EH"))
        Jplha = .5 * (besselj(l - 1, ha) - besselj(l + 1, ha))
        Jlha = besselj(l, ha)
        Kplqa = - .5 * (besselk(l - 1, qa) + besselk(l + 1, qa))
        Klqa = besselk(l, qa)
        s = l * (1/ha^2 + 1/qa^2) / (Jplha / (ha * Jlha) + Kplqa / (qa * Klqa))
        s1 = beta^2 * s / (k0 * n1)^2
        s2 = beta^2 * s / (k0 * n2)^2

        inBraket = ((1 - s) * (1 - s1) * (besselj(l-1, ha)^2 - besselj(l-2, ha) * Jlha) +
                    (1 + s) * (1 + s1) * (besselj(l+1, ha)^2 - besselj(l+2, ha) * Jlha))
        outBraket = ((1 - s) * (1 - s2) * (besselk(l-2, qa) * Klqa - besselk(l-1, qa)^2) +
                    (1 + s) * (1 + s2) * (besselk(l+2, qa) * Klqa - besselk(l+1, qa)^2))
        braketH = inBraket * (n1/h)^2 + outBraket * (Jlha/Klqa * n2/q)^2

        Ah = 2sqrt(vacuumImpedance * P / (pi * a^2 * k0 * beta * braketH))
        return Ah
    elseif modeType == "TE"
        J0ha = besselj0(ha)
        K0qa = besselk(0, qa)

        inBraket = besselj(1, ha)^2 - J0ha * besselj(2, ha)
        outBraket = K0qa * besselk(2, qa) - besselk(1, qa)^2
        braketTE = inBraket / h^2 + outBraket * (J0ha/K0qa/q)^2

        ATE = sqrt(2k0 / beta * vacuumImpedance * P / (pi * a^2)) / sqrt(braketTE) / (k0 * vacuumImpedance)
        return ATE
    elseif modeType == "TM"
        J0ha = besselj0(ha)
        K0qa = besselk(0, qa)
        J0ha_K0qa = J0ha / K0qa

        inBraket = besselj(1, ha)^2 - J0ha * besselj(2, ha)
        outBraket = K0qa * besselk(2, qa) - besselk(1, qa)^2
        braketTM = (n1/h)^2 * inBraket + (n2/q * J0ha_K0qa)^2 * outBraket

        ATM = sqrt(2vacuumImpedance * P / (pi*a^2 * beta * k0 * braketTM))
        return ATM
    else
        println("ERROR: Incorrect mode name.")
        error()
    end
end


# ############ ELECTRIC FIELDS ############
# tested

# Hybrid modes
# (A8) - (A12)
function eHybrid(r, a, n1, n2, k0, l, beta, P)
    #=
        Mode profile for the hybrid modes of a stepindex cylindrical fiber

        ...
        # Arguments
        - `r::Float`: Radial coordinate.
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `l::Int`: Azymuthal mode number.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.


        # Returns
        - 'er, ephi, ez': complex vector-column
        ...
    =#

    h = sqrt((n1*k0)^2 - beta^2)
    q = sqrt(beta^2 - (n2*k0)^2)
    ha = h * a
    qa = q * a

    Jplha = .5 * (besselj(l - 1, ha) - besselj(l + 1, ha))
    Jlha = besselj(l, ha)
    Kplqa = - .5 * (besselk(l - 1, qa) + besselk(l + 1, qa))
    Klqa = besselk(l, qa)
    s = l * (1/ha^2 + 1/qa^2) / (Jplha / (ha * Jlha) + Kplqa / (qa * Klqa))
    s1 = beta^2 * s / (k0 * n1)^2
    s2 = beta^2 * s / (k0 * n2)^2

    Jlha_Klqa = Jlha/Klqa

    inBraket = ((1 - s) * (1 - s1) * (besselj(l-1, ha)^2 - besselj(l-2, ha) * Jlha) +
                (1 + s) * (1 + s1) * (besselj(l+1, ha)^2 - besselj(l+2, ha) * Jlha))
    outBraket = ((1 - s) * (1 - s2) * (besselk(l-2, qa) * Klqa - besselk(l-1, qa)^2) +
                 (1 + s) * (1 + s2) * (besselk(l+2, qa) * Klqa - besselk(l+1, qa)^2))
    braketH = inBraket * (n1/h)^2 + outBraket * (Jlha_Klqa * n2/q)^2

    Ch = beta * sqrt(vacuumImpedance * P / (pi * a^2 * k0 * beta * braketH))
    # connection with Fam notation
    # Ah = 2/beta * Ch

    er = 0.0im
    ephi = 0.0im
    ez = 0.0im
    if r < a
        hr = h * r
        Jl_1hr = besselj(l - 1, hr)
        Jlplus1hr = besselj(l + 1, hr)

        er   =  1im * Ch/h * ((1 - s) * Jl_1hr - (1 + s) * Jlplus1hr)
        ephi = - Ch/h * ((1 - s) * Jl_1hr + (1 + s) * Jlplus1hr)
        ez   = 2Ch/beta * besselj(l, hr)
    else
        qr = q * r
        Kl_1qr = besselk(l - 1, qr)
        Klplus1qr = besselk(l + 1, qr)

        er = 1im * Ch/q * Jlha_Klqa * ((1 - s) * Kl_1qr + (1 + s) * Klplus1qr)
        ephi =  - Ch/q * Jlha_Klqa * ((1 - s) * Kl_1qr - (1 + s) * Klplus1qr)
        ez =  2Ch / beta * Jlha_Klqa * besselk(l, qr)
    end

    return [er, ephi, ez]
end

function eHybridPrime(r, a, n1, n2,
                      k0, l, beta, P)
    #=
    Derivative with respect to r

    ...
    # Arguments
    - `r::Float`: Radial coordinate.
    - `a::Float`: Fiber radius.
    - `n1::Float`: Fiber epsilon.
    - `n2::Float`: Outside epsilon.
    - `k0::Float`: Vacuum wave number k0 = omega/c.
    - `l::Int`: Azymuthal mode number.
    - `beta::Float`: Propogation constant.
    - `P::Float`: Total mode power.


    # Returns
    - 'erp, ephip, ezp': complex vector-column
    ...
    =#

    h = sqrt((n1*k0)^2 - beta^2)
    q = sqrt(beta^2 - (n2*k0)^2)
    ha = h * a
    qa = q * a

    Jplha = .5 * (besselj(l - 1, ha) - besselj(l + 1, ha))
    Jlha = besselj(l, ha)
    Kplqa = - .5 * (besselk(l - 1, qa) + besselk(l + 1, qa))
    Klqa = besselk(l, qa)
    s = l * (1/ha^2 + 1/qa^2) / (Jplha / (ha * Jlha) + Kplqa / (qa * Klqa))
    s1 = beta^2 * s / (k0 * n1)^2
    s2 = beta^2 * s / (k0 * n2)^2
    Jlha_Klqa = Jlha/Klqa

    inBraket = ((1 - s) * (1 - s1) * (besselj(l-1, ha)^2 - besselj(l-2, ha) * Jlha) +
                         (1 + s) * (1 + s1) * (besselj(l+1, ha)^2 - besselj(l+2, ha) * Jlha))
    outBraket = ((1 - s) * (1 - s2) * (besselk(l-2, qa) * Klqa - besselk(l-1, qa)^2) +
                          (1 + s) * (1 + s2) * (besselk(l+2, qa) * Klqa - besselk(l+1, qa)^2))
    braketH = inBraket * (n1/h)^2 + outBraket * (Jlha_Klqa * n2/q)^2

    Ch = beta * sqrt(vacuumImpedance * P / (pi * a^2 * k0 * beta)) / sqrt(braketH)


    erPrime   = 0.0im
    ephiPrime = 0.0im
    ezPrime   = 0.0im
    if r < a
        # not yet done
        error()
    else
        qr = q * r
        Kl_2qr = besselk(l - 2, qr)
        Klplus2qr = besselk(l + 2, qr)
        Klqr = besselk(l, qr)

        erPrime = - q/2 * 1im * Ch/q * Jlha_Klqa * (2Klqr + (1 - s) * Kl_2qr + (1 + s) * Klplus2qr)
        ephiPrime =  q/2 * Ch/q * Jlha_Klqa * ((1 - s) * Kl_2qr - (1 + s) * Klplus2qr - 2s * Klqr)
        ezPrime =  -q/2 * 2Ch / beta * Jlha_Klqa * (besselk(l - 1, qr) + besselk(l + 1, qr))
    end

    return [erPrime, ephiPrime, ezPrime]
end


# Transverse modes
function eTE(r, a, n1, n2, k0, beta, P)
    #=
        Mode profile for the TE modes with l = 0 of a stepindex cylindrical fiber

        ...
        # Arguments
        - `r::Float`: Radial coordinate.
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.


        # Returns
        - 'er, ephi, ez': complex vector-column
        ...
    =#

    h = sqrt((n1*k0)^2 - beta^2)
    q = sqrt(beta^2 - (n2*k0)^2)
    ha = h * a
    qa = q * a
    J0ha = besselj0(ha)
    K0qa = besselk(0, qa)

    inBraket = besselj(1, ha)^2 - J0ha * besselj(2, ha)
    outBraket = K0qa * besselk(2, qa) - besselk(1, qa)^2
    braketTE = inBraket / h^2 + outBraket * (J0ha/K0qa/q)^2

    CTE = 1im * sqrt(2k0 / beta * vacuumImpedance * P / (pi * a^2)) / sqrt(braketTE)
    # connection with Fam notation
    # A_TE = CTE / (k0 * vacuumImpedance)

    ephi = 0.0im
    if r < a
        ephi = CTE / h * besselj1(h * r)
    else
        ephi = - CTE / q * J0ha/K0qa * besselk(1, q*r)
    end
    return [0, ephi, 0]
end

function eTEPrime(r, a, n1, n2, k0, beta, P)
    #=
        Derivative with respect to r

        ...
        # Arguments
        - `r::Float`: Radial coordinate.
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.


        # Returns
        - 'er, ephi, er': complex vector-column
        ...
    =#

    h    = sqrt((n1*k0)^2 - beta^2)
    q    = sqrt(beta^2 - (n2*k0)^2)
    ha   = h * a
    qa   = q * a
    J0ha = besselj0(ha)
    K0qa = besselk(0, qa)

    inBraket = besselj(1, ha)^2 - J0ha * besselj(2, ha)
    outBraket = K0qa * besselk(2, qa) - besselk(1, qa)^2
    braketTE = inBraket / h^2 + outBraket * (J0ha/K0qa/q)^2

    CTE = 1im * sqrt(2k0 / beta * vacuumImpedance * P / (pi * a^2)) / sqrt(braketTE)

    ephiPrime = 0.0im
    if r < a
        ephiPrime = CTE / h * h/2 * (besselj0(h * r) - besselj(2, h * r))
    else
        ephiPrime = q/2 * CTE / q * J0ha/K0qa * (besselk(0, q*r) + besselk(2, q*r) )
    end
    return [0, ephiPrime, 0]
end


function eTM(r, a, n1, n2, k0, beta, P)
    #=
        Mode profile for the TM modes with l = 0 of a stepindex cylindrical fiber

        ...
        # Arguments
        - `r::Float`: Radial coordinate.
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.


        # Returns
        - 'er, ephi, er': complex vector-column
        ...
    =#

    h  = sqrt((n1*k0)^2 - beta^2)
    q  = sqrt(beta^2 - (n2*k0)^2)
    ha = h * a
    qa = q * a
    J0ha = besselj0(ha)
    K0qa = besselk(0, qa)
    J0ha_K0qa = J0ha / K0qa

    inBraket = besselj(1, ha)^2 - J0ha * besselj(2, ha)
    outBraket = K0qa * besselk(2, qa) - besselk(1, qa)^2
    braketTM = (n1/h)^2 * inBraket + (n2/q * J0ha_K0qa)^2 * outBraket

    CTM = sqrt(2vacuumImpedance * P / (pi*a^2 * beta * k0)) / sqrt(braketTM)
    # connection with Fam notation
    # A_TM = CTM

    er = 0.0im
    ez = 0.0im
    if r < a
        hr = h * r

        er = - 1im * beta/h * CTM * besselj1(hr)
        ez = CTM * besselj0(hr)
    else
        qr = q * r

        er = 1im * beta/q * J0ha_K0qa * CTM * besselk(1, qr)
        ez = J0ha_K0qa * CTM * besselk(0, qr)
    end
    return [er, 0, ez]
end

function eTMPrime(r, a, n1, n2, k0, beta, P)
    #=
        Derivative with the respect to r

        ...
        # Arguments
        - `r::Float`: Radial coordinate.
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.


        # Returns
        - 'er, ephi, er': complex vector-column
        ...
    =#

    h = sqrt((n1*k0)^2 - beta^2)
    q = sqrt(beta^2 - (n2*k0)^2)
    ha = h * a
    qa = q * a
    J0ha = besselj0(ha)
    K0qa = besselk(0, qa)
    J0ha_K0qa = J0ha / K0qa

    inBraket = besselj(1, ha)^2 - J0ha * besselj(2, ha)
    outBraket = K0qa * besselk(2, qa) - besselk(1, qa)^2
    braketTM = (n1/h)^2 * inBraket + (n2/q * J0ha_K0qa)^2 * outBraket

    CTM = sqrt(2vacuumImpedance * P / (pi*a^2 * beta * k0)) / sqrt(braketTM)

    erPrime = 0.0im
    ezPrime = 0.0im
    if r < a
        hr = h * r

        erPrime = - 1im * beta/h * CTM * h/2 * (besselj0(hr) - besselj(2, hr))
        ezPrime = - CTM * h * besselj1(hr)
    else
        qr = q * r

        erPrime = - q/2 * 1im * beta/q * J0ha_K0qa * CTM * (besselk(0, qr) + besselk(2, qr))
        ezPrime = - q * J0ha_K0qa * CTM * besselk(1, qr)
    end
    return [erPrime, 0, ezPrime]
end


function eFiberProfile(modeType, r, a, n1, n2,
                       k0, l, beta, P)
    #=
        Mode profile of a stepindex cylindrical fiber
        e = (er, ephi, ez)

        ...
        # Arguments
        - `modeType` : HE, EH, TE, or TM.
        - `r::Float`: Radial coordinate.
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `l::Int`: Azymuthal mode number.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.


        # Returns
        - 'er, ephi, er': complex vector-column
        ...
    =#

    if ((modeType == "HE") || (modeType == "EH"))
        return eHybrid(r, a, n1, n2, k0, l, beta, P)
    elseif modeType == "TE"
        return eTE(r, a, n1, n2, k0, beta, P)
    elseif modeType == "TM"
        return eTM(r, a, n1, n2, k0, beta, P)
    else
        println("Error! Write the proper name of mode")
        error()
        return [0.0, 0.0, 0.0im]
    end
end

function ePrimeFiberProfile(modeType, r, a,
                            n1, n2, k0, l, beta, P)
    if ((modeType == "HE") || (modeType == "EH"))
        return eHybridPrime(r, a, n1, n2, k0, l, beta, P)
    elseif modeType == "TE"
        return eTEPrime(r, a, n1, n2, k0, beta, P)
    elseif modeType == "TM"
        return eTMPrime(r, a, n1, n2, k0, beta, P)
    else
        println("Error! Write the proper name of the mode")
        error()
        return [0.0, 0.0, 0.0im]
    end
end

function eFiber(mode, R, a,
                n1, n2, k0, P)
    #=
        Mode fields of a stepindex cylindrical fiber

        ...
        # Arguments
        - `mode::customStruct` : see def above.
        - `R::array` : R = (r, phi, z).
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `l::Int`: Azymuthal mode number.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.


        # Returns
        - '[er, ephi, ez]': complex vector-column
        ...
    =#
    eProfile = [0.0, 0.0, 0.0im]
    if ((mode.modeType == "HE") || (mode.modeType == "EH"))
        eProfile = eHybrid(R[1], a, n1, n2, k0, mode.l, mode.beta, P)

        return([eProfile[1],
                mode.p * eProfile[2],
                mode.f * eProfile[3]] *
                exp(im * mode.f * mode.beta * R[3] + im * mode.p * mode.l * R[2]))
    elseif mode.modeType == "TE"
        # for TE l = 0
        eProfile = eTE(R[1], a, n1, n2, k0, mode.beta, P)

        return([eProfile[1],
                mode.p * eProfile[2],
                mode.f * eProfile[3]] *
                exp(im * mode.f * mode.beta * R[3]))
    elseif mode.modeType == "TM"
        # for TM  l =0
        eProfile = eTM(R[1], a, n1, n2, k0, mode.beta, P)

        return([eProfile[1],
                mode.p * eProfile[2],
                mode.f * eProfile[3]] *
                exp(im * mode.f * mode.beta * R[3]))
    else
        println("Error! Write the proper name of mode")
        error()
        return [0.0, 0.0, 0.0im]
    end
end


function ePrimeFiber(mode, R, a, n1, n2, k0, P)
    #=
        Derivative with respect to rho
        Mode fields of a stepindex cylindrical fiber

        ...
        # Arguments
        - `mode::customStruct` : see def above.
        - `R::array` : R = (r, phi, z).
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `l::Int`: Azymuthal mode number.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.
        ...
    =#
    eProfile = [0.0, 0.0, 0.0im]
    if ((mode.modeType == "HE") || (mode.modeType == "EH"))
        eProfile = eHybridPrime(R[1], a, n1, n2, k0, mode.l, mode.beta, P)
    elseif mode.modeType == "TE"
        eProfile = eTEPrime(R[1], a, n1, n2, k0, mode.beta, P)
    elseif mode.modeType == "TM"
        eProfile = eTMPrime(R[1], a, n1, n2, k0, mode.beta, P)
    else
        println("Error! Write the proper name of mode")
        error()
        return [0.0, 0.0, 0.0im]
    end

    return([eProfile[1],
            mode.p * eProfile[2],
            mode.f * eProfile[3]] *
           exp(im * mode.f * mode.beta * R[3] + im * mode.p * mode.l * R[2]))
end

function intensity(modeType, r, a, n1, n2,
                   k0, l, beta, P)
    e = [0.0, 0.0, 0.0im]
    if ((modeType == "HE") || (modeType == "EH"))
        e = eHybrid(r, a, n1, n2, k0, l, beta, P)
        return norm(e)^2
    elseif modeType == "TE"
        e = eTE(r, a, n1, n2, k0, beta, P)
        return norm(e)^2
    elseif modeType == "TM"
        e = eTM(r, a, n1, n2, k0, beta, P)
        return norm(e)^2
    else
        println("Error! Write the proper name of mode")
        error()
        return 0.0
    end
end


function eFiberCustomPol(S, mode, R, a,
                         n1, n2, k0, P)
    #=
        Mode fields of a stepindex cylindrical fiber
        with custom polarization defined by an arbirtary Stokes vector S.
        Custom polarization only for the hybrid modes

        ...
        # Arguments
        - `S = [S0, S1, S2, S3]` : Stokes vector.
        - `mode::customStruct` : see def above.
        - `R::array` : R = (r, phi, z).
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `l::Int`: Azymuthal mode number.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.


        # Returns
        - '[er, ephi, ez]': complex vector-column
        ...
    =#
    if ((mode.modeType == "HE") || (mode.modeType == "EH"))
        # hybrid mode are circularly polarized by default
        # orthogonal basis
        modePlus = Mode(mode.modeType,
            mode.l,
            mode.m,
            mode.f,
            +1,  # p - phase circulation direction CCW
            mode.beta)

        modeMinus = Mode(mode.modeType,
            mode.l,
            mode.m,
            mode.f,
            -1,  # p - phase circulation direction CW
            mode.beta)

        Eplus  = eFiber(modePlus,  R, a, n1, n2, k0, P)
        Eminus = eFiber(modeMinus, R, a, n1, n2, k0, P)

        Ex = 1/sqrt(2) * (Eplus + Eminus)
        Ey = -1im/sqrt(2) * (Eplus - Eminus)

        # Stokes to Jones vector
        Jv = stokes2Jones(S)

        # compound vector
        return Ex * Jv[1] + Ey * Jv[2]

    else
        # custom polarization for the TE and TM modes
        # can not be introduced
        return eFiber(mode, R, a, n1, n2, k0, P)
    end
end


function ePrimeFiberCustomPol(S, mode, R, a,
                              n1, n2, k0, P, alpha)
    #=
        Derivative with respect to alpha
        Mode fields of a stepindex cylindrical fiber
        with custom polarization defined by an arbirtary Stokes vector S.
        Custom polarization only for the hybrid modes

        ...
        # Arguments
        - `S = [S0, S1, S2, S3]` : Stokes vector.
        - `mode::customStruct` : see def above.
        - `R::array` : R = (r, phi, z).
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `l::Int`: Azymuthal mode number.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.
        - `alpha::Int`: 1 - rho, 2 - phi, 3 - z.


        # Returns
        - '[er, ephi, ez]': complex vector-column
        ...
    =#
    Eplus  = [0.0, 0.0, 0.0im]
    Eminus = [0.0, 0.0, 0.0im]
    Ex     = [0.0, 0.0, 0.0im]
    Ey     = [0.0, 0.0, 0.0im]
    Jv     = [0.0, 0.0, 0.0im]

    # \nabla_\rho
    if alpha == 1
        if ((mode.modeType == "HE") || (mode.modeType == "EH"))
            # hybrid mode are circularly polarized by default
            # orthogonal basis
            modePlus1 = Mode(mode.modeType,
                mode.l,
                mode.m,
                mode.f,
                +1,  # p - phase circulation direction
                mode.beta)

            modeMinus1 = Mode(mode.modeType,
                mode.l,
                mode.m,
                mode.f,
                -1,  # p - phase circulation direction
                mode.beta)

            Eplus  = ePrimeFiber(modePlus1,  R, a, n1, n2, k0, P)
            Eminus = ePrimeFiber(modeMinus1, R, a, n1, n2, k0, P)

            Ex = 1/sqrt(2) * (Eplus + Eminus)
            Ey = -1im/sqrt(2) * (Eplus - Eminus)

            # Stokes to Jones vector
            Jv = stokes2Jones(S)

            # compound vector
            return Ex * Jv[1] + Ey * Jv[2]
        else
            # custom polarization for the TE and TM modes
            # can not be introduced
            return ePrimeFiber(mode, R, a, n1, n2, k0, P)
        end

    # \nabla_\phi
    elseif alpha == 2
        if ((mode.modeType == "HE") || (mode.modeType == "EH"))
            e = eFiberProfile(mode.modeType, R[1], a, n1, n2, k0, mode.l, mode.beta, P)

            # Stokes to Jones vector
            Jv = stokes2Jones(S)

            prefactor = sqrt(2) * exp(im * mode.f * mode.beta * R[3]) * mode.l
            co = cos(mode.l * R[2])
            si = sin(mode.l * R[2])

            dphiEx = prefactor * [-e[1] * si, e[2] * im * co, - mode.f * e[3] * si]
            dphiEy = prefactor * [e[1] * co, im * e[2] * si, mode.f * e[3] * co]

            return 1/R[1] * (Jv[1] * dphiEx + Jv[2] * dphiEy)
        else
            # zeros for TE and TM modes
            return [0.0, 0.0, 0.0im]
        end

    # \nabla_z
    elseif alpha == 3
        if ((mode.modeType == "HE") || (mode.modeType == "EH"))
            # hybrid mode are circularly polarized by default
            # orthogonal basis
            modePlus2 = Mode(mode.modeType,
                mode.l,
                mode.m,
                mode.f,
                +1,  # p - phase circulation direction
                mode.beta)

            modeMinus2 = Mode(mode.modeType,
                mode.l,
                mode.m,
                mode.f,
                -1,  # p - phase circulation direction
                mode.beta)

            Eplus  = eFiber(modePlus2,  R, a, n1, n2, k0, P)
            Eminus = eFiber(modeMinus2, R, a, n1, n2, k0, P)

            Ex = 1/sqrt(2) * (Eplus + Eminus)
            Ey = -1im/sqrt(2) * (Eplus - Eminus)

            # Stokes to Jones vector
            Jv = stokes2Jones(S)

            # compound vector
            return im * mode.f * mode.beta * (Ex * Jv[1] + Ey * Jv[2])

        else
            # custom polarization for the TE and TM modes
            # can not be introduced
            return im * mode.f * mode.beta * eFiber(mode, R, a, n1, n2, k0, P)
        end
    else
        error()
        return [0.0, 0.0, 0.0im]
    end
end


# ############ MAGNETIC FIELDS ############
# field distributions looks OK

# Hybrid modes
function hHybrid(r, a, n1, n2, k0, l, beta, P)
    #=
        Mode profile for the hybrid modes of a stepindex cylindrical fiber
        from PHYSICAL REVIEW A 96, 023835 (2017)
        eq. (A10), (A12)

        ...
        # Arguments
        - `r::Float`: Radial coordinate.
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `l::Int`: Azymuthal mode number.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.


        # Returns
        - 'hrrr, hphi, hzzz': complex vector-column
        ...
    =#

    h = sqrt((n1*k0)^2 - beta^2)
    q = sqrt(beta^2 - (n2*k0)^2)
    ha = h * a
    qa = q * a

    Jplha = .5 * (besselj(l - 1, ha) - besselj(l + 1, ha))
    Jlha = besselj(l, ha)
    Kplqa = - .5 * (besselk(l - 1, qa) + besselk(l + 1, qa))
    Klqa = besselk(l, qa)
    s = l * (1/ha^2 + 1/qa^2) / (Jplha / (ha * Jlha) + Kplqa / (qa * Klqa))
    s1 = beta^2 * s / (k0 * n1)^2
    s2 = beta^2 * s / (k0 * n2)^2
    Jlha_Klqa = Jlha/Klqa

    inBraket = ((1 - s) * (1 - s1) * (besselj(l-1, ha)^2 - besselj(l-2, ha) * Jlha) +
                         (1 + s) * (1 + s1) * (besselj(l+1, ha)^2 - besselj(l+2, ha) * Jlha))
    outBraket = ((1 - s) * (1 - s2) * (besselk(l-2, qa) * Klqa - besselk(l-1, qa)^2) +
                          (1 + s) * (1 + s2) * (besselk(l+2, qa) * Klqa - besselk(l+1, qa)^2))
    braketH = inBraket * (n1/h)^2 + outBraket * (Jlha_Klqa * n2/q)^2

    Ah = 2sqrt(vacuumImpedance * P / (pi * a^2 * k0 * beta * braketH))


    hrrr = 0.0im
    hphi = 0.0im
    hzzz = 0.0im
    if r < a
        hr = h * r
        Jl_1hr = besselj(l - 1, hr)
        Jlplus1hr = besselj(l + 1, hr)

        hrrr = Ah * k0 * n1^2 / (2h * vacuumImpedance) * ((1 - s1) * Jl_1hr + (1 + s1) * Jlplus1hr)
        hphi = im * Ah * k0 * n1^2 / (2h * vacuumImpedance) * ((1 - s1) * Jl_1hr - (1 + s1) * Jlplus1hr)
        hzzz = im * Ah * beta * s / (k0 * vacuumImpedance) * besselj(l, hr)
    else
        qr = q * r
        Kl_1qr = besselk(l - 1, qr)
        Klplus1qr = besselk(l + 1, qr)

        hrrr = Ah * k0 * n2^2 / (2q * vacuumImpedance) * Jlha_Klqa * ((1 - s2) * Kl_1qr - (1 + s2) * Klplus1qr)
        hphi = im * Ah * k0 * n2^2 / (2q * vacuumImpedance) * Jlha_Klqa * ((1 - s2) * Kl_1qr + (1 + s2) * Klplus1qr)
        hzzz = im * Ah * beta * s / (k0 * vacuumImpedance) * Jlha_Klqa * besselk(l, qr)
    end

    return [hrrr, hphi, hzzz]
end


# Transverse modes
function hTE(r, a, n1, n2, k0, beta, P)
    #=
        Mode profile for the TE modes with l = 0 of a stepindex cylindrical fiber

        ...
        # Arguments
        - `r::Float`: Radial coordinate.
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.


        # Returns
        - '[hrrr, 0, hzzz]': complex vector-column
        ...
    =#

    h = sqrt((n1*k0)^2 - beta^2)
    q = sqrt(beta^2 - (n2*k0)^2)
    ha = h * a
    qa = q * a
    J0ha = besselj0(ha)
    K0qa = besselk(0, qa)

    inBraket = besselj(1, ha)^2 - J0ha * besselj(2, ha)
    outBraket = K0qa * besselk(2, qa) - besselk(1, qa)^2
    braketTE = inBraket / h^2 + outBraket * (J0ha/K0qa/q)^2

    ATE = sqrt(2k0 / beta * vacuumImpedance * P / (pi * a^2)) / sqrt(braketTE) / (k0 * vacuumImpedance)

    hrrr = 0.0im
    hzzz = 0.0im
    if r < a
        hr = h * r
        hrrr = - im * beta / h * ATE * besselj1(hr)
        hzzz = ATE * besselj0(hr)
    else
        qr = q * r
        hrrr = im * beta / q * ATE * J0ha/K0qa * besselk(1, qr)
        hzzz = ATE * J0ha/K0qa * besselk(0, qr)
    end
    return [hrrr, 0.0, hzzz]
end


function hTM(r, a, n1, n2, k0, beta, P)
    #=
        Mode profile for the TM modes with l = 0 of a stepindex cylindrical fiber

        ...
        # Arguments
        - `r::Float`: Radial coordinate.
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.


        # Returns
        - '[0, hphi, 0]': complex vector-column
        ...
    =#

    h = sqrt((n1*k0)^2 - beta^2)
    q = sqrt(beta^2 - (n2*k0)^2)
    ha = h * a
    qa = q * a
    J0ha = besselj0(ha)
    K0qa = besselk(0, qa)
    J0ha_K0qa = J0ha / K0qa

    inBraket = besselj(1, ha)^2 - J0ha * besselj(2, ha)
    outBraket = K0qa * besselk(2, qa) - besselk(1, qa)^2
    braketTM = (n1/h)^2 * inBraket + (n2/q * J0ha_K0qa)^2 * outBraket

    ATM = sqrt(2vacuumImpedance * P / (pi*a^2 * beta * k0)) / sqrt(braketTM)

    hphi = 0.0im
    if r < a
        hphi = - im * k0/vacuumImpedance * n1^2/h * ATM * besselj1(h*r)
    else
        hphi = im * k0/vacuumImpedance * n2^2/q * J0ha_K0qa * ATM * besselk(1, q*r)
    end
    return [0, hphi, 0]
end


function hFiber(mode, R, a, n1, n2, k0, P)
    #=
        Mode fields of a stepindex cylindrical fiber

        ...
        # Arguments
        - `mode::customStruct` : see def above.
        - `R::array` : R = (r, phi, z).
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `l::Int`: Azymuthal mode number.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.


        # Returns
        - '[hrrr, hphi, hzzz]': complex vector-column
        ...
    =#
    hProfile = [0.0, 0.0, 0.0im]

    if ((mode.modeType == "HE") || (mode.modeType == "EH"))
        hProfile = hHybrid(R[1], a, n1, n2, k0, mode.l, mode.beta, P)

        return([mode.f * mode.p * hProfile[1],
                mode.f * hProfile[2],
                mode.p * hProfile[3]] *
                exp(im * mode.f * mode.beta * R[3] + im * mode.p * mode.l * R[2]))
    elseif mode.modeType == "TE"
        # for TE l = 0
        hProfile = hTE(R[1], a, n1, n2, k0, mode.beta, P)

        return([mode.f * mode.p * hProfile[1],
                mode.f * hProfile[2],
                mode.p * hProfile[3]] *
                exp(im * mode.f * mode.beta * R[3]))
    elseif mode.modeType == "TM"
        # for TM l = 0
        hProfile = hTM(R[1], a, n1, n2, k0, mode.beta, P)

        return([mode.f * mode.p * hProfile[1],
                mode.f * hProfile[2],
                mode.p * hProfile[3]] *
                exp(im * mode.f * mode.beta * R[3]))
    else
        println("Error! Write the proper name of mode")
        error()
        return [0.0, 0.0, 0.0im]
    end
end


function hFiberCustomPol(S, mode, R,
                         a, n1, n2, k0, P)
    #=
        Mode fields of a stepindex cylindrical fiber
        with custom polarization defined by an arbirtary Stokes vector S.
        Custom polarization only for the hybrid modes

        ...
        # Arguments
        - `S = [S0, S1, S2, S3]` : Stokes vector.
        - `mode::customStruct` : see def above.
        - `R::array` : R = (r, phi, z).
        - `a::Float`: Fiber radius.
        - `n1::Float`: Fiber epsilon.
        - `n2::Float`: Outside epsilon.
        - `k0::Float`: Vacuum wave number k0 = omega/c.
        - `l::Int`: Azymuthal mode number.
        - `beta::Float`: Propogation constant.
        - `P::Float`: Total mode power.


        # Returns
        - '[hrrr, hphi, hzzz]': complex vector-column
        ...
    =#
    if ((mode.modeType == "HE") || (mode.modeType == "EH"))
        # hybrid mode are circularly polarized by default
        # orthogonal basis
        modePlus = Mode(mode.modeType,
            mode.l,
            mode.m,
            mode.f,
            +1,  # p - phase circulation direction
            mode.beta)

        modeMinus = Mode(mode.modeType,
            mode.l,
            mode.m,
            mode.f,
            -1,  # p - phase circulation direction
            mode.beta)

        Hplus  = hFiber(modePlus,  R, a, n1, n2, k0, P)
        Hminus = hFiber(modeMinus, R, a, n1, n2, k0, P)

        Hx = 1/sqrt(2) * (Hplus + Hminus)
        Hy = -1im/sqrt(2) * (Hplus - Hminus)

        # Stokes to Jones vector
        Jv = stokes2Jones(S)

        # return compound vector
        return Hx * Jv[1] + Hy * Jv[2]

    else
        # custom polarization for the TE and TM modes
        # can not be introduced
        return hFiber(mode, R, a, n1, n2, k0, P)
    end
end


# ####### Faster versions of field profiles #######

# set up parameters
function setUpPrm4Fast(mode, a, n1, n2, k0, P)
    #=
        calculates all the prefactors in the fiber profile functions

        to check see phys rev A 96, 023835 (2017), Appendix A

        # Argumnets
        - mode -- custom type wich contains all the informaiton about the mode
                        see definition above
        - a -- Fiber radius
        - n1 -- Fiber epsilon
        - n2 -- Outside epsilon
        - k0 -- vacuum wave number, k0 = \omega / c
        - P -- total mode power

        # Returns
        - Ah, ATE, ATM, Bh, BTE, BTM  -- all sets of the coefficents
                                         non-zeros for the correct type of mode
    =#
    Ah  = coefHybridElectric(0., 0., 0., 0., 0.,
                             0., 0., 0., 0., 0.)
    ATE = coefTEElectric(0., 0.)
    ATM = coefTMElectric(0., 0., 0., 0.)

    Bh  = coefHybridMagnetic(0., 0., 0., 0., 0.,
                             0., 0., 0., 0., 0.)
    BTE = coefTEMagnetic(0., 0., 0., 0.)
    BTM = coefTMMagnetic(0., 0.)



    beta = mode.beta
    l = mode.l

    h = sqrt((n1*k0)^2 - beta^2)
    q = sqrt(beta^2 - (n2*k0)^2)
    ha = h * a
    qa = q * a

    # find the A const
    if ((mode.modeType == "HE") || (mode.modeType == "EH"))
        # power to A-const
        A = P2A(mode.modeType, a, n1, n2, k0, mode.l, mode.beta, P)

        # support consts
        Jplha = .5 * (besselj(l - 1, ha) - besselj(l + 1, ha))
        Jlha = besselj(l, ha)
        Kplqa = - .5 * (besselk(l - 1, qa) + besselk(l + 1, qa))
        Klqa = besselk(l, qa)
        s = l * (1/ha^2 + 1/qa^2) / (Jplha / (ha * Jlha) + Kplqa / (qa * Klqa))
        s1 = beta^2 * s / (k0 * n1)^2
        s2 = beta^2 * s / (k0 * n2)^2

        Jlha_Klqa = Jlha/Klqa

        # electric coef
        Ahless_r1 =  im * A * beta / 2h * (1 - s)
        Ahless_r2 = -im * A * beta / 2h * (1 + s)
        Ahless_p1 = -A * beta / 2h * (1 - s)
        Ahless_p2 = -A * beta / 2h * (1 + s)
        Ahless_z1 = A

        Ahmore_r1 = im * A * beta / 2q * Jlha_Klqa * (1 - s)
        Ahmore_r2 = im * A * beta / 2q * Jlha_Klqa * (1 + s)
        Ahmore_p1 = -A * beta / 2q * Jlha_Klqa * (1 - s)
        Ahmore_p2 = A * beta / 2q * Jlha_Klqa * (1 + s)
        Ahmore_z1 = A * Jlha_Klqa

        # magnetic coef
        Bhless_r1 = A * k0/vacuumImpedance * n1^2 / 2h * (1 - s1)
        Bhless_r2 = A * k0/vacuumImpedance * n1^2 / 2h * (1 + s1)
        Bhless_p1 = im * A * k0/vacuumImpedance * n1^2 / 2h * (1 - s1)
        Bhless_p2 = -im * A * k0/vacuumImpedance * n1^2 / 2h * (1 + s1)
        Bhless_z1 = im * A * beta * s / (k0 * vacuumImpedance)

        Bhmore_r1 =  A * k0/vacuumImpedance * n2^2 / 2q * Jlha_Klqa * (1 - s2)
        Bhmore_r2 = -A * k0/vacuumImpedance * n2^2 / 2q * Jlha_Klqa * (1 + s2)
        Bhmore_p1 = im * A * k0/vacuumImpedance * n2^2 / 2q * Jlha_Klqa * (1 - s2)
        Bhmore_p2 = im * A * k0/vacuumImpedance * n2^2 / 2q * Jlha_Klqa * (1 + s2)
        Bhmore_z1 = im * A * beta * s / (k0 * vacuumImpedance) * Jlha_Klqa

        # set hybrid coefficents
        Ah  = coefHybridElectric(Ahless_r1, Ahless_r2, Ahless_p1, Ahless_p2, Ahless_z1,
                                 Ahmore_r1, Ahmore_r2, Ahmore_p1, Ahmore_p2, Ahmore_z1)
        Bh  = coefHybridMagnetic(Bhless_r1, Bhless_r2, Bhless_p1, Bhless_p2, Bhless_z1,
                                 Bhmore_r1, Bhmore_r2, Bhmore_p1, Bhmore_p2, Bhmore_z1)
        # set other const to zeros
        ATE = coefTEElectric(0., 0.)
        ATM = coefTMElectric(0., 0., 0., 0.)
        BTE = coefTEMagnetic(0., 0., 0., 0.)
        BTM = coefTMMagnetic(0., 0.)
    elseif mode.modeType == "TE"
        # power to A-const
        A = P2A(mode.modeType, a, n1, n2, k0, mode.l, mode.beta, P)

        # support consts
        J0ha_K0qa = besselj(0, ha) / besselk(0, qa)

        # electric coef
        ATEless_p1 =   im * k0 * vacuumImpedance / h * A

        ATEmore_p1 = - im * k0 * vacuumImpedance / q * J0ha_K0qa * A

        # magnetic coef
        BTEless_r1 = - im * beta / h * A
        BTEless_z1 = A

        BTEmore_r1 = im * beta / q * J0ha_K0qa * A
        BTEmore_z1 = J0ha_K0qa * A

        # set TE coefficents
        ATE = coefTEElectric(ATEless_p1, ATEmore_p1)
        BTE = coefTEMagnetic(BTEless_r1, BTEless_z1, BTEmore_r1, BTEmore_z1)
        # set other const to zeros
        Ah  = coefHybridElectric(0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)
        ATM = coefTMElectric(0., 0., 0., 0.)
        Bh  = coefHybridMagnetic(0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)
        BTM = coefTMMagnetic(0., 0.)
    elseif mode.modeType == "TM"
        # power to A-const
        A = P2A(mode.modeType, a, n1, n2, k0, mode.l, mode.beta, P)

        # support consts
        J0ha_K0qa = besselj(0, ha) / besselk(0, qa)

        # electric coef
        ATMless_r1 = - im * beta / h * A
        ATMless_z1 = A

        ATMmore_r1 = im * beta / q * J0ha_K0qa * A
        ATMmore_z1 = J0ha_K0qa * A

        # magnetic coef
        BTMless_p1 = - im * k0/vacuumImpedance * n1^2/h * A

        BTMmore_p1 = im * k0/vacuumImpedance * n2^2/q * J0ha_K0qa * A

        # set TM coefficents
        ATM = coefTMElectric(ATMless_r1, ATMless_z1, ATMmore_r1, ATMmore_z1)
        BTM = coefTMMagnetic(BTMless_p1, BTMmore_p1)
        # set other const to zeros
        Ah  = coefHybridElectric(0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)
        ATE = coefTEElectric(0., 0.)
        Bh  = coefHybridMagnetic(0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)
        BTE = coefTEMagnetic(0., 0., 0., 0.)
    else
        println("ERROR: Incorrect mode name.")
        error()
    end

    return Ah, ATE, ATM, Bh, BTE, BTM, q, h
end

# ############ ELECTRIC FIELDS ############
# faster versions

# Hybrid modes
# (A8) - (A12)
function eHybridFast(r, a, q, h, l, Ah)
    #=
        Mode profile for the hybrid modes of a stepindex cylindrical fiber

        ...
        # Arguments
        - `r::Float`: Radial coordinate.
        - 'a' -- fiber radius

        # Returns
        - 'er, ephi, ez': complex vector-column
        ...
    =#

    er = 0.0im
    ephi = 0.0im
    ez = 0.0im
    if r < a
        hr = h * r
        Jl_1hr = besselj(l - 1, hr)
        Jlplus1hr = besselj(l + 1, hr)

        er   = Ah.Ahless_r1 * Jl_1hr + Ah.Ahless_r2 * Jlplus1hr
        ephi = Ah.Ahless_p1 * Jl_1hr + Ah.Ahless_p2 * Jlplus1hr
        ez   = Ah.Ahless_z1 * besselj(l, hr)
    else
        qr = q * r
        Kl_1qr = besselk(l - 1, qr)
        Klplus1qr = besselk(l + 1, qr)

        er   = Ah.Ahmore_r1 * Kl_1qr + Ah.Ahmore_r2 * Klplus1qr
        ephi = Ah.Ahmore_p1 * Kl_1qr + Ah.Ahmore_p2 * Klplus1qr
        ez   = Ah.Ahmore_z1 * besselk(l, qr)
    end

    return [er, ephi, ez]
end


# Transverse modes
function eTEFast(r, a, q, h, ATE)
    #=
        Mode profile for the TE modes with l = 0 of a stepindex cylindrical fiber

        ...
        # Arguments
        - `r::Float`: Radial coordinate.

        # Returns
        - 'er, ephi, ez': complex vector-column
        ...
    =#
    ephi = 0.0im
    if r < a
        ephi = ATE.ATEless_p1 * besselj1(h * r)
    else
        ephi = ATE.ATEmore_p1 * besselk(1, q*r)
    end
    return [0, ephi, 0]
end


function eTMFast(r, a, q, h, ATM)
    #=
        Mode profile for the TM modes with l = 0 of a stepindex cylindrical fiber

        ...
        # Arguments
        - `r::Float`: Radial coordinate.

        # Returns
        - 'er, ephi, ez': complex vector-column
        ...
    =#


    er = 0.0im
    ez = 0.0im
    if r < a
        hr = h * r

        er = ATM.ATMless_r1 * besselj1(hr)
        ez = ATM.ATMless_z1 * besselj0(hr)
    else
        qr = q * r

        er = ATM.ATMmore_r1 * besselk(1, qr)
        ez = ATM.ATMmore_z1 * besselk(0, qr)
    end
    return [er, 0, ez]
end


function eFiberFast(mode, R, a, q, h,
                Ah, ATE, ATM)
    #=
        Mode fields of a stepindex cylindrical fiber

        ...
        # Arguments
        - `mode::customStruct` : see def above.
        - `R::array` : R = (r, phi, z).

        # Returns
        - '[er, ephi, ez]': complex vector-column
        ...
    =#
    eProfile = [0.0, 0.0, 0.0im]
    if ((mode.modeType == "HE") || (mode.modeType == "EH"))
        eProfile = eHybridFast(R[1], a, q, h, mode.l, Ah)

        return([eProfile[1],
                mode.p * eProfile[2],
                mode.f * eProfile[3]] *
                exp(im * mode.f * mode.beta * R[3] + im * mode.p * mode.l * R[2]))
    elseif mode.modeType == "TE"
        # for TE l = 0
        eProfile = eTEFast(R[1], a, q, h, ATE)

        return([eProfile[1],
                mode.p * eProfile[2],
                mode.f * eProfile[3]] *
                exp(im * mode.f * mode.beta * R[3]))
    elseif mode.modeType == "TM"
        # for TM  l =0
        eProfile = eTMFast(R[1], a, q, h, ATM)

        return([eProfile[1],
                mode.p * eProfile[2],
                mode.f * eProfile[3]] *
                exp(im * mode.f * mode.beta * R[3]))
    else
        println("Error! Write the proper name of mode")
        error()
        return [0.0, 0.0, 0.0im]
    end
end


function eFiberCustomPolFast(S, mode, R, a, q, h,
                             Ah, ATE, ATM)
    #=
        Mode fields of a stepindex cylindrical fiber
        with custom polarization defined by an arbirtary Stokes vector S.
        Custom polarization only for the hybrid modes

        ...
        # Arguments
        - `S = [S0, S1, S2, S3]` : Stokes vector.
        - `mode::customStruct` : see def above.
        - `R::array` : R = (r, phi, z).

        # Returns
        - '[er, ephi, ez]': complex vector-column
        ...
    =#
    if ((mode.modeType == "HE") || (mode.modeType == "EH"))
        # hybrid mode are circularly polarized by default
        # orthogonal basis
        modePlus = Mode(mode.modeType,
            mode.l,
            mode.m,
            mode.f,
            +1,  # p - phase circulation direction CCW
            mode.beta)

        modeMinus = Mode(mode.modeType,
            mode.l,
            mode.m,
            mode.f,
            -1,  # p - phase circulation direction CW
            mode.beta)

        Eplus  = eFiberFast(modePlus,  R, a, q, h, Ah, ATE, ATM)
        Eminus = eFiberFast(modeMinus, R, a, q, h, Ah, ATE, ATM)

        Ex = 1/sqrt(2) * (Eplus + Eminus)
        Ey = -1im/sqrt(2) * (Eplus - Eminus)

        # Stokes to Jones vector
        Jv = stokes2Jones(S)

        # compound vector
        return Ex * Jv[1] + Ey * Jv[2]

    else
        # custom polarization for the TE and TM modes
        # can not be introduced
        return eFiberFast(mode, R, a, q, h, Ah, ATE, ATM)
    end
end



# ############ MAGNETIC FIELDS ############
# faster versions

# Hybrid modes
function hHybridFast(r, a, q, h, l, Bh)
    #=
        Mode profile for the hybrid modes of a stepindex cylindrical fiber
        from PHYSICAL REVIEW A 96, 023835 (2017)
        eq. (A10), (A12)

        ...
        # Arguments
        - `r::Float`: Radial coordinate.

        # Returns
        - 'hrrr, hphi, hzzz': complex vector-column
        ...
    =#
    hrrr = 0.0im
    hphi = 0.0im
    hzzz = 0.0im
    if r < a
        hr = h * r
        Jl_1hr = besselj(l - 1, hr)
        Jlplus1hr = besselj(l + 1, hr)

        hrrr = Bh.Bhless_r1 * Jl_1hr + Bh.Bhless_r2 * Jlplus1hr
        hphi = Bh.Bhless_p1 * Jl_1hr + Bh.Bhless_p2 * Jlplus1hr
        hzzz = Bh.Bhless_z1 * besselj(l, hr)
    else
        qr = q * r
        Kl_1qr = besselk(l - 1, qr)
        Klplus1qr = besselk(l + 1, qr)

        hrrr = Bh.Bhmore_r1 * Kl_1qr + Bh.Bhmore_r2 * Klplus1qr
        hphi = Bh.Bhmore_p1 * Kl_1qr + Bh.Bhmore_p2 * Klplus1qr
        hzzz = Bh.Bhmore_z1 * besselk(l, qr)
    end

    return [hrrr, hphi, hzzz]
end


# Transverse modes
function hTEFast(r, a, q, h, BTE)
    #=
        Mode profile for the TE modes with l = 0 of a stepindex cylindrical fiber

        ...
        # Arguments
        - `r::Float`: Radial coordinate.

        # Returns
        - '[hrrr, 0, hzzz]': complex vector-column
        ...
    =#

    hrrr = 0.0im
    hzzz = 0.0im
    if r < a
        hr = h * r
        hrrr = BTE.BTEless_r1 * besselj1(hr)
        hzzz = BTE.BTEless_z1 * besselj0(hr)
    else
        qr = q * r
        hrrr = BTE.BTEmore_r1 * besselk(1, qr)
        hzzz = BTE.BTEmore_z1 * besselk(0, qr)
    end
    return [hrrr, 0.0, hzzz]
end


function hTMFast(r, a, q, h, BTM)
    #=
        Mode profile for the TM modes with l = 0 of a stepindex cylindrical fiber

        ...
        # Arguments
        - `r::Float`: Radial coordinate.

        # Returns
        - '[0, hphi, 0]': complex vector-column
        ...
    =#

    hphi = 0.0im

    if r < a
        hphi = BTM.BTMless_p1 * besselj1(h*r)
    else
        hphi = BTM.BTMmore_p1 * besselk(1, q*r)
    end

    return [0, hphi, 0]
end


function hFiberFast(mode, R, a, q, h,
                    Bh, BTE, BTM)
    #=
        Mode fields of a stepindex cylindrical fiber

        ...
        # Arguments
        - `mode::customStruct` : see def above.
        - `R::array` : R = (r, phi, z).
        - `a::Float`: Fiber radius.

        # Returns
        - '[hrrr, hphi, hzzz]': complex vector-column
        ...
    =#
    hProfile = [0.0, 0.0, 0.0im]

    if ((mode.modeType == "HE") || (mode.modeType == "EH"))
        hProfile = hHybridFast(R[1], a, q, h, mode.l, Bh)

        return([mode.f * mode.p * hProfile[1],
                mode.f * hProfile[2],
                mode.p * hProfile[3]] *
                exp(im * mode.f * mode.beta * R[3] + im * mode.p * mode.l * R[2]))
    elseif mode.modeType == "TE"
        # for TE l = 0
        hProfile = hTEFast(R[1], a, q, h, BTE)

        return([mode.f * mode.p * hProfile[1],
                mode.f * hProfile[2],
                mode.p * hProfile[3]] *
                exp(im * mode.f * mode.beta * R[3]))
    elseif mode.modeType == "TM"
        # for TM l = 0
        hProfile = hTMFast(R[1], a, q, h, BTM)

        return([mode.f * mode.p * hProfile[1],
                mode.f * hProfile[2],
                mode.p * hProfile[3]] *
                exp(im * mode.f * mode.beta * R[3]))
    else
        println("Error! Write the proper name of mode")
        error()
        return [0.0, 0.0, 0.0im]
    end
end


function hFiberCustomPolFast(S, mode, R,
                             a, q, h,
                             Bh, BTE, BTM)
    #=
        Mode fields of a stepindex cylindrical fiber
        with custom polarization defined by an arbirtary Stokes vector S.
        Custom polarization only for the hybrid modes

        ...
        # Arguments
        - `S = [S0, S1, S2, S3]` : Stokes vector.
        - `mode::customStruct` : see def above.
        - `R::array` : R = (r, phi, z).

        # Returns
        - '[hrrr, hphi, hzzz]': complex vector-column
        ...
    =#
    if ((mode.modeType == "HE") || (mode.modeType == "EH"))
        # hybrid mode are circularly polarized by default
        # orthogonal basis
        modePlus = Mode(mode.modeType,
            mode.l,
            mode.m,
            mode.f,
            +1,  # p - phase circulation direction
            mode.beta)

        modeMinus = Mode(mode.modeType,
            mode.l,
            mode.m,
            mode.f,
            -1,  # p - phase circulation direction
            mode.beta)

        Hplus  = hFiberFast(modePlus,  R, a, q, h, Bh, BTE, BTM)
        Hminus = hFiberFast(modeMinus, R, a, q, h, Bh, BTE, BTM)

        Hx = 1/sqrt(2) * (Hplus + Hminus)
        Hy = -1im/sqrt(2) * (Hplus - Hminus)

        # Stokes to Jones vector
        Jv = stokes2Jones(S)

        # return compound vector
        return Hx * Jv[1] + Hy * Jv[2]

    else
        # custom polarization for the TE and TM modes
        # can not be introduced
        return hFiberFast(mode, R, a, q, h, Bh, BTE, BTM)
    end
end
