# by Ivan Toftul
# 2018/07

using SpecialFunctions
using QuadGK
using Roots

export 
    GsphiphiWG,
    GsWGn11,
    iGs11ij,
    Gs11ijpolar,
    Gs11Polar,
    Gs11ijpolarRadModes,
    Gs11PolarRadModes
    

function delta(i, j)
    if i == j
        1
    else
        0
    end
end


function cart2pol(x, y)
    rho = sqrt(x^2 + y^2)
    phi = atan(y, x)
    # return
    rho, phi
end


# ## Analytical expressions for guided modes

# GF wg
function GsphiphiWG(vecRr, vecRs, k, beta, epsf, epsm, Rf, n)
    #=
        G = G(Rr, Rs, omega)
        Calculates Green's function of a fiber (phi phi component) only for waveguided modes
        No scattering! Only pure eigen modes!
        There is no numeric integration below, the integral is taken using complex analysis.

        Notation of subscripts:
            '1' : outside fiber
            '2' : inside fiber

        Arguments:
            vecRr = (rr, pr, zr) : position of a reciever
            vecRs = (rs, ps, zs) : position of a source
            k = omega/c : vacuum wavevector
            beta : propagation constant for n mode
            n : mode number
            epsf, epfm : permitivities of fiber and surrouding medium
            Rf : fiber radius

        Returns:
            G : complex value of GF
    =#
    dtheta = vecRr[2] - vecRs[2] 
    dz = abs(vecRr[3] - vecRs[3])

    k1 = sqrt(epsm) * k
    k2 = sqrt(epsf) * k
    krho1 = sqrt(k1*k1 - beta*beta + 0im)
    krho2 = sqrt(k2*k2 - beta*beta + 0im)

    # calc for the Rxx_num coefficients
    k1rf = krho1 * Rf
    k2rf = krho2 * Rf

    Jn1 = besselj(n, k1rf)
    Hn1 = hankelh1(n, k1rf)

    J1 = .5 * (besselj(n - 1, k1rf) - besselj(n + 1, k1rf)) / (krho1 * Jn1)
    J2 = .5 * (besselj(n - 1, k2rf) - besselj(n + 1, k2rf)) / (krho2 * besselj(n, k2rf))
    H1 = .5 * (hankelh1(n - 1, k1rf) - hankelh1(n + 1, k1rf)) / (krho1 * Hn1)

    braket = 1/krho2^2 - 1/krho1^2 

    CMM = Jn1 / Hn1 * ((braket * beta * n)^2 - 
                               Rf^2 * (J2 - J1) * (J2 * k2^2 - H1 * k1^2))
    CNM = - Jn1 / Hn1 * braket * (J1 - H1) * k1 * beta * n * Rf 
    # CMN = CNM
    CNN = Jn1 / Hn1 * ((braket*beta*n)^2 - (J2- H1) * (J2 * k2^2 - J1 * k1^2) * Rf^2) 

    # calc for the F term
    rho = vecRr[1]
    rhop = vecRs[1]
    k1r = krho1 * rho
    k1rp = krho1 * rhop

    Hnk1r = hankelh1(n, k1r)
    Hnk1rp = hankelh1(n, k1rp)
    Hnpk1r = .5 * (hankelh1(n - 1, k1r) - hankelh1(n + 1, k1r))
    Hnpk1rp = .5 * (hankelh1(n - 1, k1rp) - hankelh1(n + 1, k1rp))

    f = (CMM * Hnpk1r * Hnpk1rp +
         CNM * n * beta / (krho1 * k1) * (Hnk1r * Hnpk1rp / rho +
                                          Hnk1rp * Hnpk1r / rhop) +
         CNN * (n*n * beta*beta /(krho1*krho1 * k1*k1 * rho * rhop) * Hnk1r * Hnk1rp))
    f *= (2 - delta(0, n)) * cos(n * dtheta)

    function fooprime(func, j)
        q = 0.
        krho = 0.
        if j == 1
            q = krho1 * Rf
            krho = krho1
        elseif j == 2
            q = krho2 * Rf
            krho = krho2
        else
            error()
        end
        # C = J or H^(1)
        C = func(n, q)
        Cp = .5 * (func(n - 1, q) - func(n + 1, q))
        Cpp = .25 * (func(n - 2, q) - 2 * func(n, q) + func(n + 2, q))

        ans = beta * Rf / (krho^2 * C) * (Cp^2 / C - Cpp) + beta / (krho^3) * Cp / C
        return(ans)
    end

    A = - 4 * beta^3 * n^2 * braket * (1/krho2^4 - 1/krho1^4) - 2 * beta * n^2 * braket^2
    J2p = fooprime(besselj, 2)
    H1p = fooprime(hankelh1, 1)
    
    Dp = A + ((J2p - H1p) * (J2 * k2^2 - H1 * k1^2) + 
              (J2 - H1) * (J2p * k2^2 - H1p * k1^2)) * Rf^2 

    G = - .25 * f/Dp * exp(1im * beta * dz)
end


# total tensor
function GsWGn11(vecRr, vecRs, k, beta, epsf, epsm, Rf, n)
    #=
        G = G(Rr, Rs, omega)
        Calculates Green's tensor of a fiber only for waveguided modes for one mode n
        No scattering! Only pure eigen modes!
        There is no numeric integration below, the integral is taken using complex analysis.

        Notation of subscripts:
            '1' : outside fiber
            '2' : inside fiber

        Arguments:
            vecRr = (rr, pr, zr) : position of a reciever
            vecRs = (rs, ps, zs) : position of a source
            k = omega/c : vacuum wavevector
            beta : propagation constant for n mode
            n : azimuthal number
            epsf, epfm : permitivities of fiber and surrouding medium
            Rf : fiber radius

        Returns:
            G : complex tensor of GF
    =#
    dtheta = vecRr[2] - vecRs[2] 
    dz = abs(vecRr[3] - vecRs[3])

    k1 = sqrt(epsm) * k
    k2 = sqrt(epsf) * k
    krho1 = sqrt(k1*k1 - beta*beta + 0im)
    krho2 = sqrt(k2*k2 - beta*beta + 0im)
    
    # calc for the Rxx_num coefficients
    k1rf = krho1 * Rf
    k2rf = krho2 * Rf

    Jn1 = besselj(n, k1rf)
    Hn1 = hankelh1(n, k1rf)

    J1 = .5 * (besselj(n - 1, k1rf) - besselj(n + 1, k1rf)) / (krho1 * Jn1)
    J2 = .5 * (besselj(n - 1, k2rf) - besselj(n + 1, k2rf)) / (krho2 * besselj(n, k2rf))
    H1 = .5 * (hankelh1(n - 1, k1rf) - hankelh1(n + 1, k1rf)) / (krho1 * Hn1)

    braket = 1/krho2^2 - 1/krho1^2 

    CMM = Jn1 / Hn1 * ((braket * beta * n)^2 - 
                               Rf^2 * (J2 - J1) * (J2 * k2^2 - H1 * k1^2))
    CNM = - Jn1 / Hn1 * braket * (J1 - H1) * k1 * beta * n * Rf 
    CMN = CNM
    CNN = Jn1 / Hn1 * ((braket*beta*n)^2 - (J2- H1) * (J2 * k2^2 - J1 * k1^2) * Rf^2) 
    
    # tensor products of the vector harmonics from Danny's manuscript
    # MM = (M \otimes M) / krho1^2 without exp()
    
    rho = vecRr[1]
    rhop = vecRs[1]
    k1r = krho1 * rho
    k1rp = krho1 * rhop

    Hnk1r = hankelh1(n, k1r)
    Hnk1rp = hankelh1(n, k1rp)
    Hnpk1r = .5 * (hankelh1(n - 1, k1r) - hankelh1(n + 1, k1r))
    Hnpk1rp = .5 * (hankelh1(n - 1, k1rp) - hankelh1(n + 1, k1rp))
    
    MM11 = n^2 * Hnk1r * Hnk1rp / (krho1^2 * rho * rhop)
    MM12 = - im * n * Hnk1r * Hnpk1rp / (krho1 * rho)
    MM21 = im * n / (krho1 * rhop) * Hnpk1r * Hnk1rp
    MM22 = Hnpk1r * Hnpk1rp
    MM = [MM11 MM12 0.0;
          MM21 MM22 0.0;
           0.0  0.0 0.0]
    
    MN11 = n * beta / (k1 * krho1 * rho) * Hnk1r * Hnpk1rp
    MN12 = -im * n^2 * beta / (rho * rhop * k1 * krho1^2) * Hnk1r * Hnk1rp
    MN13 = im * n / (rho * k1) * Hnk1r * Hnk1rp
    MN21 = im * beta / k1 * Hnpk1r * Hnpk1rp
    MN22 = n * beta / (rhop * k1 * krho1) * Hnpk1r * Hnk1rp
    MN23 = - krho1 / k1 * Hnpk1r * Hnk1rp
    MN = [MN11 MN12 MN13;
          MN21 MN22 MN23;
           0.0  0.0  0.0]
    
    NM11 = n * beta / (krho1 * k1 * rhop) * Hnpk1r * Hnk1rp
    NM12 = - im * beta / k1 * Hnpk1r * Hnpk1rp
    NM21 = im * n^2 * beta / (krho1^2 * k1 * rho * rhop) * Hnk1r * Hnk1rp
    NM22 = n * beta / (krho1 * k1 * rho) * Hnk1r * Hnpk1rp
    NM31 = -im * n / (k1 * rhop) * Hnk1r * Hnk1rp
    NM32 = - krho1 / k1 * Hnk1r * Hnpk1rp
    NM = [NM11 NM12  0.0;
          NM21 NM22  0.0;
          NM31 NM32  0.0]
    
    NN11 = beta^2 / k1^2 * Hnpk1r * Hnpk1rp
    NN12 = - im * n * beta^2 / (krho1 * k1^2 * rhop) * Hnpk1r * Hnk1rp
    NN13 = im * krho1 * beta / k1^2 * Hnpk1r * Hnk1rp
    NN21 = im * n * beta^2 / (krho1 * k1^2 * rho) * Hnk1r * Hnpk1rp
    NN22 = n^2 * beta^2 / (krho1^2 * rho * rhop * k1^2) * Hnk1r * Hnk1rp
    NN23 = - n * beta / (k1^2 * rho) * Hnk1r * Hnk1rp
    NN31 = - im * krho1 * beta / k1^2 * Hnk1r * Hnpk1rp
    NN32 = - n * beta / (k1^2 * rhop) * Hnk1r * Hnk1rp
    NN33 = krho1^2 / k1^2 * Hnk1r * Hnk1rp
    NN = [NN11 NN12 NN13;
          NN21 NN22 NN23;
          NN31 NN32 NN33]
    
    # calc the derivative from the dispersion expresion
    function fooprime(func, j)
        q = 0.
        krho = 0.
        if j == 1
            q = krho1 * Rf
            krho = krho1
        elseif j == 2
            q = krho2 * Rf
            krho = krho2
        else
            error()
        end
        # C = J or H^(1)
        C = func(n, q)
        Cp = .5 * (func(n - 1, q) - func(n + 1, q))
        Cpp = .25 * (func(n - 2, q) - 2 * func(n, q) + func(n + 2, q))

        ans = beta * Rf / (krho^2 * C) * (Cp^2 / C - Cpp) + beta / (krho^3) * Cp / C
        return(ans)
    end

    A = - 4 * beta^3 * n^2 * braket * (1/krho2^4 - 1/krho1^4) - 2 * beta * n^2 * braket^2
    J2p = fooprime(besselj, 2)
    H1p = fooprime(hankelh1, 1)
    
    Dp = A + ((J2p - H1p) * (J2 * k2^2 - H1 * k1^2) + 
              (J2 - H1) * (J2p * k2^2 - H1p * k1^2)) * Rf^2 
    
    # componets of Green's tensor
    G11 = - .25 * (2 - delta(0, n)) * cos(n * dtheta) * (CMM * MM[1, 1] + CNM * NM[1, 1] + CMN * MN[1, 1] + CNN * NN[1, 1])/Dp * exp(1im * beta * dz)
    G12 = - .25 * (2 - delta(0, n)) * sin(n * dtheta) * (CMM * MM[1, 2] + CNM * NM[1, 2] + CMN * MN[1, 2] + CNN * NN[1, 2])/Dp * exp(1im * beta * dz)
    G13 = - .25 * (2 - delta(0, n)) * cos(n * dtheta) * (CMM * MM[1, 3] + CNM * NM[1, 3] + CMN * MN[1, 3] + CNN * NN[1, 3])/Dp * exp(1im * beta * dz)
    G21 = - .25 * (2 - delta(0, n)) * sin(n * dtheta) * (CMM * MM[2, 1] + CNM * NM[2, 1] + CMN * MN[2, 1] + CNN * NN[2, 1])/Dp * exp(1im * beta * dz)
    G22 = - .25 * (2 - delta(0, n)) * cos(n * dtheta) * (CMM * MM[2, 2] + CNM * NM[2, 2] + CMN * MN[2, 2] + CNN * NN[2, 2])/Dp * exp(1im * beta * dz)
    G23 = - .25 * (2 - delta(0, n)) * sin(n * dtheta) * (CMM * MM[2, 3] + CNM * NM[2, 3] + CMN * MN[2, 3] + CNN * NN[2, 3])/Dp * exp(1im * beta * dz)
    G31 = - .25 * (2 - delta(0, n)) * cos(n * dtheta) * (CMM * MM[3, 1] + CNM * NM[3, 1] + CMN * MN[3, 1] + CNN * NN[3, 1])/Dp * exp(1im * beta * dz)
    G32 = - .25 * (2 - delta(0, n)) * sin(n * dtheta) * (CMM * MM[3, 2] + CNM * NM[3, 2] + CMN * MN[3, 2] + CNN * NN[3, 2])/Dp * exp(1im * beta * dz)
    G33 = - .25 * (2 - delta(0, n)) * cos(n * dtheta) * (CMM * MM[3, 3] + CNM * NM[3, 3] + CMN * MN[3, 3] + CNN * NN[3, 3])/Dp * exp(1im * beta * dz)
    
    return [G11 G12 G13;
            G21 G22 G23;
            G31 G32 G33]
end


# GF integrand for numeric calculations
function iGs11ij(x, rr, pr, zr, rs, ps, zs, k, epsf, epsm, Rf, n, i, j)
    #=
        Fiber Green's function integrand
        Calculates the integrand of the scattered part of the Green's tensor of the fiber;

    Parameters
    ----------
        x : complex
            k_z, integration variable
        k : float
            k-vector value, 2pi/\lambda_0 = \omega/c;
        eps_out, eps_in : complex
            electric permetivity outside and inside the fiber
        rc : float
            fiber radius;
        n : int
            mode order;
        rr, pr, zr : float
            reciever coordinates;
        rs, ps, zs : float
            source coordinates
        i, j : int
            rho, phi, z tensor indeces

    Returns
    -------
        iG : complex
            one component of inetegrand iGij
    =#
    iG = 0.
    
    k1 = sqrt(epsm) * k
    k2 = sqrt(epsf) * k

    a = sqrt(k1^2 - x^2 + 0im)  # krho1
    b = sqrt(k2^2 - x^2 + 0im)  # krho2
    
    arr = a * rr
    ars = a * rs
    brc = b * Rf
    arc = a * Rf

    Hn1r = hankelh1(n, arr)
    Hn1s = hankelh1(n, ars)

    DHn1r = .5 * (hankelh1(n - 1, arr) - hankelh1(n + 1, arr))
    DHn1s = .5 * (hankelh1(n - 1, ars) - hankelh1(n + 1, ars))

    DJnb = .5 * (besselj(n - 1, brc) - besselj(n + 1, brc))
    Jnb  =  besselj(n, brc)

    DJna = .5 * (besselj(n - 1, arc) - besselj(n + 1, arc))
    Jna = besselj(n, arc)

    DHna = .5 * (hankelh1(n - 1, arc) - hankelh1(n + 1, arc))
    Hna = hankelh1(n, arc)

    braket = 1/b^2 - 1/a^2  # 1/krho2^2 - 1/krho1^2
    J1 = DJna / (a * Jna)
    J2 = DJnb / (b * Jnb)
    H1 = DHna / (a * Hna)
    
    Det = Rf^2 * (k2^2 * J2 - k1^2 * H1) * (J2 - H1) - (n * x * braket)^2
    
    preTerm = Jna / Hna / Det
    # FRESNEL COEFFICIENTS
    Rn11mm = preTerm * ((braket * n * x)^2 - (J2 - J1) * (J2 * k2^2 - H1 * k1^2) * Rf^2)
    Rn11mn = k1 * x * n * Rf * preTerm * braket * (H1 - J1)
    Rn11nm = Rn11mn
    Rn11nn = preTerm * ((braket * x * n)^2 - (J2 - H1) * (J2 * k2^2 - J1 * k1^2) * Rf^2)

    preExp = (2 - delta(0, n)) * 1im * exp(1im * x * (zr-zs)) / 8pi
    cosdp = cos(n * (pr - ps)) 
    sindp = sin(n * (pr - ps)) 

    # rr component
    if ((i == 1) && (j == 1))
        iGNrr11mm = Hn1r * Hn1s * n^2 * Rn11mm / (rr * rs * a^2)
        iGNrr11nm = DHn1r * Hn1s * n * Rn11nm * x / (k1 * rs * a)
        iGNrr11mn = DHn1s * Hn1r * n * Rn11mn * x / (k1 * rr * a)
        iGNrr11nn = DHn1r * DHn1s * Rn11nn * x^2 / k1^2
        
        iG = cosdp * preExp * (iGNrr11mm + iGNrr11nm + iGNrr11mn + iGNrr11nn)
    # rp component
    elseif ((i == 1) && (j == 2))
        iGNrp11mm = n * Hn1r * DHn1s * Rn11mm / (rr * a)
        iGNrp11nm = x * DHn1r * DHn1s * Rn11nm / k1
        iGNrp11mn = Hn1r * Hn1s * n^2 * x * Rn11mn / (k1 * rr * rs * a^2)
        iGNrp11nn = DHn1r * Hn1s * n * Rn11nn * x^2 / (k1^2 * rs * a)

        iG = sindp * preExp * (iGNrp11mm + iGNrp11nm + iGNrp11mn + iGNrp11nn)
    # rz component
    elseif ((i == 1) && (j == 3))
        iGNrz11mm = 0.
        iGNrz11nm = 0.
        iGNrz11mn = 1im * Hn1r * Hn1s * n * Rn11mn / (k1 * rr)
        iGNrz11nn = 1im * a * DHn1r * Hn1s * Rn11nn * x / k1^2

        iG = cosdp * preExp * (iGNrz11mm + iGNrz11nm + iGNrz11mn + iGNrz11nn)
    # pr component
    elseif ((i == 2) && (j == 1))
        iGNpr11mm = -DHn1r * Hn1s * n * Rn11mm / (rs * a)
        iGNpr11nm = -Hn1r * Hn1s * n^2 * Rn11nm * x / (k1 * rr * rs * a^2)
        iGNpr11mn = -DHn1r * DHn1s * Rn11mn * x / k1
        iGNpr11nn = -DHn1s * Hn1r * n * Rn11nn * x^2 / (k1^2 * rr * a)

        iG = sindp * preExp * (iGNpr11mm + iGNpr11nm + iGNpr11mn + iGNpr11nn)
    # pp component
    elseif ((i == 2) && (j == 2))
        iGNpp11mm = Rn11mm * DHn1r * DHn1s
        iGNpp11nm = n * x * Rn11nm * Hn1r * DHn1s / (k1 * rr * a)
        iGNpp11mn = n * Rn11mn * DHn1r * Hn1s * x / (k1 * rs * a)
        iGNpp11nn = n^2 * x^2 * Rn11nn * Hn1r * Hn1s / (k1^2 * rr * rs * a^2)

        iG = cosdp * preExp * (iGNpp11mm + iGNpp11nm + iGNpp11mn + iGNpp11nn)
    # pz component
    elseif ((i == 2) && (j == 3))
        iGNpz11mm = 0.
        iGNpz11nm = 0.
        iGNpz11mn = -1im * a * DHn1r * Hn1s * Rn11mn / k1
        iGNpz11nn = -1im * Hn1r * Hn1s * n * Rn11nn * x / (k1^2 * rr)

        iG = sindp * preExp * (iGNpz11mm + iGNpz11nm + iGNpz11mn + iGNpz11nn)
    # zr component
    elseif ((i == 3) && (j == 1))
        iGNzr11mm = 0.
        iGNzr11nm = -1im * Hn1r * Hn1s * n * Rn11nm / (k1 * rs)
        iGNzr11mn = 0.
        iGNzr11nn = -1im * a * DHn1s * Hn1r * Rn11nn * x / k1^2

        iG = cosdp * preExp * (iGNzr11mm + iGNzr11nm + iGNzr11mn + iGNzr11nn)
    # zp component
    elseif ((i == 3) && (j == 2))
        iGNzp11mm = 0.
        iGNzp11nm = -1im * a * DHn1s * Hn1r * Rn11nm / k1
        iGNzp11mn = 0.
        iGNzp11nn = -1im * Hn1r * Hn1s * n * Rn11nn * x / (k1^2 * rs)

        iG = sindp * preExp * (iGNzp11mm + iGNzp11nm + iGNzp11mn + iGNzp11nn)
    # zz component
    elseif ((i == 3) && (j == 3))
        iG = cosdp * preExp * Rn11nn * Hn1r * Hn1s * (a / k1)^2
    end

    return iG
end


function Gs11ijpolar(rr, pr, zr, rs, ps, zs, k, epsf, epsm, Rf, nMaxModeOrder, i, j; kzmax=-1.0, tol=1e-4, maxevals=4000)
    #=
        Fiber Green's function

    Symmetry properties (double check it!)
    -------------------

                | G11 G12 G13 |
    G(r1, r2) = | G21 G22 G23 |
                | G31 G32 G33 |

                | G11 ... G31 |
    G(r2, r1) = | ... ... G32 |
                | G13 G23 G33 |

    where '...' means something different

    Parameters
    ----------
        k : float
            k-vector value, 2pi/\lambda_0 = \omega/c;
        eps_out, eps_in : complex
            electric permetivity outside and inside the fiber;
        Rf : float
            fiber radius;
        rs, rp : float
            positions of the source and the reciever,
            r = (r, p, z) -- in polar coordinates;
        nMaxModeOrder : int
            maximal mode number in the sun;
        i, j : int
            tensor indicies over cylindrical coordinates \rho, \phi, z;
        tol : float
            relative tolerance for the sum (how many modes 'n' to consider?);
            |G^{N} - G^{N-1}|/G^{N};
        kzmax : float
            upper integration limit in kz,
            usually several k_in;

    Returns
    -------
        G : complex
            one component of Gij;
        nmodes : number of considered modes
            condition for cutting the 'n' exists
    =#
    k2 = sqrt(epsf) * k
    if kzmax == -1.0
        # stadart value of an upper limit
        kzmax = k2 * 7.0
    end

    kzReMax = kzmax
    kzImMax = k * 1e-4
    
    #=
    AREAS:

                   | Im kz
             *     |
         (2) * *   |(3)
             *   * |            (5)       Re kz
    **********-----*-----**************-----
      (1)          | *   *
                   |   * * (4)
                   |     *
                   |
    =#
    
    rel = 0.0

    # way points
    A = complex(-kzmax, 0)
    B = complex(-1.1 * k2, 0)
    C = complex(-1.1 * k2, kzImMax)
    D = complex(1.1 * k2, -kzImMax)
    E = complex(1.1 * k2, 0)
    F = complex(kzmax, 0)
    
    Gsij = 0. # Green's tensor component
    nmodes = 0 # number of considered 'n' modes
    Gnprev = 0. # previous, it sums from 0 to nmax-1; when nmax = 0 Gnprev = 0

    for n = 0:nMaxModeOrder
        Gsij += quadgk(x -> iGs11ij(x, rr, pr, zr, rs, ps, zs, 
                                   k, epsf, epsm, Rf, n, i, j), A, B, C, D, E, F; 
                                   rtol=tol, maxevals=maxevals)[1]
        rel = abs(abs(Gsij) - abs(Gnprev)) / abs(Gsij)
        if rel < tol
            break
        end
        Gnprev = Gsij
        nmodes += 1
    end

    return Gsij, nmodes
end

function Gs11Polar(vecRr, vecRs, k, epsf, epsm, Rf, nmax; kzmax=-1.0, tol=1e-4, maxevals=4000)

    Gs = zeros(ComplexF64, 3, 3)
    #=
    Returns full tensor Gs in polar coordinates

    Parameters
    ----------
        \vec{r}_s, vec{r}_r : coordinates of source and reciever 
                              in polar coordinates (r, phi, z)
    =#
    rr, pr, zr = vecRr[1], vecRr[2], vecRr[3]
    rs, ps, zs = vecRs[1], vecRs[2], vecRs[3]
    if vecRr == vecRs
        # then Gs is diagonal
        #      |rr        |
        # Gs = |    pp    |
        #      |        zz|
        for i = 1:3
            Gs[i, i] = Gs11ijpolar(rr, pr, zr, 
                            rs, ps, zs, 
                            k, epsf, epsm, Rf,
                            nmax, i, i; kzmax=kzmax, tol=tol, maxevals=maxevals)[1]
        end
    elseif zr == zs
        # if zr = zs then 
        #      |rr  rp  0 |
        # Gs = |pr  pp  0 |
        #      |0   0   zz|
        for i = 1:3
            Gs[i, i] = Gs11ijpolar(rr, pr, zr, rs, ps, zs, k, epsf, epsm, Rf, nmax, i, i;
             kzmax=kzmax, tol=tol, maxevals=maxevals)[1]
        end
        Gs[1, 2] = Gs11ijpolar(rr, pr, zr, rs, ps, zs, k, epsf, epsm, Rf,nmax, 1, 2; 
        kzmax=kzmax, tol=tol, maxevals=maxevals)[1]
        Gs[2, 1] = Gs11ijpolar(rr, pr, zr, rs, ps, zs, k, epsf, epsm, Rf,nmax, 2, 1; 
        kzmax=kzmax, tol=tol, maxevals=maxevals)[1]
        
    else
        # calc all components
        for i = 1:3, j = 1:3
            Gs[i, j] = Gs11ijpolar(rr, pr, zr, rs, ps, zs, k, epsf, epsm, Rf, nmax, i, j; 
            kzmax=kzmax, tol=tol, maxevals=maxevals)[1]
        end
    end

    return Gs
end



function Gs11ijpolarRadModes(rr, pr, zr, rs, ps, zs, k, epsf, epsm, Rf, nMaxModeOrder, i, j; tol=1e-4, maxevals=4000)
    #=
        Fiber Green's function
        Semi-analytical approach

    Area of integration:

                   | Im kz
                   |
                   |
                   |                    Re kz
    --|-----|*************|------|----------
    -k2   -k1      |      k1     k2
                   |   
                   |     
                   |

    waveguided modes considered analytically

    Parameters
    ----------
        k : float
            k-vector value, 2pi/\lambda_0 = \omega/c;
        eps_out, eps_in : complex
            electric permetivity outside and inside the fiber;
        Rf : float
            fiber radius;
        rs, rp : float
            positions of the source and the reciever,
            r = (r, p, z) -- in polar coordinates;
        nMaxModeOrder : int
            maximal mode number in the sun;
        i, j : int
            tensor indicies over cylindrical coordinates \rho, \phi, z;
        tol : float
            relative tolerance for the sum (how many modes 'n' to consider?);
            |G^{N} - G^{N-1}|/G^{N};

    Returns
    -------
        G : complex
            one component of Gij;
        nmodes : number of considered modes
            condition for cutting the 'n' exists

    =#

    # radiating modes
    k1 = sqrt(epsm) * k        
    rel = 0.0
    Gsij = 0. # Green's tensor component
    nmodes = 0 # number of considered 'n' modes
    Gnprev = 0. # previous, it sums from 0 to nmax-1; when nmax = 0 Gnprev = 0

    for n = 0:nMaxModeOrder
        Gsij += quadgk(x -> iGs11ij(x, rr, pr, zr, rs, ps, zs, 
                                   k, epsf, epsm, Rf, n, i, j), -k1, k1; rtol=tol, maxevals=maxevals)[1]
        rel = abs(abs(Gsij) - abs(Gnprev)) / abs(Gsij)
        if rel < tol
            break
        end
        Gnprev = Gsij
        nmodes += 1
    end

    return Gsij, nmodes
end


function Gs11PolarRadModes(vecRr, vecRs, k, epsf, epsm, Rf, nmax; tol=1e-4, maxevals=4000)

    Gs = zeros(ComplexF64, 3, 3)
    #=
    Returns full tensor Gs in polar coordinates

    Parameters
    ----------
        \vec{r}_s, vec{r}_r : coordinates of source and reciever 
                in polar coordinates (r, phi, z)
    =#
    rr, pr, zr = vecRr[1], vecRr[2], vecRr[3]
    rs, ps, zs = vecRs[1], vecRs[2], vecRs[3]
    if vecRr == vecRs
        # then Gs is diagonal
        #      |rr        |
        # Gs = |    pp    |
        #      |        zz|
        for i = 1:3
            Gs[i, i] = Gs11ijpolarRadModes(rr, pr, zr, 
                        rs, ps, zs, 
                        k, epsf, epsm, Rf,
                        nmax, i, i; tol=tol, maxevals=maxevals)[1]
        end
    elseif zr == zs
        # if zr = zs then 
        #      |rr  rp  0 |
        # Gs = |pr  pp  0 |
        #      |0   0   zz|
        for i = 1:3
            Gs[i, i] = Gs11ijpolarRadModes(rr, pr, zr, rs, ps, zs, k, epsf, epsm, Rf, nmax, i, i;
             tol=tol, maxevals=maxevals)[1]
        end
        Gs[1, 2] = Gs11ijpolarRadModes(rr, pr, zr, rs, ps, zs, k, epsf, epsm, Rf,nmax, 1, 2;
         tol=tol, maxevals=maxevals)[1]
        Gs[2, 1] = Gs11ijpolarRadModes(rr, pr, zr, rs, ps, zs, k, epsf, epsm, Rf,nmax, 2, 1;
         tol=tol, maxevals=maxevals)[1]
    else
        # calc all components
        for i = 1:3, j = 1:3
            Gs[i, j] = Gs11ijpolarRadModes(rr, pr, zr, rs, ps, zs, k, epsf, epsm, Rf, nmax, i, j;
             tol=tol, maxevals=maxevals)[1]
        end
    end

    return Gs
end