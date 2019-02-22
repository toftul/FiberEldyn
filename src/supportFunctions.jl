# by Ivan Toftul
# 2018 / 8 / 20

using SpecialFunctions

export
    meshgrid, 
    cart2pol,
    pol2cart,
    rotMat2D,
    oddeven,
    LegendrePnm,
    LegendrePnm_sph,
    SphericalHarmonicY,
    sphericaljn,
    sphericalyn,
    sphericalh1,
    sphericalh2,
    sphericaljnp,
    sphericalynp,
    sphericalh1p,
    sphericalh2p,
    riccatiPsi,
    riccatiPsip,
    riccatiChi,
    riccatiZeta,
    riccatiXi

# for surface plots
function meshgrid(x, y)
    XX = zeros(length(y), length(x))
    YY = zeros(length(y), length(x))

    for i = 1:length(y), j = 1:length(x)
        XX[i, j] = x[j]
        YY[i, j] = y[i]
    end

    return XX, YY
end

# Geometry
function cart2pol(Rcart)
    rho = sqrt(Rcart[1]^2 + Rcart[2]^2)
    phi = atan(Rcart[2], Rcart[1])
    return [rho, phi, Rcart[3]]
end

function pol2cart(Rpol)
    x = Rpol[1] * cos(Rpol[2])
    y = Rpol[1] * sin(Rpol[2])
    return [x, y, Rpol[3]]
end

function rotMat2D(phi)
    # rotates vector field in CCW direction
    R = [ cos(phi)  sin(phi) 0;
         -sin(phi)  cos(phi) 0;
                0         0  1]
    return R
end


function oddeven(m)
    # evaluates (-1)^m
    # yes, it is faster
    if m % 2 == 1
        return -1
    else
        return  1
    end
end


# Special funcitons
function LegendrePnm(n, m, x)
    #=
        Associated Legendre polynomials with Condom phase included

        based on Shanjie Zhang, Computation of specail functions (1996)
        section 4.4.4

        Arguments:
        x -- Argument of Pnm(x)
        n = 0, 1, 2, ...
        m = 0, 1, 2, ..., n
    =#
    if n == 0
        if m == 0
            return 1
        else
            return 0
        end
    elseif x == 1.0
        if m == 0
            return 1
        else
            return 0
        end
    elseif x == -1.0
        if m == 0
            return oddeven(n)
        else
            return 0
        end
    end
    PM = zeros(m+1, n+1)
    PM[1, 1] = 1.0

    XS = 1.0 - x*x
    XQ = sqrt(XS)
    for i = 1:m
        # eq. (4.4.27)
        PM[i+1, i+1] = - (2i - 1) * XQ * PM[i, i]
    end
    for i = 1:min(m+1, n)
        # eq. (4.4.26)
        PM[i, i+1] = (2i - 1) * x * PM[i, i]
    end
    for i = 0:m, j = i+2:n
        # eq. (4.4.11)
        PM[i+1, j+1] = ((2j - 1) * x * PM[i+1, j] - (i+j-1) * PM[i+1, j-1])/(j-i)
    end
    return PM[m+1, n+1]
end


function LegendrePnm_sph(n, m, x)
    #=
        Associated Legendre polynomials with Condom phase included
        normalized to spherical harmonics prefactor
            sqrt( (2n+1)/4pi ) * sqrt( (n-m)! / (m+n)! )

        based on Shanjie Zhang, Computation of specail functions (1996)
        section 4.4.4

        Arguments:
        x -- Argument of Pnm(x)
        n = 0, 1, 2, ...
        m = 0, 1, 2, ..., n
    =#
    if n == 0
        if m == 0
            return 1.0/sqrt(4pi)
        else
            return 0
        end
    elseif x == 1.0
        if m == 0
            return 1.0/sqrt(4pi)
        else
            return 0
        end
    elseif x == -1.0
        if m == 0
            return oddeven(n)
        else
            return 0
        end
    end
    PM = zeros(m+1, n+1)
    PM[1, 1] = 1.0/sqrt(4pi)

    XS = 1.0 - x*x
    XQ = sqrt(XS)
    for i = 1:m
        # modified eq. (4.4.27)
        PM[i+1, i+1] = - sqrt((2i + 1) / 2i) * XQ * PM[i, i]
    end
    for i = 1:min(m+1, n)
        # modified eq. (4.4.26)
        PM[i, i+1] = sqrt(2i + 1) * x * PM[i, i]
    end
    for i = 0:m, j = i+2:n
        # modified eq. (4.4.11)
        PM[i+1, j+1] = sqrt((2j+1)/(j*j - i*i)) * (sqrt(2j - 1) * x * PM[i+1, j] - sqrt((j-i-1)*(i+j-1)/(2j-3)) * PM[i+1, j-1])
    end
    return PM[m+1, n+1]
end


# spherical harmonics
# tested
function SphericalHarmonicY(n, m, theta, phi)
    #=
        based on Generalized Lorenz-Mie theories by Gerard Gouesbet
        eq. (2.92) and (2.75).

        Condom phase included in Legendre polynomials.

        Mathematica has the same definitions as in the textbook of Gerard
        (info in order to know where to check).
    =#
    if abs(m) > abs(n)
        return 0.0
    elseif n <= -1
        # from Mathematica docs for SphericalHarmonicY
        return SphericalHarmonicY(-(n+1), m, theta, phi)
    elseif m < 0
        m *= -1
        return oddeven(m) * LegendrePnm_sph(n, m, cos(theta)) * exp(- im * m * phi)
    else
        return LegendrePnm_sph(n, m, cos(theta)) * exp(im * m * phi)
    end
end


# ### adding bessel spherical funcitons

function sphericaljn(n, z)
    # https://dlmf.nist.gov/10.47
    return sqrt(pi / 2z) * besselj(n + 0.5, z)
end


function sphericalyn(n, z)
    # https://dlmf.nist.gov/10.47
    return sqrt(pi / 2z) * bessely(n + 0.5, z)
end


function sphericalh1(n, z)
    # https://dlmf.nist.gov/10.47
    return sqrt(pi / 2z) * hankelh1(n + .5, z)
end


function sphericalh2(n, z)
    # https://dlmf.nist.gov/10.47
    return sqrt(pi / 2z) * hankelh2(n + .5, z)
end


function sphericaljnp(n, z)
    # https://dlmf.nist.gov/10.51
    return (n * sphericaljn(n - 1, z) - (n + 1) * sphericaljn(n + 1, z)) / (2n + 1)
end


function sphericalynp(n, z)
    # https://dlmf.nist.gov/10.51
    return (n * sphericalyn(n - 1, z) - (n + 1) * sphericalyn(n + 1, z)) / (2n + 1)
end


function sphericalh1p(n, z)
    # https://dlmf.nist.gov/10.51
    return (n * sphericalh1(n - 1, z) - (n + 1) * sphericalh1(n + 1, z)) / (2n + 1)
end


function sphericalh2p(n, z)
    # https://dlmf.nist.gov/10.51
    return (n * sphericalh2(n - 1, z) - (n + 1) * sphericalh2(n + 1, z)) / (2n + 1)
end


# Riccati–Bessel functions
# from Gorodetski - Opticheskiie resonatory s gigantskoy dobrotnostyu (Chapter 5)
# and from http://mathworld.wolfram.com/Riccati-BesselFunctions.html
function riccatiPsi(n, z)
    # or riccatiSn as sin
    return sqrt(0.5pi * z) * besselj(n + .5, z)
end


function riccatiPsip(n, z)
    error("don't be lazy!")
end


function riccatiChi(n, z)
    # or riccatiCn as cos
    return -sqrt(0.5pi * z) * bessely(n + .5, z)
end


function riccatiZeta(n, z)
    # ζl(x) = ψl(x) − iχl(x) = xh(1)(x)
    return sqrt(0.5pi * z) * hankelh1(n + .5, z)
end


function riccatiXi(n, z)
    # ξl(x) = ψl(x) + iχl(x) = xh(2)(x)
    return sqrt(0.5pi * z) * hankelh2(n + .5, z)
end
