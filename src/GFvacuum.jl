# by Ivan Toftul
# 2018/07

using LinearAlgebra

export GFvac, GFvacPolar


function GFvac(R, k0, epsm)
    #=
    Green's function of Maxwell eq. in a free space

    Parameters
    ----------
        R = (Rx, Ry, Rz) : vec column;
        k0 : float
            vacuum k0-vector value, 2pi/\lambda_0 = \omega/c;
        epsm : permettivity of around media

    Returns
    -------
        G : complex array 3x3

    Details
    -------
        R = r1 - r2;
        r1 -- source position;
        r2 -- reciever postion;
        omega -- radiation frequency

    =#
    if norm(R) != 0
        k = sqrt(epsm) * k0
        Rmod = norm(R)
        kR = k * Rmod
        EXPikR4piR = exp(1im * kR) / (4pi * Rmod)
        CONST1 = (1 + (1im * kR - 1) / kR^2)
        CONST2 = (3 - 3im * kR - kR^2) / (kR^2 * Rmod^2)
        # return Green function
        G = EXPikR4piR * (CONST1 * I +
                          CONST2 * R*R')
        # R*R' -- tensor product
        return G
    else
        # Re G0 is divergent but Im G0 is finite
        # return only Im part
        return 1im * I * k0 * sqrt(epsm) / 6pi
    end
end


function GFvacPolar(RrPolar, RsPolar, k0, epsm)
    #=
    Green's function of Maxwell eq. in a free space
    in polar coordinates
        
        G_pol(Rr, Rs) = Q(Rr) G_car(Rr, Rs) Q^T(Rs)

    For details see https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates
    Matrix for transform of coordinates:
        e^p = Q e^c
    Vector field  A = A^c e^c = A^p e^p, so
        A^p = A^c Q^T
        
    Parameters
    ----------
        RrPolar, RsPolar : vec column;
            Rpolar = (rho, theta, z)
            Rr -- reciever postion;
            Rs -- source position;
        k0 : float
            k-vector value, 2pi/\lambda_0 = \omega/c;

    Returns
    -------
        G : complex array 3x3

    Details
    -------
        R = Rr - Rs;
        Rr -- reciever postion;
        Rs -- source position;
        omega -- radiation frequency

    =#
    RrCart = [RrPolar[1] * cos(RrPolar[2]),
              RrPolar[1] * sin(RrPolar[2]),
              RrPolar[3]]
    RsCart = [RsPolar[1] * cos(RsPolar[2]),
              RsPolar[1] * sin(RsPolar[2]),
              RsPolar[3]]
    
    G_car = GFvac(RrCart - RsCart, k0, epsm)

    # Rotation matrix
    Qr = [ cos(RrPolar[2]) sin(RrPolar[2])  0;
          -sin(RrPolar[2]) cos(RrPolar[2])  0;
                        0               0   1 ]
    
    Qs = [ cos(RsPolar[2]) sin(RsPolar[2])  0;
          -sin(RsPolar[2]) cos(RsPolar[2])  0;
                        0               0   1 ]

    return Qr * G_car * Qs'
end
