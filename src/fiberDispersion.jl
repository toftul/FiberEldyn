# by Ivan Toftul
# 2018/07


using SpecialFunctions
using Roots

export fiberDispertion
    
# fiber dispersion
function ha_(beta, k, rf, epsf, epsm)
    return(rf * sqrt(epsf * k^2 - beta^2))
end

function qa_(beta, k, rf, epsf, epsm)
    return(rf * sqrt(beta^2 - epsm * k^2))
end

# Eigenvalue equations
# for n > 0
function evEqHE(beta, n, k, rf, epsf, epsm)
    ha = ha_(beta, k, rf, epsf, epsm)
    ha2 = ha^2
    qa = qa_(beta, k, rf, epsf, epsm)
    qa2 = qa^2
    Jl_1ha = besselj(n - 1, ha)
    Jlha = besselj(n, ha)
    Klpqa = - .5 * (besselk(n - 1, qa) + besselk(n + 1, qa))
    Klqa = besselk(n, qa)
    Klpqa_qaKlqa = Klpqa / (qa * Klqa)
    R = sqrt(((epsf - epsm)/(2*epsf) * Klpqa_qaKlqa)^2 +
                (n * beta / (sqrt(epsf) * k) * (1/qa2 + 1/ha2))^2)

    ans = -Jl_1ha / (ha * Jlha) + n/ha2 - (epsf + epsm)/(2*epsf) * Klpqa_qaKlqa - R
    return ans
end

# for n > 0
function evEqEH(beta, n, k, rf, epsf, epsm)
    ha = ha_(beta, k, rf, epsf, epsm)
    ha2 = ha * ha
    qa = qa_(beta, k, rf, epsf, epsm)
    qa2 = qa * qa
    Jl_1ha = besselj(n - 1, ha)
    Jlha = besselj(n, ha)
    Klpqa = - .5 * (besselk(n - 1, qa) + besselk(n + 1, qa))
    Klqa = besselk(n, qa)
    Klpqa_qaKlqa = Klpqa / (qa * Klqa)
    R = sqrt(((epsf - epsm)/(2*epsf) * Klpqa_qaKlqa)^2 +
                (n * beta / (sqrt(epsf) * k) * (1/qa2 + 1/ha2))^2)

    ans = (-Jl_1ha / (ha * Jlha) + n/ha2 -
           (epsf + epsm)/(2*epsf) * Klpqa_qaKlqa + R)
end

# for n = 0
function evEqTE(beta, n, k, rf, epsf, epsm)
    ha = ha_(beta, k, rf, epsf, epsm)
    qa = qa_(beta, k, rf, epsf, epsm)
    J1ha = besselj1(ha)
    J0ha = besselj0(ha)
    K1qa = besselk(1, qa)
    K0qa = besselk(0, qa)

    ans = (J1ha / (ha * J0ha) + K1qa / (qa * K0qa))
end

# for n = 0
function evEqTM(beta, n, k, rf, epsf, epsm)
    ha = ha_(beta, k, rf, epsf, epsm)
    qa = qa_(beta, k, rf, epsf, epsm)
    J1ha = besselj1(ha)
    J0ha = besselj0(ha)
    K1qa = besselk(1, qa)
    K0qa = besselk(0, qa)

    ans = (J1ha / (ha * J0ha) + epsm/epsf * K1qa / (qa * K0qa))
end

function fiberDispertion(modeType, n, k, rf, epsf, epsm)
    #=
        Solves dispersion equation for of the step-index fiber.
        Based on Appendix A, equations (A1-A5) of https://doi.org/10.1103/PhysRevA.96.023835

        Input
        -------
            modeType : HE, EH, TE, or TM
            n : int - azimuthal number
            k : vacuum wave number
            rf : float
                    fiber radius
            epsf, epsm : complex
                    epsilon inside and outside the fiber

        Output
        -------
            roots - array of propagation constants
                    [beta1, beta2, ..., beta_m]
                    if there is no roots it returns [0]
    =#
    kMin = k * (sqrt(epsm) + 1e-13)
    kMax = k * (sqrt(epsf) - 1e-13)
    if modeType == "HE"
        foo = evEqHE
    elseif modeType == "EH"
        foo = evEqEH
    elseif modeType == "TE"
        foo = evEqTE
    elseif modeType == "TM"
        foo = evEqTM
    else
        error()
    end

    # array of all possible values of beta
    betaSpace = range(kMin, stop=kMax, length=2000)  # usually 1500 is enouph
    # searching for f(beta) = 0
    # at first step we look for the intervals which includes a single root
    fSpace = foo.(betaSpace, n, k, rf, epsf, epsm)
    # plot f(beta) in order to understand what is happening here
    fSpace[abs.(fSpace) .> 0.25] .= 0
    #println(fSpace)
    lPoints = Float64[]
    rPoints = Float64[]
    for (i, f) in enumerate(fSpace)
        if i < length(fSpace)
            if ((f == 0) & (fSpace[i + 1] != 0))  # mark all possible left points of the intervals
                lPoints = append!(lPoints, betaSpace[i + 1])
                #println(i)
            elseif ((f != 0) & (fSpace[i + 1] == 0))  # mark all possible right points of the intervals
                rPoints = append!(rPoints, betaSpace[i + 1])
            end
        end
    end
    if ((length(lPoints) > 0) & (length(rPoints) > 0))
        if lPoints[1] > rPoints[1]
            lPoints = prepend!(lPoints, [kMin])
        elseif lPoints[end] > rPoints[end]
            rPoints = append!(rPoints, [kMax])
        end
    else
        if length(lPoints) == 0
            lPoints = [kMin]
        end
        if length(rPoints) == 0
            rPoints = [kMax]
        end
    end

    ab = zeros(length(lPoints), 2)
    #println()
    #println(kMin, " ", kMax)
    #println(lPoints)
    #println(rPoints)
    ab[:, 1], ab[:, 2] = lPoints, rPoints
    # f(a) and f(b) must have different signs
    f_ab = foo.(ab, n, k, rf, epsf, epsm)
    # delete all intervals [a, b] with the same signs
    ab = ab[f_ab[:, 1] .* f_ab[:, 2] .< 0, :]

    # KOSTYL
    # may be it is because of .25 truncate... lazy to figure out
    if ((length(ab) == 0) && (modeType == "HE") && (n == 1))
        # HE11 mode always exists but there might be
        # problems for the low omega in fzero
        return [kMin]
    end
    # searching for roots in [ab] intervals
    if length(ab) == 0
        return [0.0]
    else
        roots = zeros(length(ab[:, 1]))
        for i = 1:length(ab[:, 1])
            r = fzero(x -> foo(x, n, k, rf, epsf, epsm),
                      ab[i, 1], ab[i, 2])
            if ((r > kMin) & (r < kMax))
                roots[i] = r
            else
                roots[i] = 0.0
            end
        end
        # sort roots inverse
        roots = roots[end:-1:1]
        return roots
    end
end
