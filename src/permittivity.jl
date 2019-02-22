# by Ivan Toftul
# 2018/08

# include general libs
using CSV
using Interpolations

export getEpsilonFuncOfWl

#=
    Structure of the data files are
    Wavelength, µm	| n 	| k
    ...             | ...   | ...

    take data from here
        https://refractiveindex.info/?shelf=main&book=Ag&page=Werner
    NOTE: download from the plot!
=#

location = "./lib/refractiveIndexData/"

AuFile = location * "AuWerner2009.csv"
AgFile = location * "AgWerner2009.csv"


function getEpsilonFuncOfWl(material, customValue=1.0)
    #=
        returns complex epsilon as a function of wavelength [m] in vacuum
    =#
    
    # small losses always present
    epsIM = 1e-5
    
    if material == "Au"
        file = AuFile
    elseif material == "Ag"
        file = AgFile
    elseif material == "Silica"
        return foo(wl) = 3.9 + epsIM
    elseif material == "Silicon"
        return bar(wl) = 11.65 + epsIM
    elseif material == "Polystyrene"
        return var(wl) = 2.5 + epsIM
    else
        return baz(wl) = customValue
    end
    
    data = CSV.read(file)
    # n -> eps
    # eps = (n + i k)^2 = n^2-k^2 + i 2nk
    epsRe = data[:, 2].^2 - data[:, 3].^2
    epsIm = 2data[:, 2] .* data[:, 3]
    
    # creating interpolated functions
    epsReInter = LinearInterpolation(data[:, 1] * 1e-6, epsRe)
    epsImInter = LinearInterpolation(data[:, 1] * 1e-6, epsIm)
    
    return qux(wl) = epsReInter(wl) + im * epsImInter(wl)
end