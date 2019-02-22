# by Ivan Toftul
# 2018/07

# export everything
export
    speedOfLight,
    mu0,
    epsilon0,
    vacuumImpedance,
    qElementary,
    mElectron,
    hbar,
    kBoltzmann,
    VnumberCRITICAL,
    mili,
    micro,
    nano,
    piko,
    femto


# Some physical consts in SI units

# speed of light
const speedOfLight = 299792458.0  # [m/s]
# vacuum permeability
const mu0 = 4pi * 1e-7  # [H / m]

# vacuum permittivity
const epsilon0 = 1.0 / (mu0 * speedOfLight^2)  # [F/m]

const vacuumImpedance = sqrt(mu0 / epsilon0)  # 376.73 [Ohm]

# elementary charge
const qElementary = 1.6e-19  # [C]
# electron mass
const mElectron = 9.1e-31  # [kg]

# Plank const
const hbar = 1.0545718 * 1e-34  # [m^2 kg / s]

# Boltzmann constant
const kBoltzmann = 1.380648528e-23  # [J/K]

# Critical V number for the signle mode spte-index waveguied
# e.g. see Chin-Lin - Chen Foundations for guided-wave optics
# p. 233
const VnumberCRITICAL = 2.4048


# multipliers
const mili   = 1e-3
const micro  = 1e-6
const nano   = 1e-9
const piko   = 1e-12
const femto  = 1e-15
