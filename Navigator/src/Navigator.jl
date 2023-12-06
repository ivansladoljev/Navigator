module Navigator
include("power_spectrum.jl")
include("mapmaker.jl")
include("analysis.jl")

using .power_spectrum
using .mapmaker
using .analysis

export get_CMBcl,linl, get_CMBdl, get_noise_spectrum, derivativeDl
export plot_CMB, vary_amplitude,vary_H0,vary_index,vary_omegab,vary_omegac,vary_tau
export plot_derivative_of_amplitude, plot_derivative_of_H0, plot_derivative_of_index, plot_derivative_of_omegab, plot_derivative_of_omegac, plot_derivative_of_tau
export get_alm, make_map, make_noise, make_noisymap,get_Cl
export get_clestimate, get_clobserved, get_dlestimate, get_dlobserved

end
