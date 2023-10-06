module Navigator
include("power_spectrum.jl")
include("mapmaker.jl")
include("analysis.jl")
include("download_data.jl")

using .power_spectrum
using .mapmaker
using .analysis

export CMBcl,linl, CMBdl, noise_spectrum,derivativeDl
export plot_CMB, vary_amplitude,vary_H0,vary_index,vary_omegab,vary_omegac,vary_tau
export derivative_of_amplitude,derivative_of_H0,derivative_of_index,derivative_of_omegab,derivative_of_omegac,derivative_of_tau
export alm,Map, Noise, NoisyMap
export clestimate,clobserved,dlestimate,dlobserved


end
