





module mapmaker

include("power_spectrum.jl")
export get_alm, make_map, make_noise, make_noisymap
using Random
using Healpix
using Plots
using Statistics
using Distributions
using .power_spectrum


function Cl(;As=3.043,ns=0.9645,H0=67.54,wb=0.02217,wc=0.1191,tau=0.0571,lmax=5000)
    cl=copy(get_CMBcl(;As,ns,H0,wb,wc,tau,lmax))
    ncl=pushfirst!(cl,0)
    nncl=pushfirst!(ncl,0)
    return nncl
end


@doc raw"""
    get_alm(;As,ns,H0,wb,wc,tau,lmax)

Calculates spherical coeficients(``a_{lm}``) from the power spectrum, given the cosmological coeficients. 
# Arguments
- `As::Float64 = 3.043`: amplitude of primordial perturbations.

- `ns::Float64 = 0.964`: spectral index of the power spectrum.

- `H0::Float64 = 67.54`: Hubble expansion rate at current time.

- `wb::Float64 = 0.02217`: the ratio of baryonic matter (``\\Omega_b *h^2``).

- `wc::Float64 = 0.1191`: the ratio of dark matter (``\\Omega_c *h^2``).

- `tau::Float64 = 0.0571`: optical depth parameter in the reionisation era.

- `lmax::Integer = 5000`: maximum multipole moment which the function calculates.
# Returns
-`Healpix.Alm{ComplexF64, Vector{ComplexF64}}`: returns the alm coeficients.
"""
function get_alm(;As=3.043,ns=0.9645,H0=67.54,wb=0.02217,wc=0.1191,tau=0.057,lmax=5000)
    synalm(get_CMBdl(;As,ns,H0,wb,wc,tau,lmax))
end


@doc raw"""
    make_map(Nside;As,ns,H0,wb,wc,tau,lmax,seed) 

Gives the map of the CMB power spectrum.
# Arguments
- `Nside::Integer, must be a power of 2``: the Nside parameter for the number of pixels in the map.

- `As::Float64 = 3.043`: amplitude of primordial perturbations.

- `ns::Float64 = 0.964`: spectral index of the power spectrum.

- `H0::Float64 = 67.54`: Hubble expansion rate at current time.

- `wb::Float64 = 0.02217`: the ratio of baryonic matter (``\\Omega_b *h^2``).

- `wc::Float64 = 0.1191`: the ratio of dark matter (``\\Omega_c *h^2``).

- `tau::Float64 = 0.0571`: optical depth parameter in the reionisation era.

- `lmax::Integer = 2*Nside`: maximum multipole moment of the power spectrum for which the map is calculated.

- `seed::Integer = 1234`: fixes the seed for the MersenneTwister random generator
# Returns
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: returns the Healpix map of the CMB in RingOrdering.

# Alternative imput
    make_map(cl,Nside)
# Alternative arguments
- `cl::Vector{Float64}`: a set of power spectrum ``C_l``s for which to calculate a map.

- `Nside::Integer, must be a power of 2 `: the Nside parameter for the number of pixels in the map.

- `seed::Integer = 1234`: fixes the seed for the MersenneTwister random generator
# Returns
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: returns the Healpix map of the CMB in RingOrdering.
"""
function make_map(Nside::Integer;As=3.043,ns=0.9645,H0=67.54,wb=0.02217,wc=0.1191,tau=0.057,lmax=2*Nside,seed=1234)
    rng = MersenneTwister(seed)
    Map=synfast(Cl(;As,ns,H0,wb,wc,tau,lmax), Nside,rng)
    return Map 
end


function make_map(cl::Vector{Float64}, Nside::Integer;seed::Integer)
    rng = MersenneTwister(seed)
    Map=synfast(cl, Nside,rng)
    return Map 
end


@doc raw"""
    make_noise(Nside;sigma)

Gives the map of the Gaussian noise with mean zero and dispersion `sigma`.
# Arguments
- `Nside::Integer, must be a power of 2 `: the Nside parameter for the number of pixels in the map.

- `sigma::Float64 = 10`: dispersion of Gaussian distribution from which the noise is drawn.

- `seed::Integer = 1234`: fixes the seed for the MersenneTwister random generator
# Returns
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: returns the Healpix map of the CMB in RingOrdering.
"""
function make_noise(Nside::Integer;sigma=10.0,seed=1234)
    rng = MersenneTwister(seed)
    gauss=Normal(0, sigma)
    Noise=rand(rng,gauss,nside2npix(Nside))
    NoiseMp=Healpix.HealpixMap{Float64, Healpix.RingOrder}(Noise)
    return NoiseMp
end


@doc raw"""
    make_noisymap(Nside;As,ns,H0,wb,wc,tau,lmax,sigma)

Gives the map of the CMB power spectrum with noise added.
# Arguments
- `Nside::Integer, must be a power of 2 = 2048`: the Nside parameter for the number of pixels in the map.

- `As::Float64 = 3.043`: amplitude of primordial perturbations.

- `ns::Float64 = 0.964`: spectral index of the power spectrum.

- `H0::Float64 = 67.54`: Hubble expansion rate at current time.

- `wb::Float64 = 0.02217`: the ratio of baryonic matter (``\\Omega_b *h^2``).

- `wc::Float64 = 0.1191`: the ratio of dark matter (``\\Omega_c *h^2``).

- `tau::Float64 = 0.0571`: optical depth parameter in the reionisation era.

- `lmax::Integer = 5000`: maximum multipole moment of the power spectrum for which the map is calculated.

- `sigma::Float64 = 10`: dispersion of Gaussian distribution from which the noise is drawn.

- `seed::Integer = 1234`: fixes the seed for the MersenneTwister random generator
# Returns
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: returns the Healpix map of the CMB in RingOrdering with npise added.

# Alternative imput
    make_noisymap(cl,Nside;sigma)
# Alternative arguments
- `cl::Vector{Float64}`: a set of power spectrum ``C_l``s for which to calculate a map.

- `Nside::Integer, must be a power of 2 `: the Nside parameter for the number of pixels in the map.

- `sigma::Float64 = 10`: dispersion of Gaussian distribution from which the noise is drawn.


- `seed::Integer = 1234`: fixes the seed for the MersenneTwister random generator
# Returns
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: returns the Healpix map of the CMB in RingOrdering with npise added.
"""
function make_noisymap(Nside::Integer;As=3.043,ns=0.9645,H0=67.54,wb=0.02217,wc=0.1191,tau=0.057,sigma=10.0,lmax=2*Nside,seed=1234)
    Mp=make_map(Nside;As,ns,H0,wb,wc,tau,lmax,seed).+make_noise(Nside;sigma,seed)
    Mapa=Healpix.HealpixMap{Float64, Healpix.RingOrder}(Mp)
    return Mapa
end


function make_noisymap(cl::Vector{Float64},Nside::Integer;sigma=10.0,seed=1234)
    Mp=make_map(cl,Nside;seed).+make_noise(Nside;sigma,seed)
    Mapa=Healpix.HealpixMap{Float64, Healpix.RingOrder}(Mp)
    return Mapa
end

end