





module mapmaker

include("power_spectrum.jl")
export alm,Map, Noise, NoisyMap
using Random
using Healpix
using Plots
using Statistics
using Distributions
using .power_spectrum
rng = MersenneTwister(1234)


function Cl(;As=3.043,ns=0.9645,H0=67.54,wb=0.02217,wc=0.1191,tau=0.0571,lmax=5000)
    cl=copy(CMBcl(;As,ns,H0,wb,wc,tau,lmax))
    ncl=pushfirst!(cl,0)
    nncl=pushfirst!(ncl,0)
    return nncl
end

@doc raw"""
    alm(;As,ns,H0,wb,wc,tau,lmax)

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
function alm(;As=3.043,ns=0.9645,H0=67.54,wb=0.02217,wc=0.1191,tau=0.057,lmax=5000)
    synalm(CMBdl(;As,ns,H0,wb,wc,tau,lmax))
end
alm()


@doc raw"""
    Map(;As,ns,H0,wb,wc,tau,lmax,Nside)

Gives the map of the CMB power spectrum.
# Arguments
- `As::Float64 = 3.043`: amplitude of primordial perturbations.
- `ns::Float64 = 0.964`: spectral index of the power spectrum.
- `H0::Float64 = 67.54`: Hubble expansion rate at current time.
- `wb::Float64 = 0.02217`: the ratio of baryonic matter (``\\Omega_b *h^2``).
- `wc::Float64 = 0.1191`: the ratio of dark matter (``\\Omega_c *h^2``).
- `tau::Float64 = 0.0571`: optical depth parameter in the reionisation era.
- `lmax::Integer = 5000`: maximum multipole moment of the power spectrum for which the map is calculated.
- `Nside::Integer, must be a power of 2 = 2048`: the Nside parameter for the number of pixels in the map.
# Returns
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: returns the Healpix map of the CMB in RingOrdering.
"""
function Map(;As=3.043,ns=0.9645,H0=67.54,wb=0.02217,wc=0.1191,tau=0.057,lmax=5000,Nside=2048)
    Map=synfast(Cl(;As,ns,H0,wb,wc,tau,lmax), Nside,lmax,rng)
    return Map 
end

@doc raw"""
    Noise(;sigma,Nside)

Gives the map of the Gaussian noise with mean zero and dispersion `sigma`.
# Arguments
- `sigma::Float64 = 10`: dispersion of Gaussian distribution from which the noise is drawn.
- `Nside::Integer, must be a power of 2 = 2048`: the Nside parameter for the number of pixels in the map.
# Returns
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: returns the Healpix map of the CMB in RingOrdering.
"""
function Noise(;sigma=10,Nside=2048)
    gauss=Normal(0, sigma)
    Noise=rand(rng,gauss,nside2npix(Nside))
    NoiseMp=Healpix.HealpixMap{Float64, Healpix.RingOrder}(Noise)
    return NoiseMp
end

@doc raw"""
    NoisyMap(;As,ns,H0,wb,wc,tau,lmax,sigma,Nside)

Gives the map of the CMB power spectrum with noise added.
# Arguments
- `As::Float64 = 3.043`: amplitude of primordial perturbations.
- `ns::Float64 = 0.964`: spectral index of the power spectrum.
- `H0::Float64 = 67.54`: Hubble expansion rate at current time.
- `wb::Float64 = 0.02217`: the ratio of baryonic matter (``\\Omega_b *h^2``).
- `wc::Float64 = 0.1191`: the ratio of dark matter (``\\Omega_c *h^2``).
- `tau::Float64 = 0.0571`: optical depth parameter in the reionisation era.
- `lmax::Integer = 5000`: maximum multipole moment of the power spectrum for which the map is calculated.
- `sigma::Float64 = 10`: dispersion of Gaussian distribution from which the noise is drawn.
- `Nside::Integer, must be a power of 2 = 2048`: the Nside parameter for the number of pixels in the map.
# Returns
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: returns the Healpix map of the CMB in RingOrdering.
"""
function NoisyMap(;As=3.043,ns=0.9645,H0=67.54,wb=0.02217,wc=0.1191,tau=0.057,lmax=5000,sigma=10,Nside=2048)
    Mp=Map(;As,ns,H0,wb,wc,tau,lmax,Nside).+Noise(;sigma,Nside)
    Mapa=Healpix.HealpixMap{Float64, Healpix.RingOrder}(Mp)
    return Mapa
end

end