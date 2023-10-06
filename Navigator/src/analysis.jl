module analysis
include("mapmaker.jl")
include("power_spectrum.jl")
export clestimate,clobserved,dlestimate,dlobserved
using Healpix
using .mapmaker
using .power_spectrum

@doc raw"""
    clobserved(map;lmax)

Gives the observed power spectrum coeficients ``C_l`` from the given map.
# Arguments
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: the map from which to infere the power spectrum. Must be in RingOrder.
- `lmax::Integer = 5000`: -`lmax::Integer = 5000`:  maximum multipole moment which the function calculates.
# Returns
- `Vector{Float64}`: an array of infered power spectrum coeficients ``C_l``.
"""
function clobserved(map;lmax=5000)
    clobs=anafast(map; lmax, mmax=nothing, niter = 1)
    popfirst!(clobs)
    popfirst!(clobs)
    return clobs
end

@doc raw"""
    dlobserved(map;lmax)

Gives the observed power spectrum coeficients ``D_l=C_l * l(l+1)/2\pi`` from the given map.
# Arguments
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: the map from which to infere the power spectrum. Must be in RingOrder.
- `lmax::Integer = 5000`: -`lmax::Integer = 5000`:  maximum multipole moment which the function calculates.
# Returns
- `Vector{Float64}`: an array of infered power spectrum coeficients ``D_l``.
"""
function dlobserved(map;lmax=5000)
    newl=linl[1:lmax-1]
    dl=clobserved(map;lmax)./((2*pi)./(newl.*(newl.+1)))
    return dl
end



@doc raw"""
   clestimate(map;lmax,sigma,Nside)

Gives the infered power spectrum coeficients ``C_l`` (with noise removed) from the given map.
# Arguments
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: the map from which to infere the power spectrum. Must be in RingOrder.
- `lmax::Integer = 5000`: -`lmax::Integer = 5000`:  maximum multipole moment which the function calculates.
-`sigma::Float64 = 10`: dispersion of Gaussian distribution from which the theoretical noise is drawn
- `Nside::Integer, must be a power of 2 = 2048`: the Nside parameter for the number of pixels in the map.
# Returns
- `Vector{Float64}`: an array of infered power spectrum coeficients ``C_l`` with noise removed.
"""
function clestimate(map;lmax=5000,sigma=10,Nside=2048)
    newl=linl[1:lmax-1]
    clest=clobserved(map,lmax).-noise_spectrum(;lmax,sigma,Nside).*(2*pi)./(newl.*(newl.+1))
    return clest
end

@doc raw"""
   dlestimate(map;lmax,sigma,Nside)

Gives the infered power spectrum coeficients ``D_l=C_l * l(l+1)/2\pi`` (with noise removed) from the given map.
# Arguments
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: the map from which to infere the power spectrum. Must be in RingOrder.
- `lmax::Integer = 5000`: -`lmax::Integer = 5000`:  maximum multipole moment which the function calculates.
-`sigma::Float64 = 10`: dispersion of Gaussian distribution from which the theoretical noise is drawn
- `Nside::Integer, must be a power of 2 = 2048`: the Nside parameter for the number of pixels in the map.
# Returns
- `Vector{Float64}`: an array of infered power spectrum coeficients ``D_l`` with noise removed.
"""
function dlestimate(map;lmax=5000,sigma=10,Nside=2048)
    dlest=dlobserved(map;lmax).-noise_spectrum(;lmax,sigma,Nside)
    return dlest
end


end