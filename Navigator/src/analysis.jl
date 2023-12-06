module analysis
include("mapmaker.jl")
include("power_spectrum.jl")
export get_clestimate, get_clobserved, get_dlestimate, get_dlobserved
using Healpix
using .mapmaker
using .power_spectrum

@doc raw"""
    get_clobserved(map;lmax)

Gives the observed power spectrum coeficients ``C_l`` from the given map.
# Arguments
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: the map from which to infere the power spectrum. Must be in RingOrder.

-`lmax::Integer = 2*Nside`:  maximum multipole moment which the function calculates. Maximum default is 2*Nside parameter of the map.
# Returns
- `Vector{Float64}`: an array of infered power spectrum coeficients ``C_l``.
"""
function get_clobserved(map::Healpix.HealpixMap{Float64, Healpix.RingOrder};lmax=2*npix2nside(size(map)[1]))
    clobs=anafast(map; lmax, mmax=nothing, niter = 1)
    popfirst!(clobs)
    popfirst!(clobs)
    return clobs
end



@doc raw"""
    get_dlobserved(map;lmax)

Gives the observed power spectrum coeficients ``D_l=C_l * l(l+1)/2\pi`` from the given map.
# Arguments
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: the map from which to infere the power spectrum. Must be in RingOrder.

-`lmax::Integer = 2*Nside`:  maximum multipole moment which the function calculates. Maximum default is 2*Nside parameter of the map.
# Returns
- `Vector{Float64}`: an array of infered power spectrum coeficients ``D_l``.
"""
function get_dlobserved(map::Healpix.HealpixMap{Float64, Healpix.RingOrder};lmax=2*npix2nside(size(map)[1]))
    newl=power_spectrum.linl[1:lmax-1]
    dl=get_clobserved(map;lmax)./((2*pi)./(newl.*(newl.+1)))
    return dl
end


@doc raw"""
    get_clestimate(map;lmax,sigma)

Gives the infered power spectrum coeficients ``C_l`` (with noise removed) from the given map.
# Arguments
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: the map from which to infere the power spectrum. Must be in RingOrder.

-`lmax::Integer = 2*Nside`:  maximum multipole moment which the function calculates. Maximum default is 2*Nside parameter of the map.

-`sigma::Float64 = 10`: dispersion of Gaussian distribution from which the theoretical noise is drawn
# Returns
- `Vector{Float64}`: an array of infered power spectrum coeficients ``C_l`` with noise removed.
"""
function get_clestimate(map;lmax=2*npix2nside(size(map)[1]),sigma=10.0)
    newl=power_spectrum.linl[1:lmax-1]
    clest=get_clobserved(map;lmax).-(power_spectrum.get_noise_spectrum(;lmax,sigma,Nside=npix2nside(size(map)[1])).*(2*pi)./(newl.*(newl.+1)))
    return clest
end



@doc raw"""
    get_dlestimate(map;lmax,sigma)

Gives the infered power spectrum coeficients ```D_l=C_l * l(l+1)/2\pi``` (with noise removed) from the given map.
# Arguments
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: the map from which to infere the power spectrum. Must be in RingOrder.

-`lmax::Integer = 2*Nside`:  maximum multipole moment which the function calculates. Maximum default is 2*Nside parameter of the map.

-`sigma::Float64 = 10.0`: dispersion of Gaussian distribution from which the theoretical noise is drawn.

# Returns
- `Vector{Float64}`: an array of infered power spectrum coeficients ``D_l`` with noise removed.
"""
function get_dlestimate(map; lmax=2*npix2nside(size(map)[1]), sigma=10.0)
    Nside=npix2nside(size(map)[1])
    dlest=get_dlobserved(map;lmax).-power_spectrum.get_noise_spectrum(;lmax,sigma,Nside)
    return dlest
end





end