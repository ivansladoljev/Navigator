using Pkg
Pkg.activate("Navigator")
include("Navigator")
using Navigator
using Random
using Healpix
using Plots
using Distributions
using Turing
using DynamicPPL
using LinearAlgebra

Random.seed!(12)

cl_test=get_CMBcl(lmax=16)
ncl=pushfirst!(cl_test,0)
nncl=pushfirst!(ncl,0)
ℓ=range(0,16)

@model function FLI_CMB_T(data::Vector{Float64}, noise::Vector{Float64})
    # setting up some variables
    nside = Healpix.npix2nside(length(data))
    lmax = nside
    mmax = lmax 
    alm = Healpix.Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
    aa = zero(alm.alm) #FIXME
    #d = data.pixels
    

    #prior on model parameters
    ln10As ~ Uniform(0.25, 0.35)
    ns     ~ Uniform(0.88, 1.06)
    h0     ~ Uniform(0.60, 0.80)
    ωb     ~ Uniform(0.1985, 0.25)
    ωc     ~ Uniform(0.08, 0.20)
    τ      ~ Normal(0.0506, 0.0086)

    # call the emulator:
    cl=get_CMBcl(As=10*ln10As, ns=ns, H0=100*h0, wb=ωb/10, wc=ωc, tau=τ, lmax=lmax)
    
    # BECAUSE HEALPIX.JL STARTS AT ALM WITH ELL=0!!!!
    cl=pushfirst!(cl,0)
    cl=pushfirst!(cl,0)
    # aa[Healpix.almIndex(alm, 0, 0)] = 0.0
    # aa[Healpix.almIndex(alm, 1, 0)] = 0.0
    # aa[Healpix.almIndex(alm, 1, 1)] = 0.0

#    alms = Healpix.synalm(cl)
    # sampling the Alms:
    for l = 2:length(cl)-1
        if l <= mmax
            maxm = l
        else
            maxm = mmax
        end
        for m = 0:maxm
            i = Healpix.almIndex(alm, l, m)
            if m==0
                aa[i] ~ Normal(0, sqrt(cl[l+1]))
                #alm.alm[i] = aa[i]
            else
                aa[i].re ~ Normal(0, sqrt(cl[l+1]/2)) #FIXME change aa to alm.alm
                aa[i].im ~ Normal(0, sqrt(cl[l+1]/2))
            end
        end
    end
    # # # make the map:
    alm.alm = aa
    m = Healpix.alm2map(alm, nside)

    # this is the likelihood:
    #data ~ MvNormal(m.pixels, diagm(noise)) 
    for i in eachindex(data)
        data[i] ~ Normal(m.pixels[i], noise[i])
    end
    #data ~ filldist(Normal(), length(data), m.pixels, noise)
    return (cl, m, alm)
    #return(alms)
end

@model function FLI_CMB_T_Cheat(data::Vector{Float64}, noise::Vector{Float64})
    # setting up some variables
    nside = Healpix.npix2nside(length(data))
    lmax = nside
    mmax = lmax 
    alm = Healpix.Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
    aa = zero(alm.alm) #FIXME
    #d = data.pixels
    

    #prior on model parameters
    ln10As ~ Uniform(0.25, 0.35)
    ns     ~ Uniform(0.88, 1.06)
    h0     ~ Uniform(0.60, 0.80)
    ωb     ~ Uniform(0.1985, 0.25)
    ωc     ~ Uniform(0.08, 0.20)
    τ      ~ Normal(0.0506, 0.0086)

    # call the emulator:
    cl=get_CMBcl(As=10*ln10As, ns=ns, H0=100*h0, wb=ωb/10, wc=ωc, tau=τ, lmax=lmax)
    
    # BECAUSE HEALPIX.JL STARTS AT ALM WITH ELL=0!!!!
    cl=pushfirst!(cl,0)
    cl=pushfirst!(cl,0)
    # aa[Healpix.almIndex(alm, 0, 0)] = 0.0
    # aa[Healpix.almIndex(alm, 1, 0)] = 0.0
    # aa[Healpix.almIndex(alm, 1, 1)] = 0.0

    alms = Healpix.synalm(cl)

    m = Healpix.alm2map(alm, nside)

    # this is the likelihood:
    #data ~ MvNormal(m.pixels, diagm(noise)) 
    for i in eachindex(data)
        data[i] ~ Normal(m.pixels[i], noise[i])
    end
    #data ~ filldist(Normal(), length(data), m.pixels, noise)
    return (cl, m, alm)
    #return(alms)
end


nside_test = 8
test_map = Navigator.make_noisymap(cl_test, nside_test; seed=1234)
random_noise = ones(length(test_map.pixels))*sqrt(10.0)

test_model_Cl = FLI_CMB_T(test_map.pixels, random_noise)
test_model_Cl_dirty = FLI_CMB_T_Cheat(test_map.pixels, random_noise)
#chain = sample(test_model_Cl, HMC(0.05, 10), 10)
chain = sample(test_model_Cl, MH(), 1000)
chain_cheat = sample(test_model_Cl_dirty, MH(), 1000)

returned_stuff = generated_quantities(test_model_Cl, chain)

anafast_test = Healpix.anafast(returned_stuff[2][2], lmax=16)
plot(anafast_test, label="anafast")
plot!(returned_stuff[2][1], label="Theory")

plot(returned_stuff[1][2])

plot(chain[["ωb"]]; colordim=:parameter, legend=true)

plot(chain[:lp])
#plot!(chain_cheat[["ωb"]]; colordim=:parameter, legend=true)
#plot(chain[[""]]; colordim=:parameter, legend=true)

chain.value.data