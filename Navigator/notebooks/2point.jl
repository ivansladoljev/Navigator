using Pkg
Pkg.activate("Navigator")
#include("Navigator")
using Navigator
using Healpix
using Random
using LinearAlgebra
using Turing
using Serialization
using PairPlots

cl_true=get_CMBcl(lmax=64)
pushfirst!(cl_true,0)
pushfirst!(cl_true,0) ## getting the fiducial cl from l=0 to lmax

alm=Healpix.synalm(cl_true,64,MersenneTwister(111))
Map=Healpix.alm2map(alm,64) ##making a test Map of the fiducial cls (reason why I went first alms then map is for reproducibility-only way for the seed to listen)
cl_est=get_clestimate(Map,lmax=64,sigma=10.)
pushfirst!(cl_est,0)
pushfirst!(cl_est,0) ## get the estimated cls from map [0,lmax]

n=2:64
clnoise=get_noise_spectrum(lmax=64,sigma=10.,Nside=64).*((2*pi)./(n.*(n.+1)))
pushfirst!(clnoise,0)
pushfirst!(clnoise,0) ##get the theoretical noise spectrum

function make_sigmasq(cl::Vector{Float64},noise::Vector{Float64})
    l=size(cl)[1]
    sigma=zeros(l)
    for i in range(1,l)
    sigma[i]=(1/(2i+1)).*(cl[i].^2 .+ noise[i].^2)
    end
return sigma
end
noise=make_sigmasq(cl_true,clnoise) ##get the noise = cl^2 + N^2
#N=Diagonal(noise)

using Plots
plot(0:64,cl_est) #testing what we got
plot!(0:64,cl_true)

@model function CMB_model(data::Vector{Float64},noise::Vector{Float64}) ##turing model - you can remove the addition parameters if you want
    lmax=size(data)[1]-1
    ln10As ~ Uniform(0.25, 0.35)
    ns     ~ Uniform(0.88, 1.06)
    h0     ~ Uniform(0.60, 0.80)
    ωb     ~ Uniform(0.1985, 0.25)
    ωc     ~ Uniform(0.08, 0.20)
    τ      ~ Normal(0.0506, 0.0086)
  #  yₚ     ~ Normal(1.0, 0.0025)

   # θ = [10*ln10As, ns]#, 100*h0, ωb/10, ωc, τ]
    function c(A, B,C,D,E,F,G)
            t=get_Cl(As=A, ns=B,H0=C,wb=D,wc=E,tau=F,lmax=G)
            return t
    end

    cl=c(10*ln10As,ns,100*h0, ωb/10, ωc, τ,lmax)
# H0=C,wb=D,wc=E,,100*h0, ωb/10, ωc,
    #compute theoretical prediction
    #pred = iΓ * cl # ./(yₚ^2)
    #pred = theory_planck(θ) ./(yₚ^2)
    #compute likelihood
     #data ~ Normal(pred, I)
    for i in 1:size(data)[1]
         if i <=2
            data[i]=0
        else
            data[i] ~ Normal(cl[i], noise[i])
        end
    end

    return cl,data
end


chain=sample(CMB_model(cl_true,noise), NUTS(0.65), 50000) #mcmc 

#next is ploting the chains
#Pkg.add("StatsPlots")
#Pkg.add("CairoMakie")
#using CairoMakie, PairPlots

#pairplot(chain)

###if you want to save the chains

#f = open("INSERTNAME.jls", "w")
#serialize(f, chain)
#close(f)

###or open them again
#f = open("NAMEYOUGAVE.jls", "r")
#chain = deserialize(f)
#close(f)
