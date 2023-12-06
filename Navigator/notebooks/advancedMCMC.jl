

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
using Distributions
using DynamicPPL
using AdvancedHMC
using Zygote
using LogDensityProblems
using ForwardDiff


cl_true=get_Cl(lmax=64)


@model function CMB(data::Vector{Float64})
    lmax = size(data)[1]-1

    
   
    ln10As ~ Uniform(0.25, 0.35)
    ns     ~ Uniform(0.88, 1.06)

    function c(A, B,C)
        t=get_Cl(As=A, ns=B,lmax=C)
        return t
    end

    cl=c(10*ln10As, ns,lmax)
 return cl
   
    
end


### transforming turing cl model into logdensityproblem
CMB(cl_true)
logdens=LogDensityFunction(CMB(cl_true))

###sample the advanced mcmc  ---works but not sure of the results???
D=2
n_samples, n_adapts = 2_000, 1_00
metric = DiagEuclideanMetric(D)
hamiltonian = Hamiltonian(metric, logdens, ForwardDiff)
integrator = Leapfrog(0.005)
kernel = HMCKernel(Trajectory{MultinomialTS}(integrator, GeneralisedNoUTurn()))
adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))
samples, stats = sample(hamiltonian, kernel, initial_θ, n_samples, adaptor, n_adapts; progress=true)
samples

LogDensityProblems.dimension(logdens)


#### analytical gradients

###dcl/dtheta
x=[3.043,0.9645,67.54,0.02217,0.1191,0.0571]
function cl(x)  
    get_Cl(As=x[1],ns=x[2],H0=x[3],wb=x[4],wc=x[5],tau=x[6])
end

derivative=ForwardDiff.jacobian(cl, x)
dCl=derivative[1:65,:] #### derivatives of cl with respect to parameters, if you want to change lmax change 65 to new_lmax+1

####dpsi/dcl

function derivative_cl(alm::Alm{ComplexF64, Vector{ComplexF64}}, cl::Vector{Float64})
    nul=zeros(size(cl)[1])
    newcl=alm2cl(alm)
    pushfirst!(newcl,0)
    pushfirst!(newcl,0)
    for i in 2:size(cl)[1]-1
        nul[i+1]=(1/2 + i)*cl[i+1] - 0.5*newcl[i+1]/(cl[i+1]^2)
    end
    
    return nul
end

alm=synalm(cl_true) #### testing
derivative_cl(alm,cl_true)


####dpsi/dalm 

nside_test = 64
test_map = Navigator.make_noisymap(cl_test, nside_test; seed=1234) ## map
alm=synalm(cl_true)

####returns a vector, each of the elements is a derivative mesured at specific alm (ordering same as healpix)
function derivative_alm(map::HealpixMap{Float64, Healpix.RingOrder},alm::Alm{ComplexF64, Vector{ComplexF64}},cl::Vector{Float64})
    nside=npix2nside(size(map)[1])
    lmax = nside 
    mmax = lmax 
    inverse_noise=ones(length(map.pixels))*(sqrt(10.0)^-1)

    firstpart=map2alm(Healpix.HealpixMap{Float64, Healpix.RingOrder}(inverse_noise.*(map.-alm2map(alm,nside))),lmax=nside)

    newalm = Healpix.Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
    for l = 0:lmax
        if l <= mmax
            maxm = l
        else
            maxm = mmax
        end
        for m = 0:maxm
            i = Healpix.almIndex(alm, l, m)
            if l<= 1
                newalm.alm[i]=0.0
            else
                if m==0
                    newalm.alm[i] = alm.alm[i]/cl[l+1]
                else
                    newalm.alm[i] = 2*alm.alm[i]/cl[l+1]
    
                end
            end
        end
    end  
    final=newalm+firstpart
    return final.alm

    
end


derivative_alm(test_map,alm,cl_true)



###trying AdvancedHMC with our model

@model function CMB_model(data::Vector{Float64}, noise::Vector{Float64})
    nside = Healpix.npix2nside(length(data))
    lmax = nside
    mmax = lmax 
    alm = Healpix.Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
   

    ln10As ~ Uniform(0.25, 0.35)
    ns     ~ Uniform(0.88, 1.06)
    h0     ~ Uniform(0.60, 0.80)
    ωb     ~ Uniform(0.1985, 0.25)
    ωc     ~ Uniform(0.08, 0.20)
    τ      ~ Normal(0.0506, 0.0086)

    function c(A, B,C,D,E,F,G)
        t=get_Cl(As=A, ns=B,H0=C,wb=D,wc=E,tau=F,lmax=G)
    return t
    end
    
    cl=c(10*ln10As,ns,100*h0, ωb/10, ωc, τ,lmax)


    for l = 0:lmax
        if l <= mmax
            maxm = l
        else
            maxm = mmax
        end
        for m = 0:maxm
            i = Healpix.almIndex(alm, l, m)
            if l<= 1
                alm.alm[i]=0.0
            else
                if m==0
                    alm.alm[i] ~ Normal(0, sqrt(cl[l+1]))
                
                else
                    alm.alm[i].re ~ Normal(0, sqrt(cl[l+1]/2))  
                    alm.alm[i].im ~ Normal(0, sqrt(cl[l+1]/2))
            
                end
            end
        end
    end  
    m = Healpix.alm2map(alm, nside)
    for i in 1:size(data)[1]
         data[i] ~ Normal(m.pixels[i], noise[i])
    end
    
    return cl,alm,m
   
    
end

cl_test=Navigator.get_Cl(lmax=64) ##fiducial cls
nside_test = 64
test_map = Navigator.make_noisymap(cl_test, nside_test; seed=1234) ## map
random_noise = ones(length(test_map.pixels))*sqrt(10.0) 

cmb_model=CMB_model(test_map.pixels,random_noise)
lp=LogDensityFunction(cmb_model)


#### look at here how to cunstruct a gradient https://turinglang.org/AdvancedHMC.jl/stable/#Convenience-Constructors
