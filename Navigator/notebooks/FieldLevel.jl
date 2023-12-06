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
 
cl_test=Navigator.get_Cl(lmax=16) ##fiducial cls
ℓ=range(0,16)


nside_test = 16
test_map = Navigator.make_noisymap(cl_test, nside_test; seed=1234) ## map
random_noise = ones(length(test_map.pixels))*sqrt(10.0)  ## noise variance


@model function CMB_model(data::Vector{Float64}, noise::Vector{Float64}) ###model
    nside = Healpix.npix2nside(length(data))
    lmax = nside
    mmax = lmax 
    alm = Healpix.Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
   
    ln10As ~ Uniform(0.25, 0.35)
    ns     ~ Uniform(0.88, 1.06)

    function c(A, B,C)
        t=get_Cl(As=A, ns=B,lmax=C)
        return t
    end

    cl=c(10*ln10As, ns,lmax)

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

###sampling
cmb=CMB_model(test_map.pixels, random_noise)
sample(cmb,MH(),1000)
gibs=Gibbs(HMC(0.05,10,:ln10As, :ns,:cl), MH(:alm,:m))
chn = sample(cmb, gibs,5000;save_state=true)

###saving the chains
f = open("NAME.jls", "w")
serialize(f, chn)
close(f)


###returned cl,alms and map
stuff=generated_quantities(cmb,chn)
stuff[1][1] 

## to see what are our sampled parameters 
names(chn, :parameters)

##parameter values
n=get(chn, :ns)
a=get(chn, :ln10As)
plot(a[1][:],n[1][:],markersize=1.5)

## separate chains 
A=group(chn,  "ln10As")
N=group(chn,  "ns")


## this should be the cov matrix for every alm (right cl for every  l) but isnt working
cova=zeros(153)
s=0
for i in 0:16
     for j in i+1:17
        cova[s+j-i]=cl_test[j]

    end
    s+=17-i
end
covMH=Diagonal(cova)








