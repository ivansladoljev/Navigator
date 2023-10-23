using Pkg
Pkg.activate("Navigator")
include("Navigator")
using Navigator
using Random
using Healpix
using Plots
using Distributions
using Turing


rng = MersenneTwister(1234)

function generator(cl::Vector{Float64};seed::Integer=1234) 
    rng = MersenneTwister(seed)
    cl_size = length(cl)
    lmax = cl_size-1
    mmax = cl_size-1
    alm = Healpix.Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
    for l = 0:lmax
        if l <= mmax
            maxm = l
        else
            maxm = mmax
        end
        for m = 0:maxm
            i = Healpix.almIndex(alm, l, m)
            if m==0
                gauss=Normal(0, sqrt(cl[l+1]))
                alm.alm[i]=rand(rng,gauss)
            else
                gauss=Normal(0, sqrt(cl[l+1]/2))
                real=rand(rng,gauss)
                comp=rand(rng,gauss)
                alm.alm[i]=complex(real,comp)
            end
        end
    end
    alm
end


#c=get_CMBcl(lmax=16)

#ncl=pushfirst!(c,0)
#nncl=pushfirst!(ncl,0)

#proba=generator(nncl)
#sizeof(proba)




#newm = Healpix.alm2map(proba, 16)
#plot(newm)

#a=get_dlobserved(newm)
#x=2:4096

##d=get_CMBdl(lmax=4096)

#plot(x,a)
#plot!(x,d)




cl=get_CMBcl(lmax=4)
ncl=pushfirst!(cl,0)
nncl=pushfirst!(ncl,0)

function sample_alm(cl::Vector{Float64}) 
    alm=zeros(length(cl)-1)
    for i in 1:length(cl)-1
        alm[i]= rand(rng,Normal(0, sqrt(cl[i])))
    end
    return alm
end


@model function pr(cl::Vector{Float64}) 
    al=zeros(length(cl))
    for i in 1:length(cl)
        al[i] ~ Normal(0, sqrt(cl[i]))
    end
    return al
end

generator(nncl)

a=pr(nncl)

zz=rand(rng,a)


@model function pr1(cl::Vector{Float64}) 
    al=zeros(length(cl))
    for i in 1:length(cl)
        real ~ Normal(0, sqrt(cl[i]))
        imag ~ Normal(0, sqrt(cl[i])) 
        al[i]=complex(real,imag)
    end
    return al
end

b=pr1(nncl)

zz=rand(rng,b)



r=Healpix.Alm{ComplexF64, Vector{ComplexF64}}
r
map = Healpix.alm2map(r, 4)

function alm_to_Heal(alm)
    size = 3
    lmax = size-1
    mmax = size-1
    Alm = Healpix.Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
    for l = 0:lmax
        if l <= mmax
            maxm = l
        else
            maxm = mmax
        end
        for m = 0:maxm
            i = Healpix.almIndex(alm, l, m)
            Alm[i]=alm[i]
        end
    end
    Alm
end

alm_to_Heal(r)

#@model function sample_cl()
