using Pkg
Pkg.activate("Navigator")
include("Navigator")
using .Navigator
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




cl=get_CMBcl(lmax=8)
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




@model  function alm_syn(cl::Vector{Float64}) 
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
                alm.alm[i]~Normal(0, sqrt(cl[l+1]))
            else
                alm.alm[i].re ~Normal(0, sqrt(cl[l+1]/2))
                alm.alm[i].im ~Normal(0, sqrt(cl[l+1]/2))
            end
        end
    end
    alm
end

c=alm_syn( nncl)
rng = MersenneTwister(1234)
g=rand(rng,c)


function alm_to_Heal(g)
    lmax=8
    mmax=8
    Alm = Healpix.Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
    for l = 0:lmax
        if l <= mmax
            maxm = l
        else
            maxm = mmax
        end
        for m = 0:maxm
            i = Healpix.almIndex(Alm, l, m)
            Alm[i]=g[i]
        end
    end
    return Alm
end

alm_to_Heal(g)




function test(cl)
    rng = MersenneTwister(1234)
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
            if l<= 1
                alm.alm[i]=0.0
            else
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
    end
    alm
end

molim=test(get_Cl(lmax=8))
newm = Healpix.alm2map(molim, 8)
plot(newm)


myvector=[0 0.1 * i for i in 1:3]

function foo(x)
    a = sin(x)
    b = 0.2 + a
    c = asin(b)
    return c
end


Zygote.gradient(foo, 3)







using Pkg
Pkg.activate("Navigator")
Pkg.precompile()
Pkg.instantiate()
Pkg.resolve()
include("Navigator")
using Navigator
using Random
using Healpix
using Plots
using Distributions
using Turing
using DynamicPPL
using LinearAlgebra
using StatsPlots
using Optim
using NPZ
using Serialization
#Pkg.add("AdvancedMH")
#using AdvancedMH


Turing.setadbackend(:forwarddiff)
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
   
    #function fakecurve(x,A,B)
    #    f= A*exp(-x*B)
      #  return f
    #end
    #prior on model parameters
    ln10As ~ Uniform(0.25, 0.35)
    ns     ~ Uniform(0.88, 1.06)
    h0     ~ Uniform(0.60, 0.80)
    ωb     ~ Uniform(0.1985, 0.25)
    ωc     ~ Uniform(0.08, 0.20)
    τ      ~ Normal(0.0506, 0.0086)
    
   

   # x=0:lmax
  #  test=fakecurve.(x,ln10As,ns)
    
   
    cl=get_Cl(As=10*ln10As, ns=ns,lmax=lmax)
    

    for l = 0:lmax-1
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
                    #alm.alm[i] ~ Normal(0, test[l])
                    #alm.alm[i] = aa[i]
                else
                    alm.alm[i].re ~ Normal(0, sqrt(cl[l+1]/2))  #FIXME change aa to alm.alm
                    alm.alm[i].im ~ Normal(0, sqrt(cl[l+1]/2))
                  # alm.alm[i].re ~ Normal(0, test[l]/2)
                   #alm.alm[i].im ~ Normal(0, test[l]/2)
                    #real~ Normal(0, sqrt(cl[l+1]/2))
                   # comp ~ Normal(0, sqrt(cl[l+1]/2))
                   # alm.alm[i]=complex(real,comp)
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


nside_test = 8
test_map = Navigator.make_noisymap(cl_test, nside_test; seed=1234)
random_noise = ones(length(test_map.pixels))*sqrt(10.0)

test_model_Cl = FLI_CMB_T(test_map.pixels, random_noise)


chain = sample(test_model_Cl, HMC(0.05,10), 1000)
chain2=sample(test_model_Cl, MH(), 100)


returned_stuff = generated_quantities(test_model_Cl, chain2)
returned_stuff[2][1]
#rr=generated_quantities(CMB1(test_map.pixels, random_noise),s)


g=Gibbs(NUTS(100,0.65,:ln10As, :ns,:cl), MH( :alm,:m))
chn = sample(FLI_CMB_T(test_map.pixels, random_noise), g,1000)

gg=Gibbs(MH(:ln10As,:ns),MH(:cl,:alm,:m) )
w=sample(FLI_CMB_T(test_map.pixels, random_noise),gg,100)

:cl=> AdvancedMH.StaticProposal(HMC(0.05,10,:cl))


@model function CMB1(data::Vector{Float64}, noise::Vector{Float64})
    # setting up some variables
    nside = Healpix.npix2nside(length(data))
    lmax = nside
    mmax = lmax 
   
   
    function fakecurve(x,A,B)
        f= A*exp(-x*B)
        return f
    end
    #prior on model parameters
    ln10As ~ Uniform(0.25, 0.35)
    ns     ~ Uniform(0.88, 1.06)
  #  h0     ~ Uniform(0.60, 0.80)
  #  ωb     ~ Uniform(0.1985, 0.25)
  #  ωc     ~ Uniform(0.08, 0.20)
   # τ      ~ Normal(0.0506, 0.0086)
    
   

    x=0:lmax
    test=fakecurve.(x,ln10As,ns)
    
    
   # cl=get_Cl(As=10*ln10As, ns=ns,lmax=lmax)
    
 return test
end

@model function CMB2(data::Vector{Float64}, noise::Vector{Float64})
    # setting up some variables
    nside = Healpix.npix2nside(length(data))
    lmax = nside
    mmax = lmax 
    alm = Healpix.Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
   

   
   # function fakecurve(x,A,B)
    #    f= A*exp(-x*B)
       # return f
    #end
    #prior on model parameters
    ln10As ~ Uniform(0.25, 0.35)
    ns     ~ Uniform(0.88, 1.06)
   # h0     ~ Uniform(0.60, 0.80)
   # ωb     ~ Uniform(0.1985, 0.25)
   # ωc     ~ Uniform(0.08, 0.20)
    #τ      ~ Normal(0.0506, 0.0086)
    
   

   # x=0:lmax
   # test=fakecurve.(x,ln10As,ns)
    
    function c(A, B,C)
        t=get_Cl(As=A, ns=B,lmax=C)
        return t
    end

    cl=c(10*ln10As, ns,lmax)

    for l = 0:lmax-1
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
                    #alm.alm[i] ~ Normal(0, test[l])
                    #alm.alm[i] = aa[i]
                else
                    alm.alm[i].re ~ Normal(0, sqrt(cl[l+1]/2))  #FIXME change aa to alm.alm
                    alm.alm[i].im ~ Normal(0, sqrt(cl[l+1]/2))
                  # alm.alm[i].re ~ Normal(0, test[l]/2)
                   #alm.alm[i].im ~ Normal(0, test[l]/2)
                    #real~ Normal(0, sqrt(cl[l+1]/2))
                   # comp ~ Normal(0, sqrt(cl[l+1]/2))
                   # alm.alm[i]=complex(real,comp)
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

s=sample(CMB1(test_map.pixels, random_noise),HMC(0.005,10) , 100)
s2=sample(CMB2(test_map.pixels, random_noise),HMC(0.005,10) , 100)

params = collect.(eachrow(s.value[1:100, 1, 1]))
#s2=sample(CMB1(test_map.pixels, random_noise),MH() , 100)


#stuff = generated_quantities(CMB2(test_map.pixels, random_noise), s2)[1]



g2=Gibbs(NUTS(100,0.65,:ln10As, :ns,:cl), MH( :alm,:m))
chn2 = sample(CMB2(test_map.pixels, random_noise), g2,1000;save_state=true)

stuff = generated_quantities(CMB2(test_map.pixels, random_noise), chn2)[1][3]
#plot(stuff)

plot(chn2[["ln10As"]]; colordim=:parameter, legend=true)

plot(chn2[:lp])


bestfit_Planck = optimize(CMB2(test_map.pixels, random_noise), MAP(), Optim.Options(iterations=1000, allow_f_increases=true))

#npzwrite("chains_CMB.jls", chn2)
serialize("chain-file.jls", chn2)



@model function simpleCMB(data::Vector{Float64}, noise::Vector{Float64})
    # setting up some variables
    nside = Healpix.npix2nside(length(data))
    lmax = nside
    mmax = lmax 
    alm = Healpix.Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
   
    function fakecurve(x,A,B)
        f= A*exp(-x*B)
        return f
    end
    #prior on model parameters
    ln10As ~ Uniform(0.25, 0.35)
    ns     ~ Uniform(0.88, 1.06)
   # h0     ~ Uniform(0.60, 0.80)
   # ωb     ~ Uniform(0.1985, 0.25)
    #ωc     ~ Uniform(0.08, 0.20)
   # τ      ~ Normal(0.0506, 0.0086)
    
   

    x=0:lmax
    test=fakecurve.(x,ln10As,ns)
    
   
   # cl=get_Cl(As=10*ln10As, ns=ns,lmax=lmax)
    

    for l = 0:lmax-1
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
                    #alm.alm[i] ~ Normal(0, sqrt(cl[l+1]))
                    alm.alm[i] ~ Normal(0, test[l])
                    #alm.alm[i] = aa[i]
                else
                    #alm.alm[i].re ~ Normal(0, sqrt(cl[l+1]/2))  #FIXME change aa to alm.alm
                    #alm.alm[i].im ~ Normal(0, sqrt(cl[l+1]/2))
                   alm.alm[i].re ~ Normal(0, test[l]/2)
                   alm.alm[i].im ~ Normal(0, test[l]/2)
                    #real~ Normal(0, sqrt(cl[l+1]/2))
                   # comp ~ Normal(0, sqrt(cl[l+1]/2))
                   # alm.alm[i]=complex(real,comp)
                end
            end
        end
    end  
    
    m = Healpix.alm2map(alm, nside)
    for i in 1:size(data)[1]
         data[i] ~ Normal(m.pixels[i], noise[i])
    end
    
 return test,alm,m
 

   
  

end


sas=sample(simpleCMB(test_map.pixels, random_noise), MH(), 100)

stuff = generated_quantities(simpleCMB(test_map.pixels, random_noise), sas)[1]
collect.(eachrow(sas.value[1, :, 1]))


using Serialization
chn2 = deserialize("chain-file.jls")

plot(chn2[:lp])
