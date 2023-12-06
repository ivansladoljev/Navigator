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
using Zygote
using Healpix.ChainRulesCore

Turing.setadbackend(:zygote)
Random.seed!(12)

cl_test=get_CMBcl(lmax=16)
ncl=pushfirst!(cl_test,0)
nncl=pushfirst!(ncl,0)
ℓ=range(0,16)


function test(;As=3.043,ns=0.9645,H0=67.54,wb=0.02217,wc=0.1191,tau=0.0571,lmax=2)
    a=zeros(lmax+1)
   i=0
    for i in lmax+1
        if i <=2
            a[i]=0.0
        else
            a[i]=get_CMBcl(As=As, ns=ns, H0=H0, wb=wb, wc=wc, tau=tau, lmax=lmax)[i-2]
        end
        print(a[i])
    end
    return a
end   
test()

#function fakecurve(l,As,ns,tau)
   # f=As+ tau*exp( ns*(-l))
    #return f
#end

@model function FLI_CMB_T(data::Vector{Float64}, noise::Vector{Float64})
    # setting up some variables
    nside = Healpix.npix2nside(length(data))
    lmax = nside
    mmax = lmax 
    alm = Healpix.Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
    #alm = Healpix.Alm{Any, Any}(lmax, mmax)
    #aa = zero(alm.alm) #FIXME
    #d = data.pixels
    

    function fakecurve(x::Any,A::Any,B::Any)
        f= A*exp(-x*B)
        return f
    end
    #prior on model parameters
    #ln10As ~ Uniform(0.25, 0.35)
   # ns     ~ Uniform(0.88, 1.06)
  #  h0     ~ Uniform(0.60, 0.80)
  #  ωb     ~ Uniform(0.1985, 0.25)
  #  ωc     ~ Uniform(0.08, 0.20)
   # τ      ~ Normal(0.0506, 0.0086)
    
   
   
   x=2:lmax
   A~Uniform(1,100)
   B~Uniform(1,100)



    # call the emulator:
   # c=get_CMBcl(As=10*ln10As, ns=ns, H0=100*h0, wb=ωb/10, wc=ωc, tau=τ, lmax=lmax)
    test=fakecurve.(x,A,B)
   # cl = get_Cl()
  #  for i in lmax+1
    #    cl.at[i].set[cl[i]]
   # end
    
   # D~MvNormal(c,I)
 
    # BECAUSE HEALPIX.JL STARTS AT ALM WITH ELL=0!!!!
    #cl=pushfirst!(cl,0)
    #cl=pushfirst!(cl,0)
    # aa[Healpix.almIndex(alm, 0, 0)] = 0.0
    # aa[Healpix.almIndex(alm, 1, 0)] = 0.0
    # aa[Healpix.almIndex(alm, 1, 1)] = 0.0
    
    #alms = Healpix.synalm(test)
    # sampling the Alms:
    #
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
                    #alm.alm[i] ~ Normal(0, sqrt(test[l+1])
                    alm.alm[i] ~ Normal(0, test[l])
                    #alm.alm[i] = aa[i]
                else
                    #alm.alm[i].re ~ Normal(0, sqrt(test[l+1]/2)) #FIXME change aa to alm.alm
                   # alm.alm[i].im ~ Normal(0, sqrt(test[l+1]/2))
                   alm.alm[i].re ~ Normal(0, test[l]/2)
                   alm.alm[i].im ~ Normal(0, test[l]/2)
                    #real~ Normal(0, sqrt(cl[l+1]/2))
                   # comp ~ Normal(0, sqrt(cl[l+1]/2))
                   # alm.alm[i]=complex(real,comp)
                end
            end
        end
    end  
    # # # make the map:
   # alm.alm = aa
   #
    m = Healpix.alm2map(alm, nside)

    # this is the likelihood:
    #data ~ MvNormal(m.pixels, diagm(noise)) 
  
    for i in 1:size(data)[1]
       data[i] ~ Normal(m.pixels[i], noise[i])
    end
    #
    #data ~ filldist(Normal(), length(data), m.pixels, noise)


    #d=zeros(size(data)[1])
    #for i in size(d[1])
       #d[i] ~ Normal(m.pixels[i], noise[i])
   #end

    return test,alm,data
    #return (cl, m, alm)
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
   
    m = Healpix.alm2map(alms, nside)
   
       # this is the likelihood:
       #data ~ MvNormal(m.pixels, diagm(noise)) 
    for i in eachindex(data)
        data[i] ~ Normal(m.pixels[i], noise[i])
    end
       #data ~ filldist(Normal(), length(data), m.pixels, noise)
    return (cl, m, alms)
end


nside_test = 8
test_map = Navigator.make_noisymap(cl_test, nside_test; seed=1234)
random_noise = ones(length(test_map.pixels))*sqrt(10.0)

test_model_Cl = FLI_CMB_T(test_map.pixels, random_noise)
test_model_Cl_dirty = FLI_CMB_T_Cheat(test_map.pixels, random_noise)


chain = sample(test_model_Cl, HMC(0.005,10), 10)
chain2=sample(test_model_Cl, MH(), 100)
chain_cheat = sample(test_model_Cl_dirty, MH(), 1000)

returned_stuff = generated_quantities(test_model_Cl, chain2)
returned_stuff[5][2]

map=Healpix.alm2map(returned_stuff[5][2],8)
plot(map)
#anafast_test = Healpix.anafast(returned_stuff[5][2], lmax=8


#plot(anafast_test, label="anafast"
plot!(returned_stuff[2][1], label="Theory")

plot(returned_stuff[2][2])

plot(chain[["ln10As"]]; colordim=:parameter, legend=true)

plot(chain[:lp])
plot(chain_cheat[["ωb"]]; colordim=:parameter, legend=true)
#plot(chain[[""]]; colordim=:parameter, legend=true)

chain.value.data

newsampler=Gibbs( HMC(0.05, 10 , :ln10As, :ns, :h0, :ωb, :ωc, :τ ), MH(:aa), HMC(0.05, 10 ,:data)  )

newchain=sample(test_model_Cl,newsampler,100)      





function cc()
    cl = get_Cl()
    for i in 5000+1
     cl[i].set[cl[i]]
    end
    return cl
end
cc()



using Turing
using Capse
using Statistics
using SimpleChains
using Turing
using Optim
using JSON
using NPZ
using BenchmarkTools
using Healpix
using Plots
using LaTeXStrings
using ForwardDiff





@model function proba(data::Vector{Float64}, noise::Vector{Float64})
# setting up some variables
    nside = Healpix.npix2nside(length(data))
    lmax = nside
    mmax = lmax 
    alm = Healpix.Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
#aa = zero(alm.alm) #FIXME
#d = data.pixels


#prior on model parameters
    ln10As ~ Uniform(2.5, 3.5)
    ns     ~ Uniform(0.88, 1.06)
    h0     ~ Uniform(60, 80)
    ωb     ~ Uniform(0.01985, 0.025)
    ωc     ~ Uniform(0.08, 0.20)
    τ      ~ Normal(0.0506, 0.0086)


    NN_dict = JSON.parsefile("Navigator/src/nn_setup.json")


    weights_folder = "Navigator/src/DATA/weights/weights_cosmopowerspace_10000/"

    weights_TT = npzread(weights_folder*"weights_TT_lcdm.npy")
    ℓ = npzread(weights_folder*"l.npy")
    trained_emu_TT = Capse.init_emulator(NN_dict, weights_TT, Capse.SimpleChainsEmulator)
    CℓTT_emu = Capse.CℓEmulator(TrainedEmulator = trained_emu_TT, ℓgrid = ℓ,
    InMinMax = npzread(weights_folder*"inMinMax_lcdm.npy"),
    OutMinMax = npzread(weights_folder*"outMinMaxCℓTT_lcdm.npy"));

    cl=Capse.get_Cℓ([ln10As,ns,h0,ωb,ωc,τ], CℓTT_emu)[1:lmax-1]
    return cl

end


t= proba(test_map.pixels, random_noise)
chaiz2=sample(t, MH(), 100)

chainz = sample(t, HMC(0.005,10), 10)



for l = 0:length(a[i])-1
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
                alm.alm[i] ~ Normal(0, sqrt(a[l+1]))
                #alm.alm[i] = aa[i]
            else
                alm.alm[i].re ~ Normal(0, sqrt(a[l+1]/2)) #FIXME change aa to alm.alm
                alm.alm[i].im ~ Normal(0, sqrt(a[l+1]/2))
                #real~ Normal(0, sqrt(cl[l+1]/2))
               # comp ~ Normal(0, sqrt(cl[l+1]/2))
               # alm.alm[i]=complex(real,comp)
            end
        end
    end
end    


 

#function fakecurve(x::Any,A::Float64,B::Float64)
    #f= A*exp(-x*B)
    #return f
#end    

x=[1,2,3]
fakecurve.(x,3.0,3.0)
Zygote.gradient(fakecurve,x[1],3,3)
#Zygote.jacobian((x,A,B)->A.*exp.(-x.*B),[1,1,1],3,3)


function der() 
    rng = MersenneTwister(1234)
    lmax=8
   
   x=2:lmax
   A=rand(rng,Uniform(1,10))
   B=rand(rng,Uniform(0.0000001,1))

   function fakecurve(x::Any,A::Float64,B::Float64)
        f= A*exp(-x*B)
        return f
    end    
    
    test=fakecurve.(x,A,B)
    function alms(test)
        lmax = 8
        mmax = 8
        alm = Healpix.Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
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
                        gauss=Normal(0, sqrt(test[l]))
                        alm.alm[i]=rand(rng,gauss)
                    else
                        gauss=Normal(0, sqrt(test[l]/2))
                        real=rand(rng,gauss)
                        comp=rand(rng,gauss)
                        alm.alm[i]=complex(real,comp)
                    end
                end
            end
        end 
       return  alm
    end
    
    #return Zygote.jakobian(test(A,B),[A,B])[1]


    dercl=zeros(lmax-1)
    deralm=zeros(lmax-1)
    for i in 1:lmax-1
         dercl[i]= Zygote.gradient(fakecurve,x[i],A,B)[1]
    end 

   for i in 1:lmax-1
        deralm[i]=Zygote.gradient(alms,test[i])
    end 
    return test,alms(test),dercl,deralm
end

der()[2]


map=Healpix.alm2map(der()[2],8)
plot(map)



function foo(x,y)
    a = sin(x)+y
    b = 0.2+ a
    c = asin(b)
    return c
end
zz=3
yy=0

Zygote.gradient(foo, zz,yy)

function fakecurve(x,A,B)
    f= A*exp(-x*B)
    return f
end 
Zygote.gradient(fakecurve,3,3,3)

aaa=zeros(3)
for i in 1:3
   aaa[i]= Zygote.gradient(fakecurve,x[i],1,3)[1]
end
Zygote.gradient(fakecurve,x[1],1,3)[1]
aaa



function der2() 
    rng = MersenneTwister(1234)
    lmax=8
   
   x=2:lmax
   A=rand(rng,Uniform(1,10))
   B=rand(rng,Uniform(0.0000001,1))

   function fakecurve(x::Any,A::Float64,B::Float64)
        f= A*exp(-x*B)
        return f
    end    
    
    test=fakecurve.(x,A,B)
    alms=synalm(test)
    #return Zygote.jakobian(test(A,B),[A,B])[1]


    dercl=zeros(lmax-1)
    deralm=zeros(lmax-1)
    for i in 1:lmax-1
         dercl[i]= Zygote.gradient(fakecurve,x[i],A,B)[1]
    end 


    Zygote.gradient(synalm,test)



    return test,alms,dercl,deralm
end

der2()[2]


function rr(x)
    x.+1
end

xxx=[1,1,1]
Zygote.jacobian(rr,xxx)[1]
