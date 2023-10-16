module power_spectrum
export get_CMBcl,linl, get_CMBdl, get_noise_spectrum, derivativeDl
export plot_CMB, vary_amplitude,vary_H0,vary_index,vary_omegab,vary_omegac,vary_tau
export plot_derivative_of_amplitude, plot_derivative_of_H0, plot_derivative_of_index, plot_derivative_of_omegab, plot_derivative_of_omegac,plot_derivative_of_tau
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



NN_dict = JSON.parsefile("Navigator/src/nn_setup.json")


weights_folder = "Navigator/src/DATA/weights/weights_cosmopowerspace_10000/"

weights_TT = npzread(weights_folder*"weights_TT_lcdm.npy")
ℓ = npzread(weights_folder*"l.npy")
trained_emu_TT = Capse.init_emulator(NN_dict, weights_TT, Capse.SimpleChainsEmulator)
CℓTT_emu = Capse.CℓEmulator(TrainedEmulator = trained_emu_TT, ℓgrid = ℓ,
InMinMax = npzread(weights_folder*"inMinMax_lcdm.npy"),
OutMinMax = npzread(weights_folder*"outMinMaxCℓTT_lcdm.npy"));


linl = 2:5000

@doc raw"""
    get_CMBdl(;As,ns,H0,wb,wc,tau,lmax)

Calculates the normalised CMB power spectrum (``D_l=C_l * l(l+1)/2\pi``) given the cosmological parameters. 
# Arguments
- `As::Float64 = 3.043`: amplitude of primordial perturbations.
- `ns::Float64 = 0.964`: spectral index of the power spectrum.
- `H0::Float64 = 67.54`: Hubble expansion rate at current time.
- `wb::Float64 = 0.02217`: the ratio of baryonic matter (``\\Omega_b *h^2``).
- `wc::Float64 = 0.1191`: the ratio of dark matter (``\\Omega_c *h^2``).
- `tau::Float64 = 0.0571`: optical depth parameter in the reionisation era.
- `lmax::Integer = 5000`: maximum multipole moment which the function calculates.

# Returns
-`Vector{Float64}`: an array of normalised correlation coefitients [given in ``\mu K^2``] starting from l=2 to lmax.
"""
function get_CMBdl(;As=3.043,ns=0.9645,H0=67.54,wb=0.02217,wc=0.1191,tau=0.0571,lmax=5000)
    return Capse.get_Cℓ([As,ns,H0,wb,wc,tau], CℓTT_emu)[1:lmax-1]
end



@doc raw"""
    get_CMBcl(;As,ns,H0,wb,wc,tau,lmax)

Calculates CMB power spectrum (``C_l``) given the cosmological parameters. 
# Arguments
- `As::Float64 = 3.043`: amplitude of primordial perturbations.
- `ns::Float64 = 0.964`: spectral index of the power spectrum.
- `H0::Float64 = 67.54`: Hubble expansion rate at current time.
- `wb::Float64 = 0.02217`: the ratio of baryonic matter (``\\Omega_b``).
- `wc::Float64 = 0.1191`: the ratio of dark matter (``\\Omega_c``).
- `tau::Float64 = 0.0571`: optical depth parameter in the reionisation era.
- `lmax::Integer = 5000`: maximum multipole moment which the function calculates.

# Returns
-`Vector{Float64}`: an array of correlation coefitients[given in ``\mu K^2``] starting from l=2 to lmax.
"""
function get_CMBcl(;As=3.043,ns=0.9645,H0=67.54,wb=0.02217,wc=0.1191,tau=0.0571,lmax=5000)
    newl=linl[1:lmax-1]
    a= get_CMBdl(;As,ns,H0,wb,wc,tau,lmax).*(2*pi)./(newl.*(newl.+1))
    return a 
end



function plot_CMB(;As=3.043,ns=0.9645,H0=67.54,wb=0.02217,wc=0.1191,tau=0.0571,lmax=5000)
    linl = 2:lmax
    plot(linl, CMBdl(;As,ns,H0,wb,wc,tau,lmax), label="CMB power spectrum")
    xlabel!(L"l")
    ylabel!(L"$D_l$ [$\mu\mathrm{K}^2$]")
end

function vary_amplitude()
    p1=plot()
    a=palette(:YlOrRd_9, 11)
    A_s_space= 2.5:0.1:3.5
    for i in A_s_space
        linl = 2:5000
        var_amp=Capse.get_Cℓ([i,0.9645,67.54,0.02217,0.1191,0.0571], CℓTT_emu)
        plot!(p1,linl,var_amp,label="$i",palette=a)
        #plot!(xscale=:log10)
    end
    xlabel!(L"l")
    ylabel!(L"$D_l$ [$\mu\mathrm{K}^2$]") 
    p1
end


function vary_index()
    p2=plot()
    a=palette(:YlOrRd_9, 21)
    testns= 0.88:0.01:1.06
    for i in testns
        linl = 2:5000
        varns=Capse.get_Cℓ([3.043,i,67.54,0.02217,0.1191,0.0571], CℓTT_emu)
        plot!(p2,linl,varns,label="$i",palette=a)
    end
    xlabel!(L"l")
    ylabel!(L"$D_l$ [$\mu\mathrm{K}^2$]")
    p2
end

function vary_H0()
    p3=plot()
    testH0= 60:1:80
    a=palette(:YlOrRd_9, 21)
    for i in testH0
        linl = 2:5000
        varH0=Capse.get_Cℓ([3.043,0.9645,i,0.02217,0.1191,0.0571], CℓTT_emu)
        plot!(p3,linl,varH0,label="$i",palette=a)
    end
    xlabel!(L"l")
    ylabel!(L"$D_l$ [$\mu\mathrm{K}^2$]")
    p3
end

function vary_omegab()
    p4=plot()
    testomb= 0.019:0.001:0.03
    a=palette(:YlOrRd_9, 14)
    for i in testomb
        linl = 2:5000
        varomb=Capse.get_Cℓ([3.043,0.9645,67.54,i,0.1191,0.0571], CℓTT_emu)
        plot!(p4,linl,varomb,label="$i",palette=a)
    end
    xlabel!(L"l")
    ylabel!(L"$D_l$ [$\mu\mathrm{K}^2$]")
    p4

end

function vary_omegac()
    p5=plot()
    testomc= 0.08:0.01:0.2
    a=palette(:YlOrRd_9, 14)
    for i in testomc
        linl = 2:5000
        varomb=Capse.get_Cℓ([3.043,0.9645,67.54,0.02217,i,0.0571], CℓTT_emu)
        plot!(p5,linl,varomb,label="$i",palette=a)
    end
    xlabel!(L"l")
    ylabel!(L"$D_l$ [$\mu\mathrm{K}^2$]")
p5
end

function vary_tau()
    p6=plot()
    testtau= 0.05:0.003:0.086
    a=palette(:YlOrRd_9,14)
    for i in testtau
    linl = 2:5000
    vartau=Capse.get_Cℓ([3.043,0.9645,67.54,0.02217,0.1191,i], CℓTT_emu)
    plot!(p6,linl,vartau,label="$i",palette=a)
    end
    xlabel!(L"l")
    ylabel!(L"$D_l$ [$\mu\mathrm{K}^2$]")
    p6
end

theory(θ) = Capse.get_Cℓ(θ, CℓTT_emu)
x=[3.043,0.9645,67.54,0.02217,0.1191,0.0571]
derivativeDl=ForwardDiff.jacobian(theory, x)

function plot_derivative_of_amplitude()
    plot(linl,derivativeDl[:,1],label=L"derivative of $A_s$")
    xlabel!(L"l")
    ylabel!(L"dD_l / d\theta \, [\mu\mathrm{K}^2]")    
end

function plot_derivative_of_index()
    plot(linl,derivativeDl[:,2],label=L"derivative of $n_s$")
    xlabel!(L"l")
    ylabel!(L"dD_l / d\theta \, [\mu\mathrm{K}^2]")    
end

function plot_derivative_of_H0()
plot(linl,derivativeDl[:,3],label=L"derivative of $H_0$")
xlabel!(L"l")
ylabel!(L"dD_l / d\theta \, [\mu\mathrm{K}^2]")    
end

function plot_derivative_of_omegab()
    plot(linl,derivativeDl[:,4],label=L"derivative of  $\omega_b$")
    xlabel!(L"l")
    ylabel!(L"dD_l / d\theta \, [\mu\mathrm{K}^2]")    
end

function plot_derivative_of_omegac()
plot(linl,derivativeDl[:,5],label=L"derivative of $ omega_b$")
xlabel!(L"l")
ylabel!(L"dD_l / d\theta \, [\mu\mathrm{K}^2]")    
end

function plot_derivative_of_tau()
plot(linl,derivativeDl[:,6],label=L"derivative of $\tau$")
xlabel!(L"l")
ylabel!(L"dD_l / d\theta \, [\mu\mathrm{K}^2]")    
end


@doc raw"""
    get_noise_spectrum(;lmax,sigma,Nside)

Gives the theoretical normalized noise power spectrum for gaussian noise with mean 0 and variance `sigma^2`.
#Arguments
-`lmax::Integer = 5000`:  maximum multipole moment which the function calculates.
-`sigma::Float64 = 10`: dispersion of Gaussian distribution from which the noise is drawn
-`Nside::Integer, must be a power of 2, default=2048`: the Nside parameter related to a number of pixels of a Healpix map
#Returns
-`Vector{Float64}`: an array of the theoretical noise power spectrum [given in ``\mu K^2``] starting from l=2 to lmax.
"""
function get_noise_spectrum(;lmax::Integer=5000,sigma::Float64=10.0,Nside::Integer=2048)
    noisearray=zeros(lmax-1)
    for i in 1:lmax-1
        noisearray[i]=sigma^2*(4*pi)/nside2npix(Nside)
    end
    newl=linl[1:lmax-1]
    noise=noisearray./((2*pi)./(newl.*(newl.+1)))
    return noise
end

end
