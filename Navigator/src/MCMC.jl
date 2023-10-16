include("power_spectrum.jl")
using .power_spectrum

function make_sigmasq(cl::Vector{Float64},noise::Vector{Float64})
    l=size(cl)[1]
    sigma=zeros(l)
    for i in range(1,l)
    sigma[i]=(1/(2i+1)).*(cl[i].^2 .+ noise[i].^2)
    end
return sigma
end

