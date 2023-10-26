



# Navigator
Creates and reads CMB maps. The program runs with underlying Healpix and Capse's emulator frameworks.


# Installation

Step 1: Download the package using `git clone https://github.com/ivansladoljev/Navigator.git` in the directory you are using.

Step 2: Then create a folder `DATA` in the `Navigator/src/` folder (this can be done by going to the directory where you dowloaded the Navigator file, typing `cd Navigator/src` in the terminal to go into the right folder and creating DATA using  `mkdir DATA` command to create the folder).

Step 3: In there you will have to download the chain weights for Capse emulator from `https://zenodo.org/record/8187935/files/chains_weights.zip?download=1`, unzip them and add them to the `DATA` folder
(this can be done by `wget https://zenodo.org/record/8187935/files/chains_weights.zip` command and unziping them by writing `unzip chains_weights.zip` in the terminal). 

Step 4: Now you successfully downloaded everything you need to operate the code. If you want to use it, just write  `using Pkg, Pkg.activate("Navigator"), include("Navigator"), using Navigator` and you are good to go.


# Introduction 

With this program you can create your custom CMB temperature maps and analyse them. 

## CMB power spectrum
To generate the CMB power spectrum use the **get_CMBcl()** or **get_CMBdl()** functions which uses Capse's trained neural network to generate a set of CMB power spectra, given the inital cosmological parameters. 

### get_CMBcl(;As,ns,H0,wb,wc,tau,lmax)

Calculates CMB power spectrum ($`C_l`$) given the cosmological parameters. 

Arguments
- `As::Float64 = 3.043`: amplitude of primordial perturbations.
- `ns::Float64 = 0.964`: spectral index of the power spectrum.
- `H0::Float64 = 67.54`: Hubble expansion rate at current time.
- `wb::Float64 = 0.02217`: the ratio of baryonic matter ($`\Omega_b`$).
- `wc::Float64 = 0.1191`: the ratio of dark matter ($`\Omega_c`$).
- `tau::Float64 = 0.0571`: optical depth parameter in the reionisation era.
- `lmax::Integer = 5000`: maximum multipole moment which the function calculates.

 Returns
-`Vector{Float64}`: an array of correlation coefitients[given in $`\mu K^2`$] starting from l=2 to lmax.

### get_CMBdl(;As,ns,H0,wb,wc,tau,lmax)

Calculates the normalised CMB power spectrum ($` D_l= C_l * l(l+1)/2\pi`$) given the cosmological parameters. 
 
Arguments
- `As::Float64 = 3.043`: amplitude of primordial perturbations.
- `ns::Float64 = 0.964`: spectral index of the power spectrum.
- `H0::Float64 = 67.54`: Hubble expansion rate at current time.
- `wb::Float64 = 0.02217`: the ratio of baryonic matter ($`\Omega_b *h^2`$).
- `wc::Float64 = 0.1191`: the ratio of dark matter ($`\Omega_c *h^2`$).
- `tau::Float64 = 0.0571`: optical depth parameter in the reionisation era.
- `lmax::Integer = 5000`: maximum multipole moment which the function calculates.

Returns
-`Vector{Float64}`: an array of normalised correlation coefitients [given in $` \mu K^2`$] starting from l=2 to lmax.


### get_noise_spectrum(;lmax,sigma,Nside)

Gives the theoretical normalized noise power spectrum for gaussian noise with mean 0 and variance $`sigma^2`$.

Arguments
-`lmax::Integer = 5000`:  maximum multipole moment which the function calculates.
-`sigma::Float64 = 10`: dispersion of Gaussian distribution from which the noise is drawn
-`Nside::Integer, must be a power of 2, default=2048`: the Nside parameter related to a number of pixels of a Healpix map

Returns
-`Vector{Float64}`: an array of the theoretical noise power spectrum [given in $`\mu K^2`$] starting from l=2 to lmax.

## Making maps
From the generated set of $`C_l`$s or $`D_l`$s you can make your own custom map of the CMB, with differing levels of pixelization, using the **make_map()** function. You can also make the noise map with **make_noise()**, or add the gaussian noise to you existing CMB map with **make_noisymap()** function.

###  make_map(Nside;As,ns,H0,wb,wc,tau,lmax,seed) 

Gives the map of the CMB power spectrum.

Arguments
- `Nside::Integer, must be a power of 2``: the Nside parameter for the number of pixels in the map.

- `As::Float64 = 3.043`: amplitude of primordial perturbations.

- `ns::Float64 = 0.964`: spectral index of the power spectrum.

- `H0::Float64 = 67.54`: Hubble expansion rate at current time.

- `wb::Float64 = 0.02217`: the ratio of baryonic matter ($`\Omega_b *h^2`$).

- `wc::Float64 = 0.1191`: the ratio of dark matter ($`\Omega_c *h^2`$).

- `tau::Float64 = 0.0571`: optical depth parameter in the reionisation era.

- `lmax::Integer = 2*Nside`: maximum multipole moment of the power spectrum for which the map is calculated.

- `seed::Integer = 1234`: fixes the seed for the MersenneTwister random generator

Returns
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: returns the Healpix map of the CMB in RingOrdering.

Alternative imput
    make_map(cl,Nside)

Alternative arguments
- `cl::Vector{Float64}`: a set of power spectrum $`C_l`$s for which to calculate a map.

- `Nside::Integer, must be a power of 2 `: the Nside parameter for the number of pixels in the map.

- `seed::Integer = 1234`: fixes the seed for the MersenneTwister random generator

Returns
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: returns the Healpix map of the CMB in RingOrdering.

###  make_noise(Nside;sigma)

Gives the map of the Gaussian noise with mean zero and dispersion $`sigma`$.

Arguments
- `Nside::Integer, must be a power of 2 `: the Nside parameter for the number of pixels in the map.

- `sigma::Float64 = 10`: dispersion of Gaussian distribution from which the noise is drawn.

- `seed::Integer = 1234`: fixes the seed for the MersenneTwister random generator

Returns
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: returns the Healpix map of the CMB in RingOrdering.

### make_noisymap(Nside;As,ns,H0,wb,wc,tau,lmax,sigma)

Gives the map of the CMB power spectrum with noise added.

Arguments
- `Nside::Integer, must be a power of 2 = 2048`: the Nside parameter for the number of pixels in the map.

- `As::Float64 = 3.043`: amplitude of primordial perturbations.

- `ns::Float64 = 0.964`: spectral index of the power spectrum.

- `H0::Float64 = 67.54`: Hubble expansion rate at current time.

- `wb::Float64 = 0.02217`: the ratio of baryonic matter ($`\Omega_b *h^2`$).

- `wc::Float64 = 0.1191`: the ratio of dark matter ($`\Omega_c *h^2`$).

- `tau::Float64 = 0.0571`: optical depth parameter in the reionisation era.

- `lmax::Integer = 5000`: maximum multipole moment of the power spectrum for which the map is calculated.

- `sigma::Float64 = 10`: dispersion of Gaussian distribution from which the noise is drawn.

- `seed::Integer = 1234`: fixes the seed for the MersenneTwister random generator

Returns
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: returns the Healpix map of the CMB in RingOrdering with npise added.

Alternative imput
    make_noisymap(cl,Nside;sigma)

Alternative arguments
- `cl::Vector{Float64}`: a set of power spectrum $`C_l`$s for which to calculate a map.

- `Nside::Integer, must be a power of 2 `: the Nside parameter for the number of pixels in the map.

- `sigma::Float64 = 10`: dispersion of Gaussian distribution from which the noise is drawn.


- `seed::Integer = 1234`: fixes the seed for the MersenneTwister random generator

Returns
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: returns the Healpix map of the CMB in RingOrdering with npise added.


## Analysing maps

If you already have a Healpix map, you can infere the $`C_l`$s or the $`D_l`$s from it with **get_clobserved()**/**get_dlobserved()** functions. If you want to remove the theoretical noise from your obsevations and get the pure power spectrum, use the **get_clestimate()**/**get_dlestimate()** functions instead.

### get_clobserved(map;lmax)

Gives the observed power spectrum coeficients ``C_l`` from the given map.

Arguments
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: the map from which to infere the power spectrum. Must be in RingOrder.

-`lmax::Integer = 2*Nside`:  maximum multipole moment which the function calculates. Maximum default is 2*Nside parameter of the map.

Returns
- `Vector{Float64}`: an array of infered power spectrum coeficients $`C_l`$.

### get_dlobserved(map;lmax)

Gives the observed power spectrum coeficients $`D_l=C_l * l(l+1)/2\pi`$ from the given map.

Arguments
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: the map from which to infere the power spectrum. Must be in RingOrder.

-`lmax::Integer = 2*Nside`:  maximum multipole moment which the function calculates. Maximum default is 2*Nside parameter of the map.

Returns
- `Vector{Float64}`: an array of infered power spectrum coeficients $`D_l`$.

### get_clestimate(map;lmax,sigma)

Gives the infered power spectrum coeficients $`C_l`$ (with noise removed) from the given map.

Arguments
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: the map from which to infere the power spectrum. Must be in RingOrder.

-`lmax::Integer = 2*Nside`:  maximum multipole moment which the function calculates. Maximum default is 2*Nside parameter of the map.

-`sigma::Float64 = 10`: dispersion of Gaussian distribution from which the theoretical noise is drawn

Returns
- `Vector{Float64}`: an array of infered power spectrum coeficients $`C_l`$ with noise removed.

### get_dlestimate(map;lmax,sigma)

Gives the infered power spectrum coeficients $`D_l=C_l * l(l+1)/2\pi`$ (with noise removed) from the given map.

Arguments
-`Healpix.HealpixMap{Float64, Healpix.RingOrder}`: the map from which to infere the power spectrum. Must be in RingOrder.

-`lmax::Integer = 2*Nside`:  maximum multipole moment which the function calculates. Maximum default is 2*Nside parameter of the map.

-`sigma::Float64 = 10.0`: dispersion of Gaussian distribution from which the theoretical noise is drawn.

Returns
- `Vector{Float64}`: an array of infered power spectrum coeficients $`D_l`$ with noise removed.
