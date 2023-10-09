# Navigator
Creates and reads CMB maps.


# Installation

Step 1: Download the package using `git clone https://github.com/ivansladoljev/Navigator.git` in the directory you are using.

Step2: Then create a folder `DATA` in the `Navigator/src/` folder. In there you will have to downoald the chain weights for Capse emulator from `https://zenodo.org/record/8187935/files/chains_weights.zip?download=1` and add them to the `DATA` folder. Everything should be up and running then.



# Introduction 

With this program you can create your custom CMB temperature maps. The program runs with underlying Healpix and Capse's emulator frameworks.
To generate the CMB power spectrum use the CMBcl() or CMBdl() functions which uses Capse's trained neural network to generate a set of CMB power spectra, given the inital cosmological parameters. 
From the generated set of $`C_l`$s or $`D_l`$s you can make your own custom map of the CMB, with differing levels of pixelization, using the Map() function. You can also make the noise mape with Noise(), or add the gaussian noise to you existing CMB map with NoisyMap() function.
If you already have a Healpix map, you can infere the $`C_l`$s or the $`D_l`$s from it with clobserved()/dlobserved() functions. If you want to remove the theoretical noise from your obsevations and get the pure power spectrum, use the clestimate()/dlestimate functions instead.
