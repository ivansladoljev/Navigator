# Navigator
Creates and reads CMB maps. The program runs with underlying Healpix and Capse's emulator frameworks.


# Installation

Step 1: Download the package using `git clone https://github.com/ivansladoljev/Navigator.git` in the directory you are using.

Step 2: Then create a folder `DATA` in the `Navigator/src/` folder (this can be done by going to the directory where you dowloaded the Navigator file, typing `cd Navigator/src` in the terminal to go into the right folder and creating DATA using  `mkdir DATA` command to create the folder).

Step 4: In there you will have to download the chain weights for Capse emulator from `https://zenodo.org/record/8187935/files/chains_weights.zip?download=1`, unzip them and add them to the `DATA` folder
(this can be done by `wget https://zenodo.org/record/8187935/files/chains_weights.zip` command and unziping them by writing `unzip chains_weights.zip` in the terminal). 


Step 3: Now you successfully downloaded everything you need to operate the code. If you want to use it, just write  `using Pkg, Pkg.activate("Navigator"), include("Navigator"), using Navigator` and you are good to go.


# Introduction 

With this program you can create your custom CMB temperature maps. 


To generate the CMB power spectrum use the get_CMBcl() or get_CMBdl() functions which uses Capse's trained neural network to generate a set of CMB power spectra, given the inital cosmological parameters. 
From the generated set of $`C_l`$s or $`D_l`$s you can make your own custom map of the CMB, with differing levels of pixelization, using the make_map() function. You can also make the noise map with make_noise(), or add the gaussian noise to you existing CMB map with make_noisymap() function.
If you already have a Healpix map, you can infere the $`C_l`$s or the $`D_l`$s from it with get_clobserved()/get_dlobserved() functions. If you want to remove the theoretical noise from your obsevations and get the pure power spectrum, use the get_clestimate()/get_dlestimate functions instead.
