# elnes_fitter
OOP python code to fit elnes simulation to expt

run in directory with .broadspec/.msa files.

initialize the various simluation spectra needed: 

eg $LiF  = Spectrum("LiF_hole_sc")

Then setup the fitter with the expt data, and load the sims:

$y = Elnes_fitter("mix_deconv")
$y.broadspec_loader(LiF)


Then fit and plot the spectra with:

$y.spectrum_fitter()
$y.final_spectrum_plotter()

As you might expect, this code was made by a grad student and I bear no responsibility for what you do with it, nor offer any gaurentees that it will work. 
