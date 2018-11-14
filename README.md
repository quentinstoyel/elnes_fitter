# elnes_fitter
OOP python code to fit elnes simulation to expt

run in directory with .broadspec/.msa files.

initialize the various simluation spectra needed ($ represents a python console): 

eg $LiF  = Spectrum("LiF_hole_sc")

Then setup the fitter with the expt data, and load the sims:

$expt_data = Elnes_fitter("mix_deconv")
$expt_data.broadspec_loader(LiF)


Then fit and plot the spectra with:

$y.spectrum_fitter()
$y.final_spectrum_plotter()


For the example files, there are three parts to the experimental spectrum, Li, LiF, and Li2CO3. 
To fit these: 

from main import * 
LiF = Spectrum("LiF_hole_sc") 
Li = Spectrum("Li_partial_hole_sc")
LCO = Spectrum("Li2CO3")

data_fitter = Elnes_fitter("mix_deconv")

data_fitter.broadspec_loader(LiF)
data_fitter.broadspec_loader(Li)
data_fitter.broadspec_loader(LCO) 
data_fitter.spectrum_fitter()
data_fitter.final_spectrum_plotter()


As you might expect, this code was made by a grad student and I bear no responsibility for what you do with it, nor offer any gaurentees that it will work. 
