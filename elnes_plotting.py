# coding: utf-8
from main import Spectrum, Elnes_fitte
from main import Spectrum, Elnes_fitter
x = Spectrum("Li2CO3")
y=Elnes_fitter()
y.broadspec_loader(x)
y.Spectrum_plotter(x)
