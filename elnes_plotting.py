# coding: utf-8
from main import Elnes_fitter, Spectrum
LiF  = Spectrum("LiF_hole_sc")
Li = Spectrum("Li_partial_hole_sc")
y = Elnes_fitter("mix_deconv")
y.broadspec_loader(LiF)
y.point_getter(LiF)
y.point_getter(LiF)
y.peak_align(LiF)
from scipy.optimize import minimize
y.point_getter(Li)
y.point_getter(Li)
y.spectrum_plotter(Li)
y.point_getter(Li)
y.spectrum_plotter(Li)
y.plot_maker()
y.broadspec_loader(Li)
y.point_getter(Li)
y.point_getter(Li)
y.peak_align(Li)
fit = minimize(y.minimizer, [LiF.scale_factor, Li.scale_factor,y.expt_data.v_offset],method='Powell')
fit
LiF.scale_factor = fit.x[0]
Li.scale_factor = fit.x[1]
y.expt_data.v_offset = fit.x[2]
y.spectrum_plotter(Li)
y.spectrum_plotter(LiF)
y.spectrum_plotter(y.expt_data)
y.plot_maker
