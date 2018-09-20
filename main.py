import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd
from scipy.integrate import simps
from scipy.optimize import minimize


class Elnes_fitter:
    #class to manipulate various elnes spc
    def __init__(self, expt_data):
        # Elnes fitter needs to be called with a string corresponding to the
        # .msa file
        self.element_edge = "Lithium K Edge, 55eV"
        self.edge_location = 55
        self.plot_region = [55, 70]
        self.smoother_factor = 40  # how much to smooth the experimental data
        self.vector_length = 2000  # how long we want the vectors
        # initializing figures
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.cid = []
        expt_data_spectrum = Spectrum(str(expt_data))
        # loads the experimental data right off the bat
        self.expt_loader(expt_data_spectrum)
        self.simulation_list = []
        return

    def __call__(self, spectrum):
        # self.broadspec_loader(spectrum)
        # self.point_getter(spectrum)
        return

    def broadspec_loader(self, spectrum):
        # Function loads a wien2k broadspec file, given case name
        # pulls up two plots comparing it to the expt data_type
        # allows scissor offset to be set
        data = pd.read_csv(spectrum.case_name + ".broadspec",
                           header=6, delim_whitespace=True)
        spectrum.energy_values = np.array(data['#'])
        spectrum.intensities = np.array(data['Energy'])
        spectrum.data_type = "simulation"
        spectrum.scissor_shift = self.edge_location
        spectrum.scale_factor = max(
            self.expt_data.intensities) / max(spectrum.intensities)  # normalizes the spectra to the expt data
        # adds the spectra to the sim list
        self.simulation_list.append(spectrum)
        self.point_getter(spectrum)
        self.point_getter(spectrum)  # still need to call this twice
        self.peak_align(spectrum)
        return

    def expt_loader(self, spectrum):
        # function that loads an msa file and saves the spectrum to expt data
        # also sets the length of energy/intensity vectors to correct value
        # also runs rolling average of the data to minimize noise
        data = pd.read_csv(spectrum.case_name + '.msa',
                           header=16, skipfooter=1, engine='python')
        spectrum.energy_values = np.array(data['#SPECTRUM'])
        spectrum.intensities = np.array(
            data['    : Spectral Data Starts Here'])
        new_energy_array = np.linspace(spectrum.energy_values[self.nearest_index_value(spectrum.energy_values, self.plot_region[
                                       0])], spectrum.energy_values[self.nearest_index_value(spectrum.energy_values, self.plot_region[1])], self.vector_length)
        # interpolating to fill all the data points
        spectrum.intensities = scipy.interpolate.spline(
            spectrum.energy_values, spectrum.intensities, new_energy_array) - np.average(spectrum.intensities[0:10])
        spectrum.energy_values = new_energy_array
        spectrum.data_type = "experimental"
        self.expt_data = spectrum
        # smooths the data, as desired
        self.data_smoother(self.expt_data, self.smoother_factor)
        self.expt_data.legend_label = "Experimental Data"

        return

    def data_smoother(self, spectrum, divisor):
        # function that adds adjacent data points to smooth noise in
        # experimental data, and then interpolates the data back to its
        # original length

        #spectrum.energy_values = spectrum.energy_values[::divisor]
        y_values = [sum(spectrum.intensities[
            i:i + divisor]) / divisor for i in range(0, len(spectrum.intensities), divisor)]

        x_values = spectrum.energy_values[::divisor]

        # generates a function/mapping that returns the interpolated value
        splined_mapping = scipy.interpolate.UnivariateSpline(
            x_values, y_values)
        # splined_mapping.set_smoothing_factor(0.5)
        # calling function with dataset of appropriate length sets the
        # intensities length correctly
        spectrum.intensities = splined_mapping(spectrum.energy_values)

        return

    def onclick(self, event):
        # see point_getter
        ix, iy = event.xdata, event.ydata
        self.coords.append((ix, iy))
        if len(self.coords) == 2:
            self.fig.canvas.mpl_disconnect(self.cid)
            plt.close("all")

        return

    def point_getter(self, spectrum):
        # in combination with onclick, function takes a spectrum, plots it, and
        # waits for it to be clicked twice, and then returns coordinates of
        # clicks into coords
        self.coords = []
        self.spectrum_plotter(self.expt_data)
        self.spectrum_plotter(spectrum)
        self.plot_maker()
        self.cid = self.ax.figure.canvas.mpl_connect(
            'button_press_event', self.onclick)
        # python 2 madness for python3, use input
        raw_input(
            "click graph twice (simulation peak, then corrseponding expt peak) and press enter to continue")
        return

    def peak_align(self, spectrum):
        # aligns peak to appropriate location, sets its data length correctly
        # requires point getter to have been run and coords to have non zero length
        # as well
        self.peak_finder(spectrum, self.coords[0][0])
        self.peak_finder(self.expt_data, self.coords[1][0])
        spectrum.scissor_shift = (self.expt_data.energy_values[
                                  self.expt_data.peak_index] - spectrum.energy_values[spectrum.peak_index])
        new_energy_array = np.linspace(spectrum.energy_values[self.nearest_index_value(spectrum.energy_values + spectrum.scissor_shift, self.plot_region[
                                       0])], spectrum.energy_values[self.nearest_index_value(spectrum.energy_values + spectrum.scissor_shift, self.plot_region[1])], self.vector_length)
        spectrum.intensities = scipy.interpolate.spline(
            spectrum.energy_values, spectrum.intensities, new_energy_array)
        spectrum.energy_values = new_energy_array

        self.spectrum_plotter(self.expt_data)
        self.spectrum_plotter(spectrum)
        self.plot_maker()
        return

    def peak_finder(self, spectrum, xcoord):
        # function takes a spectrum and a peak coordinate (float) as input and
        # sets spectrum.peak_index to the index of the maximum value within 5
        # ev of xcoord
        search_range = [xcoord - 2.5, xcoord + 2.5]
        search_range_index = [self.nearest_index_value(spectrum.energy_values + spectrum.scissor_shift, search_range[
            0]), self.nearest_index_value(spectrum.energy_values + spectrum.scissor_shift, search_range[1])]
        spectrum.peak_index = np.argmax(
            spectrum.intensities[search_range_index[0]:search_range_index[1]]) + search_range_index[0]
        return

    def nearest_index_value(self, ordered_list, search_value):
        # given an ordered list (all increasing) and a value, this helper
        # function returns the index of the nearst value in the list to the
        # search value
        niv = np.argmin(np.abs(ordered_list - search_value))
        return niv

    def final_spectrum_plotter(self):
        # once all spectras are aligned and fitted, this function plots the
        # final aligned sum, the expt data and each individual spectrum
        final_intensities = sum(
            sim.intensities * sim.scale_factor for sim in self.simulation_list)
        self.spectrum_plotter(self.expt_data)
        for sim in self.simulation_list:
            self.spectrum_plotter(sim)
        plt.plot(self.expt_data.energy_values,
                 final_intensities, label="Fitted Spectra")
        self.plot_maker()

    def spectrum_fitter(self):
        # function fits all spectra loaded into the simulation list to the expt
        # data, by integrating the area under the cures
        guess = np.array([0])  # sets virtical offset guess to zero
        for sim in self.simulation_list:
            guess = np.append(guess, sim.scale_factor)

        fit = minimize(self.minimizer, guess, method='Powell')
        self.expt_data.intensities = self.expt_data.intensities + \
            fit.x[0]  # setting virtical v_offset
        scale_factor_index = 1
        for sim in self.simulation_list:
            sim.scale_factor = fit.x[scale_factor_index]
            scale_factor_index += 1

    def minimizer(self, scale_factor_array):
        # function to be minimized, input is a vector, satrting with the virtical offset, and followed by the scale
        # factors of each component
        # returns the absolute value of the difference between the sum of all
        # simulated curves and the experimental data
        expt_offset = scale_factor_array[0]
        scale_factor_array = scale_factor_array[1:]
        #*args will be the fit parameters, function to be minimized
        x_int_values = self.expt_data.energy_values
        y_int_values = self.expt_data.intensities + expt_offset - \
            np.dot(np.transpose(np.array(
                [sim.intensities for sim in self.simulation_list])), scale_factor_array)

        integral_value = simps(np.abs(y_int_values), x_int_values)
        return integral_value


    def plot_maker(self):
        #sets up all of the aesthetics of plots and generates the legend
        sefl.ax.title(self.element_edge)
        self.ax.set_xlim(self.plot_region)
        self.ax.set_xlabel("Energy loss (eV)")
        self.ax.set_ylabel("Intensity (arb units)")
        self.ax.legend()
        self.fig.show()

        self.fig = plt.figure()  # need these to reset the figure
        self.ax = self.fig.add_subplot(111)

        return

    def spectrum_plotter(self, spectrum):
        self.ax.plot(spectrum.energy_values + spectrum.scissor_shift,
                     spectrum.scale_factor * spectrum.intensities + spectrum.v_offset, label=spectrum.legend_label)


class Spectrum:
    #class for the various elnes spectra
    def __init__(self, case_name):
        self.case_name = str(case_name)
        self.scissor_shift = 0
        self.v_offset = 0
        self.scale_factor = 1
        self.legend_label = str(case_name)
        self.energy_values = np.array([])
        self.intensities = np.array([])
        self.data_type = ""
