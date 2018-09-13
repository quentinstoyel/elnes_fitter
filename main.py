import numpy as np
import scipy
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import simps


class Elnes_fitter:

    # hard variables here?

    def __init__(self,expt_data):
        self.element_edge = "Lithium K Edge, 55eV"
        self.edge_location = 55
        self.plot_region = [55, 80]
        self.smoother_factor = 5 #how much to smooth the experimental data
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.cid = []
        expt_data_spectrum = Spectrum(str(expt_data))
        self.expt_loader(expt_data_spectrum) #loads the experimental data right off the bat
        return

    def __call__(self, spectrum):
        #self.broadspec_loader(spectrum)
        #self.point_getter(spectrum)
        return

    def broadspec_loader(self, spectrum):
        # Function to load a wien2k broadspec file, given case name
        data = pd.read_csv(spectrum.case_name + ".broadspec",
                           header=6, delim_whitespace=True)
        spectrum.energy_values = np.array(data['#'])
        spectrum.intensities = np.array(data['Energy'])
        spectrum.data_type = "simulation"
        spectrum.scissor_shift = self.edge_location
        spectrum.scale_factor = max(self.expt_data.intensities)/max(spectrum.intensities)

        return



    def expt_loader(self, spectrum):
        # function that loads an msa file
        data = pd.read_csv(spectrum.case_name + '.msa',
                           header=16, skipfooter=1,engine='python')
        spectrum.energy_values = np.array(data['#SPECTRUM'])
        spectrum.intensities = np.array(
            data['    : Spectral Data Starts Here'])
        spectrum.data_type = "experimental"
        self.expt_data = spectrum
        self.data_smoother(self.expt_data,self.smoother_factor)
        self.expt_data.legend_label = "Experimental Data"

        return

    def data_smoother(self, spectrum, divisor):
        # function that adds adjacent data points to smooth noise in
        # experimental data
        spectrum.energy_values = spectrum.energy_values[::divisor]
        spectrum.intensities = [sum(spectrum.intensities[
                                    i:i + divisor]) / divisor for i in range(0, len(spectrum.intensities), divisor)]

        return

    def onclick(self, event):
        #see point_getter
        ix, iy = event.xdata, event.ydata
        self.coords.append((ix, iy))
        if len(self.coords) == 2:
            self.fig.canvas.mpl_disconnect(self.cid)
            print(self.coords)
            plt.close(self.fig)

        return

    def point_getter(self, spectrum):
        #in combination with onclick, function takes a spectrum, plots it and
        self.coords = []
        self.spectrum_plotter(self.expt_data)
        self.spectrum_plotter(spectrum)
        self.plot_maker()
        self.cid = self.ax.figure.canvas.mpl_connect(
            'button_press_event', self.onclick)

        return


    def peak_align(self,spectrum):
        self.peak_finder(spectrum,self.coords[0][0])
        self.peak_finder(self.expt_data,self.coords[1][0])
        spectrum.scissor_shift = (self.expt_data.energy_values[self.expt_data.peak_index] - spectrum.energy_values[spectrum.peak_index])


        self.spectrum_plotter(self.expt_data)
        self.spectrum_plotter(spectrum)
        self.plot_maker()
        return


    def peak_finder(self, spectrum, xcoord):
        #function takes a spectrum and a peak coordinate (float) as input and sets spectrum.peak_index to the index of the maximum value within 5 ev of xcoord
        search_range = [xcoord - 2.5, xcoord + 2.5]
        search_range_index = [np.argmin(np.abs(spectrum.energy_values + spectrum.scissor_shift - search_range[
                                        0])), np.argmin(np.abs(spectrum.energy_values + spectrum.scissor_shift - search_range[1]))]
        spectrum.peak_index = np.argmax(
            spectrum.intensities[search_range_index[0]:search_range_index[1]]) + search_range_index[0]
        return

    def integrator(self, spectrum):
        return

    def plot_maker(self):
        self.ax.set_xlim(self.plot_region)
        self.ax.set_xlabel("Energy loss (eV)")
        self.ax.set_ylabel("Intensity (arb units)")
        self.ax.legend()
        self.fig.show()

        self.fig = plt.figure() #need these to reset the figure
        self.ax = self.fig.add_subplot(111)


        return

    def spectrum_plotter(self, spectrum):
        print spectrum.scissor_shift
        self.ax.plot(spectrum.energy_values + spectrum.scissor_shift,
                     spectrum.scale_factor * spectrum.intensities, label=spectrum.legend_label)


class Spectrum:

    def __init__(self, case_name):
        self.case_name = str(case_name)
        self.scissor_shift = 0
        self.scale_factor = 1
        self.legend_label = str(case_name)
        self.energy_values = np.array([])
        self.intensities = np.array([])
        self.data_type = ""
