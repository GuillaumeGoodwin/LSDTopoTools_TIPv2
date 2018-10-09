"""
    This package processes topographic data in order to extract marsh platforms
    Guillaume C.H. Goodwin
    Released unedr GPL3
"""

# Load useful Python packages
import os
import sys
import numpy as np

#from osgeo import gdal, osr, gdalconst
import matplotlib.pyplot as plt
#from osgeo.gdalconst import *
import cPickle as pickle

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D

import timeit


import copy

import pandas as bb


import LSDMarshPlatform_functions as fct

##########################################################################################################
##########################################################################################################
class Land_surface (np.ndarray):
    def __new__ (Land_surface, x_length, y_length):
        print 'In __new__ with class %s' % Land_surface
        return np.ndarray.__new__(Land_surface, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length
        self[np.isnan(self)] = 0
        self = 0 * self

    # We have no copy method because the ndarray.copy method is good enough

    def set_attribute (self, select_array, select_value, attribute_array, Nodata_value, classification = False):
        new_array = self.copy()
        Select = np.where (select_array >= select_value)
        Discard = np.where (np.logical_and(select_array < select_value, select_array != Nodata_value))
        Nodata = np.where (attribute_array == Nodata_value)
        new_array[Discard] = 0
        new_array[Nodata] = Nodata_value
        if classification == False:
            new_array[Select] = attribute_array [Select]
        else:
            new_array[Select] = select_value
        return new_array


    def label_connected (self, Nodata_value):
        new_array = self.copy()
        Ref = self.copy()
        import scipy.ndimage as scim
        array, numfeat = scim.label(self)
        for value in np.arange(1, np.amax(array)+1, 1):
            line = np.where(array == value)
            new_array[line] = value
            new_array[Ref == Nodata_value] = Nodata_value
        return new_array


    def extract_outlines (self):
        print 'extracting outlines'
        M_labels = range (int(np.amin(self[self>0])), int(np.amax(self[self>0]))+1, 1)
        Outlines = Polyline()
        for M in M_labels:
            print '... Label ', M
            M_outlines = fct.Pandas_outline(self,M,1)
            if len(M_outlines) > 0:
                if M_outlines.size > 0:
                    print '.... Appending label'
                    Outlines.append(M_outlines)
        return Outlines



    def relative_relief (self, Nodata_value):
        # We calculate the relative relief of the DEM to have values of elevation between 0 and 1
        Relief = self-np.amin(self[self > Nodata_value])
        Rel_relief = Relief/np.amax(Relief)
        Rel_relief[self <= Nodata_value] = Nodata_value

        return Rel_relief





    def define_search_space (self, search_arr1, search_arr2, Nodata_value, opt):
        """
       This function defines a search space (Search_space) within a 2-D array, based on the combined values of 2 2-D arrays (DEM and Slope) of the same dimensions. It defines the threshold for the selection of the search space according to a threshold value (opt). It is set to ignore elements with a specific value (Nodata_value).
        Args:
            DEM (2D numpy array): a 2-D array (here a DEM) used as a first condition for the definition of the search space
            Slope (2D numpy array): a 2-D array (here a DEM) used as a second condition for the definition of the search space
            Nodata_value (float): The value for ignored elements
            opt (float): the value of the threshold for the selection of the search space

        Returns:
            Search_space (2D numpy array): The resulting search space array. Search_space has a value of 0 for non-selected elements and 1 for selected elements.
            Crossover (2D numpy array): The array resulting of the multiplication of relative slope and relative relief.
            bins (1D array): the value bins for the Crossover array
            hist (1D array): the value hist for the Crossover array
            Inflecion_point(float): the value of the threshold for the search space selection.

        Author: GCHG
        """

        print 'Generating a search space ...'
        SS_object = Search_space(self.shape[0], self.shape[1]); SS_object = 0 * SS_object

        Rr1 = search_arr1.relative_relief(Nodata_value)
        Rr2 = search_arr2.relative_relief(Nodata_value)

        # This is one type of option for search spaces. We can probs make more later
        Product = Rr1 * Rr2
        Product[self <= Nodata_value] = Nodata_value

        bins, hist, step = Product.PDF(Nodata_value, 100)

        Inflexion_point = 0

        # We now find the slope of that curve
        hist_der = fct.Derivate(bins, hist)

        # If the slope gets above the -1 threshold, now that we have hit the closest point to the origin.
        # We call it the inflexion point even though it's not really an inflexion point.
        for j in range(1, len(hist)-1, 1):
            if hist_der[j] < opt and hist_der[j+1] >= opt:
                Inflexion_point = bins[j]

        # Points within the search space should have a Crossover value above the inflexion point
        Search = np.where(Product > Inflexion_point)
        SS_object[Search] = 1

        # We get rid of the borders of the DEM because otherwise it will be difficult to work with the smaller slope array
        SS_object[0,:] = 0; SS_object[-1,:] = 0; SS_object[:,0] = 0; SS_object[:,-1] = 0

        # And update the search locations for the shaved edges
        Search = np.where(SS_object == 1)

        # If this happens, your landscape is weird
        if np.amax(SS_object) == 0:
            print; print " ... Your search space is empty! Are you sure there's a marsh platform here?";  print
            quit()

        return SS_object




    def PDF (self, Nodata_value, numbins):
        Data1D = self.ravel()

        Max_distribution = max(Data1D)
        if len(Data1D[Data1D>Nodata_value]) == 0:
            Min_distribution = -1
        else:
            Min_distribution = min(Data1D[Data1D>Nodata_value])

        bin_size = (Max_distribution - Min_distribution) / numbins

        X_values = np.arange(Min_distribution, Max_distribution, bin_size)

        hist, bins = np.histogram (Data1D, X_values, density=True)
        hist=hist/sum(hist)
        bins=bins[:-1]

        return bins, hist, bin_size






    def plot_map (self, background, save_dir, figname, title, Nodata_value):
        print ' Plotplotplotplot'
        twin  = self.copy()

        fig_height = min(np.floor(twin.shape[1])/5, 50)
        fig_width = min(np.floor(twin.shape[1])/5, 50)

        fig=plt.figure(title, facecolor='White',figsize=[fig_height,fig_width])

        ax1 = plt.subplot2grid((2,1),(0,0),colspan=1, rowspan=1)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')

        Vmin = min(np.amin(twin[twin!=Nodata_value])*0.95, np.amin(twin[twin!=Nodata_value])*1.05)
        Vmax = max(np.amax(twin)*0.95, np.amax(twin)*1.05)

        Map = ax1.imshow(background, interpolation='None', cmap=plt.cm.binary, vmin=np.amin(background[background!=Nodata_value]), vmax=np.amax(background), alpha = 0.6)
        Map = ax1.imshow(twin, interpolation='None', cmap=plt.cm.gist_earth, vmin=Vmin, vmax=Vmax, alpha = 0.6)
        ax2 = fig.add_axes([0.1, 0.98, 0.85, 0.02])
        scheme = plt.cm.gist_earth; norm = colors.Normalize(vmin=Vmin, vmax=Vmax)
        cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)


        #Now plot map properties
        ax2 = plt.subplot2grid((2,1),(1,0),colspan=1, rowspan=1)
        ax2.tick_params(axis='map value', colors='black')
        ax2.tick_params(axis='pdf', colors='black')

        bins, hist, step = self.PDF(Nodata_value, 100)

        hist[0] = 0

        ax2.plot(bins, hist, 'k', linewidth = 3)


        plt.savefig(save_dir+figname+'.png')





    def label_connected (self, Nodata_value):
        new_array = self.copy()
        Ref = self.copy()
        import scipy.ndimage as scim
        array, numfeat = scim.label(self)
        for value in np.arange(1, np.amax(array)+1, 1):
            line = np.where(array == value)
            new_array[line] = value
            new_array[Ref == Nodata_value] = Nodata_value
        return new_array






##########################################################################################################
##########################################################################################################
class Search_space (Land_surface):
    def __new__ (Search_space, x_length, y_length):
        print 'In __new__ with class %s' % Search_space
        return np.ndarray.__new__(Search_space, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length
        self[np.isnan(self)] = 0
        self = 0 * self


    def extract_full_scarps (self, DEM_arr, slope_arr, Nodata_value, opt, clean = True):
        #1 Get peaks
        Peaks, slopes_nomax = self.Local_maxima (slope_arr, Nodata_value)
        #2 Get Initiators
        Scarps, slopes_nomax = self.Initiate_ridge (slopes_nomax, Peaks, Nodata_value)
        #3 Continue scarps until they end.
        scarp_order = 3
        for order in range(scarp_order, 50):
            Scarps, slopes_nomax = self.Continue_ridge (slopes_nomax, Scarps, Nodata_value, order)

        if clean == True:
            Scarps =Scarps.Clean_scarps (DEM_arr, Nodata_value, opt)

        return Scarps


    def Local_maxima (self, value_arr, Nodata_value):
        """
        This function is the first stage of a routing process used to identify lines of maximum slopes.
        This function identifies multiple local maxima in an array (Slope), within a predefined search space (Search_space). The identified maxima are given a value of Order.

        Args:
            Slope (2D numpy array): the input 2-D array, here issued from a slope raster.
            Search_space (2D numpy array): the search space array in which to look for local maxima.
            Order (int): the value given to the local maxima points.

        Returns:
            Peaks (2D numpy array): a 2-D array where the local maxima have a value of Order and other elements are null.
            Slope_copy (2D numpy array): a copy of the input array where the value of the selected local maxima has been set to 0.

        Author: GCHG
        """
        print 'Finding local maxima ...'

        values_nomax = Land_surface(value_arr.shape[0], value_arr.shape[1])
        values_nomax = values_nomax.set_attribute (value_arr, Nodata_value, value_arr, Nodata_value, classification = False)

        Search = np.where(self == 1) # the searched locations

        Peaks = Scarps_arr(value_arr.shape[0], value_arr.shape[1]); Peaks = 0*Peaks

        for i in range(len(Search[0])):
            x=Search[0][i]; y=Search[1][i] # coordinates of the kernel's centre
            Kernel_value = fct.kernel (value_arr, 3, x, y)
            Kernel_search = fct.kernel(self, 3, x, y)

            # if the centre of the kernel is its maximum and is not an isolated point
            if Kernel_value[1,1] == np.amax(Kernel_value) and np.amax(Kernel_search[Kernel_search<=Kernel_search[1,1]] > 0):
                Peaks[x,y] = 1 # The kernel centre becomes a local peak
                values_nomax[x,y] = 0 # The slope of the modified data array drops to 0

        return Peaks, values_nomax





    def Initiate_ridge (self, value_arr, peaks_arr, Nodata_value):
        """
        This function is the second stage of a routing process used to identify lines of maximum slopes.
        This function identifies multiple duplets of elements in an array (Slope), within a predefined search space (Search_space) and within the neighbourhood of the local maxima identified in a second input array (Peaks). The identified elements are given a value of Order. To make this function work, the input array Slope should be the output array Slope_copy of the function peak_flag.

        Args:
            Slope (2D numpy array): the input 2-D array, here issued from a slope raster where the local maximal values have been replaced by 0.
            Search_space (2D numpy array): the search space array.
            Peaks (2D numpy array): A 2-D array containing elements with a value of 1. These elements have the same indices as the elements with a value of 0 in Slope.
            Order (int): the value given to the identified elements. it should be superior by 1 to the value of Order in the function peak_flag.

        Returns:
            Ridges (2D numpy array): a 2-D array where the identified elements have a value of Order. This array is modified from the Peaks array and therefore also contains elements of a value equal to the Order in the function peak_flag.
            Slope_copy (2D numpy array): a copy of the input array where the value of the selected elements has been set to 0.

        Author: GCHG
        """

        print ' ... Starting ridges ...'
        values_nomax = Land_surface(value_arr.shape[0], value_arr.shape[1]); values_nomax = 0*values_nomax
        values_nomax = values_nomax.set_attribute (value_arr, Nodata_value, value_arr, Nodata_value, classification = False)

        Search = np.where(self == 1) # the searched locations
        Search_peaks = np.where(peaks_arr == 1) # the searched locations where the peaks are

        ridge_arr = Scarps_arr(peaks_arr.shape[0], peaks_arr.shape[1]); ridge_arr = 0*ridge_arr
        ridge_arr = ridge_arr.set_attribute (peaks_arr, Nodata_value, peaks_arr, Nodata_value, classification = False)


        # Define Kernels
        for i in range(len(Search_peaks[0])):
            x=Search_peaks[0][i]; y=Search_peaks[1][i] # coordinates of the kernel's centre
            Kernel_values = fct.kernel (value_arr, 3, x, y)
            Kernel_values_nomax = fct.kernel (values_nomax, 3, x, y)
            Kernel_ridges = fct.kernel (ridge_arr, 3, x, y)
            Kernel_search = fct.kernel (self, 3, x, y)

            # 1/ If there are no other peaks, we have two ridge starters
            if np.count_nonzero(Kernel_ridges) == 1:
                Ridge_starter1 = np.where (Kernel_values_nomax == np.amax (Kernel_values_nomax))
                X1=Ridge_starter1[0][0]; Y1=Ridge_starter1[1][0]

                # if it is within the initial search space
                if self[x+X1-1, y+Y1-1] != 0:
                    ridge_arr[x+X1-1, y+Y1-1] = 2
                    values_nomax[x+X1-1, y+Y1-1] = 0

                    # Look for a second ridge starter
                    Ridge_starter2 = np.where (Kernel_values_nomax == np.amax (Kernel_values_nomax))
                    X2=Ridge_starter2[0][0]; Y2=Ridge_starter2[1][0]
                    Distance = np.sqrt((X2-X1)**2+(Y2-Y1)**2)

                    # if it is within the initial search space AND not next to the first ridge starter
                    if self[x+X2-1, y+Y2-1] != 0 and Distance > np.sqrt(2):
                        ridge_arr[x+X2-1, y+Y2-1] = 2
                        values_nomax[x+X2-1, y+Y2-1] = 0

                    # Otherwise, look for second ridge starter elsewhere in the kernel
                    elif self[x+X2-1, y+Y2-1] != 0 and Distance <= np.sqrt(2):
                        for j in np.arange(0,9,1):
                            Kernel_values_nomax[X2, Y2] = 0

                            Ridge_starter2 = np.where (Kernel_values_nomax == np.amax (Kernel_values_nomax))
                            X2=Ridge_starter2[0][0]; Y2=Ridge_starter2[1][0]
                            Distance = np.sqrt((X2-X1)**2+(Y2-Y1)**2)

                            if self[x+X2-1, y+Y2-1] != 0 and Distance > np.sqrt(2):
                                ridge_arr[x+X2-1, y+Y2-1] = 2
                                values_nomax[x+X2-1, y+Y2-1] = 0
                                break


            # 2/ If there are two peaks, we have one ridge starter
            elif np.count_nonzero(Kernel_ridges) == 2:
                Ridge_starter1 = np.where (Kernel_values_nomax == np.amax (Kernel_values_nomax))
                X1=Ridge_starter1[0][0]; Y1=Ridge_starter1[1][0]

                # if it is within the initial search space
                if self[x+X1-1, y+Y1-1] != 0:
                    ridge_arr[x+X1-1, y+Y1-1] = 2
                    values_nomax[x+X1-1, y+Y1-1] = 0

        return ridge_arr, values_nomax




    def Continue_ridge (self, value_arr, peaks_arr, Nodata_value, Order):
        """
        This function is the third and final stage of a routing process used to identify lines of maximum slopes.
        IMPORTANT: this function is meant to be run several times! It requires the incrementation of the Order value with each iteration.
        This function identifies multiple elements in an array (Slope), within a predefined search space (Search_space) and within the neighbourhood of the local maxima identified in a second input array (Peaks).  The identified elements are given a value of Order. To make this function work, the input array Slope should be the output array Slope_copy of the function initiate_ridge.

        Args:
            Slope (2D numpy array): the input 2-D array, here issued from a slope raster where the elements selected in the initiate_ridge function have been replaced by 0.
            Search_space (2D numpy array): the search space array.
            Peaks (2D numpy array): A 2-D array containing elements with a value of 1. These elements have the same indices as the elements with a value of 0 in Slope.
            Order (int): the value given to the identified elements. On the first iteration it should be superior by 1 to the value of Order in the function initiate_ridge. the value of Order then needs to be incremented with every iteration.

        Returns:
            Ridges (2D numpy array): a 2-D array where the identified elements have a value of Order. This array is modified from the Peaks array and therefore also contains elements of a value equal to the Order in the functions peak_flag and initiate_ridge.
            Slope_copy (2D numpy array): a copy of the input array where the value of the selected elements has been set to 0.

        Author: GCHG
        """

        print ' ... Prolongating ridges ...'



        values_nomax = Land_surface(value_arr.shape[0], value_arr.shape[1]); values_nomax = 0*values_nomax
        values_nomax = values_nomax.set_attribute (value_arr, Nodata_value, value_arr, Nodata_value, classification = False)

        Search = np.where(self == 1) # the searched locations
        Search_peaks = np.where(peaks_arr == Order-1) # the searched locations where the peaks are

        ridge_arr = Scarps_arr(peaks_arr.shape[0], peaks_arr.shape[1]); ridge_arr = 0*ridge_arr
        ridge_arr = ridge_arr.set_attribute (peaks_arr, Nodata_value, peaks_arr, Nodata_value, classification = False)


        # Define Kernels
        for i in range(len(Search_peaks[0])):
            x=Search_peaks[0][i]; y=Search_peaks[1][i] # coordinates of the kernel's centre

            Kernel_values = fct.kernel (value_arr, 3, x, y)
            Kernel_values_nomax = fct.kernel (values_nomax, 3, x, y)
            Kernel_ridges = fct.kernel (ridge_arr, 3, x, y)
            Kernel_search = fct.kernel (self, 3, x, y)

            # Count the number of nonzero points in the kernel of the ridge array
            Ridge_count = np.count_nonzero(Kernel_ridges)

            # If there are only the 2 previous ridge points, draw a third point that is far enough from the previous point
            if Ridge_count == 2:
                New_point = np.where (Kernel_values_nomax == np.amax (Kernel_values_nomax))
                X=New_point[0][0]; Y=New_point[1][0]
                Grandad_point = np.where (Kernel_ridges == Order-2)
                Xgd=Grandad_point[0][0]; Ygd=Grandad_point[1][0]
                Distance = np.sqrt((X-Xgd)**2+(Y-Ygd)**2)

                if self[x+X-1, y+Y-1] != 0 and Distance > np.sqrt(2):
                    ridge_arr[x+X-1, y+Y-1] = Order
                    values_nomax[x+X-1, y+Y-1] = 0

                elif self[x+X-1, y+Y-1] != 0 and Distance <= np.sqrt(2):
                    for j in np.arange(0,9,1):
                        Kernel_values_nomax[X, Y] = 0

                        New_point = np.where (Kernel_values_nomax == np.amax (Kernel_values_nomax))
                        X=New_point[0][0]; Y=New_point[1][0]
                        Distance = np.sqrt((X-Xgd)**2+(Y-Ygd)**2)

                        if self[x+X-1, y+Y-1] != 0 and Distance > np.sqrt(2):
                            ridge_arr[x+X-1, y+Y-1] = Order
                            values_nomax[x+X-1, y+Y-1] = 0
                            break

        return ridge_arr, values_nomax








##########################################################################################################
##########################################################################################################
class Scarps_arr (Land_surface):
    def __new__ (Scarps_arr, x_length, y_length):
        print 'In __new__ with class %s' % Scarps_arr
        return np.ndarray.__new__(Scarps_arr, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length
        self[np.isnan(self)] = 0
        self = 0 * self


    def Clean_scarps (self, DEM, Nodata_value, opt):
        """
        This function eliminates some of the ridges (Peaks) identified by the trio of functions (peak_flag, initiate_ridge and continue_ridge). The elimination process depends on local relief, which uses a DEM (DEM) and a threshold value (opt). It is set to ignore elements with a value of  Nodata_value.

        Args:
            Peaks (2D numpy array): the input 2-D arraym which is the output of the ridge identification process.
            DEM (2D numpy array): the DEM array used as a base for the elimination of unnecessary ridges.
            Nodata_value (float): The value for ignored elements.
            opt (float): The value of the threshold to eliminate unnecessary ridges.

        Returns:
            Peaks (2D numpy array): a 2-D array much like the input Peaks array, but the unnecessary elemets have been reset to 0.

        Author: GCHG
        """

        print "Cleaning up ridges ..."
        DEM_copy = np.copy(DEM)
        DEM_copy[DEM_copy==Nodata_value] = 0
        Search_ridge = np.where (self != 0)

        Cutoff = np.percentile(DEM_copy,75)
        Threshold = np.amax(DEM_copy[DEM_copy<Cutoff])
        DEM_copy[DEM_copy>Threshold]=Threshold

        for i in range(len(Search_ridge[0])):
            x=Search_ridge[0][i]; y=Search_ridge[1][i] # coordinates of the kernel's centre
            Kernel_DEM = fct.kernel (DEM_copy, 9, x, y)
            Kernel_DEM[Kernel_DEM==Nodata_value]=0

            if np.amax(Kernel_DEM)/Threshold < opt:
                self[x,y] = 0

        Search_ridge = np.where (self != 0)
        for i in range(len(Search_ridge[0])):
            x=Search_ridge[0][i]; y=Search_ridge[1][i] # coordinates of the kernel's centre
            Kernel_ridges = fct.kernel (self, 9, x, y)
            # If there aren't at least 8 ridge points in the neighbourhood of 10 by 10
            if np.count_nonzero(Kernel_ridges) < 8:
                self[x,y] = 0

        return self


    def clean_scarps_from_marsh(self, Marsh):
        Search_false_scarp = np.where (self > 0)
        for i in range(len(Search_false_scarp[0])):
            x = Search_false_scarp[0][i]; y = Search_false_scarp[1][i]
            Kernel_marsh = fct.kernel (Marsh, 3, x, y)
            if np.count_nonzero (Kernel_marsh) == 0:
                self[x, y] = 0

        # We get rid of the sticky-outy bits
        Search_ridge = np.where (self > 0)
        for i in range(len(Search_ridge[0])):
            x=Search_ridge[0][i]; y=Search_ridge[1][i]
            Kernel_ridges = fct.kernel (self, 9, x, y)
            if np.count_nonzero(Kernel_ridges) < 8:
                self[x,y] = 0

        return self


    def extract_platforms (self, DEM, SS_0, Nodata_value, opt):
        """
        This function builds a marsh platform array by using the Peaks array as a starting point. It uses the DEM array to establish conditions on the elements to select. the opt parameter sets a threshold value to eliminate superfluous elements. It is set to ignore elements with a value of Nodata_value.

        Args:
            DEM (2D numpy array): the DEM array.
            Peaks (2D numpy array): the 2-D array of ridge elements, which is the output of the ridge identification and cleaning process.
            Nodata_value (float): The value for ignored elements.
            opt (float): The value of the threshold to eliminate unnecessary elements.

        Returns:
            Marsh (2D numpy array): a 2-D array where the marsh platform elements are identified by strictly positive values. Other elements have a valuof 0 or Nodata_value.

        Author: GCHG
        """

        DEM_copy = Land_surface(DEM.shape[0], DEM.shape[1]); DEM_copy = 0*DEM_copy
        DEM_copy = DEM_copy.set_attribute (DEM, Nodata_value, DEM, Nodata_value, classification = False)
        Marsh = Marsh_platform(DEM.shape[0], DEM.shape[1]); Marsh = 0*Marsh

        print "Initiate platform ..."
        Counter = 1
        Marsh = Marsh.initiate_from_scarps (self, DEM, Nodata_value, Counter)


        print ' ... Build platform ...'
        while Counter < 100:
            Counter = Counter+1
            Marsh, DEM_copy = Marsh.prolong_platform (self, DEM, Nodata_value, Counter)




        print ' ... defining the elimination of low platforms ...'
        Platform = Marsh_platform(Marsh.shape[0], Marsh.shape[1]); Platform = 0*Platform
        Platform = Platform.set_attribute (Marsh, Nodata_value, Marsh, Nodata_value, classification = False)

        Cutoff_Index, Cutoff_Z = Platform.define_cutoff_Z(DEM, SS_0, Nodata_value, opt)


        print 'Cutoff elevation =', Cutoff_Z



        Marsh[Platform<Cutoff_Z] = 0




        print " ... Fill high areas left blank ..."
        Marsh = Marsh.fill_high_gaps(DEM, Cutoff_Index, Nodata_value, 3)


        print ' ... Fill the interior of pools ...'
        Marsh = Marsh.fill_pools(self, DEM)








        # Reapply the cutoff because the straight line thing is ugly
        """NB: The method changes here because we redefine Cutoff_Z instead of using the same.
        This may change accuracy... the old method is below just in case"""
        Platform = Marsh_platform(Marsh.shape[0], Marsh.shape[1]); Platform = 0*Platform
        Platform = Platform.set_attribute (Marsh, Nodata_value, Marsh, Nodata_value, classification = False)

        Cutoff_Index, Cutoff_Z = Platform.define_cutoff_Z(DEM, SS_0, Nodata_value, opt)
        Marsh[Platform<Cutoff_Z] = 0

        #Platform = np.copy(Marsh)
        #Platform[Platform > 0] = DEM [Platform > 0]
        #Marsh[Platform<Cutoff_Z] = 0



        # We fill in the wee holes
        Marsh = Marsh.fill_wee_holes (self, 105)


        print ' ... Adding the ridges'
        # We get rid of scarps that do not have a marsh next to them
        self  = self.clean_scarps_from_marsh(Marsh)

        # We put the scarps in the platform
        Search_scarps = np.where (self > 0)
        Marsh[Search_scarps] = 110


        print " ... eliminate patches of empty elements ..."
        Marsh = Marsh.fill_high_gaps(DEM, Cutoff_Index, Nodata_value, 3)

        print ' ... Fill the interior of pools ...'
        Marsh = Marsh.fill_pools(self, DEM)

        print ' ... defining the elimination of low platforms ...'
        """NB: The method changes here because we redefine Cutoff_Z instead of using the same.
        This may change accuracy... the old method is below just in case"""
        Platform = Marsh_platform(Marsh.shape[0], Marsh.shape[1]); Platform = 0*Platform
        Platform = Platform.set_attribute (Marsh, Nodata_value, Marsh, Nodata_value, classification = False)

        Cutoff_Index, Cutoff_Z = Platform.define_cutoff_Z(DEM, SS_0, Nodata_value, opt)
        Marsh[Platform<Cutoff_Z] = 0

        #Platform = np.copy(Marsh)
        #Platform[Platform > 0] = DEM [Platform > 0]
        #Marsh[Platform<Cutoff_Z] = 0

        Marsh[DEM == Nodata_value] = Nodata_value


        return self, Marsh



##########################################################################################################
##########################################################################################################
class Marsh_platform (Land_surface):
    def __new__ (Marsh_platform, x_length, y_length):
        print 'In __new__ with class %s' % Marsh_platform
        return np.ndarray.__new__(Marsh_platform, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length
        self[np.isnan(self)] = 0
        self = 0 * self


    def calc_labelled_area (self, label_value):
        Selected = np.where(self == label_value)
        Labelled_area = len(Selected[0])

        return Labelled_area


    def initiate_from_scarps (self, Scarps, DEM, Nodata_value, Counter):
        Counter = 1
        Search_ridges = np.where (Scarps > 0)
        for i in range(len(Search_ridges[0])):
            x=Search_ridges[0][i]; y=Search_ridges[1][i]
            Kernel_ridges = fct.kernel (Scarps, 3, x, y)
            Kernel_DEM = fct.kernel (DEM, 3, x, y)

            Marsh_point = np.where (np.logical_and (Kernel_DEM >= Kernel_DEM[1,1], Kernel_ridges == 0))
            for j in range(len(Marsh_point[0])):
                X=Marsh_point[0][j]; Y=Marsh_point[1][j]
                self[x+X-1, y+Y-1] = Counter

        Search_marsh_start = np.where (self == 1)
        for i in range(len(Search_marsh_start[0])):
            x=Search_marsh_start[0][i]; y=Search_marsh_start[1][i]
            Kernel_marsh = fct.kernel (self, 3, x, y)
            Kernel_ridges = fct.kernel (Scarps, 3, x, y)
            if np.count_nonzero(Kernel_marsh) <=2:
                self[x,y] = 0

        return self




    def prolong_platform (self, Scarps, DEM, Nodata_value, Counter):

        DEM_copy = Land_surface(DEM.shape[0], DEM.shape[1]); DEM_copy = 0*DEM_copy
        DEM_copy = DEM_copy.set_attribute (DEM, Nodata_value, DEM, Nodata_value, classification = False)

        Search_marsh = np.where (self == Counter-1)
        for i in range(len(Search_marsh[0])):
            x = Search_marsh[0][i]; y = Search_marsh[1][i]
            Kernel_DEM = fct.kernel (DEM, 3, x, y)
            Kernel_DEM_copy = fct.kernel (DEM_copy, 3, x, y)
            Kernel_ridges = fct.kernel (Scarps, 3, x, y)
            Kernel_marsh = fct.kernel (self, 3, x, y)
            Big_Kernel_DEM = fct.kernel (DEM, 11, x, y)
            Big_Kernel_DEM_copy = fct.kernel (DEM_copy, 11, x, y)


            Conditions = np.zeros((len(Kernel_DEM), len(Kernel_DEM[0,:])), dtype = np.float)
            # 1: free space
            Condition_1 = np.where (np.logical_and(Kernel_ridges == 0, Kernel_marsh == 0)); Conditions[Condition_1] = 1
            # 2: not topped
            Condition_2 = np.where (np.logical_and(Kernel_DEM_copy > np.amax(Big_Kernel_DEM_copy)-0.20, Conditions == 1)); Conditions[Condition_2] = 2


            #This is a distance thing to make sure you don't cross the ridges agin
            Here_be_ridges = np.where (Kernel_ridges != 0)
            Here_be_parents = np.where (Kernel_marsh == Counter-1)

            for j in range(len(Condition_2[0])):
                X=Condition_2[0][j]; Y=Condition_2[1][j]
                Distance_to_ridges = []
                Distance_to_parents = []

                for k in range(len(Here_be_ridges[0])):
                    Xr=Here_be_ridges[0][k]; Yr=Here_be_ridges[1][k]
                    Distance = np.sqrt((X-Xr)**2+(Y-Yr)**2)
                    Distance_to_ridges.append(Distance)

                for k in range(len(Here_be_parents[0])):
                    Xp=Here_be_parents[0][k]; Yp=Here_be_parents[1][k]
                    Distance = np.sqrt((X-Xp)**2+(Y-Yp)**2)
                    Distance_to_parents.append(Distance)

                if len(Distance_to_ridges)>0:
                    if len(Distance_to_parents)>0:
                        if min(Distance_to_ridges) > min(Distance_to_parents):
                            self[x+X-1, y+Y-1] = Counter
                else:
                    self[x+X-1, y+Y-1] = Counter
                    DEM_copy[x+X-1, y+Y-1] = 0

        return self, DEM_copy



    def define_cutoff_Z (self, DEM, SS_0, Nodata_value, opt):
        self[self > 0] = DEM [self > 0]
        Platform_bins, Platform_hist, Platform_step = self.PDF(Nodata_value, 100)

        #1. Find the highest and biggest local maximum of frequency distribution
        # Initialize Index
        Index = len(Platform_hist)-1
        # Initiate Cutoff_Z value
        Cutoff_Z = 0

        for j in range(1,len(Platform_hist)-1):
            if Platform_hist[j]>0.9*max(Platform_hist[1:]) and Platform_hist[j]>Platform_hist[j-1] and Platform_hist[j]>Platform_hist[j+1]:
                Index  = j

        #2. Now run a loop from there toward lower elevations.
        Counter = 0
        for j in range(Index,0,-1):
            # See if you cross the mean value of frequency. Count for how many indices you are under.
            if Platform_hist[j] < np.mean(Platform_hist):
                Counter = Counter + 1
            # Reset the counter value if you go above average again
            else:
                Counter = 0

            #If you stay long enough under (10 is arbitrary for now), initiate cutoff and stop the search
            if Counter > opt:
                Cutoff = j
                Cutoff_Z = Platform_bins[Cutoff]
                break

        # If you stay under for more than 5, set a Cutoff_Z value but keep searching
        if Counter > opt/2:
            Cutoff = j
            Cutoff_Z = Platform_bins[Cutoff]


        # Make sure the Cutoff_Z is higher than the lowest initial search space elevation
        """THIS IS UNPUBLISHED STUFF. IT MIGHT CHANGE THE ACCURACY. TEST IT A COUPLE OF TIMES"""
        SS_1 = SS_0*DEM
        bins, hist, step = SS_1.PDF(Nodata_value, 100)
        zero_point = np.where(bins == 0.)[0]; hist[zero_point] = 0

        low_proportion = 0.01 * max(hist)
        low_proportion_index = np.where(hist > low_proportion)[0][0]

        SS_0_cutoff_Z = bins[low_proportion_index]

        if SS_0_cutoff_Z > Cutoff_Z:
            Cutoff_Z = SS_0_cutoff_Z

        return Index, Cutoff_Z



    def fill_high_gaps(self, DEM, Index, Nodata_value, value):
        Platform_bins, Platform_hist, Platform_step = self.PDF(Nodata_value, 100)
        Search_marsh_condition = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)
        Search_marsh = np.where (DEM >= Platform_bins[Index])
        Search_marsh_condition [Search_marsh] = 1
        Search_marsh_2 = np.where (np.logical_and(self == 0, Search_marsh_condition == 1))
        self[Search_marsh_2] = value

        return self




    def fill_pools(self, Scarps, DEM):
        for Iteration in np.arange(0,10,1):
            Counter = 100
            while Counter > 2:
                print Counter
                Counter = Counter-1
                Search_marsh = np.where (self == Counter+1)
                Non_filled = 0
                for i in range(len(Search_marsh[0])):
                    x = Search_marsh[0][i]; y = Search_marsh[1][i]
                    Kernel_DEM = fct.kernel (DEM, 3, x, y)
                    Kernel_ridges = fct.kernel (Scarps, 3, x, y)
                    Kernel_marsh = fct.kernel (self, 3, x, y)

                    if Non_filled <len(Search_marsh[0]):
                        if np.count_nonzero(Kernel_marsh) > 6:
                            Condition = np.where (np.logical_and(Kernel_marsh == 0, Kernel_ridges == 0))
                            for j in range(len(Condition[0])):
                                X=Condition[0][j]; Y=Condition[1][j]
                                self[x+X-1, y+Y-1] = Counter
                        else:
                            Non_filled = Non_filled + 1

        return self



    def fill_wee_holes (self, Scarps, value):
        Search_marsh = np.where (np.logical_and(self == 0, Scarps == 0))
        for i in range(len(Search_marsh[0])):
            x = Search_marsh[0][i]; y = Search_marsh[1][i]
            Kernel_marsh = fct.kernel (self, 3, x, y)
            if np.count_nonzero(Kernel_marsh) == 8:
                self[x,y] = value

        return self



    def edges_to_shp(self, dir, name, envidata, Nodata_value):

        self_labels = self.label_connected (Nodata_value)

        if not os.path.isfile(dir+name+'_Out.pkl'):
            Outlines = self_labels.extract_outlines()
            Outlines.Save_as_pickle(dir,name+'_Out.pkl')

        else:
            print 'Loading outlines from Pickle'
            Outlines = pickle.load( open(dir+name+'_Out.pkl', "rb" ) )

        print Outlines
        Outlines.save_to_shp (envidata, dir, name, multiple = False)





        quit()

        """if not os.path.isfile(dir+file[:-4]+'.pkl'):
            Outlines = Marsh_labels.extract_outlines()
            Outlines.Save_as_pickle(Output_dir+'Marsh_metrics/',str(site)+'_Out.pkl')

        else:
            print 'Loading outlines from Pickle'
            Outlines = pickle.load( open(file[:-4]+'.pkl', "rb" ) )
        Outlines.save_to_shp (envidata_DEM, DEM, Output_dir+'Shapefiles/', site, multiple = False)



        +'Marsh_metrics/'+str(site)+'_Out.pkl'):
            print 'Extracting outlines from array'
            Outlines = Marsh_labels.extract_outlines()
            Outlines.Save_as_pickle(Output_dir+'Marsh_metrics/',str(site)+'_Out.pkl')
        else:
            print 'Loading outlines from Pickle'
            Outlines = pickle.load( open( Output_dir+'Marsh_metrics/'+str(site)+'_Out.pkl', "rb" ) )
        Outlines.save_to_shp (envidata_DEM, DEM, Output_dir+'Shapefiles/', site, multiple = False)

        """




##########################################################################################################
##########################################################################################################
class Marsh_outline (Land_surface):
    def __new__ (Marsh_outline, x_length, y_length):
        print 'In __new__ with class %s' % Marsh_outline
        return np.ndarray.__new__(Marsh_outline, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length
        self[np.isnan(self)] = 0
        self = 0 * self

    def complete_outline_from_array (self, array, Nodata_value):
        print '\nExtracting outline...'
        Start = timeit.default_timer()
        new_array = self.copy()

        new_array, Outline_value = fct.Surface_outline (new_array, array, 2, Nodata_value)

        Stop = timeit.default_timer()
        print '  Extraction runtime : ', Stop - Start , 's'

        return new_array, Outline_value

    def tightrope_outline_from_array (self, array):
        print '\nExtracting outline...'
        Start = timeit.default_timer()
        new_array = self.copy()
        new_array, Outline_value = fct.Tightrope_outline (new_array, array, 2)
        Stop = timeit.default_timer()
        print '  Extraction runtime : ', Stop - Start , 's'

        return new_array, Outline_value



    def reduce_to_marsh_labels (self, Marsh_labels, Nodata_value):
        # This one relies on the fact that the longest continuous line should be the marsh outline for each label
        print '\nReducing outline...'
        Start = timeit.default_timer()

        new_array = 0 * self.copy()

        # Find out how many marsh labels you have
        M_Labels = range (int(np.amin(Marsh_labels[Marsh_labels>0])), int(np.amax(Marsh_labels[Marsh_labels>0]))+1, 1)

        # Find out how many outline labels you have
        L_Labels = range (int(np.amin(self[self>0])), int(np.amax(self[self>0]))+1, 1)

        # Make a list counting the elements for each label. The 0 index is for the label value, the 1 index is for the number of elements
        Num_Elements = [[],[]]
        for label in L_Labels:
            Elements = np.where(self == label)
            num_elements = len (Elements[0])
            Num_Elements[0].append(label); Num_Elements[1].append(num_elements)

        for i in range(len(M_Labels)):
            # If the label has a non-negligible area
            Labelled_area = Marsh_labels.calc_labelled_area (M_Labels[i])

            if Labelled_area >1 :
                Most_pop_index = np.where(np.asarray(Num_Elements[1])==max(np.asarray(Num_Elements[1])))[0][0]

                Most_pop_label = Num_Elements[0][Most_pop_index]
                new_array[self == Most_pop_label] = Most_pop_label

                # Mow remove it.
                del Num_Elements[0][Most_pop_index]; del Num_Elements[1][Most_pop_index]

        Stop = timeit.default_timer()
        print '  Reduction runtime : ', Stop - Start , 's'

        return new_array



    def trim_to_main_stem (self, Marsh_array, Nodata_value):
        print '\nTrimming the messy outlines...'
        Start = timeit.default_timer()
        new_array = self.copy()
        # Find out how many labels you have
        Labels = range (int(np.amin(self[self>0])), int(np.amax(self[self>0]))+1, 1)
        #And loop through the labels
        for lab in range(len(Labels)):
            print '  This is the label: ', lab+1, '/', len(Labels), ' (', Labels[lab], ')'
            new_array = fct.Billhook (new_array, Marsh_array, Labels[lab])
        Stop = timeit.default_timer()
        print '  Gardening runtime : ', Stop - Start , 's'

        return new_array



    def calc_outline_length (self, Scale, Short = False):
        print '\nCalculating outline length...'
        Start = timeit.default_timer()
        #Setup the environmnents
        Length_array = self.copy()
        Labels_array = self.copy()
        #Initiate the vector
        Polylines = Polyline()
        # Find out how many labels you have
        Labels = range (int(np.amin(self[self>0])), int(np.amax(self[self>0]))+1, 1)
        #And loop through thev labels
        for lab in range(len(Labels)):
            print '\nThis is the label: ', lab+1, '/', len(Labels), ' (', Labels[lab], ')'
            # Measure the length of the stitched line for this label value
            This_polyline, Length_array, Code_array = fct.Measure_all_lines (Length_array, Labels[lab], Scale)
            #Stitch the diverging starts
            print
            # Rehabilitate this in a different manner
            This_polyline, Length_array, Code_array = fct.Stitch_diverging_starts (Length_array, Labels_array, Labels[lab], This_polyline, Code_array, Scale)
            print
            #Stitch outward going branches
            #This_polyline, Length_array = fct.Graft_diverging_branch (Length_array, Labels_array, Labels[lab], This_polyline, Code_array, Scale)
            print

            #Stitch inward going branches, but not yet.
            #new_array, Line_row, Line_col, Line_
            #Line_dist, Line_code = fct.Graft_converging_branch (new_array, Labels_array, Labels[lab], Line_row, Line_col, Line_dist, Line_code, Code_array, Scale)

            Polylines.append(This_polyline)

        Stop = timeit.default_timer()
        print '  Surveying runtime : ', Stop - Start , 's'

        return Length_array, Polylines

#######################################################################################################
##########################################################################################################
class Point (tuple):
    """With this class we make a point that has row and column attributes."""

    def __new__(self, x, y):
        return tuple.__new__(Point, (x, y))

    def row (self):
        row = self[0]
        return row

    def col (self):
        col = self[1]
        return col

    def Dist_to_point (self, other):
        return np.sqrt((self.row()-other.row())**2 + (self.col()-other.col())**2)


##########################################################################################################
##########################################################################################################
class Line (bb.DataFrame):
    """The Line class is a Pandas DataFrame.
    It is super important"""
    @property
    def _constructor(self):
        return Line


    def set_first_point(self, M_code, L_code, row, col, scale = 1):
        """You must always initiate a line with a first point"""
        self['M_code'] = [M_code]
        self['L_code'] = [L_code]
        self['rowcol'] = [Point(row,col)]
        self['length'] = [0.001*scale]


    def add_element (self, position, M_code, L_code, row, col, scale = 1):
        """position must be positive"""
        length = scale * Point(row,col).Dist_to_point(self['rowcol'].iloc[position-1]) + self['length'].iloc[position-1]
        element = bb.DataFrame({'M_code': [M_code],'L_code': [L_code], 'rowcol':[Point(row,col)], 'length':[length]})
        top = self[0:position]; bottom = self[position:]
        self=bb.concat((top,element,bottom))
        self = self.reset_index(drop=True)
        return self


    def remove_element (self,position):
        self = self.drop(self.index[position])
        self = self.reset_index(drop=True)
        return self


    def start_point (self):
        return self['rowcol'].iloc[0]


    def end_point (self):
        return self['rowcol'].iloc[-1]


    def add_attribute_list (self, attr_name, attr_list):
        """Only use after there are no more points to add"""
        if type(attr_list) is bool:
            self[str(attr_name)] = bb.Series(attr_list, index=self.index)
        elif type(attr_list) is list:
            if len(attr_list) == len(self['rowcol']):
                self[str(attr_name)] = bb.Series(attr_list, index=self.index)
        elif type(attr_list) is np.float64:
            self[str(attr_name)] = bb.Series(attr_list, index=self.index)
        return self


    def remove_attribute (self, attr_name):
        """Same here, do not use if you are still modifying the line geometry"""
        self.drop(columns=[str(attr_name)])


    def save_line_to_shp (self, Envidata, Enviarray, save_dir, file_name):
        fct.Line_to_shp(self, Envidata, Enviarray, save_dir, file_name)


    def Pandaline_to_shp (self, Envidata, save_dir, site_name):
        value_range = range(min(self['L_code']),max(self['L_code'])+1)
        for L in value_range:
            To_save = self.loc[self['L_code'] == L]
            M_code = To_save['M_code'][0]
            L_code = To_save['L_code'][0]
            file_name = str(site_name)+'_'+str(M_code)+'_'+str(L_code)
            fct.Line_to_shp(To_save, Envidata, save_dir, file_name)


    def prepare_for_plotting(self,colour,opaque = True):
        Line_row = []; Line_col = []
        for i in range(len(self['rowcol'])):
            Line_row.append(self['rowcol'].iloc[i].row()); Line_col.append(self['rowcol'].iloc[i].col())
        if opaque == False:
            Line = Line2D(Line_col, Line_row, color = colour, alpha = 0.0)
        else:
            Line = Line2D(Line_col, Line_row, color = colour, alpha = 1, linewidth = 0.7)
        return Line


    def extract_values_from_basemap (self, basemap, attr_name, Nodata_value):
        prop_list = []
        for i in range(len(self['rowcol'])):
            point = self['rowcol'].iloc[i]
            if int(point.row()) < basemap.shape[0] and int(point.col()) < basemap.shape[1] and int(point.row()) > 0 and int(point.col()) > 0:
                value = basemap[int(point.row()), int(point.col())]
            else:
                value = Nodata_value
            prop_list.append(value)
        self.add_attribute_list(attr_name,prop_list)
        return self


    def recalc_length (self,scale):
        self['length'].iloc[0] = 0.001 * scale
        for i in range(1,len(self['length'])):
            self['length'].iloc[i] = self['rowcol'].iloc[i].Dist_to_point(self['rowcol'].iloc[i-1])+ self['length'].iloc[i-1]
        return self


    def Line_transects (self, spacing, length, Envidata, Enviarray, save_dir, site_name):
        Transects = Transect()
        value_range = range(min(self['L_code']),max(self['L_code'])+1)
        for L in value_range:
            To_save = self.loc[self['L_code'] == L]
            M_code = To_save['M_code'][0]
            L_code = To_save['L_code'][0]
            in_name = str(site_name)+'_'+str(M_code)+'_'+str(L_code)+'.shp'
            out_name = str(site_name)+'_'+str(M_code)+'_'+str(L_code)+'_Tr.shp'

            fct.Make_transects(save_dir+'Shapefiles/'+in_name, save_dir+'Shapefiles/'+out_name, spacing, length)

            os.system('rm '+save_dir+'Shapefiles/'+str(site_name)+'_'+str(M_code)+'_'+str(L_code)+'.*')

            if not os.path.isfile(save_dir+'Shapefiles/'+str(site_name)+'_'+str(M_code)+'_'+str(L_code)+'_Tr.dbf'):
                print 'Missed a transect in the saving procedure'
            else:
                Trans = fct.Shp_to_transects (save_dir+"Shapefiles/", out_name, M_code, L_code, Envidata, Enviarray)
                if Trans.size>0:
                    Transects = Transects.append(Trans)

            os.system('rm '+save_dir+'Shapefiles/'+str(site_name)+'_'+str(M_code)+'_'+str(L_code)+'_Tr.*')

        return Transects


##########################################################################################################
##########################################################################################################
class Transect (Line):
    """The Transect class is a Pandas DataFrame.
    It is also super important"""
    @property
    def _constructor(self):
        return Transect


    def Save_as_pickle (self,filepath,filename):
        with open(filepath+filename, 'wb') as handle:
            pickle.dump(self,handle)


    def assign_transect_code (self):
        Tr_codes = []
        for i in range(len(self['rowcol'])//2):
            Tr_codes.append(int(i+1)); Tr_codes.append(int(i+1))
        self = self.add_attribute_list('T_code',Tr_codes)
        return self


    def add_transect_element  (self, position, M_code, L_code, T_code, row, col, scale = 1):
        """position must be positive"""
        length = scale * Point(row,col).Dist_to_point(self['rowcol'].iloc[position-1]) + self['length'].iloc[position-1]
        element = bb.DataFrame({'M_code': [M_code],'L_code': [L_code],'T_code': [T_code], 'rowcol':[Point(row,col)], 'length':[length]})
        top = self[0:position]; bottom = self[position:]
        self=bb.concat((top,element,bottom))
        self = self.reset_index(drop=True)
        return self


    def slope_1D (self, scale = 1):
        if self.size > 0:
            dZ_list = []
            step = max(self.index)+1
            for i in range(0,len(self['rowcol']),step):

                Start = i; End = i+step
                dZ_list.append(0)
                for j in range(1,step):
                    dZ_list.append((self['Z'].iloc[Start+j] - self['Z'].iloc[Start+j-1])/scale)
            self = self.add_attribute_list ('dZ', dZ_list)
        return self


    def get_bearing (self):
        if self.size > 0:
            #print self
            bear_list = []
            step = max(self.index)+1
            for i in range(0,len(self['rowcol']),step):
                Start = i; End = i+step-1
                bear_list.append(0)
                row = self['rowcol'].iloc[End][0]-self['rowcol'].iloc[Start][0]
                col = self['rowcol'].iloc[End][1]-self['rowcol'].iloc[Start][1]
                for j in range(1,step):
                    bear = np.pi/2 + np.arctan2(row, col)# * 180 / np.pi
                    if bear < 0:
                        bear = bear + 2*np.pi + 360
                    bear_list.append(bear)
            self = self.add_attribute_list ('bear', bear_list)
        return self


    def single_subdivide (self,sub_number):
        if len(self['rowcol']) == 2:
            sub_row = []; sub_col = []
            for sub in range(1,sub_number):
                point_1 = self['rowcol'].iloc[0]; point_2 = self['rowcol'].iloc[-1]
                sub_row.append(point_1.row() + sub * float(point_2.row()-point_1.row())/sub_number)
                sub_col.append(point_1.col() + sub * float(point_2.col()-point_1.col())/sub_number)
            for sub in range(0,sub_number-1):
                position = len(self['rowcol'])-1
                self = self.add_transect_element (position, self['M_code'].iloc[0], self['L_code'].iloc[0], self['T_code'].iloc[0], sub_row[sub], sub_col[sub])
        return self


    def multiple_subdivide (self,sub_number):
        """Only works if the initial line only has two points."""
        new_self = Transect(); new_self.set_first_point(0, 0, 0, 0)
        new_self = new_self.add_attribute_list ('T_code', [0])
        L_range = range(min(self['L_code']),max(self['L_code'])+1)
        for L in L_range:
            L_self = self.loc[self['L_code'] == L]
            T_range = range(min(self['T_code']),max(self['T_code'])+1)
            for T in T_range:
                T_self = L_self.loc[L_self['T_code'] == T]
                T_self = T_self.single_subdivide(sub_number)
                new_self = new_self.append(T_self)
        new_self = new_self.iloc[1:]
        return new_self


    def orient_seaward (self):
        """Leave the Nodata_values as they are for now."""
        """Only use if you aready have a Z attribute"""
        if self.size > 0:
            step = max(self.index)+1
            for i in range(0,len(self['rowcol']),step):
                Start = i; End = i+step
                if self['Z'].iloc[Start] < self['Z'].iloc[End-1]:
                    New = self.iloc[Start:End].sort_index(ascending = False)
                    New = New.reset_index(drop = True)
                    self.iloc[Start:End] = New
        return self


    def select_transect (self, Nodata_value):
        """Only use if you already have Marsh and Z attribute"""
        if self.size > 0:
            Select_list = []
            step = max(self.index)+1
            for i in range(0,len(self['rowcol']),step):
                Start = i; End = i+step
                Zeros = np.where(np.asarray(self['Marsh'].iloc[Start:End]) < 0.5)
                Ones = np.where(np.asarray(self['Marsh'].iloc[Start:End]) > 0.5)
                #print Ones
                #print Zeros
                if len(Zeros[0]) == 0 or len(Ones[0]) == 0:
                    for j in range(step):
                        Select_list.append(False)
                elif min(Zeros[0])<max(Ones[0]):
                    for j in range(step):
                        Select_list.append(False)
                elif min(self['Z'].iloc[Start:End]) <= Nodata_value:
                    for j in range(step):
                        Select_list.append(False)
                else:
                    for j in range(step):
                        Select_list.append(True)

                #print Select_list
                #quit()
            self.add_attribute_list ('select', Select_list)
        return self


    def get_statistics (self):
        """ Objects are are for Mean, Stdev, """
        if self.size > 0:
            step = max(self.index)+1
            Mean = self[:step].copy(deep=True); Stdev = self[:step].copy(deep=True)
            for i in range(0,step):
                if len (self.loc[self['select'] == True]) > 0:
                    Selected = self.loc[self['select'] == True].loc[i]
                    for attr_name in ['Z','dZ']:
                        Mean[attr_name].loc[i] = np.mean(Selected[attr_name])
                        Stdev[attr_name].loc[i] = np.std(Selected[attr_name])
                    Mean['length'].loc[i] = i; Mean['select'].loc[i] = True
                    Stdev['length'].loc[i] = i; Stdev['select'].loc[i] = True
        else:
            Mean = self.copy(deep=True); Stdev = self.copy(deep=True)
        return Mean, Stdev



    def squish(self, profile_length):

        transect_len = profile_length
        # Make a new transect object
        Squished_self = Transect()

        # Give it the columns you want
        Squished_self['T_code'] = [1]
        Squished_self['rowcol'] = [Point(1,1)]

        content_ini = (1.,)
        for x in range(transect_len):
            content = content_ini + (x,)
            content_ini = content
        Squished_self['bearing'] = [1]

        i = 1
        for t in range(min(self['T_code']), max(self['T_code'])+1):
            this_transect = self[self['T_code']==t]

            code = this_transect['T_code'].values[0]

            coord_start, coord_end = this_transect['rowcol'].values[0], this_transect['rowcol'].values[-1]
            Z_list = tuple(this_transect['Z'].values)
            dZ_list = tuple(this_transect['dZ'].values)
            bearing = this_transect['bear'].values[1]
            select = this_transect['select'].values[0]

            element = bb.DataFrame({'T_code': [code,code], 'rowcol':[coord_start,coord_end], 'Z_dZ':[Z_list,dZ_list], 'bearing':[bearing,bearing], 'select':[select,select]})

            top = Squished_self[0:i]; bottom = Squished_self[i:]
            Squished_self=bb.concat((top,element,bottom))
            Squished_self = Squished_self.reset_index(drop=True)

            i+=2

        Squished_self = Squished_self[1:]

        return Squished_self



    def save_to_fullshp(self, envidata, save_file):
        fct.Transect_to_fullshp(self, envidata, save_file)







##########################################################################################################
##########################################################################################################
class Polyline (list):
    """A polyline object has Line or Polyline objects inside. Let's say for now it just has Lines."""
    def __init__(self):
        list.__init__(self)


    def Save_as_pickle (self,filepath,filename):
        with open(filepath+filename, 'wb') as handle:
            pickle.dump(self,handle)


    def save_to_shp (self, Envidata, save_dir, site_name, multiple = True):
        if multiple is True:
            for i in range(len(self)):
                Outline = self[i]
                Outline.Pandaline_to_shp(Envidata, save_dir, site_name)
        else:
            Outline = self
            file_name = str(site_name)+'_MarshOutline'
            fct.Polyline_to_shp (Outline, Envidata, save_dir, file_name)


    def Transects_to_fullshp (self, envidata, save_file):
        for i in range(len(self)):
            self[i].save_to_fullshp(envidata, save_file)



    def Polyline_transects (self, length, Envidata, Enviarray, save_dir, site_name, needs_saving = False):
        #if you haven't done so yet, save your polyline to a .shp
        if needs_saving is True:
            self.save_to_shp (Envidata, Enviarray, save_dir, site_name)

        in_name = str(site_name)+'_MarshOutline.shp'
        out_name = str(site_name)+'_OutlineProfiles.shp'

        # Produce transects for the outline shapefiles
        fct.Make_transects(save_dir+'Shapefiles/'+in_name, save_dir+'Shapefiles/'+out_name, spacing, length)


    def squish_properties(self, profile_length):
        # This one condenses the properties of a long pandas into a shorter pandas
        for i in range(len(self)):
            self[i] = self[i].squish(profile_length)
        return self


    def get_attribute_from_basemap (self, refinement, basemap, attr_name, Nodata_value):
        for i in range(len(self)):
            if self[i].size > 0:
                self[i] = self[i].multiple_subdivide(refinement)
                self[i] = self[i].extract_values_from_basemap (basemap, attr_name, Nodata_value)
                if attr_name == 'Z':
                    self[i] = self[i].orient_seaward()
                    self[i] = self[i].slope_1D()
                    self[i] = self[i].get_bearing()
                if attr_name == 'Marsh':
                    self[i] = self[i].select_transect(Nodata_value)
        return self


    def Polyline_stats (self):
        Mean_polyline = Polyline()
        Stdev_polyline = Polyline()
        for i in range(len(self)):
            Transects = self[i]
            if i == 0:
                Bigtransect = Transects.copy(deep = True)
            else:
                Bigtransect = bb.concat([Bigtransect,Transects])
            Transects_mean,Transects_stdev = Transects.get_statistics()
            Mean_polyline.append(Transects_mean)
            Stdev_polyline.append(Transects_stdev)
        Big_mean, Big_stdev = Bigtransect.get_statistics()
        return Mean_polyline, Stdev_polyline, Big_mean, Big_stdev, Bigtransect




    def plot_transects_basemap_and_profiles (self, basemap, mask, bg, save_dir, figname, Nodata_value):
        twin  = basemap.copy()
        twinmask  = mask.copy()

        #Make the canvas
        fig_height = min(np.floor(twin.shape[1])/5, 25)
        fig_width = min(np.floor(twin.shape[1])/5, 10)
        fig=plt.figure(figname, facecolor='White',figsize=[fig_height,fig_width])

        ax1 = plt.subplot2grid((1,2),(0,0),colspan=1, rowspan=1)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')
        plt.xlabel('x (m)', fontsize=18)
        plt.ylabel('y (m)', fontsize=18)

        ax1.annotate('a.', fontsize = 14,
             xy=(0.95, 0.05),
             xycoords='axes fraction',color='white')

        # configure the basemap
        Vmin = min(np.amin(twin[twin>Nodata_value])*0.95, np.amin(twin[twin>Nodata_value])*1.05)
        Vmax = max(np.amax(twin)*0.95, np.amax(twin)*1.05)
        twinmask = np.ma.masked_where(twinmask == 1, twinmask)

        Map = ax1.imshow(bg, interpolation='None', cmap=plt.cm.gray, vmin=0, vmax=210, alpha = 1.0)
        Map = ax1.imshow(twin, interpolation='None', cmap=plt.cm.gist_earth, vmin=Vmin, vmax=Vmax, alpha = 0.6)
        Map = ax1.imshow(twinmask, interpolation='None', cmap= plt.cm.gray, vmin=Vmin, vmax=Vmax, alpha = 0.4)

        ax2 = fig.add_axes([0.13, 0.98, 0.345, 0.02])
        scheme = plt.cm.gist_earth; norm = colors.Normalize(vmin=Vmin, vmax=Vmax)
        cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)
        cb1.set_label('Elevation (m a.m.s.l.)', fontsize=18)


        # Draw the lines, panda style
        colour = 0
        for i in range(len(self)):
            Pandaline = self[i]
            if Pandaline.size > 0:
                colour+=1
                L_labels = range (1,max(Pandaline['L_code'])+1)
                for L in L_labels:
                    Pandaline_slice = Pandaline.loc[Pandaline['L_code'] == L]
                    #print Pandaline_slice
                    if Pandaline_slice.size > 0:
                        T_labels = range (1,max(Pandaline_slice['T_code'])+1)
                        for T in T_labels:
                            Pandaline_zest = Pandaline_slice.loc[Pandaline_slice['T_code'] == T]
                            To_draw = Pandaline_zest.prepare_for_plotting(plt.cm.jet(colour*50),opaque = Pandaline_zest['select'].iloc[0])
                            ax1.add_line(To_draw)



        # Make the other plot
        ax3 = plt.subplot2grid((1,2),(0,1),colspan=1, rowspan=1)
        ax3.tick_params(axis='x', colors='black')
        ax3.tick_params(axis='y', colors='black')
        plt.ylabel('Elevation (m a.m.s.l.)', fontsize=18)
        plt.xlabel('Transect length (m)', fontsize=18)

        ax3.annotate('b.', fontsize = 14,
             xy=(0.05, 0.05),
             xycoords='axes fraction')

        ax3.set_xlim(0,10)
        #ax3.set_ylim(ymin=-2, ymax = 10)

        colour = 0
        for i in range(len(self)):
            Pandaline = self[i]
            if Pandaline.size > 0:
                colour+=1
                L_labels = range (1,max(Pandaline['L_code'])+1)
                for L in L_labels:
                    Pandaline_slice = Pandaline.loc[Pandaline['L_code'] == L]
                    #print Pandaline_slice
                    if Pandaline_slice.size > 0:
                        T_labels = range (1,max(Pandaline_slice['T_code'])+1)
                        for T in T_labels:
                            Pandaline_zest = Pandaline_slice.loc[Pandaline_slice['T_code'] == T]
                            if Pandaline_zest['select'].iloc[0] == True:
                                ax3.plot(Pandaline_zest['Z'], color = plt.cm.jet(colour*50))

        ax3.axvline(5, color='black', lw=1.0, alpha=0.8)

        ax4 = fig.add_axes([0.55, 0.98, 0.345, 0.02])
        scheme = plt.cm.jet; norm = colors.Normalize(vmin=0, vmax=colour)
        cb2 = matplotlib.colorbar.ColorbarBase(ax4, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)
        cb2.set_label('Platform number', fontsize=18)

        plt.savefig(save_dir+figname+'_Tr.png', bbox_inches='tight')








    def plot_on_basemap(self,basemap, mask, save_dir, figname, Nodata_value):
        'Plotting on basemap'
        twin  = basemap.copy()
        #Make the canvas
        fig_height = min(np.floor(twin.shape[1])/5, 50); fig_width = min(np.floor(twin.shape[1])/5, 50)
        fig=plt.figure(figname, facecolor='White',figsize=[fig_height,fig_width])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')

        # Make the basemap
        Vmin = min(np.amin(twin[twin!=Nodata_value])*0.95, np.amin(twin[twin!=Nodata_value])*1.05)
        Vmax = max(np.amax(twin)*0.95, np.amax(twin)*1.05)

        Map = ax1.imshow(twin, interpolation='None', cmap=plt.cm.gist_earth, vmin=Vmin, vmax=Vmax, alpha = 0.6)
        ax2 = fig.add_axes([0.1, 0.98, 0.85, 0.02])
        scheme = plt.cm.gist_earth; norm = colors.Normalize(vmin=Vmin, vmax=Vmax)
        cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)

        # Draw the lines, panda style
        for i in range(len(self)):
            Pandaline = self[i]
            L_labels = range (1,max(Pandaline['L_code'])+1)
            for j in L_labels:
                Pandaline_slice = Pandaline.loc[Pandaline['L_code'] == j]
                To_draw = Pandaline_slice.prepare_for_plotting(i,0.5)
                ax1.add_line(To_draw)
                ax1.scatter(Pandaline_slice['rowcol'].iloc[0].col(),Pandaline_slice['rowcol'].iloc[0].row())

        plt.savefig(save_dir+figname+'.png')



    def plot_transects_on_basemap(self,basemap, save_dir, figname, Nodata_value):
        'Plotting on basemap'
        twin  = basemap.copy()
        #Make the canvas
        fig_height = min(np.floor(twin.shape[1])/5, 50); fig_width = min(np.floor(twin.shape[1])/5, 50)
        fig=plt.figure(figname, facecolor='White',figsize=[fig_height,fig_width])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')

        # Make the basemap
        Vmin = min(np.amin(twin[twin!=Nodata_value])*0.95, np.amin(twin[twin!=Nodata_value])*1.05)
        Vmax = max(np.amax(twin)*0.95, np.amax(twin)*1.05)

        Map = ax1.imshow(twin, interpolation='None', cmap=plt.cm.gist_earth, vmin=Vmin, vmax=Vmax, alpha = 0.6)
        ax2 = fig.add_axes([0.1, 0.98, 0.85, 0.02])
        scheme = plt.cm.gist_earth; norm = colors.Normalize(vmin=Vmin, vmax=Vmax)
        cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)

        # Draw the lines, panda style
        for i in range(len(self)):
            Pandaline = self[i]
            if Pandaline.size > 0:
                L_labels = range (1,max(Pandaline['L_code'])+1)
                for L in L_labels:
                    Pandaline_slice = Pandaline.loc[Pandaline['L_code'] == L]
                    #print Pandaline_slice
                    if Pandaline_slice.size > 0:
                        T_labels = range (1,max(Pandaline_slice['T_code'])+1)
                        for T in T_labels:
                            Pandaline_zest = Pandaline_slice.loc[Pandaline_slice['T_code'] == T]
                            To_draw = Pandaline_zest.prepare_for_plotting(i,opaque = Pandaline_zest['select'].iloc[0])
                            ax1.add_line(To_draw)
                            ax1.scatter(Pandaline_slice['rowcol'].iloc[0].col(),Pandaline_slice['rowcol'].iloc[0].row())

        plt.savefig(save_dir+figname+'_Tr.png')




    def plot_property(self, save_dir, figname, draw):
        fig=plt.figure(figname, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')


        # Draw the lines.
        Structure_list = self.structure()
        print Structure_list
        print len(Structure_list)

        if len(Structure_list) == 4 and Structure_list[-1] is Point:
            for i in range(len(self)):
                for j in range(len(self[i])):
                    Line1 = self[i][j]
                    #print Line1
                    if Line1[-1] == True:
                        # This is what you plot
                        ax1.plot(Line1[draw])#, color = plt.cm.jet(i*20), alpha = 0.5+j/100)

        elif len(Structure_list) == 5 and Structure_list[-1] is Point:
            for h in range(len(self)):
                    for i in range(len(self[h])):
                        for j in range(len(self[h][i])):
                            Line1 = self[h][i][j]

                            if Line1[-1] == True:
                                # This is what you plot
                                ax1.plot(Line1[draw])#, color = plt.cm.jet(i*20), alpha = 0.5+j/100)

        plt.savefig(save_dir+figname+'.png')


    def plot_property_stats(self, save_dir, figname, draw):
        fig=plt.figure(figname, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')

        # Draw the lines.
        Structure_list = self.structure()

        if len(Structure_list) == 4 and Structure_list[-1] is Point:
            for i in range(len(self)):
                ax1.plot(self[i][-2], color = plt.cm.jet(i*20), alpha = 0.5/100)

        elif len(Structure_list) == 5 and Structure_list[-1] is Point:
            for h in range(len(self)):
                    for i in range(len(self[h])):
                        #print self[h][i][-2]
                        #ax1.plot(self[h][i][-2])#, color = plt.cm.jet(i*20), alpha = 1)
                        if max( self[h][i][-2]) > 0:
                            ax1.fill_between(range(0,21), self[h][i][-2]-self[h][i][-1], self[h][i][-2]+self[h][i][-1], alpha = 0.3)#, color = plt.cm.jet(i*20), alpha = 0.3)

        plt.savefig(save_dir+figname+'.png')



    """def select_few_longest (self):
        New_polyline = fct.Select_few_longest(self)
        return New_polyline"""



    """def select_few_longest (self, Nodata_value, num):
        # num is the number of lines you select
        empty = np.zeros(self.shape, dtype = np.float)

        twin = np.copy(self)

        values = range (np.amin(self[self>0]), np.amax(self), 1)
        line_lengths = []
        for value in values:
            line_lengths.append(len(np.where(self == value)[0]))

        line_lengths = np.asarray(line_lengths)
        Longest = np.where(line_lengths == np.amax(line_lengths))
        print values[Longest[0][0]], line_lengths[Longest[0][0]]
        array_2[array == values[Longest[0][0]]] = values[Longest[0][0]]

        if num > 0:
            for i in range(num):
                line_lengths[Longest[0][0]] = 0
                Longest = np.where(line_lengths == np.amax(line_lengths))
                print Longest[0][0]
                self[twin == values[Longest[0][0]]] = values[Longest[0][0]]

        return self"""



    """def select_few_longest (array, Nodata_value, num):
    # num is the number of lines you select
    array_2 = np.zeros(array.shape, dtype = np.float)

    values = range (np.amin(array[array>0]), np.amax(array), 1)
    line_lengths = []
    for value in values:
        line_lengths.append(len(np.where(array == value)[0]))

    line_lengths = np.asarray(line_lengths)
    Longest = np.where(line_lengths == np.amax(line_lengths))
    print values[Longest[0][0]], line_lengths[Longest[0][0]]
    array_2[array == values[Longest[0][0]]] = values[Longest[0][0]]

    if num > 0:
        for i in range(num):
            line_lengths[Longest[0][0]] = 0
            Longest = np.where(line_lengths == np.amax(line_lengths))
            print Longest[0][0]
            array_2[array == values[Longest[0][0]]] = values[Longest[0][0]]

    return array_2"""

    def vectorise (self, Nodata_value):
        from matplotlib.lines import Line2D

        Outlines = []; Outlines_row = []; Outlines_col = []

        Labels = range (int(np.amin(self[self>0])), int(np.amax(self[self>0]))+1, 1)

        for lab in range(len(Labels)):
            Select = np.where (self == lab)

            Line_row = Select[0]; Line_col = Select[1]
            Line = Line2D(Line_col, Line_row)

            Outlines_row.append(Line_row); Outlines_col.append(Line_col)
            Outlines.append(Line)



        twin  = self.copy()

        fig=plt.figure('title', facecolor='White',figsize=[np.floor(twin.shape[1])/5,np.floor(twin.shape[0])/5])

        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')

        Vmin = min(np.amin(twin[twin!=Nodata_value])*0.95, np.amin(twin[twin!=Nodata_value])*1.05)
        Vmax = max(np.amax(twin)*0.95, np.amax(twin)*1.05)

        Map = ax1.imshow(twin, interpolation='None', cmap=plt.cm.gist_earth, vmin=Vmin, vmax=Vmax, alpha = 0.6)
        ax2 = fig.add_axes([0.1, 0.98, 0.85, 0.02])
        scheme = plt.cm.gist_earth; norm = colors.Normalize(vmin=Vmin, vmax=Vmax)
        cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)


        for Line in Outlines:
            ax1.add_line(Line)

        plt.savefig('/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshAnalysis/Example_Data/Output/Figures/'+'TEST'+'.png')


        return Outlines, Outlines_row, Outlines_col
