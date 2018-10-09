"""
Example_Plots.py


    This contains wrapper functions that simplify plotting raster
    and vector data for publication-ready figures.

    Guillaume CH Goodwin, Simon Mudd and Fiona Clubb, June 2017
    Released under GPL3

"""


#------------------------------------------------------------------
#0. Set up display environment in putty if you are working on a terminal with no graphical interface.
import matplotlib
matplotlib.use('Agg')

#----------------------------------------------------------------
#1. Load useful Python packages
import os
import sys

import numpy as np
import functools
import math as mt
import cmath
import scipy as sp
import scipy.stats as stats
from datetime import datetime
import cPickle
from pylab import *
import functools
import itertools as itt
from osgeo import gdal, osr
from osgeo import gdal, gdalconst
from osgeo.gdalconst import *
from copy import copy
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.patches import Rectangle
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import timeit



import sphviewer as sph
#------------------------------------------------------------------
# Import the marsh-finding functions
from LSDMarshPlatform_functions import ENVI_raster_binary_to_2d_array
from LSDMarshPlatform_functions import ENVI_raster_binary_from_2d_array
from LSDMarshPlatform_functions import kernel
from LSDMarshPlatform_functions import Outline
#from LSDMarshPlatform_functions import Distribution


#------------------------------------------------------------------
#2. Set up the important variables

# Name your data input directory
#Input_dir = "//csce.datastore.ed.ac.uk/csce/geos/users/s1563094/Software/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/"
# Name your results output directory
#Output_dir = "//csce.datastore.ed.ac.uk/csce/geos/users/s1563094/Software/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/"
# NB: When naming your work directories, make sure the syntax is compatible with the OS you are using. Above is an example in the Windows command shell



#Select site names. Simply add a site in the list to analyse multiple sites simultaneously.
Sites = ["FEL"]

# Set the value for empty DEM cells
Nodata_value = -9999


################################################################################
################################################################################
def list_all_directories (directory, dir_list):
    """
    This function looks into a directory and checks:
    - whether there are subdirectories
    - what these subdirectories are

    args:
    - basedir (String): the directory

    Returns:
    - dir_list (List): a list of the directories in the base directory, in order of finding

    NB: this is a recursive function

    """

    for root, dirs, files in os.walk(directory, topdown=True):
        if len(dirs) > 0:
            keep_going = True
            for dir in dirs:
                full_dir =  directory + dir + '/'
                if os.path.isdir(full_dir):
                    dir_list.append(full_dir)
                    list_all_directories(full_dir, dir_list)
        else:
            keep_going = False

    return dir_list



################################################################################
################################################################################
def Figure1 (Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This plots the extracted marsh platform on a hillshade

    Args:
        Input_dir (str): Name your data input directory
        Output_dir (str): Name your results output directory
        Sites (str list): A list of strings. The file names are modified based on these sites

    Author: GCHG
    """
    #Plot 1: Draw the platform on a DEM, superimposed on a hillshade
    for site in Sites:
        fig=plt.figure(1, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)

        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_hs_clip.bil" % (site))
        DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (Input_dir+"%s_DEM_clip.bil" % (site))
        Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array (Output_dir+"Marsh.bil")

        # Make a mask!
        Platform_mask = np.ma.masked_where(Platform <=0, Platform)
        Platform_mask[Platform_mask>0] = DEM[Platform_mask>0]
        DEM_mask = np.ma.masked_where(DEM <= Nodata_value, DEM)
        HS_mask = np.ma.masked_where(HS <= Nodata_value, HS)

        # Make a map!
        Map_HS = ax1.imshow(HS_mask, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(HS_mask), vmax = np.amax(HS_mask))
        Map_DEM = ax1.imshow(DEM_mask, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM_mask[DEM_mask!=Nodata_value]), vmax = np.amax(DEM_mask), alpha = 0.5)
        Map_Marsh = ax1.imshow(Platform_mask, interpolation='None', cmap=plt.cm.gist_earth, vmin=np.amin(DEM[DEM!=Nodata_value]), vmax=np.amax(DEM), alpha = 0.5)

        plt.savefig(Output_dir[:-1]+'_Marsh_DEM.png')



################################################################################
################################################################################
def Figure2 (Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This plots the extracted marsh platform on a hillshade

    Args:
        Input_dir (str): Name your data input directory
        Output_dir (str): Name your results output directory
        Sites (str list): A list of strings. The file names are modified based on these sites

    Author: GCHG
    """

    #Plot 1: Draw the platform on a DEM, superimposed on a hillshade
    for site in Sites:
        site_dir = Input_dir[:-5]
        dirlist = list_all_directories (site_dir, [])

        lastyear = dirlist[-1]
        last_DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (lastyear+"%s_DEM_clip.bil" % (site))
        last_HS, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (lastyear+"%s_hs_clip.bil" % (site))
        last_Marsh, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (lastyear+"Marsh.bil")
        last_Marsh[last_Marsh > 0] = 1

        for dir in dirlist[:-1]:
            if os.path.isfile(dir+"Marsh.bil"):
                first_DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (dir+"%s_DEM_clip.bil" % (site))
                first_Marsh, post_DEM, envidata_DEM  = ENVI_raster_binary_to_2d_array (dir+"Marsh.bil")
                first_Marsh[first_Marsh > 0] = 1
                firstyear = dir
                break

        Marsh_diff = last_Marsh - first_Marsh
        Marsh_diff[last_DEM == Nodata_value] = Nodata_value
        Marsh_diff[first_DEM == Nodata_value] = Nodata_value

        Constant_marsh = np.where(np.logical_and(last_Marsh == 1, Marsh_diff == 0))
        Marsh_diff[Constant_marsh] = 2

        DEM_diff = last_DEM - first_DEM
        DEM_diff[last_DEM == Nodata_value] = Nodata_value
        DEM_diff[first_DEM == Nodata_value] = Nodata_value

        DEM_diff [DEM_diff>100] = Nodata_value
        DEM_diff [DEM_diff<-100] = Nodata_value

        DEM_diff_erosion = 0 * np.ndarray((last_DEM.shape[0], last_DEM.shape[1]), dtype =np.float) + Nodata_value
        DEM_diff_erosion [Marsh_diff == -1] = DEM_diff[Marsh_diff == -1]


        DEM_diff_accretion = 0 * np.ndarray((last_DEM.shape[0], last_DEM.shape[1]), dtype =np.float) + Nodata_value
        DEM_diff_accretion [Marsh_diff == 1] = DEM_diff[Marsh_diff == 1]


        DEM_diff_stable = 0 * np.ndarray((last_DEM.shape[0], last_DEM.shape[1]), dtype =np.float) + Nodata_value
        DEM_diff_stable [Marsh_diff == 2] = DEM_diff[Marsh_diff == 2]



        fig=plt.figure(1, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((10,10),(0,0),colspan=6, rowspan=6)
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        ax2 = plt.subplot2grid((10,10),(0,6),colspan=4, rowspan=10)
        ax2.set_xlabel('state (m)', fontsize = 12)
        ax2.set_ylabel('dz (m)', fontsize = 12)

        ax3 = plt.subplot2grid((10,10),(6,0),colspan=6, rowspan=4)
        ax3.set_xlabel('z (m)', fontsize = 12)
        ax3.set_ylabel('dz (m)', fontsize = 12)

        # Make a mask!
        DEM_diff_erosion_mask = np.ma.masked_where(DEM_diff_erosion == Nodata_value, DEM_diff_erosion)
        DEM_diff_accretion_mask = np.ma.masked_where(DEM_diff_accretion == Nodata_value, DEM_diff_accretion)
        DEM_diff_stable_mask = np.ma.masked_where(DEM_diff_stable == Nodata_value, DEM_diff_stable)


        E = [abs(np.amin(DEM_diff_erosion_mask)), np.amax(DEM_diff_erosion_mask)]
        A = [abs(np.amin(DEM_diff_accretion_mask)), np.amax(DEM_diff_accretion_mask)]
        S = [abs(np.amin(DEM_diff_stable_mask)), np.amax(DEM_diff_stable_mask)]


        # Make a map!
        #Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
        #Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
        Map_Marsh = ax1.imshow(-DEM_diff_erosion_mask, interpolation='None', cmap=plt.cm.seismic, vmin=-max(E), vmax=max(E), alpha = 0.8)
        Map_Marsh = ax1.imshow(DEM_diff_accretion_mask, interpolation='None', cmap=plt.cm.PiYG, vmin=-max(A), vmax=max(A), alpha = 0.8)
        Map_Marsh = ax1.imshow(DEM_diff_stable_mask, interpolation='None', cmap=plt.cm.RdGy, vmin=-max(S), vmax=max(S), alpha = 0.5)



        # Second plot
        re = ax2.violinplot(DEM_diff_erosion[DEM_diff_erosion!=Nodata_value], [1], points=60, widths=0.3, showmeans=True, showextrema=False, showmedians=False)
        re['cmeans'].set_color('r')
        for vp in re['bodies']:
            vp.set_facecolor('r')
            vp.set_alpha(0.8)

        ra = ax2.violinplot(DEM_diff_accretion[DEM_diff_accretion!=Nodata_value].ravel(), [-1], points=60, widths=0.3, showmeans=True, showextrema=False, showmedians=False)
        ra['cmeans'].set_color('g')
        for vp in ra['bodies']:
            vp.set_facecolor('g')
            vp.set_alpha(0.8)


        rs = ax2.violinplot(DEM_diff_stable[DEM_diff_stable!=Nodata_value].ravel(), [0], points=60, widths=0.3, showmeans=True, showextrema=False, showmedians=False)
        rs['cmeans'].set_color('k')
        for vp in rs['bodies']:
            vp.set_facecolor('k')
            vp.set_alpha(0.8)


    # Third plot
    X = last_DEM[DEM_diff_erosion!=Nodata_value]
    Y = DEM_diff_erosion[DEM_diff_erosion!=Nodata_value]
    #ax3.scatter(X.ravel(), Y.ravel(), c='r', alpha=0.1, linewidth = 0)

    X = last_DEM[DEM_diff_accretion!=Nodata_value]
    Y = DEM_diff_accretion[DEM_diff_accretion!=Nodata_value]
    #ax3.scatter(X.ravel(), Y.ravel(), c='g', alpha=0.1, linewidth = 0)

    X = last_DEM[DEM_diff_stable!=Nodata_value]
    Y = DEM_diff_stable[DEM_diff_stable!=Nodata_value]
    #ax3.scatter(X.ravel(), Y.ravel(), c='k', alpha=0.1, linewidth = 0)


    where = np.where(np.logical_and(Marsh_diff>=-1, Marsh_diff != 0))
    x = last_DEM[where].ravel()#; x = x[x!=Nodata_value]
    y = DEM_diff[where].ravel()#; y = y[y!=Nodata_value]

    heatmap_32, extent_32 = myplot(x,y, nb=32)



    map = ax3.imshow(heatmap_32, extent=extent_32, origin='lower', aspect='auto')
    #ax3.set_title("Smoothing over 32 neighbors")




    plt.savefig(Output_dir[:-5]+'diff_%s_%s_Marsh_DEM.png' % (firstyear[-5:-1], lastyear[-5:-1]))





################################################################################
################################################################################
def Plot_platform_on_hillshade(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This plots the extracted marsh platform on a hillshade

    Args:
        Input_dir (str): Name your data input directory
        Output_dir (str): Name your results output directory
        Sites (str list): A list of strings. The file names are modified based on these sites

    Author: GCHG
    """
    #Plot 1: Draw the platform on a DEM, superimposed on a hillshade
    for site in Sites:
        fig=plt.figure(1, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)

        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_hs.bil" % (site), site)
        DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (Input_dir+"%s.bil" % (site), site)
        Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)


        # Make a mask!
        Platform_mask = np.ma.masked_where(Platform <=0, Platform)
        Platform_mask[Platform_mask>0] = DEM[Platform_mask>0]

        # Make a map!
        #Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
        Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
        #Map_Marsh = ax1.imshow(Platform_mask, interpolation='None', cmap=plt.cm.gist_earth, vmin=np.amin(DEM[DEM!=Nodata_value]), vmax=np.amax(DEM), alpha = 0.5)

        plt.savefig(Output_dir+'Platform_DEM_%s.png' % (site))





def Plot_Elevation_PDF(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This plots the extracted marsh platform on a hillshade

    Args:
        Input_dir (str): Name your data input directory
        Output_dir (str): Name your results output directory
        Sites (str list): A list of strings. The file names are modified based on these sites

    Author: GCHG
    """
    #Plot 1: Draw the platform on a DEM, superimposed on a hillshade
    for site in Sites:
        fig=plt.figure(1, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)

        # Name the axes
        ax1.set_xlabel('Elevation (m)', fontsize = 12)
        ax1.set_ylabel('Probability Distribution (m)', fontsize = 12)

        # Load the relevant data
        DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (Input_dir+"%s.bil" % (site), site)
        Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)
        Platform [Platform>0] = DEM [Platform>0]

        bins_z, hist_z = Distribution (DEM, Nodata_value)
        bins_M, hist_M = Distribution (Platform, Nodata_value)

        Elevation_range_z = np.arange(min(bins_z[bins_z!=Nodata_value]), max(bins_z), 0.1)
        Elevation_range_M = np.arange(min(bins_z[bins_z!=Nodata_value]), max(bins_M), 0.1)
        Ratio = (max(hist_z[hist_z < 0.2])/max(hist_M[hist_M < 0.2]))
        hist_z_copy = hist_z / Ratio
        hist_M[0] = 0


        ax1.fill_between( bins_z, -5, hist_z_copy, color=plt.cm.gray(0), alpha = 0.5, linewidth = 0.0)
        ax1.plot( bins_M, hist_M, '-r', linewidth = 2.0)


        # Set the ticks
        A = 0.01
        #for x in range(len(hist_M)-1):
            #if hist_M[x]==0 and hist_M[x+1]>0:
                #A = bins_M[x]
                #break
        #xmin = max(-5,A)
        ymax = max(hist_M[hist_M<0.2])

        #ax1.set_xlim (xmin = xmin)
        ax1.set_ylim (ymin = 0, ymax = ymax*1.05)

        #majorLocator1 = MultipleLocator(np.floor(100*ymax)/200)
        #ax1.yaxis.set_major_locator(majorLocator1)
        #majorLocator2 = MultipleLocator(1)
        #ax1.xaxis.set_major_locator(majorLocator2)


        plt.savefig(Output_dir+'Elevation_PDF_%s.png' % (site))






def Plot_marsh_outline_on_hillshade(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This draws the marsh outline on a DEM, superimposed on a hillshade.

    Args:
        Input_dir (str): Name your data input directory.
        Output_dir (str): Name your results output directory.
        Sites (str list): A list of strings. The file names are modified based on these sites.

    Author: GCHG
    """

    for site in Sites:
        fig=plt.figure(2, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)

        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_DEM_hs.bil" % (site), site)
        DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (Input_dir+"%s.bil" % (site), site)
        Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)

        # Outline the marsh
        Platform[Platform > 0] = 1
        Marsh_outline = Outline (Platform, 2, Nodata_value)


        # Make a mask!
        Outline_mask = np.ma.masked_where(Marsh_outline <=1, Marsh_outline)

        # Make a map!
        Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
        Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
        Map_Marsh = ax1.imshow(Outline_mask, interpolation='None', cmap=plt.cm.Reds, vmin = 0, vmax = 2, alpha = 1)


    plt.savefig(Output_dir+'Platform_outline_%s.png' % (site))


def Plot_marsh_reference_on_hillshade(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This draws the marsh reference on a hillshade

    Args:
        Input_dir (str): Name your data input directory.
        Output_dir (str): Name your results output directory.
        Sites (str list): A list of strings. The file names are modified based on these sites.

    Author: GCHG
    """
    #Plot 3: Draw the marsh map and reference outline, superimposed on a hillshade
    for site in Sites:
        fig=plt.figure(3, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)


        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_hs.bil" % (site), site)
        DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (Input_dir+"%s.bil" % (site), site)
        Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)
        Reference, post_Reference, envidata_Reference = ENVI_raster_binary_to_2d_array (Input_dir+"%s_ref.bil" % (site), site)

        # Outline the reference
        Reference[Reference > 0] = 1
        Ref_outline = Outline (Reference, 2, Nodata_value)


        # Make a mask!
        Outline_mask = np.ma.masked_where(Ref_outline <=1, Ref_outline)


        # Make a map!
        Platform_mask = np.ma.masked_where(Platform <=0, Platform)
        Platform_mask[Platform_mask>0] = DEM[Platform_mask>0]

        Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
        Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
        Map_Marsh = ax1.imshow(Platform_mask, interpolation='None', cmap=plt.cm.gist_earth, vmin=np.amin(DEM[DEM!=Nodata_value]), vmax=np.amax(DEM), alpha = 0.5)

        Map_Marsh = ax1.imshow(Outline_mask, interpolation='None', cmap=plt.cm.Reds, vmin = 0, vmax = 2, alpha = 1)

    plt.savefig(Output_dir+'Reference_map_%s.png' % (site))



def Plot_confusion_map_on_hillshade(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This draws the marsh reference on a hillshade

    Args:
        Input_dir (str): Name your data input directory.
        Output_dir (str): Name your results output directory.
        Sites (str list): A list of strings. The file names are modified based on these sites.

    Author: GCHG
    """
    #Plot 4: Draw the confusion map, superimposed on a hillshade
    for site in Sites:
        fig=plt.figure(4, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)


        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_DEM_hs.bil" % (site), site)
        Confusion, post_Confusion, envidata_Confusion = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Confusion.bil" % (site), site)


        # Make a mask!
        Confusion_mask = np.ma.masked_where(Confusion == Nodata_value, Confusion)


        # Make a map!
        Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
        Map_DEM = ax1.imshow(Confusion_mask, interpolation='None', cmap=plt.cm.Spectral, vmin = -2, vmax = 2, alpha = 0.5)


    plt.savefig(Output_dir+'Confusion_%s.png' % (site))
