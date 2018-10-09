"""
LSDMarshPlatform_Marsh_ID.py

This is your driver file to run the marsh platform extraction.

Please read the README and the instructions in this script before you run it.

Authors: Guillaume CH Goodwin and Simon Marius Mudd

"""

#------------------------------------------------------------------
#0. Set up display environment if you are working on a terminal with no GUI.
import matplotlib
matplotlib.use('Agg')

#------------------------------------------------------------------

# Useful Python packages
import numpy as np
import cPickle as pickle
import timeit
import os

# A very useful package
from LSDMarshPlatform_functions import ENVI_raster_binary_to_2d_array
from LSDMarshPlatform_functions import ENVI_raster_binary_from_2d_array

# The main functions for the marsh identification
from LSDMarshPlatform_functions import *
from LSDMarshPlatform_classes import *

# Retained directories from Guillaume
# "//csce.datastore.ed.ac.uk/csce/geos/users/s1563094/Software/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/"


def MarshID(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL_DEM_clip"], opt1 = -2.0, opt2 = 0.85, opt3 = 8.0, compare_with_digitised_marsh = False):
    """
    This function wraps all the marsh ID scripts in one location

    Args:
        Input_dir (str): Name your data input directory
        Output_dir (str): Name your results output directory
        Sites (str list): A list of strings. The file names are modified based on these sites
        opt1 (flt): first optimisation
        opt2 (flt): 2nd optimisation
        opt3 (flt): 3rd optimisation
        compare_with_digitised_marsh (bool): If true, this will compare the data with a digitised marsh platform

    Author:
        GCHG, Modified by SMM 02/10/2017
    """
    #------------------------------------------------------------------
    print("Welcome to the marsh ID program!")
    print("I am opening the file: "+Input_dir)

    # Set the value for empty DEM cells
    Nodata_value = -9999

    # Timing
    Start = timeit.default_timer()


    for site in Sites:

        # If the analysis has not already been done:
        if not os.path.isfile(Input_dir+site+"_Full_profile.shp"):

            print("Loading input data from site: "+site)
            # When loading input data, please make sure the naming convention shown here is respected.
            print(" Loading DEM")
            DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array (Input_dir+"%s_DEM_clip.bil" % (site))
            #HS, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array (Input_dir+"%s_hs_clip.bil" % (site))
            #Make a proper object for the DEM
            DEM_object = Land_surface(DEM.shape[0], DEM.shape[1])
            DEM_object = DEM_object.set_attribute (DEM, Nodata_value, DEM, Nodata_value, classification = False)

            print " Loading Slopes"
            # check to get the correct slope raster
            slope_fname = site+"_slope_clip.bil"
            if not os.path.isfile(Input_dir+slope_fname):
                slope_fname = site+"_SLOPE_clip.bil"
            Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array (Input_dir+slope_fname)
            #Make a proper object for the DEM
            Slope_object = Land_surface(DEM.shape[0], DEM.shape[1])
            Slope_object = Slope_object.set_attribute (Slope, Nodata_value, Slope, Nodata_value, classification = False)

            # This is the search space
            SS_0 = DEM_object.define_search_space(DEM_object, Slope_object, Nodata_value, -2.0)


            # If there is no corrected scarps file:
            if not os.path.isfile(Input_dir+site+"_Scarps_C.shp"):
                print "Extracting scarps"
                # Potential for terraces
                SS_1 = SS_0 * DEM_object
                Scarps = SS_0.extract_full_scarps (DEM_object, Slope_object, Nodata_value, 0.85)

            else:
                #run the scarps extraction
                print "Loading scarps"
                # Convert Scarps shp to raster (with GDAL).
                # Make sure it's in the same place and has the same dimensions
                xmin = str(envidata_DEM [0][0])
                ymax = str(envidata_DEM [0][3])
                xmax = str(envidata_DEM [0][0]+DEM.shape[1])
                ymin = str(envidata_DEM [0][3]-DEM.shape[0])

                os.system("gdal_rasterize -of ENVI -burn 1 -tr 1 1 -te " + xmin + " " + ymin +" " + xmax +" " + ymax + " " + Input_dir+site +"_Scarps_C.shp "+ Input_dir+site +"_Scarps_C.bil")

                Scarps_array, post_sc, envidata_sc =  ENVI_raster_binary_to_2d_array (Input_dir+"%s_Scarps_C.bil" % (site))

                Scarps = Scarps_arr(Scarps_array.shape[0], Scarps_array.shape[1])
                Scarps = Scarps.set_attribute (Scarps_array, 1, DEM, Nodata_value, classification = True)


            # Then make the platforms and simplified scarps
            Scarps, Marsh = Scarps.extract_platforms(DEM_object, SS_0, Nodata_value, 8)
            Marsh_object = Marsh_platform(Marsh.shape[0], Marsh.shape[1])
            Marsh_object = Marsh_object.set_attribute (Marsh, 1, DEM, Nodata_value, classification = True)
            ENVI_raster_binary_from_2d_array(envidata_DEM, Input_dir+site+"_Marsh.bil", post_DEM, Marsh)
            print "saving the platform edges"
            Marsh_object.edges_to_shp(Input_dir, site, envidata_DEM, Nodata_value)

            #os.system("rm " + Input_dir + site + "_Scarps_C.*")

            #quit()

            # Final step, get rid of corrected scarps?





















        if os.path.isfile(Input_dir+"Digitised.shp"):
            # This is the complete ridges
            #Ridges_binary = Land_surface(DEM.shape[0], DEM.shape[1])
            #Ridges_binary = Ridges_binary.set_attribute (Scarps, Nodata_value, Scarps, Nodata_value, classification = False)
            #Ridges_binary[Ridges_binary > 0] = 1
            # Potential for terraces
            #Ridges_DEM = Ridges_binary * DEM
            CRS = 'EPSG:27700' # This is WGS84
            print 'Comparing to digitised marsh'
            os.system("gdal_rasterize -a id -tr 1 1 -a_srs " + CRS + " " + Input_dir+"Digitised.shp " + Input_dir+"Digitised.tif")
            print 'clipping'
            src = Input_dir+"Digitised.tif"
            dst = Input_dir+"Digitised.bil"
            domain = Input_dir+"domain.shp"
            #os.system("gdalwarp -overwrite -of ENVI -t_srs " + CRS + " -cutline " + domain + " -crop_to_cutline " + src + " " + dst)

            os.system("gdalwarp -overwrite -of ENVI -t_srs " + CRS + " -cutline " + domain + " -crop_to_cutline " + src + " " + dst)
            # When loading input data, please make sure the naming convention shown here is respected.
            Platform_work, post_Platform, envidata_Platform =  ENVI_raster_binary_to_2d_array (Input_dir+"Marsh.bil")
            Reference, post_Reference, envidata_Reference =  ENVI_raster_binary_to_2d_array (Input_dir+"Digitised.bil")
            Confusion_matrix, Performance, Metrix = Confusion (Platform_work, Reference, Nodata_value)
            ENVI_raster_binary_from_2d_array (envidata_Platform, Input_dir+"Confusion.bil" , post_Platform, Confusion_matrix)

            #cPickle.dump(Performance,open(Output_dir+"Performance.pkl", "wb"))
            #cPickle.dump(Metrix,open(Output_dir+"Metrix.pkl", "wb"))




    # Comment these 2 lines if you don't want to know how long the script run for.
    Stop = timeit.default_timer()
    print 'Runtime = ', Stop - Start , 's'


    quit()



    Marsh, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array (Input_dir+"Marsh.bil")
    Marsh_object = Marsh_platform(Marsh.shape[0], Marsh.shape[1])
    Marsh_object = Marsh_object.set_attribute (Marsh, Nodata_value, Marsh, Nodata_value, classification = False)



    Marsh_binary = Marsh_platform(DEM.shape[0], DEM.shape[1])
    Marsh_binary = Marsh_binary.set_attribute (Marsh_object, Nodata_value, Marsh_object, Nodata_value, classification = False)
    Marsh_binary[Marsh_binary > 0] = 1
    Marsh_binary[Marsh_binary <= 0] = 0

    # Potential for terraces
    Marsh_DEM = Marsh_binary * DEM_object

    print np.amax(Marsh_binary), np.amin(Marsh_binary)
    print np.amax(DEM_object), np.amin(DEM_object)
    print np.amax(Marsh_DEM), np.amin(Marsh_DEM)
    print Marsh_binary[0,0], DEM_object[0,0], Marsh_DEM[0,0]


    Marsh_DEM[Marsh_DEM <= 0] = Nodata_value

    Marsh_DEM.plot_map (HS, Input_dir[:-5], Input_dir[-5:-1]+site+'Marsh', 'NoTitle', Nodata_value)



    #Marsh_DEM[Scarps > 0] = Nodata_value

    #Marsh_DEM.plot_map (Input_dir[:-5], Input_dir[-5:-1]+site+'Marsh_noscarp', 'NoTitle', Nodata_value)

    DEM_work = np.copy(DEM)
    Search_space, Scarps, Platform = MARSH_ID(DEM, Slope, Nodata_value, opt1, opt2, opt3)
    Platform_work = np.copy(Platform)
    Scarps[Scarps == 0] = Nodata_value

    # Here is where you save your output files for use in a GIS software
    print "Saving marsh features"
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, Output_dir+"%s_Search_space.bil" % (site), post_DEM, Search_space)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, Output_dir+"%s_Scarps.bil" % (site), post_DEM, Scarps)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, Output_dir+"%s_Marsh.bil" % (site), post_DEM, Platform)


    # Disable the following section if you do not wish to compare your results to a reference marsh
    if compare_with_digitised_marsh:
        # When loading input data, please make sure the naming convention shown here is respected.
        print " Loading detected Marsh"
        Platform_work, post_Platform, envidata_Platform =  ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)
        print "Loading reference marsh"
        Reference, post_Reference, envidata_Reference =  ENVI_raster_binary_to_2d_array (Input_dir+"%s_ref.bil" % (site), site)
        print "Evaluating the performance of the detection"
        Confusion_matrix, Performance, Metrix = Confusion (Platform_work, Reference, Nodata_value)
        new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_Platform, Output_dir+"%s_Confusion.bil" % (site),                                                                                post_Platform, Confusion_matrix)
        cPickle.dump(Performance,open(Output_dir+"%s_Performance.pkl" % (site), "wb"))
        cPickle.dump(Metrix,open(Output_dir+"%s_Metrix.pkl" % (site), "wb"))



    Scarps, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array (Input_dir+"Scarps.bil")
    Scarp_object = Scarps_arr(Scarps.shape[0], Scarps.shape[1])
    Scarp_object = Scarp_object.set_attribute (Scarps, Nodata_value, Scarps, Nodata_value, classification = False)
