"""
    This package processes topographic data in order to extract marsh platforms
    Guillaume C.H. Goodwin
    Released unedr GPL3
"""

# Load useful Python packages
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal, osr, gdalconst
from osgeo.gdalconst import *
import cPickle

import pandas as bb

from osgeo import ogr, osr
from shapely.geometry import LineString

import fiona


import LSDMarshPlatform_classes as cl
#---------------------------------------------------------------
def ENVI_raster_binary_to_2d_array(file_name):
    """
    This function transforms a raster into a numpy array.

    Args:
        file_name (ENVI raster): the raster you want to work on.
        gauge (string): a name for your file

    Returns:
        image_array (2-D numpy array): the array corresponding to the raster you loaded
        pixelWidth (geotransform, inDs) (float): the size of the pixel corresponding to an element in the output array.

    Source: http://chris35wills.github.io/python-gdal-raster-io/
    """


    driver = gdal.GetDriverByName('ENVI')

    driver.Register()

    inDs = gdal.Open(file_name, GA_ReadOnly)

    if inDs is None:
        print "Couldn't open this file: " + file_name
        print "Perhaps you need an ENVI .hdr file? "
        sys.exit("Try again!")
    else:
        print "%s opened successfully" %file_name

        #print '~~~~~~~~~~~~~~'
        #print 'Get image size'
        #print '~~~~~~~~~~~~~~'
        cols = inDs.RasterXSize
        rows = inDs.RasterYSize
        bands = inDs.RasterCount

        #print "columns: %i" %cols
        #print "rows: %i" %rows
        #print "bands: %i" %bands

        #print '~~~~~~~~~~~~~~'
        #print 'Get georeference information'
        #print '~~~~~~~~~~~~~~'
        geotransform = inDs.GetGeoTransform()
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]

        #print "origin x: %i" %originX
        #print "origin y: %i" %originY
        #print "width: %2.2f" %pixelWidth
        #print "height: %2.2f" %pixelHeight

        # Set pixel offset.....
        print '~~~~~~~~~~~~~~'
        print 'Convert image to 2D array'
        print '~~~~~~~~~~~~~~'
        band = inDs.GetRasterBand(1)
        print band
        image_array = band.ReadAsArray(0, 0, cols, rows)
        image_array_name = file_name
        print type(image_array)
        print image_array.shape

        return image_array, pixelWidth, (geotransform, inDs)




#---------------------------------------------------------------------
def ENVI_raster_binary_from_2d_array(envidata, file_out, post, image_array):
    """
    This function transforms a numpy array into a raster.

    Args:
        envidata: the geospatial data needed to create your raster
        file_out (string): the name of the output file
        post: coordinates for the goegraphical transformation
        image_array (2-D numpy array): the input raster

    Returns:
        new_geotransform
        new_projection: the projection in which the raster
        file_out (ENVI raster): the raster you wanted

    Source: http://chris35wills.github.io/python-gdal-raster-io/
    """
    driver = gdal.GetDriverByName('ENVI')

    original_geotransform, inDs = envidata

    rows, cols = image_array.shape
    bands = 1

    # Creates a new raster data source
    outDs = driver.Create(file_out, cols, rows, bands, gdal.GDT_Float32)

    # Write metadata
    originX = original_geotransform[0]
    originY = original_geotransform[3]

    outDs.SetGeoTransform([originX, post, 0.0, originY, 0.0, -post])
    outDs.SetProjection(inDs.GetProjection())

    #Write raster datasets
    outBand = outDs.GetRasterBand(1)
    outBand.WriteArray(image_array)

    new_geotransform = outDs.GetGeoTransform()
    new_projection = outDs.GetProjection()

    print "Output binary saved: ", file_out

    return new_geotransform,new_projection,file_out



##########################################################################################################
##########################################################################################################
def Pandas_outline (Surface_array, M_code, scale):
    """
    This magic function from the internet (https://stackoverflow.com/questions/24539296/outline-a-region-in-a-graph) takes an array of positives v. negatives and outlines the border between the two. It gives you a nice Pandas Dataframe at the end.

    Args:

    Returns:

    Author: GCHG
    """

    image = Surface_array
    maskimg = np.zeros(Surface_array.shape, dtype='int')
    maskimg[image == M_code] = 3

    x0 = 0; x1 = Surface_array.shape[1]
    y0 = 0; y1 = Surface_array.shape[0]

    # our image with the numbers 1-3 is in array maskimg
    # create a boolean image map which has trues only where maskimg[x,y] == 3
    mapimg = (maskimg == 3)


    # a vertical line segment is needed, when the pixels next to each other horizontally
    #   belong to diffferent groups (one is part of the mask, the other isn't)
    # after this ver_seg has two arrays, one for row coordinates, the other for column coordinates
    ver_seg = np.where(mapimg[:,1:] != mapimg[:,:-1])
    # the same is repeated for horizontal segments
    hor_seg = np.where(mapimg[1:,:] != mapimg[:-1,:])

    # if we have a horizontal segment at 7,2, it means that it must be drawn between pixels
    #   (2,7) and (2,8), i.e. from (2,8)..(3,8)
    # in order to draw a discountinuous line, we add Nones in between segments
    l = []
    for p in zip(*hor_seg):
        l.append((p[1], p[0]+1))
        l.append((p[1]+1, p[0]+1))
        l.append((np.nan,np.nan))

    # and the same for vertical segments
    for p in zip(*ver_seg):
        l.append((p[1]+1, p[0]))
        l.append((p[1]+1, p[0]+1))
        l.append((np.nan, np.nan))

    # now we transform the list into a numpy array of Nx2 shape
    segments = np.array(l)
    # now we need to know something about the image which is shown
    #   at this point let's assume it has extents (x0, y0)..(x1,y1) on the axis
    #   drawn with origin='lower'
    # with this information we can rescale our points



    if len(segments) > 0:

        segments[:,0] = (x0 + (x1-x0) * segments[:,0] / mapimg.shape[1]) - scale/2.
        segments[:,1] = (y0 + (y1-y0) * segments[:,1] / mapimg.shape[0]) - scale/2.

        #This is me now ^^
        M_outline = stitch_segments(segments,M_code)

    else: M_outline = []


    return M_outline


##########################################################################################################
##########################################################################################################
def stitch_segments (segments, M_code, reset_length = False, select_longest = True):
    """This function stitches tiny segments into something that makes sense"""

    #Initiate the objects
    segments[np.isnan(segments)] = 0

    Original_seg_length = len(segments)


    L_code = 1
    # This is the list of outlines
    M_outlines = cl.Polyline()

    # loop until all the segments are stitched
    while np.amax(segments) != 0:
    #while len(segments) >= 1:
        #print 'These be the segments'
        #print segments [:12]
        #print len(segments)
        #print
        # Setup the pandas outline
        M_outline = cl.Line()


        #Start = timeit.default_timer()

        #Find the first non-null element in the array
        nonzero_x = np.where(segments[:,0]!=0)[0][0]
        nonzero_y = np.where(segments[:,1]!=0)[0][0]

        # choose a this segment to be the first element in the Pandas
        M_outline.set_first_point(M_code, L_code, segments[nonzero_y,1], segments[nonzero_x,0], 1)
        M_outline = M_outline.add_element(1,M_code, L_code, segments[nonzero_y+1,1], segments[nonzero_x+1,0])

        segments = segments[3:]

        # Now stitch the connected segments together in both directions. Because we're 2D, bruv!
        for d in [-1,0]:
            for i in range(3,Original_seg_length,3):
                # Find a segment that starts at the end of the previous segment

                #print M_outline

                last_point = M_outline['rowcol'].iloc[d]
                find_next = np.where(np.logical_and(segments[:,1] == last_point.row(), segments[:,0] == last_point.col()))[0]

                #print find_next

                if len(find_next) > 0:

                    # Add the segment. Invert it if need be.
                    for j in range(len(find_next)):
                        if find_next[j]/3. == find_next[j]//3:
                            M_outline = M_outline.add_element(-d*len(M_outline['rowcol']),M_code, L_code, segments[find_next[j]+1,1], segments[find_next[j]+1,0])

                            #print segments[find_next[j]:find_next[j]+3,:]

                            #segments = np.delete(segments, np.s_[find_next[j]:find_next[j]+3], 0)
                            segments[find_next[j]:find_next[j]+3,:] = 0

                        elif (find_next[j]-1)/3. == (find_next[j]-1)//3:
                            M_outline = M_outline.add_element(-d*len(M_outline['rowcol']),M_code, L_code, segments[find_next[j]-1,1], segments[find_next[j]-1,0])

                            #print segments[find_next[j]-1:find_next[j]+2,:]

                            #segments = np.delete(segments, np.s_[find_next[j]-1:find_next[j]+2], 0)
                            segments[find_next[j]-1:find_next[j]+2,:] = 0

                else:
                    break

            #print
            #print M_outline

            #Stop = timeit.default_timer()
            #time = Stop - Start
            #print 'runtime:', time
            #quit()



        #Recalculate the distances to make it nicer. Takes a loooong time
        if reset_length is True:
            M_outline = M_outline.recalc_length(1)

        #add this outline to the list of outlines
        M_outlines.append(M_outline)

        # update the L_code counter
        print 'Processed Line:', L_code

        #print 'These be the new segments'
        #print segments
        #print len(segments)
        #print

        #print M_outline
        #Stop = timeit.default_timer()
        #time = Stop - Start
        #print 'runtime:', time
        #quit()


        L_code+=1

    #This part is to group everything in one Pandas instead of a silly list
    Final_outlines = M_outlines[0]
    for i in range(1,len(M_outlines)):
        Final_outlines = bb.concat([Final_outlines,M_outlines[i]])

    if select_longest is True:
        Lengths = []
        value_range = [min(Final_outlines['L_code']),max(Final_outlines['L_code'])]
        for L in value_range:
            Lengths.append(len(Final_outlines.loc[Final_outlines['L_code'] == L]))
        Longest = max(Lengths)
        for L in value_range:
            if len(Final_outlines.loc[Final_outlines['L_code'] == L]) < 0.5 * Longest:
                Final_outlines = Final_outlines[Final_outlines.L_code != L]

    return Final_outlines


##########################################################################################################
##########################################################################################################
def Polyline_to_shp (line, Envidata, save_dir, file_name):
    """
    This function takes all these masses of wriggly lines and turns each of them into a shapefile. This is going to take a lot of space, so we should probably think of deleting the shapefiles when we're done with them.

    Args:
        Surface_array (2D numpy array): a 2-D array containing the surface to outline with the value 1. Undesirable elements have the value 0 or Nodata_value.
        Outline_array (2D numpy array): a 2-D array destined to store the outline.
        Outline_value (float): The value to be given to outline cells
        Nodata_value (float): The value for empty cells

    Returns:
        Nothing. It just saves a bunch of shapefiles

    Author: GCHG
    """

    X_origin = Envidata[0][0]; X_cell_width = Envidata[0][1]
    Y_origin = Envidata[0][3]; Y_cell_width = Envidata[0][5]



    spatialReference = osr.SpatialReference() #will create a spatial reference locally to tell the system what the reference will be
    spatialReference.ImportFromProj4('+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=375,-111,431,0,0,0,0 +units=m +no_defs') #here we define this reference to be utm Zone 48N with wgs84...

    # Now convert it to a shapefile with OGR
    driver = ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource('%s/%s.shp' % (save_dir,file_name))
    layer = ds.CreateLayer('', spatialReference, ogr.wkbLineString)
    # Add one attribute
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    defn = layer.GetLayerDefn()


    ## If there are multiple geometries, put the "for" loop here
    #Coord_list = []
    for subline in line:
        x_offset,y_offset = 0,0

        value_range = range(min(subline['L_code']),max(subline['L_code'])+1)
        for L in value_range:
            To_save = subline.loc[subline['L_code'] == L]

            #Initiate list of coordinates
            Coordinates = []
            for i in range(len(To_save['rowcol'])):
                x = To_save['rowcol'].iloc[i].col()
                y = To_save['rowcol'].iloc[i].row()
                coord_pair = (X_cell_width * x + X_origin + x_offset, Y_cell_width * y + Y_origin + y_offset)
                Coordinates.append(coord_pair)

            # Make a Shapely geometry
            poly = LineString(Coordinates)

            # Create a new feature (attribute and geometry)
            feat = ogr.Feature(defn)
            feat.SetField('id', 123)
            # Make a geometry, from Shapely object
            geom = ogr.CreateGeometryFromWkb(poly.wkb)
            feat.SetGeometry(geom)
            layer.CreateFeature(feat)
    #feat = geom = None  # destroy these
    # Save and close everything
    ds = layer = feat = geom = None

    print "file saved:" '%s/%s.shp' % (save_dir,file_name)


##########################################################################################################
##########################################################################################################
def Derivate (x, y):
    x_step = (max(x)-min(x)) / len(x)
    dy = np.zeros(len(y), dtype = np.float)
    for j in range(1, len(y), 1):
        dy[j] = (y[j]-y[j-1])/x_step
    return dy



#-----------------------------------------------------------------------------------------------------------


def Outline (Raster, Outline_value, Nodata_value):
    """
    This simple function takes a 2-D array (Raster) and attributes a specific value (Outline value) to elements at the limit of a bloc of elements with identical values. Effectively, it draws an outline around a group of elements with the same value. It is set to ignore elements with a specific value (Nodata_value).

    Args:
        Raster (2D numpy array): the 2-D array
        Outline_value (float): The value associated to the outline. Be smart and select a different value from those already in your 2-D array.
        Nodata_value (float): The value for ignored elements

    Returns:
        Raster (2D numpy array): the 2-D array, with the outlines given their own value.

    Author: GCHG
    """
    P1 = np.where(Raster[:,1:] != Raster[:,:-1])
    Raster[P1] = Outline_value

    P2 = np.where(Raster[1:,:] != Raster[:-1,:])
    Raster[P2] = Outline_value

    for i in range(len(Raster)):
        for j in range(len(Raster[0,:])):
            if Raster[i,j] == Outline_value:
                K = kernel (Raster, 3, i, j)
                if np.mean(K) < 0:
                    Raster[i,j] = 0

    return Raster





#-----------------------------------------------------------------------------------------------------
def kernel (array, kernel_size, x_centre, y_centre):
    """
    This function defines a square kernel within an array (array), centred on (x_centre, y_centre). The is of a width of kernel_size.
    Args:
        array (2D numpy array): a 2-D array.
        kernel_size (float): the width of the square defining the size of the kernel. kernel_size MUST be an ODD number to account for the central element.
        x_centre (int): The index of the element in the 1st dimension.
        y_centre (int): The index of the element in the 2nd dimension.

    Returns:
        kernel (2D numpy array): The kernel of selected elements.

    Author: GCHG
    """

    if (-1)**kernel_size < 0:
        X_to_0 = x_centre
        X_to_End = len(array)-x_centre
        Y_to_0 = y_centre
        Y_to_End = len(array[0,:])-y_centre

        Lim_left = x_centre - min(np.floor(kernel_size/2), X_to_0)
        Lim_right = x_centre + min(np.floor(kernel_size/2)+1, X_to_End)
        Lim_top = y_centre - min(np.floor(kernel_size/2), Y_to_0)
        Lim_bottom = y_centre + min(np.floor(kernel_size/2)+1, Y_to_End)

        kernel = array [int(Lim_left):int(Lim_right), int(Lim_top):int(Lim_bottom)]

    else:
        print
        print " ... WARNING: you need to choose an odd kernel size, buddy"
        print
        pass

    return kernel




#---------------------------------------------------------------
def MARSH_ID (DEM, Slope, Nodata_value, opt1, opt2, opt3):
    """
    This is the master function for marsh identification. It defines in which order the functions define_search_space, peak_flag, initiate_ridge, Continue_ridge, Clean_ridges, Fill_marsh are executed. It is set to repeat the iteration of the Continue_ridge function 50 times.

    Args:
        DEM (2D numpy array): the input DEM array.
        Slope (2D numpy array): the input Slope array.
        Nodata_value (float): The value for ignored elements.
        opt1 (float): The value of the threshold used in the define_search_space function.
        opt2 (float): The value of the threshold used in the Clean_ridges function.
        opt3 (float): The value of the threshold used in the Fill_marsh function.

    Returns:
        Search_space (2D numpy array): The output search space of the define_search_space function.
        Ridge (2D numpy array): The output ridges of the peak_flag, initiate_ridge, Continue_ridge, Clean_ridges functions.
        Marsh (2D numpy array): The output marsh platform of the Fill_marsh function.

    Author: GCHG
    """

    DEM_work = np.copy(DEM); Slope_work = np.copy(Slope);

    Platform = np.copy(DEM_work)
    Ridge = np.copy(DEM_work)
    Marsh = np.copy(DEM_work)

    Platform[Platform != Nodata_value] = 0
    Summit = np.where (Platform==np.amax(Platform))
    Platform[Summit] = 1


    Search_space, Crossover, bins, hist, Inflexion_point = define_search_space (DEM_work, Slope_work, Nodata_value,opt1)

    Order = 1
    Ridge, Slope_temp = peak_flag (Slope_work, Search_space, Order)

    Order = Order+1
    Ridge, Slope_temp = initiate_ridge (Slope_temp, Search_space, Ridge, Order)

    while Order < 50:
        Order = Order+1
        Ridge, Slope_temp = Continue_ridge (Slope_temp, Search_space, Ridge, Order)

    Ridge = Clean_ridges (Ridge, DEM_work, Nodata_value, opt2)

    Marsh = Fill_marsh (DEM_work, Ridge, Nodata_value, opt3)


    print "My hovercraft is full of eels!"
    print


    return Search_space, Ridge, Marsh




#-----------------------------------------------------------------------------------------------------
def Confusion (Subject, Reference, Nodata_value):
    """
    This function compares a Subject 2-D array to a Reference 2-D array and returns an array of differences, which we call a confusion array or confusion map if it look like a map. It then calculates a number of metrics relative to the adequation between the subject and the reference. It is set to ignore elements with a value of Nodata_value.

    To learn more about confusion matrices and their associated metrics, please visit the Wikipedia page: https://en.wikipedia.org/wiki/Confusion_matrix

    Args:
        Subject (2D numpy array): the input array. This is the one you want to test
        Reference (2D numpy array): the reference array. This one is supposed to contain correct information
        Nodata_value (float): The value for ignored elements.

    Returns:
        Confusion_matrix (2D numpy array): an array containing the values 1 (True Positive), 2 (True Negative), -1 (False Positive) and -2 (False Negative).
        Performance (1D numpy array): the number of (respectively) True Positives, True Negatives, False Positives and False Negatives in Confusion_matrix.
        Metrix (1D numpy array): The values of (respectively) Accuracy, Reliability, Sensitivity, F1 derived from the Performance array.

    Author: GCHG
    """

    Height = len(Subject[:,0]); Width = len(Subject[0,:])
    Height_R = len(Reference[:,0]); Width_R = len(Reference[0,:])

    print Height, Width
    print Height_R, Width_R

    H = min (Height, Height_R)
    W = min (Width, Width_R)

    Confusion_matrix = Nodata_value*np.ones((Height, Width), dtype = np.float)

    Subject_marsh = np.where (np.logical_and(Subject != 0, Subject != Nodata_value))
    Reference_marsh = np.where (np.logical_and(Reference != 0, Reference != Nodata_value))

    Subject[Subject_marsh] = 1.
    Reference[Reference_marsh] = 1.

    for i in range (H):
        for j in range (W):
            if Subject[i,j] == 1 and Reference[i,j] == 1: # TRUE POSITIVE
                Confusion_matrix[i,j] = 1
            elif Subject[i,j] == 0 and Reference[i,j] == 0: # TRUE NEGATIVE
                Confusion_matrix[i,j] = 2
            elif Subject[i,j] == 1 and Reference[i,j] == 0: # FALSE POSITIVE
                Confusion_matrix[i,j] = -1
            elif Subject[i,j] == 0 and Reference[i,j] == 1: # FALSE NEGATIVE
                Confusion_matrix[i,j] = -2

    True_positive = np.sum(Confusion_matrix[Confusion_matrix == 1])
    True_negative = np.sum(Confusion_matrix[Confusion_matrix == 2])/2
    False_positive = -np.sum(Confusion_matrix[Confusion_matrix == -1])
    False_negative = -np.sum(Confusion_matrix[Confusion_matrix == -2])/2

    Reliability = True_positive / (True_positive+False_positive)
    Sensitivity = True_positive / (True_positive+False_negative)
    Accuracy = (True_positive+True_negative) / (True_positive+True_negative+False_positive+False_negative)
    F1 = 2*True_positive/(2*True_positive+False_positive+False_negative)

    Performance = np.array([True_positive,True_negative,False_positive,False_negative])
    Metrix = np.array([Accuracy, Reliability, Sensitivity, F1])


    return Confusion_matrix, Performance, Metrix
