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

#import LSDMOA_functions as fct

import copy

import pandas as bb

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
            if M_outlines.size > 0:
                print '.... Appending label'
                Outlines.append(M_outlines)
        return Outlines









    def plot_map (self, save_dir, figname, title, Nodata_value):
        print ' Plotplotplotplot'
        twin  = self.copy()

        fig_height = min(np.floor(twin.shape[1])/5, 50)
        fig_width = min(np.floor(twin.shape[1])/5, 50)

        fig=plt.figure(title, facecolor='White',figsize=[fig_height,fig_width])

        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')

        Vmin = min(np.amin(twin[twin!=Nodata_value])*0.95, np.amin(twin[twin!=Nodata_value])*1.05)
        Vmax = max(np.amax(twin)*0.95, np.amax(twin)*1.05)

        Map = ax1.imshow(twin, interpolation='None', cmap=plt.cm.gist_earth, vmin=Vmin, vmax=Vmax, alpha = 0.6)
        ax2 = fig.add_axes([0.1, 0.98, 0.85, 0.02])
        scheme = plt.cm.gist_earth; norm = colors.Normalize(vmin=Vmin, vmax=Vmax)
        cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)

        plt.savefig(save_dir+figname+'.png')



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

##########################################################################################################
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


    def Pandaline_to_shp (self, Envidata, Enviarray, save_dir, site_name):
        value_range = range(min(self['L_code']),max(self['L_code'])+1)
        for L in value_range:
            To_save = self.loc[self['L_code'] == L]
            M_code = To_save['M_code'][0]
            L_code = To_save['L_code'][0]
            file_name = str(site_name)+'_'+str(M_code)+'_'+str(L_code)
            fct.Line_to_shp(To_save, Envidata, Enviarray, save_dir, file_name)


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

            print in_name

            fct.Make_transects(save_dir+'Shapefiles/'+in_name, save_dir+'Shapefiles/'+out_name, spacing, length)

            print out_name


            print



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
        """Only use if you already have a Marsh attribute"""
        if self.size > 0:
            Select_list = []
            step = max(self.index)+1
            for i in range(0,len(self['rowcol']),step):
                Start = i; End = i+step
                Zeros = np.where(np.asarray(self['Marsh'].iloc[Start:End]) < 0.5)
                Ones = np.where(np.asarray(self['Marsh'].iloc[Start:End]) > 0.5)
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



##########################################################################################################
##########################################################################################################
class Polyline (list):
    """A polyline object has Line or Polyline objects inside. Let's say for now it just has Lines."""
    def __init__(self):
        list.__init__(self)


    def Save_as_pickle (self,filepath,filename):
        with open(filepath+filename, 'wb') as handle:
            pickle.dump(self,handle)


    def save_to_shp (self, Envidata, Enviarray, save_dir, site_name):
        for i in range(len(self)):
            Outline = self[i]
            Outline.Pandaline_to_shp(Envidata, Enviarray, save_dir, site_name)


    def Polyline_transects (self, spacing, length, Envidata, Enviarray, save_dir, site_name, needs_saving = False):
        #if you haven't done so yet, save your polyline to a .shp
        if needs_saving is True:
            self.save_to_shp (Envidata, Enviarray, save_dir, site_name)
        #Make a Polyline to store all those transects
        All_transects = Polyline()
        # For each Line in the Polyline
        for i in range(len(self)):
            Transects = self[i].Line_transects(spacing, length, Envidata, Enviarray, save_dir, site_name)
            All_transects.append(Transects)
        return All_transects


    def get_attribute_from_basemap (self, refinement, basemap, attr_name, Nodata_value):
        for i in range(len(self)):
            if self[i].size > 0:
                self[i] = self[i].multiple_subdivide(refinement)
                self[i] = self[i].extract_values_from_basemap (basemap, attr_name, Nodata_value)
                if attr_name == 'Z':
                    self[i] = self[i].orient_seaward()
                    self[i] = self[i].slope_1D()
                    #print self[i]
                    self[i] = self[i].get_bearing()
                    #print self[i]
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

        ax3.set_xlim(0,20)
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

        ax3.axvline(10, color='black', lw=1.0, alpha=0.8)

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
