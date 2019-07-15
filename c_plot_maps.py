#-------------------------------------------------------------------------------
# Name:        plot_maps
# Purpose:     different plots based on maps of BaWue, Bayern, Hessen and
#              Rheinland-Pfalz
#
# Author:      mueller
#
# Created:     26.03.2013
# Copyright:   (c) mueller 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from matplotlib import pyplot as plt
from matplotlib import cm
import mpl_toolkits.basemap.pyproj as pyproj
from mpl_toolkits.basemap import Basemap, shiftgrid
from mpl_toolkits.basemap import addcyclic
import matplotlib as mpl
from scipy.ndimage.filters import minimum_filter, maximum_filter
from matplotlib.mlab import griddata
import numpy as np
import os
import pandas as pd
import a_functions


# ###################################################
# http://matplotlib.org/basemap/users/examples.html #
# ###################################################

class Plot_Map(object):
    r"""
    General plotting tool

    Init(): Generates an object including basemap and matplotlib objects
    ---------
        Keywords: border, res, title, **kwargs
            border: optional
                defines borders and projection of the screen:
                ['bw', 'by', 'he', 'rp','ni', 'all', 'euro_atl', 'south_germany']
            res: string, optional
                resolution coarse to fine: 'c','l','i','h','f', default = 'c'
            title: string, optional
                title of figure, default = None

    Call(): Shortcut for saving/showing figure
    ---------
        Keywords: savefig
            savefig: string, optional
                defines path, filename and file extension of figure
                if not given, figure is shown

    Other Functions()
    ---------
        plot_scatter()
        plot_raster()
        plot_text()
        plot_contours()

    Remark:
    -------
    gk-coordinates can be transformed into lons/lats using
    plot_map.gk_to_lonlat


    **kwargs
    --------
    key word arguments controlling the map and basemap appearance:

    boarder_line: bolean
        draw contour of states, default = True
    paral: bolean or dict used by mpl_toolkits.basemap.drawmeridians
        draw parallels, default = True
    merid: bolean or dict for mpl_toolkits.basemap.drawmeridians
        draw meridians, default = True
    backg: bolean
        draw altitude map as background, default = True
    """

    def __init__(self,
                 border = 'all',
                 res = 'c',
                 figsize = (10,10),
                 title = None,
                 basepath = r'P:\2017_SYNOPSE_II',
                 **kwargs):

        #get NiedSim basepath
        self.basepath = basepath
        self.fig = plt.figure(figsize = figsize)

        #changing this, changes the space of the basemap plot inside the figure
        if title:
            self.orig_ax_position = [.05,.05,.95,.9]
            self.ax = plt.axes(self.orig_ax_position)
            self.ax.set_title(title, fontsize = 20)
        else:
            self.orig_ax_position = [.1,.01,.9,.99]
            self.ax = plt.axes(self.orig_ax_position)

        #specify area and projection
        lons_b, lats_b, lon_0, lat_0, projection, gk = self.get_borders(border)


        # llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
        # are the lat/lon values of the lower left and upper right corners
        # of the map.
        # resolution = 'i' means use intermediate resolution coastlines.
        # lon_0, lat_0 are the central longitude and latitude of the projection.

        self.m = Basemap(llcrnrlon = lons_b[0],
                         llcrnrlat = lats_b[0],
                         urcrnrlon = lons_b[1],
                         urcrnrlat = lats_b[1],
                         resolution = res,
                         projection = projection,
                         lon_0 = lon_0,
                         lat_0 = lat_0,
                         ax = self.ax)
        self._manipulate_map(border, gk, **kwargs)





    def __call__(self, savefig = None):
        if savefig!=None:
            self.fig.savefig(savefig)
        else:
            self.fig.show()
        ##plt.close('all')


    def _manipulate_map(self,border,gk, **kwargs):
        r"""
        Manipulate basemap

        Parameters
        ----------
        m: basemap object
            object containing basemap information
        boarder_line: bolean
            draw contour of states, default = True
        paral: bolean or dict used by mpl_toolkits.basemap.drawmeridians
            draw parallels, default = True
        merid: bolean or dict for mpl_toolkits.basemap.drawmeridians
            draw meridians, default = True
        backg: bolean, default = True
            draw altitude background map


        Returns
        -------
        None

        """


        if kwargs.get('backg', True):
            grid = os.path.join(self.basepath, '02_Import',
                                '03_QGIS_project_germany', 'xyz_germany_utm.dat')
            coord = pd.read_csv(grid, sep = ';')
            lons, lats = a_functions.utm32n_to_lon_lat(coord.ix[:,0].values,
                                                       coord.ix[:,1].values)
            alpha = kwargs.get('alpha', 1.)
            x, y = self.m(lons,lats)
            z = coord.ix[:,2].values
            color=cm.terrain
            scatter = self.m.scatter(x,y,c=z,
                                     s=30, marker = 's', cmap=color,
                                     vmin = 0, vmax = 2000, alpha = alpha)

            if kwargs.get('cb', False):
                # ColorbarBase derives from ScalarMappable and puts a colorbar
                # in a specified axes, so it has everything needed for a
                # standalone colorbar.  There are many more kwargs, but the
                # following gives a basic continuous colorbar with ticks
                # and labels.
                box = self.ax.get_position()
#                self.ax.set_position([box.x0, box.y0,box.width*0.83,
#                                      box.height])
#                self.ax2 = self.fig.add_axes([0.85, 0.1, 0.05, 0.8])
                self.ax.set_position([box.x0, box.y0+0.05,box.width*0.83,
                                      box.height*0.9])
             #   self.ax2 = self.fig.add_axes([0.85, 0.0949, 0.05, 0.82])
                norm = mpl.colors.Normalize(vmin=0, vmax=2000)
                cb1 = mpl.colorbar.ColorbarBase(self.ax2,
                                                cmap = color,
                                                norm = norm)
                cb1.set_label('elevation (m)', fontsize=20)
                for t in cb1.ax.get_yticklabels():
                    t.set_fontsize(15)

            plt.setp(scatter, edgecolors='None')

        if kwargs.get('boarder_line', True):
            # Shapefile with missing north of Germany
            #shapefile_path =os.path.join(self.basepath,
             #              '03_Import','04_borders','states','region_geog')
            # Shapefile with north of Germany
            shapefile_path =os.path.join(self.basepath,
                           '03_Import','07_GIS','GRASS_Niedsim','raw_data',
                           'vg2500_geo84','vg2500_bld')
            shapefile_path = os.path.join(self.basepath,'02_Import','countries')
            self.m.readshapefile(shapefile_path, 'Bayern')



        # draw default parallels and meridians if keywords are not given
        # if keyword is false, drawnothing
        # if keywords is dictionary use it to describe the meridians and parallels
        if kwargs.get('paral', True):
            # label on left and bottom of map.
            circles = np.arange(-90.,90,2.)
            self.m.drawparallels(circles,labels=[1,0,0,1], fontsize=20)
        elif kwargs.get('paral', True) == False:
            pass
        else:
            self.m.drawparallels(**kwargs['paral'])
        if kwargs.get('merid', True):
            # label on left and bottom of map.
            meridians = np.arange(0.,360.,2.)
            self.m.drawmeridians(meridians,labels=[1,0,0,1], fontsize=20)
        elif kwargs.get('merid', True) == False:
            pass
        else:
            self.m.drawmeridians(**kwargs['merid'])



    def plot_scatter(self, lons, lats, *args, **kwargs):
        r"""
        Creates a scatter plot to be used e.g. stations

        Parameters
        ----------
        lons : 1-D np.array of floats
            Longitudes in degrees
        lats : 1-D np.array of floats
            Latitudes in degrees

        *args, **kwargs [from plt.scatter()]

        Returns
        -------
        None

        """
        cb       = kwargs.pop('cb', False)
        cb_label = kwargs.pop('cb_label', False)

        # convert to map projection coords.
        x, y = self.m(lons,lats)

        scatter = self.m.scatter(x,y,*args, **kwargs)

        if cb == True:
            x0, y0, width, height = self.orig_ax_position

            if cb_label:
                self.ax.set_position([x0, y0, width*0.79, height*0.95])
                self.ax2 = self.fig.add_axes([0.80, y0, .05, height*0.95] )
            else:
                self.ax.set_position([x0, y0, width*0.9, height*0.95])
                self.ax2 = self.fig.add_axes([0.85, y0, 0.05, height*0.95])
            # Set the colormap and norm to correspond to the data for which
            # the colorbar will be used.
            vmin = kwargs.get('vmin', np.nanmin(kwargs.get('c')))
            vmax = kwargs.get('vmax', np.nanmax(kwargs.get('c')))
            norm = mpl.colors.Normalize(vmin = vmin,
                                        vmax = vmax)

            # ColorbarBase derives from ScalarMappable and puts a colorbar
            # in a specified axes, so it has everything needed for a
            # standalone colorbar.  There are many more kwargs, but the
            # following gives a basic continuous colorbar with ticks
            # and labels.
            cb1 = mpl.colorbar.ColorbarBase(
                                    self.ax2,
                                    cmap = kwargs.get('cmap',None),
                                    norm = norm)

            for t in cb1.ax.get_yticklabels():
                t.set_fontsize(20)
            if cb_label:
                cb1.set_label(cb_label, fontsize=20)
            #remove the edgecolor
            plt.setp(scatter, edgecolors='None')


    def plot_raster(self, data, lons, lats, *args, **kwargs):

        r"""
        Creates a scatter plot to be used e.g. stations

        Parameters
        ----------
        data : 2-D np.array of floats
            Raster data
        lons : 2-D np.array of floats
            Longitudes in degrees of same dimension as data
        lats : 2-D np.array of floats
            Latitudes in degrees of same dimension as data
        *args, **kwargs [from plt.pcolor()]

        Returns
        -------
        None

        """
        cb = kwargs.pop('cb', False)

        #masked nan data
        mdata=np.ma.masked_array(data, mask=np.isnan(data))
        plot = self.m.pcolormesh(lons, lats, mdata, latlon=True,*args, **kwargs)

        # ColorbarBase derives from ScalarMappable and puts a colorbar
        # in a specified axes, so it has everything needed for a
        # standalone colorbar.  There are many more kwargs, but the
        # following gives a basic continuous colorbar with ticks
        # and labels.
        if cb == True:
            box = self.ax.get_position()
            self.ax.set_position([box.x0, box.y0,box.width*0.83,
                                  box.height])
            self.ax2 = self.fig.add_axes([0.85, 0.1, 0.05, 0.8] )
            norm = mpl.colors.Normalize(vmin=0, vmax=2000)
            cb1 = mpl.colorbar.ColorbarBase(self.ax2,
                                            cmap = kwargs.get('cmap', None),
                                            norm = norm)
            for t in cb1.ax.get_yticklabels():
                t.set_fontsize(20)


    def plot_text(self, lons, lats, text, *args, **kwargs):
        r"""
        Creates a text string at the given locatioin

        Parameters
        ----------
        lons : float
            Longitude in degrees
        lats : float
            Latitude in degrees
        text : string
            Text to be plotted
        *args, **kwargs [from plt.scatter()]

        Returns
        -------
        None

        """
        x, y = self.m(lons,lats)
        if np.isscalar(text):
            plt.text(x, y, text, *args, **kwargs)
        else:
            for ii in range(text.shape[0]):
                plt.text(x[ii], y[ii], text[ii], *args, **kwargs)

    def plot_slp_contour(self, data, lons, lats, label_extr=True, dstep=1.):
        r"""
        Creates a contour plot based on regular gridded data
        Plot designed for sea level pressure anomalies
        Red isobars for positive, blue isobars for negative pressure anomalies

        Parameters
        ----------
        data : 2-D np.array of floats
            Raster data
        lons : 1-D np.array of floats
            Longitudes in degrees of same dimension as data
        lats : 1-D np.array of floats
            Latitudes in degrees of same dimension as data
        dstep: float, optional
            Distance between the contour lines, default: 1 unit (1hPa)
        label_extr: logical, optional
            If true: label high and low pressure centers


        Returns
        -------
        None

        """

        # transform to nx x ny regularly spaced 5km native projection grid
        projection_grid=5000.
        nx = int((self.m.xmax-self.m.xmin)/projection_grid)+1
        ny = int((self.m.ymax-self.m.ymin)/projection_grid)+1
        topodat = self.m.transform_scalar(data,lons,lats,nx,ny,order=0)

        #calculate contour levels
        minval=int(data.min()/dstep)*dstep
        maxval=int((data.max()+dstep)/dstep)*dstep
        clevels=np.arange(minval,maxval,dstep)
        clevelh=clevels[clevels>=0.]
        clevell=clevels[clevels<0]

        # plot contour lines
        lons, lats = np.meshgrid(lons, lats)
        x, y = self.m(lons, lats)
        if len(clevelh>0):
            cs = self.m.contour(x,y,data,clevelh,colors='r',linewidths=1.)
        if len(clevell>0):
            cs = self.m.contour(x,y,data,clevell,colors='b',linewidths=1.)

        #label extreme values
        if label_extr:
            #limitl=data.mean()+-2*dstep
            #limith=data.mean()+2*dstep
            limitl=-dstep
            limith=+dstep

            local_min,local_max=extrema(data)
            lowvals=data[local_min]
            highvals=data[local_max]

            xlows = x[local_min]; xhighs = x[local_max]
            ylows = y[local_min]; yhighs = y[local_max]

            if lowvals.min()<limitl:
                for x,y,p in zip(xlows,ylows,lowvals):
                    if p <limitl:
                        plt.text(x,y,'L',fontsize=24,fontweight='bold',
                        ha='center',va='center',color='b')

            if highvals.max()>limith:
                for x,y,p in zip(xhighs,yhighs,highvals):
                    if p >limith:
                        plt.text(x,y,'H',fontsize=24,fontweight='bold',
                        ha='center',va='center',color='r')

    def extrema(mat,mode='wrap',window=10):
        """find the indices of local extrema (min and max)
        in the input array."""
        mn = minimum_filter(mat, size=window, mode=mode)
        mx = maximum_filter(mat, size=window, mode=mode)
        # (mat == mx) true if pixel is equal to the local max
        # (mat == mn) true if pixel is equal to the local in
        # Return the indices of the maxima, minima
        return np.nonzero(mat == mn), np.nonzero(mat == mx)


    def get_borders(self, border):
        """ Border areas in grid projection specifications"""

        projection = 'tmerc'
        lon_0=10
        lat_0=47

        if border=='bw':
            lats_b = [47,50]
            lons_b = [7,11]
            gk = 'gk3'
        elif border=='by':
            lats_b = [47.2,50.6]
            lons_b = [9,14]
            gk = 'gk4'
        elif border=='he':
            lats_b = [49,52]
            lons_b = [7,11]
            gk = 'gk3'
        elif border=='rp':
            lats_b = [48.5,51]
            lons_b = [6,8.5]
            gk = 'gk2'
        elif border=='ni':
            lats_b = [51.2,56]
            lons_b = [6.5,11.7]
            gk = 'gk3'
        elif border=='all':
            lats_b = [47,55.6]
            lons_b = [6, 16]
            gk = 'gk4'
        elif border=='1km':
            lats_b = [47,55.6]
            lons_b = [6, 16]
            gk = 'gk4'
        elif border=='south_germany':
            lats_b = [47,50]
            lons_b = [7.5,14]
            gk = 'gk4'
        elif border=='euro_atl':
            lats_b = [23,70]
            lons_b = [-40,60]
            projection = 'lcc'
            gk = 'gk3'
            lon_0=0
            lat_0=25

        return lons_b, lats_b,lon_0,lat_0,projection,gk

def plot_stats(HDF5, ids, gk, Map = False, id_labels = False,
               name_labels = False, **kwargs):
    """get and plot cordinates

    Parameters
    ----------
    HDF5: object
        HDF5 object created by b_get_data
    ids: 'string' or [list of strings]
        ids of the stations
    gk: string
        gk system of the input data i.e. 'gk2','gk3' or 'gk4'
    Map: object, optional
        Basemap object, if not given, it is created by
        c_plot_maps.Plot_Map()
    color: str, optional
        color of the marker, default = 'r'
    id_labels: bool, optional
        If True, plot the station ID, default = False
    name_labels: bool, optional
        If True, plot the station Name, default = False
    """
    border = kwargs.pop('border', 'all')
    backg  = kwargs.pop('backg', True)
    if Map == False:
        Map = Plot_Map(border = border, backg = backg, **kwargs)
    coordinates = HDF5.get_coordinates(ids)
    metadata = HDF5.get_metadata(ids)
    x = 'x{:s}'.format(gk[-1])
    y = 'y{:s}'.format(gk[-1])
    lons, lats = a_functions.gk_to_lonlat(coordinates[x],coordinates[y], gk)
    Map.plot_scatter(lons, lats, **kwargs)
    if id_labels == True:
        for ii, istat in enumerate(metadata['id']):
            Map.plot_text(lons[ii], lats[ii], istat, ha='right',va='top')
    if name_labels == True:
        for ii, istat in enumerate(metadata['name']):
            Map.plot_text(lons[ii], lats[ii], istat, ha='right',va='top')
    return Map


if __name__ == '__main__':
    pass
##    lons_r=np.arange(100)+6
##    lats_r=np.arange(100)+47
##    x,y = np.meshgrid(lons_r,lats_r)
##    data_r=np.arange(10000).reshape(100,100)
    ##data_r[-1,0]=100

##    lons_s=np.arange(4)+4
##    lats_s=np.arange(4)+47




    Map = Plot_Map(border = 'by')
##    Map.plot_raster(x,y,data_r,cmap = mpl.cm.Blues)



##    3494110	5631120 Kirchhain-Klaeranlage

##
##    lon,lat = gk_to_lonlat(xcoord_borders,
##                           ycoord_borders,
##                           'gk3')
####    Map.plot_scatter(lon,lat, color = 'r')
##    Map.plot_scatter(lon, lat, marker = 'x', s=200, color='r', linewidths=3)
####    Map.plot_text(lon, lat, np.array(['Stuttgart',
####                                      'Hohenstaufen']),
####                                     size=20)
##    Map()
    plt.show()

