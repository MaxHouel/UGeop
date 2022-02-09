#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 16:19:43 2021

@author: houel
"""

import os
import time
import numpy as np
import osgeo
from osgeo import gdal
import osgeo.ogr
import osgeo.osr
import utm
import mgrs
import pyproj


class Geo(): 
    def __init__(self) :
        None
    def inside(self, inside, extent):
        """
        
        Parameters
        ----------
        inside : str
            Json format.
        extent : str
            Json format.

        Returns
        -------
        bool
            Return if 'inside' is fully inside 'extent'.

        """
        json1 = """{}""".format(inside)
        json2 = """{}""".format(extent)
        geom1 = osgeo.ogr.CreateGeometryFromJson(json1)
        geom2 = osgeo.ogr.CreateGeometryFromJson(json2)
        if geom1.Within(geom2):
            return True
        else :
            return False
       
    def intersect(self, inside, extent):
        """
        
        Parameters
        ----------
        inside : str
            Json format.
        extent : str
            Json format.

        Returns
        -------
        bool
            Return if 'inside' intersect 'extent'.

        """
        json1 = """{}""".format(inside)
        json2 = """{}""".format(extent)
        geom1 = osgeo.ogr.CreateGeometryFromJson(json1)
        geom2 = osgeo.ogr.CreateGeometryFromJson(json2)
        inter = geom1.Intersection(geom2)
        wkt = inter.ExportToWkt()
        if wkt!='POLYGON EMPTY' :
            return True
        else :
            return False

class CoordConv():
    def getDegree(self, m):
        """

        Parameters
        ----------
        m : float
            meter.

        Returns
        -------
        int
            Convert meter in degree.

        """
        return m/100000
    	
    def getMeter(self, d):
        """
        

        Parameters
        ----------
        m : float
            degree.

        Returns
        -------
        int
            Convert meter in degree.

        """
        return d*100000
    
    def getPixelCoord(self, y, x, img) :
        """
        

        Parameters
        ----------
        y : float
            y loc.
        x : float
            x loc.
        img : str
            path to img.

        Returns
        -------
        lon : float
            return longitude of a y location of a matrix.
        lat : float
            return latitude of a x location of a matrix.

        """
        ds = gdal.Open(img)
        prj = osgeo.osr.SpatialReference(wkt=ds.GetProjection())
        epsg = prj.GetAttrValue('AUTHORITY', 1)
        target = osgeo.osr.SpatialReference()
        target.ImportFromEPSG(int(epsg))
        transform = osgeo.osr.CoordinateTransformation(prj, target)
        geo_matrix = ds.GetGeoTransform()     
        ul_x = geo_matrix[0]
        ul_y = geo_matrix[3]
        x_dist = geo_matrix[1]
        y_dist = geo_matrix[5]
        _x = x * x_dist + ul_x
        _y = y * y_dist + ul_y           
        point= osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
        point.AddPoint(_x, _y)
        point.Transform(transform)
        lon, lat = point.GetX(), point.GetY()
        return (lon, lat)

    def LatLontoXY(self, lat, lon, prj):
        """
        

        Parameters
        ----------
        lat : float
            latitude.
        lon : float
            longitude.
        prj : str
            wkt projection.

        Returns
        -------
        x : float
            x loc of lat loc.
        y : float
            y loc of long loc.

        """
        converter = osgeo.osr.SpatialReference()
        converter.ImportFromWkt(prj)
        pyprj = converter.ExportToProj4()
        p = pyproj.Proj(pyprj, preserve_units=True)
        x, y = p(lat, lon)
        return x, y
    
    def XYtoLatLon(self, x, y, prj):
        """
        

        Parameters
        ----------
        y : float
            y loc.
        x : float
            x loc.
        img : str
            path to img.

        Returns
        -------
        lon : float
            return longitude of a y location of a matrix.
        lat : float
            return latitude of a x location of a matrix.

        """
        converter = osgeo.osr.SpatialReference()
        converter.ImportFromWkt(prj)
        pyprj = converter.ExportToProj4()
        p = pyproj.Proj(pyprj, preserve_units=True)
        lat, lon = p(x, y, inverse=True)
        return lon, lat            
    
    def getTile(self,lat, long) :
        """

        Parameters
        ----------
        lat : int
            Latitude to convert.
        long : int
            Longitude to convert.

        Returns
        -------
        tile : str
            Give corresponding MGRS tile from lat / long input.

        """
        m = mgrs.MGRS()
        c=m.toMGRS(lat,long)
        c = c.decode('utf-8')
        tile = 'T' + c[:5]
        return tile
    def wgsTOutm(self,lat,long):
        """

        Parameters
        ----------
        lat : int
            Latitude to convert.
        long : int
            Longitude to convert.

        Returns
        -------
        zone : str
            UTM conversion corresponding to lat / long input.

        """
        utm_coords = utm.from_latlon(lat,long)
        return utm_coords