#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 16:19:44 2021

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
import tqdm
import geopandas as gpd

class Vector():
    def __init__(self) :
        None
    def getGeoInfo(self,path) :
        """
        
        Parameters
        ----------
        path : str
            path to the vector file (format .shp / .geojson ...).

        Returns
        -------
        infos : dict
            Dictionary with Extent / fields + values of each feature.

        """
        file = osgeo.ogr.Open(path, 1)
        layer = file.GetLayer()
        layerDfn = layer.GetLayerDefn()
        fields = []
        for i in range(layerDfn.GetFieldCount()):
            fields.append(layerDfn.GetFieldDefn(i).GetName())
        infos = {}
        f = 0 
        for feature in layer :
            geom = feature.GetGeometryRef()
            wkt = geom.ExportToWkt()
            geoJson = geom.ExportToJson()
            tmp={'WKT':wkt,
                 'JSON' : geoJson}
            for field in fields :
                tmp.update({field:feature.GetField(field)})
            infos.update({f:tmp})
            f+=1
        return infos
        
    def addField(self, shp, fld, val):
        data = osgeo.ogr.Open(shp, 1)
        layer = data.GetLayer()
        if type(val)==int :
            T = osgeo.ogr.OFTInteger
        if type(val)==float :
            T = osgeo.ogr.OFTReal
        if type(val)==str :
            T = osgeo.ogr.OFTString
        fldDfn = osgeo.ogr.FieldDefn(fld, T)
        layer.CreateField(fldDfn)
        for feature in layer :
            feature.SetField(fld, val)
            layer.SetFeature(feature)
            feature = None
        data = None
        
    def sepFeatures(self, shp, dst, fieldRef, prefix='Feature'):
        if not os.path.isdir(dst) :
            os.mkdir(dst)
        file = osgeo.ogr.Open(shp, 0)
        layer = file.GetLayer()
        layerDfn = layer.GetLayerDefn()
        for feature in layer :
            srs = osgeo.osr.SpatialReference()
            srs.ImportFromEPSG(4326)
            geom = feature.GetGeometryRef()
            geoType = geom.GetGeometryName()
            if geoType == 'MULTIPOLYGON' :
                wkb = osgeo.ogr.wkbMultiPolygon
            elif geoType == 'POINT' :
                wkb = osgeo.ogr.wkbPoint
            elif geoType == 'POLYGON' :
                wkb = osgeo.ogr.wkbPolygon
            elif geoType == 'LINESTRING':
                wkb = osgeo.ogr.wkbLineString
            name = feature.GetField(fieldRef)
            if type(name) == float :
                name = int(name)
            outJSON = dst + '/{}_{}.geojson'.format(prefix ,name)
            outDriver = osgeo.ogr.GetDriverByName("GeoJSON") 
            if os.path.exists(outJSON):
                outDriver.DeleteDataSource(outJSON)
            outDataSource = outDriver.CreateDataSource(outJSON)
            outLayer = outDataSource.CreateLayer("{}".format(name), srs,  geom_type=wkb)
            for i in range(0, layerDfn.GetFieldCount()):
                fieldDefn = layerDfn.GetFieldDefn(i)
                outLayer.CreateField(fieldDefn)            
            outDefn = outLayer.GetLayerDefn()
            outFeature = osgeo.ogr.Feature(outDefn)
            for i in range(0, outDefn.GetFieldCount()):
                outFeature.SetField(outDefn.GetFieldDefn(i).GetNameRef(), feature.GetField(i))
            outFeature.SetGeometry(geom)
            outLayer.CreateFeature(outFeature)
            outFeature = None
            feature = None
            outDataSource = None
            print(outJSON, "created")
        file = None

    def getLayerExtent(self, shp, output = None):
        dataSource = osgeo.ogr.Open(shp, 1)
        layer = dataSource.GetLayer()
        extent = layer.GetExtent()
        dataSource = None
        ring = osgeo.ogr.Geometry(osgeo.ogr.wkbLinearRing)
        ring.AddPoint(extent[0],extent[2])
        ring.AddPoint(extent[1], extent[2])
        ring.AddPoint(extent[1], extent[3])
        ring.AddPoint(extent[0], extent[3])
        ring.AddPoint(extent[0],extent[2])
        poly = osgeo.ogr.Geometry(osgeo.ogr.wkbPolygon)
        poly.AddGeometry(ring)
        if output == 1 :
            dst = shp.split(".")[0] + "_extent.geojson"    
            outDriver = osgeo.ogr.GetDriverByName("GeoJSON")
            if os.path.exists(dst):
                outDriver.DeleteDataSource(dst)
            outDataSource = outDriver.CreateDataSource(dst)
            outLayer = outDataSource.CreateLayer("extent", geom_type=osgeo.ogr.wkbPolygon)
            idField = osgeo.ogr.FieldDefn("id", osgeo.ogr.OFTInteger)
            outLayer.CreateField(idField)
            featureDefn = outLayer.GetLayerDefn()
            feature = osgeo.ogr.Feature(featureDefn)
            feature.SetGeometry(poly)
            feature.SetField("id", 1)
            outLayer.CreateFeature(feature)
            feature = None
            outDataSource = None
            dataSource = None
            return print(dst, "computed")
        else :
            jsonExtent = poly.ExportToJson()
            return jsonExtent

    def rasterize(self, shp, ref, nClass, dst = 'MEM.tif') :
        if type(shp) == str :
            dtype = shp.split('.')[-1:][0]
            if dtype == 'shp' :
                driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
            elif dtype =='geojson' :
                driver = osgeo.ogr.GetDriverByName('GeoJSON')
            DS = driver.Open(shp, 0)
            Layer = DS.GetLayer()
        else :
            Layer = shp
        dataset = gdal.Open(ref, gdal.GA_ReadOnly)
        geotransform = dataset.GetGeoTransform()
        projection = dataset.GetProjection()
        driver = gdal.GetDriverByName('GTiff') # In memory dataset
        target_ds = driver.Create(dst, dataset.RasterXSize, dataset.RasterYSize, 1, gdal.GDT_Byte)
        target_ds.SetGeoTransform( geotransform )
        target_ds.SetProjection( str(projection) )
        gdal.RasterizeLayer( target_ds, [1], Layer, burn_values=[0], options=["ATTRIBUTE={}".format(nClass)]) # CHANGE EN FONCTION DU SHAPE
        labeled_pixels = target_ds.GetRasterBand(1).ReadAsArray()
        return labeled_pixels
    
            
    def buffer(self, shp, dst, bufSize, style = 1) :
        dataSource = osgeo.ogr.Open(shp, 1)
        layer = dataSource.GetLayer()
        defUnit = layer[0].GetGeometryRef().GetSpatialReference().ExportToWkt().split('UNIT[')[-1]
        if defUnit.startswith('"degree"'):
            bufSize = bufSize / 100000
        data = gpd.read_file(shp)
        buffer = data.copy()
        buffer.geometry = buffer['geometry'].buffer(bufSize, cap_style = style)
        buffer.to_file(dst, driver='GeoJSON')

            
            
    def clip(self, input, output, overlay):
        driver = osgeo.ogr.GetDriverByName('GeoJSON')
        if type(input) is str :
            DS = osgeo.ogr.Open(input, 0)
            inlayer = DS.GetLayer()
        else : 
            inlayer = input.GetLayer()
        feature = inlayer[0]
        geom = feature.GetGeometryRef()
        geoType = geom.GetGeometryName()
        if geoType == 'MULTIPOLYGON' :
            wkb = osgeo.ogr.wkbMultiPolygon
        elif geoType == 'POINT' :
            wkb = osgeo.ogr.wkbPoint
        elif geoType == 'POLYGON' :
            wkb = osgeo.ogr.wkbPolygon
        elif geoType == 'LINESTRING':
            wkb = osgeo.ogr.wkbLineString
        overlayDS = osgeo.ogr.Open(overlay, 0)
        overlayer = overlayDS.GetLayer()
        if os.path.exists(output):
            driver.DeleteDataSource(output)
        outDataSource = driver.CreateDataSource(output)
        outLayer = outDataSource.CreateLayer('', geom_type=wkb)
        osgeo.ogr.Layer.Clip(inlayer, overlayer, outLayer)
        DS = None
        overlayDS = None
        outDataSource = None
