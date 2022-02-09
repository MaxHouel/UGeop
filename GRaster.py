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


class Raster():
    def __init__(self) :
        None
    def getGeoInfo(self,path) :
        """
        
        Parameters
        ----------
        path : str
            DESCRIPTION.

        Returns
        -------
        dict
            Dictionnary with raster info : projection, Full Geometry, 
            number of Bands, x and y Resolution, Width and Height of the array.

        """
        data = gdal.Open(path)
        prj = data.GetProjection()
        geo = data.GetGeoTransform()
        nBands = data.RasterCount    
        xres = geo[1]
        yres = geo[5]
        shape = data.ReadAsArray().shape
        upx, xres, xskew, upy, yskew, yres = data.GetGeoTransform()
        cols = data.RasterXSize
        rows = data.RasterYSize
        ulx = upx + 0*xres + 0*xskew
        uly = upy + 0*yskew + 0*yres
        llx = upx + 0*xres + rows*xskew
        lly = upy + 0*yskew + rows*yres
        lrx = upx + cols*xres + rows*xskew
        lry = upy + cols*yskew + rows*yres
        urx = upx + cols*xres + 0*xskew
        ury = upy + cols*yskew + 0*yres
        STRWKT= 'POLYGON (({} {},{} {},{} {},{} {},{} {}))'.format(llx, lly, lrx, lry, urx, ury, ulx, uly, llx, lly)
        bbox = [ulx, uly, lrx, lry]
        geom = osgeo.ogr.CreateGeometryFromWkt(STRWKT)
        wkt = geom.ExportToWkt()
        json = geom.ExportToJson()
        return {'Projection' : prj,
                'Geometry' : geo,
                'Bbox' : bbox,
                'nBands' : nBands,
                'xRes' : xres,
                'yRes' : yres,
                'Width' : shape[-1:][0],
                'Height' : shape[-2:][0],
                'WKT' : wkt,
                'JSON' : json}      
    
    def asArray(self,path, band = None):
            """
            Parameters
            ----------
            path : str
                path to the raster file to open
            band : int, optional
                Integer to open a specific band. The default is None.
    
            Returns
            -------
            TYPE array
                Return the array of the raster file opened.
    
            """
            data = gdal.Open(path,gdal.GA_ReadOnly)
            if data is not None :
                num_bands = data.RasterCount
                if band is None :
                        array = data.ReadAsArray()
                else :
                    b = data.GetRasterBand(band)
                    array = b.ReadAsArray()
                return array
            
            else :
                return print('No data in {}'.format(path))       

    def buildVRT(self,list_files, dst, noData = None) :
        """

        Parameters
        ----------
        list_files : list of str
            List of file to compile.
        dst : str
            Path to the .vrt file to save.

        Returns
        -------
        Saving at dst.

        """
        gdal.BuildVRT(destName = dst,
                      srcDSOrSrcDSTab = list_files,
                      resolution = 'average',
                      allowProjectionDifference = True,
                      srcNodata = noData
                      )
        
    def vrtToTiff(self,vrtFile, dst, rsmple_alg = None) :
        """

        Parameters
        ----------
        vrtFile : str
            path to a .vrt file.
        dst : str
            Path to the .tif file to save.

        Returns
        -------
        Saving at dst.

        """
        gdal.Warp(destNameOrDestDS  = dst,
                       srcDSOrSrcDSTab  = vrtFile,
                       format = 'GTiff',
                       resampleAlg=rsmple_alg)
        
    def merge(self,list_files, dst, rsmple_alg = None, srcNodata = None) :
        """

        Parameters
        ----------
        list_files : list of str
            List of files to combine in one.
        dst : str
            Path to the destionation tif.

        Returns
        -------
        Saving .tif file combining files in the list
        Delete temporary .vrt file.
        """
        tmp = 'tmp.vrt'
        Raster.buildVRT(self, list_files, tmp, noData = srcNodata)
        Raster.vrtToTiff(self,tmp, dst, rsmple_alg = rsmple_alg )
        os.remove(tmp)
        
    def tiling(self,src, dst, size_w, step = None, dtype = 'Float32', noData = None) :
        """
        
        Parameters
        ----------
        src : str
            Path to the file to tile.
        dst : str
            Path to the tiles. (Automatically the script add
                                additional number of tile
            ex : dst = tile -> output = tile_{n}.tif
        size_w : int
            width size of the tile.
        size_h : int
            height size of the tile.
        dtype : str, optional
            Radiometric resolution : Float32 / Float64 / UInt16 / Byte
                . The default is 'Float32'.
        noData : int, optional
            Set a specific value for the nan data
            . The default is None.

        Returns
        -------
        Saved tiles at dst.

        """
        tile = 0
        if step == None :
            step = size_w
        im = Raster.asArray(self, src)
        if dtype == 'Float32':
            outputType = gdal.GDT_Float32
        if dtype == 'Float64':
            outputType = gdal.GDT_Float64
        if dtype == 'UInt16':
            outputType = gdal.GDT_UInt16 
        if dtype == 'Byte':
            outputType = gdal.GDT_Byte
        if len(im.shape) > 2 :
            bands, h, w =  im.shape
        else :
            h, w = im.shape
        for i in range(0,w, step):
            if (i+size_w) > w :
                i = w - size_w
            for j in range(0, h, step):
                if (j+size_w) > h :
                    j = h - size_w
                tile += 1
                destName = dst + '_{}.tif'.format(tile)
                if not os.path.isfile(destName):
                    gdal.Translate(destName = destName,
                                   srcDS = src,
                                   format = 'GTiff',
                                   outputType = outputType,
                                   srcWin = [i, j, size_w, size_w],
                                   noData = noData)

    
    def resizing(self,src, dst, width = None, height = None, xRes = None, yRes = None, dtype = 'Float32', resampleAlg = 0) :
        """
        Resizing process can be choose with shape size or pixel resolution
        
        Parameters
        ----------
        src : str
            Path to the file to resize.
        dst : str
            Output path.
        width : int, optional
            Width size for resizing
            . The default is None.
        height : int, optional
            Height size for resizing
            . The default is None.
        xRes : int, optional
            x resolution for resizing
            . The default is None.
        yRes : int, optional
            y resolution for resizing
            . The default is None.
        dtype : str, optional
            Radiometric resolution : Float32 / Float64 / UInt16 / Byte
            . The default is 'Float32'.
        resampleAlg : int, optional
            Resampling algorithm : 0 near / 1 bilinear / 2 cubic / 3 cubicspline
            4 Lanczos / 5 Average / 6 mode
            . The default is 0 for near
        Returns
        -------
        Resized file at dst.

        """
        if dtype == 'Float32':
            outputType = gdal.GDT_Float32
        if dtype == 'Float64':
            outputType = gdal.GDT_Float64
        if dtype == 'UInt16':
            outputType = gdal.GDT_UInt16 
        if dtype == 'Byte':
            outputType = gdal.GDT_Byte
        if xRes is None :
            if yRes is None :
                gdal.Warp(dst,
                          src,
                          width = width,
                          height = height,
                          outputType = outputType,
                          resampleAlg=resampleAlg)
        if width is None :
            if height is None :
                gdal.Warp(dst,
                          src,
                          xRes = xRes,
                          yRes = yRes,
                          outputType = outputType,
                          resampleAlg=resampleAlg)
    
    
    def crop(self,src, dst, vec, dtype = 'Float32') :
        """

        Parameters
        ----------
        src : str
            Path to file to crop.
        dst : str
            Output path.
        vec : str
            Path to the vector file for cropping.
        dtype : str, optional
            Radiometric resolution : Float32 / Float64 / UInt16 / Byte
            . The default is 'Float32'.

        Returns
        -------
        Crop file at dst

        """
        if dtype == 'Float32':
            outputType = gdal.GDT_Float32
        if dtype == 'Float64':
            outputType = gdal.GDT_Float64
        if dtype == 'UInt16':
            outputType = gdal.GDT_UInt16 
        if dtype == 'Byte':
            outputType = gdal.GDT_Byte
        gdal.Warp(dst,
                  src,
                  format = 'GTiff',
                  cutlineDSName = vec,
                  outputType = outputType,
                  cropToCutline = True)
        
    def save_tiff(self,dst_filename, nparray, prj, geom,
                  dtype = 'Float32', nodata_value =-9999.9): 
        """

        Parameters
        ----------
        dst_filename : str
            Path for the saving output.
        nparray : narray
            Array to save.
        ref_proj : str
            Path to the referenced file for projection.
        ref_geom : str
            Path to the referenced file for geometry.
        dtype : str, optional
            Radiometric resolution : Float32 / Float64 / UInt16 / Byte
            . The default is 'Float32'.
        nodata_value : int, optional
            Set a specific value for the nan data
            . The default is -9999.

        Returns
        -------
        Saved array at dst_filename.

        """
        if prj.endswith('.tif') :
            proj_ref = gdal.Open(prj, gdal.GA_ReadOnly)
            prj = proj_ref.GetProjectionRef()
        if type(geom)!=tuple:
            geom_ref = gdal.Open(geom, gdal.GA_ReadOnly)
            geom = geom_ref.GetGeoTransform()
        if len(nparray.shape) == 2 :
            [cols, rows] = nparray.shape
            bands = 1
        else :
            [bands, cols, rows] = nparray.shape
        if dtype == 'Float32':
            outputType = gdal.GDT_Float32
        if dtype == 'Float64':
            outputType = gdal.GDT_Float64
        if dtype == 'UInt16':
            outputType = gdal.GDT_UInt16 
        if dtype == 'Byte':
            outputType = gdal.GDT_Byte
        outdata = gdal.GetDriverByName("GTiff").Create(str(dst_filename),rows, cols, bands, outputType)
        outdata.SetGeoTransform(geom)
        outdata.SetProjection(prj)
        if bands > 1 :
            for band in range(1, bands + 1) :
                outdata.GetRasterBand(band).SetNoDataValue(nodata_value)
                outdata.GetRasterBand(band).WriteArray(nparray[:,:,band-1])
        else :
                outdata.GetRasterBand(bands).SetNoDataValue(nodata_value)
                outdata.GetRasterBand(bands).WriteArray(nparray)
        outdata.FlushCache()     

    def reproject(self, ref, img, dst) :
        ds = gdal.Open(img)
        geoTransform = ds.GetGeoTransform()
        old_prj = osgeo.osr.SpatialReference()
        old_prj.ImportFromWkt(ds.GetProjection())
        refds = gdal.Open(ref)
        new_prj = osgeo.osr.SpatialReference()
        new_prj.ImportFromWkt(refds.GetProjection())
        gdal.Warp(dst,
                  img,
                  srcSRS = old_prj,
                  dstSRS = new_prj)
        
    def clip(self, img, dst, bbox, prj, dtype = 'Float32', noData = None):
        if dtype == 'Float32':
            outputType = gdal.GDT_Float32
        if dtype == 'Float64':
            outputType = gdal.GDT_Float64
        if dtype == 'UInt16':
            outputType = gdal.GDT_UInt16 
        if dtype == 'Byte':
            outputType = gdal.GDT_Byte
        gdal.Translate(destName = dst,
                       srcDS = img,
                       format = 'GTiff',
                       outputType = outputType,
                       projWin = bbox,
                       projWinSRS = prj,
                       noData = noData)
