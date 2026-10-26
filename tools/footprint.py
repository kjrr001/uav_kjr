#!/usr/bin/env python3
"""Outline of a DSM's data (largest connected area, holes filled) as a polygon.
footprint(dsm_tif, epsg=None) -> (n, 2) array in the DSM's CRS, or in EPSG:<epsg> if given."""
import json
import numpy as np
from osgeo import gdal, ogr, osr
from scipy import ndimage

gdal.UseExceptions()


def footprint(dsm, epsg=None, res=2.0, simplify=3.0):
    m = gdal.Translate("", dsm, format="MEM", xRes=res, yRes=res, resampleAlg="near")
    z = m.GetRasterBand(1).ReadAsArray()
    v = np.isfinite(z) & (z > -1000)
    v = ndimage.binary_fill_holes(ndimage.binary_opening(v, iterations=2))
    lab, n = ndimage.label(v)
    v = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
    mk = gdal.GetDriverByName("MEM").Create("", v.shape[1], v.shape[0], 1, gdal.GDT_Byte)
    mk.SetGeoTransform(m.GetGeoTransform()); mk.SetProjection(m.GetProjection())
    mk.GetRasterBand(1).WriteArray(v.astype(np.uint8))
    sr = osr.SpatialReference(wkt=m.GetProjection())
    ds = ogr.GetDriverByName("Memory").CreateDataSource("")
    ly = ds.CreateLayer("p", sr); ly.CreateField(ogr.FieldDefn("v", ogr.OFTInteger))
    gdal.Polygonize(mk.GetRasterBand(1), mk.GetRasterBand(1), ly, 0)
    g = max([f.GetGeometryRef().Clone() for f in ly], key=lambda g: g.GetArea()).Simplify(simplify)
    if epsg:
        t = osr.SpatialReference(); t.ImportFromEPSG(epsg)
        t.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        g.Transform(osr.CoordinateTransformation(sr, t))
    return np.array(json.loads(g.ExportToJson())["coordinates"][0])
