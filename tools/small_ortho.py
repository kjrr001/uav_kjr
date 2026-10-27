#!/usr/bin/env python3
"""Small ortho JPG at the size of another JPG (the small colour hillshade), white outside the data.
Usage: small_ortho.py <ortho.tif (RGBA)> <size from .jpg> <out.jpg>"""
import sys

import numpy as np
from osgeo import gdal
from PIL import Image

gdal.UseExceptions()
src, like, dst = sys.argv[1:4]
w, h = Image.open(like).size
a = gdal.Translate("", src, format="MEM", width=w, height=h, resampleAlg="average").ReadAsArray()
rgb = np.where(a[3:4] > 127, a[:3], 255).transpose(1, 2, 0).astype(np.uint8)
Image.fromarray(rgb).save(dst, quality=90)
