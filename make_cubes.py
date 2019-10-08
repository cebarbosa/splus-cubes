# -*- coding: utf-8 -*-
""" 

Created on 08/10/19

Author : Carlos Eduardo Barbosa

"""
from __future__ import print_function, division

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D

import os

def get_zps_dr1(tile, field="ZP"):
    """ Read the table containing the zero points for a given tile and given
    bands. """
    zpfile = os.path.join(os.getcwd(), "tables/ZPfiles_Feb2019",
                          "{}_ZP.cat".format(tile))
    zpdata = Table.read(zpfile, format="ascii")
    zps = dict([(t["FILTER"], t[field]) for t in zpdata])
    return zps

def make_stamps_stripe82(catalogs=None, redo=False, magmin=15, magmax=18):
    bands = ['U', 'F378', 'F395', 'F410', 'F430', 'G', 'F515', 'R', 'F660',
             'I', 'F861', 'Z']
    header_keys = ["OBJECT", "FILTER", "EXPTIME", "GAIN", "TELESCOP",
                   "INSTRUME"]
    ps = 0.55 # pixel scale [arcsec / pixel]
    catalogs = os.listdir(cat_dir) if catalogs is None else catalogs
    for catalog in catalogs:
        tile = catalog.split("_")[1]
        zps = get_zps_dr1(tile)
        cat = Table.read(os.path.join(cat_dir, catalog), format="ascii")
        ########################################################################
        # Filtering objects to be used from the catalog
        idx = np.where((cat["r_auto"] > magmin) & (cat["r_auto"] < magmax))[0]
        cat = cat[idx]
        ########################################################################
        coords = SkyCoord(cat["RA"].data * u.degree, cat["Dec"].data * u.degree)
        tmp_dir = os.path.join(outdir, tile)
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)
        # Producing all stamps of a given image
        for band in bands:
            imgfile = os.path.join(tiles_dir, tile,
                               "{}_{}_swp.fits".format(tile, band))
            data = fits.getdata(imgfile)
            # fnu: spectral flux density
            fnu = data * np.power(10, -0.4 * zps[band])
            # Surface brightness using fnu
            S = fnu / ps**2 # [erg / (s cm^2 Hz arcsec^2)]
            h = fits.getheader(imgfile)
            wcs = WCS(h)
            for i, source in enumerate(cat):
                size = 5 * source["A"] * source["KrRadDet"] * u.pixel
                cutout = Cutout2D(S, coords[i], size=size, wcs=wcs)
                out = os.path.join(tmp_dir, "{}.fits".format(
                                   source["ID"].replace("griz", band)))
                if os.path.exists(out) and not redo:
                    continue
                hdu = fits.ImageHDU(cutout.data)
                for key in header_keys:
                    hdu.header[key] = h[key]
                hdu.header["EXTNAME"] = band
                hdu.header.update(cutout.wcs.to_header())
                hdu.writeto(out, overwrite=True)
        # Joining stamps in a cube
        for source in cat:
            stamps = [os.path.join(tmp_dir, "{}.fits".format(
                      source["ID"].replace("griz", band))) for band in bands]
            if not all([os.path.exists(_) for _ in stamps]):
                continue
            cubename = os.path.join(outdir, source["ID"].replace("griz",
                                                                 "fits"))
            data = np.array([fits.getdata(stamp) for stamp in stamps])
            hdu = fits.PrimaryHDU(data)
            hdu.writeto(cubename, overwrite=True)
            # Deleting stamps
            for stamp in stamps:
                os.remove(stamp)
        os.rmdir(tmp_dir)

if __name__ == "__main__":
    data_dir = "/home/kadu/Dropbox/splus-ifusci/data/splus"
    cat_dir = os.path.join(data_dir, "catalogs_dr1")
    tiles_dir = os.path.join(data_dir, "tiles_dr1")
    outdir = os.path.join(data_dir, "cubes_dr1")
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    cat_dir = os.path.join(data_dir, "catalogs_dr1")
    catalogs = os.listdir(cat_dir)[:1]
    make_stamps_stripe82(catalogs, redo=False, magmin=16, magmax=16.03)