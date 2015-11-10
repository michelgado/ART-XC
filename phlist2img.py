# -*- coding: utf-8 -*-
"""
Make a psf map from a photon list made by ray-tracing procedure.

v0.02, Hart
Should to be rewritten in vector form.
"""
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt

#...oooOOO000OOOooo...,,,...oooOOO000OOOooo...
def det2pix(ph_coord):
    "Map DET coordinates to pixels without any rotation"
    pix_out = img_size/2. + np.round(ph_coord/pix_size)
    return pix_out

def det2pix_rot(ph_x, ph_y, phi):
    u = ph_x*np.cos(phi)-ph_y*np.sin(phi)
    v = ph_x*np.sin(phi)+ph_y*np.cos(phi)     
    x_out = img_size/2. + np.round(u/pix_size)
    y_out = img_size/2. + np.round(v/pix_size)
    return x_out, y_out
#...oooOOO000OOOooo...,,,...oooOOO000OOOooo...

#Read input fits
infile_path = 'crab_photons_1_out.fits'
input_data = pf.open(infile_path)
photon_list = input_data[1].data
input_data.close()

#Initialize image
pix_size = 0.027
img_size = 700
img = np.zeros((img_size,img_size))

# Select all good photons
xs, ys, phis = [],[],[]
for ph in photon_list[:]:
    ph_stat, ph_xdet, ph_ydet = float(ph[6]),\
        float(ph[-3]),float(ph[-2])
    phi = float(ph[2])
    #Select all good photons
    if ph_stat==0:
        xs.append(ph_xdet)
        ys.append(ph_ydet)
        phis.append(phi)

# Calculate offsets to make image centered
x_offset = np.mean(xs)
y_offset = np.mean(ys)
# Map photons to image, with phi rotation and two-axis offset
for (ph_xdet, ph_ydet, phi) in zip(xs, ys, phis):
    pix_x, pix_y = det2pix_rot(ph_xdet-x_offset, ph_ydet-y_offset, phi)
    if ( pix_x<(img_size) and pix_x>=0 and\
         pix_y<(img_size) and pix_y>=0):
         img[pix_x, pix_y]+=1             

psfheader = pf.Header()
psfheader.append(('CRVAL1', '0.00000', 'Value at ref. pixel on axis 1'))
psfheader.append(('CRVAL2', '0.00000', 'Value at ref. pixel on axis 2'))
psfheader.append(('CRPIX1', str(img_size/2.), 'Reference pixel on axis 1'))
psfheader.append(('CRPIX2', str(img_size/2.), 'Reference pixel on axis 2'))
psfheader.append(('CRPIX1', str(pix_size), 'Coordinate increment on axis 1'))
psfheader.append(('CRPIX2', str(pix_size), 'Coordinate increment on axis 2'))
psfheader.append(('CTYPE1', 'DETX', 'name of data axis 1'))
psfheader.append(('CTYPE2', 'DETY', 'name of data axis 2'))
psfheader.append(('CUNIT1', 'mm', 'DETX coordinate units'))
psfheader.append(('CUNIT2', 'mm', 'DETX coordinate units'))
psfheader.append(('WCSNAME', 'DET', 'Default WCS'))
psffile = pf.ImageHDU(header=psfheader, data = img/np.sum(np.ravel(img)))
psffile.writeto('psf.fits', clobber=True)