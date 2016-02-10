# -*- coding: utf-8 -*-
"""
Make a psf map from a photon list made by ray-tracing procedure.

v0.03, Hart
Should to be rewritten in vector form.
"""
import pyfits as pf
import numpy as np

#...oooOOO000OOOooo...,,,...oooOOO000OOOooo...
#Define some basic parameters

pixel_phys_size = 0.0135 #mm, size of CCD pixel
binning = 2 
defocus = 7 #mm, private communications from Alex T.
focal_lenght = 2765 #mm,from Roman K.
eff_focal_lenght = focal_lenght - defocus #to use defocus in future
pix_size = np.arctan(binning*pixel_phys_size/eff_focal_lenght) # radians
pix_size = np.rad2deg(pix_size) #now in degrees
img_size = 700

#...oooOOO000OOOooo...,,,...oooOOO000OOOooo...

#Read input fits
infile_path = '../line8keV_photons_out.fits'
input_data = pf.open(infile_path)
photon_list = input_data[1].data
input_data.close()


# Select all good photons
xs, ys, phis = [],[],[]
for ph in photon_list[:10000]:
    ph_stat, ph_xdet, ph_ydet = float(ph[6]),\
        float(ph[-2]) ,float(ph[-1])
    phi = float(ph[2])
    #Select all good photons
    if ph_stat==0:
        xs.append(ph_xdet)
        ys.append(ph_ydet)
        phis.append(phi)
#Print some info, just to entertain user
print 'There is '+str(len(phis))+' photons to add!'

# Calculate offsets to make image centered
x_offset = np.mean(xs)
y_offset = np.mean(ys)

phlist = np.array([zip(xs, ys, phis)])

def det2pix_rot(ph_x, ph_y, phi):
    "Map DET coordinates to pixels with rotation applied"
    u = ph_x*np.cos(phi)-ph_y*np.sin(phi) #apply 
    v = ph_x*np.sin(phi)+ph_y*np.cos(phi) #      rotation
    u_deg = np.rad2deg(np.arctan(u/eff_focal_lenght)) # transform to 
    v_deg = np.rad2deg(np.arctan(v/eff_focal_lenght)) #   angular values, deg
    x_out = img_size/2. + np.round(u_deg/pix_size)
    y_out = img_size/2. + np.round(v_deg/pix_size)
    return x_out, y_out


def worker(pix_list):
    #Initialize image
    img = np.zeros((img_size,img_size))
    # Map photons to image, with phi rotation and two-axis offset
    for (ph_xdet, ph_ydet, phi) in pix_list:
        pix_x, pix_y = det2pix_rot(ph_xdet-x_offset, ph_ydet-y_offset, phi)
        if ( pix_x<(img_size) and pix_x>=0 and\
            pix_y<(img_size) and pix_y>=0):
            img[pix_x, pix_y]+=1             
    return img

import multiprocessing
#Set appropriate number of cores to use
N_proc = 2
pool = multiprocessing.Pool(processes = N_proc)
data = pool.map(worker, phlist)

img = np.zeros((img_size,img_size))
for part in data:
    img = img+part

psfheader = pf.Header()
psfheader.append(('CRVAL1', '0.00000', 'Value at ref. pixel on axis 1'))
psfheader.append(('CRVAL2', '0.00000', 'Value at ref. pixel on axis 2'))
psfheader.append(('CRPIX1', str(img_size/2.), 'Reference pixel on axis 1'))
psfheader.append(('CRPIX2', str(img_size/2.), 'Reference pixel on axis 2'))
psfheader.append(('CRPIX1', str(pix_size), 'Coordinate increment on axis 1'))
psfheader.append(('CRPIX2', str(pix_size), 'Coordinate increment on axis 2'))
psfheader.append(('CTYPE1', 'DETX', 'name of data axis 1'))
psfheader.append(('CTYPE2', 'DETY', 'name of data axis 2'))
psfheader.append(('CUNIT1', 'deg', 'DETX coordinate units'))
psfheader.append(('CUNIT2', 'deg', 'DETX coordinate units'))
psfheader.append(('WCSNAME', 'DET_ANG', 'Default WCS'))
psffile = pf.ImageHDU(header=psfheader, data = img/np.sum(np.ravel(img)))
psffile.writeto('psfnew.fits', clobber=True)
