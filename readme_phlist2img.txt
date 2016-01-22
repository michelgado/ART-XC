phlist2img.py
Utility to convert photon list files produced by ray-tracing into fits images.
Current PSF image agreement:
700x700 pixels, 0.027 mm pixel size, center of image (350;350) correspond to center-of-mass.

Rotation is applied, value is taken from 'PHI' column of incoming file for each photon, therefore lists with photons
from several sources results may be incorrect.

v0.03alpha 22 Jan 2016 Changed to angular coordinates. Version is untested.
v0.02      10 Nov 2015 Fixed FITS headers.
v0.01       9 Nov 2015
