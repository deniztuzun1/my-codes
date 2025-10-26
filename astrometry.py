#plate soln uncertainty from astrometry.net's corr.fits file
from astropy.io import fits
import numpy 
from math import*

tbl = fits.open("corr.fits")[1].data
                
dtype=(numpy.record, [('field_x', '>f8'), ('field_y','>f8'), ('field_ra', '>f8'), ('field_dec', '>f8'),
('index_x', '>f8'), ('index_y', '>f8'), ('index_ra','>f8'), ('index_dec', '>f8'), ('index_id', '>i4'),
('field_id', '>i4'), ('match_weight', '>f8'), ('FLUX','>f4'), ('BACKGROUND', '>f4')])

i=0
diff = 0.

while (i < int(tbl.shape[0])):
     diffsquare = (tbl.field_ra[i] - tbl.index_ra[i])**2
     i += 1
rms_ra = sqrt(diffsquare/i)*3600
print(f'uncertainty ra= {rms_ra}')

i=0
diff = 0

while (i < int(tbl.shape[0])):
     diffsquare = (tbl.field_dec[i] - tbl.index_dec[i])**2
     i += 1
rms_dec = sqrt(diffsquare/i)*3600
print(f'uncertainty dec= {rms_dec}')