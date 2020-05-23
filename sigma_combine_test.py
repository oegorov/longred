from auxiliary import combine_sigma
import numpy as np
from astropy.io import fits
import os
import glob

w_dir = "MyDir"

file_root = "prefix"
file_suffix = "suffix.fits"

file_out="final.fits"

treshold = 5


os.chdir(w_dir)

texp=0
file_list = glob.glob("{}*{}".format(file_root,file_suffix))
for i, f in enumerate(file_list):
    with fits.open(f) as hdu:
        if i == 0:
            images=np.ndarray((len(file_list),hdu[0].data.shape[0],hdu[0].data.shape[1]),dtype=np.float32)
            header_ref=hdu[0].header
        texp += hdu[0].header["EXPTIME"]
        images[i,:,:] = hdu[0].data
header_ref["EXPTIME"] = texp
header_ref.add_history("LONGRED: {} files were combined".format(len(file_list)))
header_ref.add_history("LONGRED: Sigma-clipping used with TR={}".format(treshold))
total = combine_sigma(images,cr_treshold=treshold)/len(file_list)

fits.writeto(file_out,total, header_ref, output_verify='fix', overwrite=True)

