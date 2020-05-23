from auxiliary import do_distortion
from auxiliary import derive_vecshift
import os
import numpy as np
from astropy.io import fits,ascii
import matplotlib.pyplot as plt
import time
from scipy.signal import find_peaks
from reduction import gaussian


# xx=np.arange(40)
# yy=gaussian(xx,[1,20,5,0.1])
# yy=np.pad(yy,30,mode='wrap')[40:-20]
# yy=np.pad(yy,15,mode='wrap')
# plt.plot(yy)
# # plt.plot(xx,yy)
# plt.show()

# os.chdir("/Users/mors/Science/test3/")
#
#
# with fits.open("obj_s1_VPHG1200@540_1x4_00_lin.fits") as hdu:
#     vecref=np.median(hdu[0].data[350:400,2400:2750],axis=0)
# with fits.open("sunsky_s1_VPHG1200@540_1x4_lin.fits") as hdu:
#     vec=np.median(hdu[0].data[350:400,2400:2750],axis=0)
# dx=derive_vecshift(vec,vecref,plot=True)
# print(dx)





tab=ascii.read("/Users/mors/Science/PyWorkspace/longred/gratings/TDS/calibration_lines_ini.txt",header_start=1)

line_lam=np.array(tab['Lam'])
line_intens=np.array(tab["Intens"])
fwhm=5.
# lam=np.linspace(3550,5900,2400)
lam=np.linspace(4200,5500,2400)
intens_lim=15
spec=np.zeros_like(lam)
for ind in range(len(line_lam)):
    if line_intens[ind]>=intens_lim:
        spec+=gaussian(lam,[line_intens[ind],line_lam[ind],fwhm/2.35428,0])

with fits.open("/Users/mors/Science/TDS-red/Manga_0104/neon_s1_G_1x1_warp_test.fits") as hdu:
    crow=int(round(hdu[0].data.shape[0]/2))
    obs_spec=np.median(hdu[0].data[crow-100:crow+101,:]+10,axis=0)

plt.figure(figsize=(20,10))
plt.subplot(211)
plt.plot(lam,np.log10(spec+10))

for ind in range(len(line_lam)):
    if line_intens[ind] >= intens_lim:
        plt.annotate("{}".format(int(round(line_lam[ind]))), (line_lam[ind],np.log10(line_intens[ind]+10)),
                     (line_lam[ind]-10, np.log10(line_intens[ind] + 10)+0.4),
                     rotation=90,arrowprops={'arrowstyle':'->'})

plt.subplot(212)
xx=np.arange(obs_spec.shape[0])
plt.plot(xx,np.log10(obs_spec))

peaks,_=find_peaks(np.log10(obs_spec),prominence=(.01,900),distance=10)
print(peaks)
for p in peaks:
    plt.annotate("{}".format(p), (p,np.log10(obs_spec[p])),
                     (p-10,np.log10(obs_spec[p])+0.1),
                     rotation=90,arrowprops={'arrowstyle':'->'})

plt.tight_layout()
plt.savefig("/Users/mors/Science/test1/TDS_green.pdf",dpi=150)

plt.show()

