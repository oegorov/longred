"""
All the main steps of the longslit data reduction pipeline are here
==============================================================
"""

import os
from misc import message
from auxiliary import *
import numpy as np
from astropy.io import fits, ascii
#from ccdproc import cosmicray_lacosmic as lacosmic
from lacosmic import lacosmic
from matplotlib import pyplot as plt
from scipy.ndimage.filters import uniform_filter
from scipy import signal
from astropy.modeling import models, fitting
from astropy.time import Time
from matplotlib import gridspec
import scipy.interpolate as interpolate
import re


def gaussian(x, param):
    #param = [a, x0, sigma, c]
    a,x0,sigma,c=param
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + c



def window_stdev(arr, nfiles, radius):
    c1 = uniform_filter(arr, [nfiles,radius,radius], mode='reflect')
    c2 = uniform_filter(arr*arr, [nfiles,radius,radius], mode='reflect')
    return ((c2 - c1*c1)**.5)[int(nfiles/2),:,:]



def cre_bias_or_dark(filein=None, fileout=('meanbias_b.fits','meanbias_r.fits'), w_dir=os.getcwd(), raw_dir='',
                    overscan=None, mode="TDS", filebias = None, filenoise=None, memory_limit=800):
    """
        Creates meanbias or dark (if filebias is not None) frame in w_dir from available bias or dark frames (if any)

        :param filein: list of names of input raw bias frames [if mode="TDS" then should be in format (['blue_channel'], ['red_channel'])]
        :param fileout: name of output file [if mode="TDS" then should be in format ('blue_chan', 'red_chan')]
        :param filebias: name of bias files to be subtracted from initial dark frames [if mode="TDS" then should be in format ('blue_chan', 'red_chan')]
        :param w_dir: output directory
        :param raw_dir: directory with data files
        :param mode: "TDS" or "SAO" [if mode="TDS" then should be in format (['blue_channel'], ['red_channel'])
        :param memory_limit: maximum amount of memory (in Mb) to be used during file combination.
                If total size of all images extends it - the combination will proceed in chunks
        :return: True or False
    """

    if not filein:
        if not filebias:
            message("No any raw bias frames was found; Meanbias is not created.")
        else:
            message("No any raw dark frames was found; Median dark frame is not created.")
        return False


    if mode == "TDS":
        nchannels=2
        if filebias is None:
            red_type='bias'
            filebias=(None,None)
        else:
            red_type='dark'
    else:
        nchannels=1
        filein=(filein,)
        fileout=(fileout,)
        if filebias is None:
            red_type='bias'
        else:
            red_type='dark'
        filebias=(filebias,)
        if filenoise:
            filenoise = (filenoise,)


    for chan in range(nchannels):
        nimages = len(filein[chan])
        nmaxpix = int(np.round(memory_limit*1e6 / (nimages * 64)))

        with fits.open(os.path.join(raw_dir, filein[chan][0])) as hdu:
            head_out = hdu[0].header
            if 'FILENAME' in head_out.keys():
                del head_out['FILENAME']
            ny, nx = (head_out["NAXIS2"], head_out["NAXIS1"])
            image_out = np.ndarray(shape=(ny, nx), dtype=float)
            if filebias[chan] is None:
                if filenoise:
                    image_error_out = np.ndarray(shape=(ny, nx), dtype=float)
                else:
                    image_error_out =[]
            nchunk = (nx * ny // nmaxpix) + 1

        no_bias_exist = True
        if filebias[chan] is not None:
            cur_bias_file = os.path.join(w_dir, filebias[chan])
            if os.path.isfile(cur_bias_file):
                with fits.open(cur_bias_file) as hdu:
                    bias_for_subtract = hdu[0].data
                    no_bias_exist = False
            else:
                bias_for_subtract=np.zeros(shape=(ny,nx), dtype=float)
        else:
            bias_for_subtract=np.zeros(shape=(ny,nx), dtype=float)

        nx_chunk=int(np.ceil(nx / nchunk))

        if nchannels > 1:
            add_mess="Channel {}: ".format(chan+1)
        else:
            add_mess=""

        message("{}Found {} of {} frames. Processing using {} chunks".format(add_mess,nimages, red_type,nchunk))

        for xstart in range(0, nx, nx_chunk):
            if (nx-xstart) < nx_chunk:
                nx_chunk=(nx-xstart)
            image_frames = np.ndarray(shape=(nimages, ny, nx_chunk), dtype=float)
            for i, f in enumerate(filein[chan]):
                with fits.open(os.path.join(raw_dir, f)) as hdu:
                    if red_type=='dark':
                        norm = hdu[0].header.get("exptime")
                        if not norm:
                            norm=hdu[0].header.get("exposure")
                        norm=float(norm)
                    else:
                        norm=1.
                    image_frames[i, :, :] = (hdu[0].data[:,xstart:(xstart+nx_chunk)] - bias_for_subtract[:,xstart:(xstart+nx_chunk)])/norm
            if nimages >=5:
                image_out[:, xstart:(xstart + nx_chunk)]=combine_sigma(image_frames,cr_treshold=2)/nimages
            else:
                image_out[:,xstart:(xstart+nx_chunk)] = combine_pairs(image_frames,cr_treshold=15)#np.median(image_frames, axis=0)
            if filebias[chan] is None:
                if filenoise:
                    image_error_out[:, xstart:(xstart + nx_chunk)] = window_stdev(image_frames, nimages, 20)
                else:
                    image_error_out.append(np.median(window_stdev(image_frames, nimages, 20)))


        if red_type == 'bias':
            head_out.add_history("LONGRED:  meanbias frame obtained")
            head_out.add_history("LONGRED: {} bias frames were used".format(nimages))
            # image_out=check_overscan(image_out,mode=mode)
            if filenoise:
                fits.writeto(os.path.join(w_dir, filenoise[chan]), np.float32(image_error_out), head_out,
                             output_verify='fix', overwrite=True)
            else:
                if len(image_error_out) > 1:
                    head_out["RNOISE"] = np.median(image_error_out)
                else:
                    head_out["RNOISE"] = np.median(image_error_out[0])
        else:
            head_out["BUNIT"]="counts/s"
            head_out.add_history("LONGRED:  median dark frame obtained")
            head_out.add_history("LONGRED: {} dark frames were used".format(nimages))
            cur_mode=mode
            if mode == "TDS":
                if chan == 0: cur_mode='TDS_b'
                else: cur_mode='TDS_r'

            force_overscan=None
            if overscan is not None and any(overscan):
                force_overscan=[0,head_out["NAXIS1"],0,head_out["NAXIS2"]]
                for j,v in enumerate(overscan):
                    applied=0
                    if v is not None and v>=2:
                        applied+=1
                        if j==0 or j==2:
                            force_overscan[j]=v//2
                        else:
                            force_overscan[j] -= v//2
                    if applied == 0:
                        force_overscan=None
            image_out = check_overscan(image_out, mode=cur_mode, force_overscan=force_overscan)
            image_out.ravel()[image_out.ravel() < 0] = 0
            if not no_bias_exist:
                head_out.add_history("LONGRED: Bias {} was subtracted".format(filebias[chan]))
        fits.writeto(os.path.join(w_dir,fileout[chan]), np.float32(image_out), head_out, output_verify='fix', overwrite=True)
        if red_type == 'bias':
            message("Meanbias {} is created".format(fileout[chan]))
        else:
            if not no_bias_exist:
                add_mess = "Bias frame {} was subtracted".format(filebias[chan])
            else:
                add_mess = "Bias frame was not found and not subtracted"
            message("Median dark frame {} is created; {}".format(fileout[chan], add_mess))

    return True





def cre_calibration(filein=None, fileout='flat.fits', w_dir=os.getcwd(), filebias=None, filedark=None,  raw_dir='', mode="TDS",
             cr_treshold=3, cr_mode="sigma", file_error_out = None, gain = None, imgtype="flat", readnoise=None,
                    cut=None, overscan=None):
    """
    На вход подаем список исходных калибровочных файлов (flat, sunsky, neon или 13dots; лежат в директории raw_dir),
    на выходе - собранную из них соответствующую калибровку.
    Вычитается bias, dark (если они заданы), производится сложение с использованием sigma-clipping (только если число кадров >=5)
    или pair - попарного сравнения
    Если задан fileerror_out, то считается файл с ошибками
    Если imgtype = 'neon' или 'sunsky' - файлы совмещаются по оси X.
    Если imgtype = '13dots' - файлы совмещаются по оси Y.
    Если imgtype = 'neon' или '13dots' - изображения нормируются на экспозицию перед сложением
    Если imgtype = 'flat' или 'sunsky' - используется нормировка на медианное значение интенсивности в каждом файле
    """
    if not filein:
        message("Nothing selected for processing => {} will not be created.".format(imgtype))
        return False

    do_xshift = False
    do_yshift = False
    if imgtype in ['neon', 'sunsky']:
        do_xshift = True
    if imgtype in ['13dots']:
        do_yshift = True
    normtype="median"
    if imgtype in ['neon', '13dots']:
        normtype = "exposure"

    nfiles = len(filein)

    message("######==== Start processing of {} data.....".format(imgtype))
    norm=np.ones(nfiles,dtype=float)
    for i,f in enumerate(filein):
        with fits.open(os.path.join(raw_dir,f)) as hdu:
            texp=hdu[0].header.get("EXPTIME")
            if not texp:
                texp = hdu[0].header.get("EXPOSURE")
            if i == 0:
                if gain is None:
                    gain = hdu[0].header.get("gain")
                head_out = hdu[0].header
                head_out['BZERO'] = 0
                if 'FILENAME' in head_out.keys():
                    del head_out['FILENAME']
                image_frames = np.ndarray(shape=(nfiles, head_out["NAXIS2"], head_out["NAXIS1"]), dtype=float)
                if filebias:
                    with fits.open(os.path.join(w_dir,filebias)) as hdu_bias:
                        bias = hdu_bias[0].data
                        if readnoise is None:
                            bias_squared_error = hdu_bias[0].header.get("RNOISE")**2
                        else:
                            bias_squared_error = (readnoise/gain) ** 2
                    head_out.add_history("LONGRED: bias subtracted: {}".format(filebias))
                else:
                    if readnoise is not None:
                        bias_squared_error = (readnoise/gain)**2
                    else: bias_squared_error=0
                    bias = 0
                if filedark:
                    with fits.open(os.path.join(w_dir,filedark)) as hdu_dark:
                        dark = hdu_dark[0].data
                    head_out.add_history("LONGRED: dark frame subtracted: {}".format(filedark))
                else:
                    dark = np.zeros((head_out["NAXIS2"], head_out["NAXIS1"]),dtype=float)
                if file_error_out:
                    error_squared = np.ndarray(shape=(nfiles,head_out["NAXIS2"], head_out["NAXIS1"]), dtype=float)
            force_overscan = None
            if cut is not None and any(cut) and overscan:
                force_overscan = [0, head_out["NAXIS1"], 0, head_out["NAXIS2"]]
                applied = 0
                for j, v in enumerate(cut):
                    if v is not None and v >= 2:
                        applied += 1
                        if j == 0 or j == 2:
                            force_overscan[j] = v // 2
                        else:
                            force_overscan[j] -= v // 2
                if applied == 0:
                    force_overscan = None
            if cut is not None and any(cut):
                force_cut = [0, head_out["NAXIS1"], 0, head_out["NAXIS2"]]
                applied = 0
                for j, v in enumerate(cut):
                    if v is not None and v > 0:
                        applied += 1
                        if j == 0 or j == 2:
                            force_cut[j] = v
                        else:
                            force_cut[j] -= v
                if applied == 0:
                    force_cut = None
            else:
                force_cut = None

            image_frames[i, :, :] = check_overscan(hdu[0].data-bias, mode=mode, force_overscan=force_overscan)-dark*texp
            if normtype == "exposure":
                image_frames[i, :, :] = image_frames[i, :, :]/texp
                norm[i]=texp
            if file_error_out:
                error_squared[i,:,:] = (np.abs(dark)*texp)/gain+bias_squared_error
                # if normtype == "exposure":
                #     error_squared[i, :, :] = error_squared[i, :, :] / texp**2

    if normtype =='median':
        y_medlimits=(int(np.round(0.3*head_out["NAXIS2"])),int(np.round(0.7*head_out["NAXIS2"])))
        x_medlimits=(int(np.round(0.3*head_out["NAXIS1"])),int(np.round(0.7*head_out["NAXIS1"])))
        m = np.median(image_frames[:,y_medlimits[0]:y_medlimits[1],
                      x_medlimits[0]:x_medlimits[1]])
        for i in range(nfiles):
            mcur=np.median(image_frames[i, y_medlimits[0]:y_medlimits[1],
                      x_medlimits[0]:x_medlimits[1]])
            image_frames[i, :, :] = image_frames[i, :, :] * m / float(mcur)
            norm[i] = float(mcur)/m
            # if file_error_out:
            #     error_squared[i, :, :] = error_squared[i, :, :] * m**2/mcur**2

    if (do_xshift or do_yshift) and (image_frames.shape[0] > 1):
        dx,dy=derive_shift(image_frames, ref_frame=0, bin=20, max_ampl=(10.,10.), do_shift = (do_xshift,do_yshift))
        image_frames = shift_images(image_frames, -dx,-dy)
        if file_error_out:
            error_squared = shift_images(error_squared, -dx,-dy)
        add_row="\n".join(["   File {} shifted by dx = {:.2f} and dy = {:.2f} pixels".format(filein[i],-dx[i],-dy[i])
                 for i in range(nfiles)])
        message("Several images of {} type were shifted along X or Y:\n{}".format(imgtype,add_row))

    if nfiles != 1:
        if (cr_mode == 'sigma') and nfiles >= 5:
            image_out=combine_sigma(image_frames, cr_treshold = cr_treshold)
        else:
            image_out=combine_pairs(image_frames, cr_treshold = cr_treshold*10, ncycles = 2)
        if file_error_out:
            error_image = np.sqrt(np.nansum(error_squared/norm[:,None,None]**2,axis=0)+
                                  np.abs(np.nansum(image_frames/norm[:,None,None],axis=0)/gain))
    else:
        image_out = image_frames[0,:,:]
        if file_error_out:
            error_image= np.sqrt(error_squared[0,:,:]/norm[0]**2+np.abs(image_out)/norm[0]/gain)


    image_out=remove_overscan(image_out,mode=mode,cut=force_cut)
    if mode=='TDS_r':
        image_out=np.fliplr(image_out)
    if file_error_out:
        error_image=remove_overscan(error_image,mode=mode,cut=force_cut)
        if mode == 'TDS_r':
            error_image = np.fliplr(error_image)
    if imgtype in ['sunsky','flat']:
        norm = np.median(image_out)
    else:
        norm = np.median(image_out[image_out>np.nanstd(image_out)])

    if imgtype == 'flat':
        head_out.add_history("LONGRED: Superflat is created; {} initial frames were used".format(nfiles))
    else:
        head_out.add_history("LONGRED: {} created; {} initial frames were used".format(imgtype.capitalize(),nfiles))
    fits.writeto(os.path.join(w_dir,fileout), np.float32(image_out / norm), head_out, output_verify='fix', overwrite=True)
    if file_error_out:
        fits.writeto(os.path.join(w_dir, file_error_out), np.float32(error_image / norm), head_out, output_verify='fix',
                     overwrite=True)

    message("{} frame {} is created using {} initial files".format(imgtype.capitalize(), os.path.join(w_dir, fileout),nfiles))
    return True







def cre_ini(filein=None, fileout_prefix='obj', fileout_suffix='', w_dir=os.getcwd(), raw_dir='',filebias=None,
            filedark=None,  mode="TDS", errorout_suffix=None, gain=None, readnoise=None, cr_clean=True, cut=None,overscan=None):
    """
    Prepare initial frames for object and star: subtract bias and dark
    """

    if not filein:
        return (False, None)

    nfiles = len(filein)
    all_prepared_files=[]
    if filebias:
        with fits.open(os.path.join(w_dir, filebias)) as hdu_bias:
            bias = hdu_bias[0].data
            if readnoise is None:
                bias_squared_error = hdu_bias[0].header.get("RNOISE") ** 2
            else:
                if gain is None:
                    gain_use=hdu_bias[0].header.get("GAIN")
                else:
                    gain_use=float(gain)
                bias_squared_error = (readnoise/gain_use)** 2
    else:
        if readnoise is None:
            bias_squared_error = 0
        else:
            if gain is None:
                with fits.open(os.path.join(w_dir, filein[0])) as hdu:
                    gain_use = hdu[0].header.get("GAIN")
            else:
                gain_use = float(gain)
            bias_squared_error = (readnoise / gain_use) ** 2
        bias = 0
    if filedark:
        with fits.open(os.path.join(w_dir, filedark)) as hdu_dark:
            dark = hdu_dark[0].data
    else:
        dark = 0

    for i,f in enumerate(filein):
        message(" .......Process: {} of {}".format(i + 1, nfiles),noheader=True)
        with fits.open(os.path.join(raw_dir,f)) as hdu:
            texp = hdu[0].header.get("EXPTIME")
            if not texp:
                if 'EXPOSURE' in hdu[0].header.keys():
                    texp = hdu[0].header.get("EXPOSURE")
                    del hdu[0].header["EXPOSURE"]
                    hdu[0].header["EXPTIME"]=texp
                else:
                    texp=0

            if not gain:
                gain_use=hdu[0].header.get("GAIN")
            else:
                gain_use=gain
            hdu[0].header['BZERO'] = 0
            hdu[0].header['BSCALE'] = 1.0
            hdu[0].header["BUNITS"] = "counts"
            if 'FILENAME' in hdu[0].header.keys():
                del hdu[0].header['FILENAME']
            force_overscan = None
            if cut is not None and any(cut) and overscan:
                force_overscan = [0,hdu[0].header["NAXIS1"], 0,  hdu[0].header["NAXIS2"]]
                applied=0
                for j, v in enumerate(cut):
                    if v is not None and v > 2:
                        applied+=1
                        if j == 0 or j == 2:
                            force_overscan[j] = v // 2
                        else:
                            force_overscan[j] -= v // 2
                if applied == 0:
                    force_overscan=None
            if cut is not None and any(cut):
                force_cut = [0, hdu[0].header["NAXIS1"], 0, hdu[0].header["NAXIS2"]]
                applied = 0
                for j, v in enumerate(cut):
                    if v is not None and v > 0:
                        applied += 1
                        if j == 0 or j == 2:
                            force_cut[j] = v
                        else:
                            force_cut[j] -= v
                if applied == 0:
                    force_cut = None
            else:
                force_cut = None

            image=check_overscan(hdu[0].data-bias,mode=mode,force_overscan=force_overscan)-dark*texp
            if errorout_suffix:
                error_image=remove_overscan(np.sqrt(np.abs(image+dark*texp)/gain_use+bias_squared_error),mode=mode,cut=force_cut)
            image = remove_overscan(image, mode=mode,cut=force_cut)
            if filebias:
                hdu[0].header.add_history("LONGRED: bias subtracted: {}".format(filebias))
            if filedark:
                hdu[0].header.add_history("LONGRED: dark frame subtracted: {}".format(filedark))
            if cr_clean:
                if errorout_suffix:
                    errors=error_image
                else:
                    errors=None
                image,crmask=lacosmic(image, 1.5,5.,3., error=errors)

                if i==0:
                    cr_mask_out = np.ndarray((nfiles,image.shape[0],image.shape[1]),dtype=bool)
                cr_mask_out[i,:,:]=crmask
                hdu[0].header.add_history("LONGRED: CR cleaned using lacosmic")

            if mode == 'TDS_r':
                image=np.fliplr(image)

            fits.writeto(os.path.join(w_dir, "{}_{:02}{}.fits".format(fileout_prefix, i,fileout_suffix)), np.float32(image),
                         header=hdu[0].header, overwrite=True)#, output_verify='fix')
            all_prepared_files.append("{}_{:02}{}.fits".format(fileout_prefix, i,fileout_suffix))
            if errorout_suffix:
                if mode == 'TDS_r':
                    error_image = np.fliplr(error_image)
                fits.writeto(os.path.join(w_dir, "{}_{:02}{}.fits".format(fileout_prefix, i, errorout_suffix)), np.float32(error_image),
                         header=hdu[0].header, overwrite=True, output_verify='fix')
                all_prepared_files.append("{}_{:02}{}.fits".format(fileout_prefix, i, errorout_suffix))
    f=os.path.join(w_dir, "{}_{}.fits".format(fileout_prefix, 'crmask'))
    if os.path.isfile(f):
        os.remove(f)
    if cr_clean:
        all_prepared_files.append("{}_{}.fits".format(fileout_prefix, 'crmask'))
        fits.writeto(f, np.float32(cr_mask_out),overwrite=True, output_verify='fix')
    if len(all_prepared_files)==0:
        all_prepared_files=None
    return (True, all_prepared_files)



def get_specslope(filein, transform_matrix=None, do_plot=True, w_dir=os.getcwd(),oversample=50):
    """
    Trace most prominent spectrum on obj or star frames to get slope of the spectrum
    """
    slope=0
    fit_g=fitting.LevMarLSQFitter()
    if type(transform_matrix) == str:
        with fits.open(os.path.join(w_dir,transform_matrix)) as hdu:
            transform_matrix=hdu[0].data
    if do_plot:
        plt.figure(figsize=(5,5))
        plt.subplot(111)
    for ind, f in enumerate(filein):
        with fits.open(os.path.join(w_dir,f)) as hdu:
            im=np.array(hdu[0].data,dtype=float)
            if transform_matrix is not None:
                im=do_distortion(im,transform_matrix)

        if ind ==0:
            s = im.shape
            win_analyse = s[1]/3
            win = 5
            ncols = 15
            yborders = 100
            win_around_peak=15
            bin_over_peak=20
            cols = np.linspace(s[1] / 2 - win_analyse, s[1] / 2 + win_analyse + 1, ncols, dtype=int)
            peaks = np.ndarray((len(filein),ncols),dtype=float)
            col_reper=cols[int(round((ncols-1)/2))]

        for cind, col in enumerate(cols):
            cur_vec=np.nanmedian(im[yborders:-yborders,col-win:col+win+1],axis=1)#,oversample,1)
            peaks_cur, _ = signal.find_peaks(cur_vec, prominence=[(np.nanmax(cur_vec)-np.nanmin(cur_vec))*0.6,None])
            if len(peaks_cur) == 1:
                cur_vec = cur_vec[peaks_cur[0] - win_around_peak: peaks_cur[0] + win_around_peak + 1]
                # peaks[ind,cind]=float(peaks_cur[0])/oversample+yborders
                yy = np.linspace(peaks_cur[0] - win_around_peak, peaks_cur[0] + win_around_peak, (2 * win_around_peak + 1) * bin_over_peak)
                g_init = models.Gaussian1D(amplitude=np.max(signal.resample_poly(cur_vec,oversample,1)), mean=peaks_cur[0],
                                           stddev=win_around_peak / 5.) + models.Const1D(np.min(cur_vec))
                #
                g = fit_g(g_init, yy + 0.5, signal.resample_poly(cur_vec, bin_over_peak, 1))
                peaks[ind,cind] = np.float(g[0].mean.value)
            else:
                peaks[ind,cind]=np.nan
            if col == col_reper:
                peaks_reper=peaks[ind,cind]
        peaks[ind,:]-=peaks_reper

        if do_plot:
            plt.scatter(cols,peaks[ind,:])

    x = np.array(np.tile(cols, len(filein)), dtype=float)
    p,z=goodpoly(x,peaks.flatten(),None,1)
    # z = np.polyfit(x, peaks.flatten(), 1)
    # p = np.poly1d(z)

    if do_plot:
        plt.plot([0,s[1]],p([0,s[1]]),'--')
        plt.xlim(0,s[1])
        plt.ylim(-1.1,1.1)
        plt.draw()
        plt.pause(0.001)

    return np.flip(z)









def calc_geometry(filein, fileerr= None, fileout_transform='transform_matrix.fits', w_dir=os.getcwd(),do_plot=True, dotmask=False,
                  nstripes = 50, win=5, snr=100, peaks_mindist=50,prominence=[None,None],
                  poly_order=2, oversample=50, trace_star=None, save_test=None):
    """
    Here we calculate lines curvature based on the calibration (neon etc.) images.
    Алгоритм такой:
    Вся щель разбивается на nstripes полос. В каждой полосе, шириной 2*win+1 пикселей, интегрируется спектр
    В нем ищем пики с отношением S/N выше заданного в snr,
    с минимальным расстоянием между ними в peaks_mindist и prominance
    Далее пока так: выбираем только те строки, где число найденных пиков равно медианному значению среди всех строк
    Считаем эти пики истинными. Уточняем положение каждого из этих пиков путем
    вписывания гауссианы в кусок спектра с оверсемплингом bin=50 и шириной вокруг пика ±8 пикселей
    Фиттируем полиномом вдоль щели каждую найденную линию - эти полиномы и задают геометрическую модель
    """

    with fits.open(os.path.join(w_dir,filein)) as hdu:
        image=hdu[0].data
        header=hdu[0].header
        if save_test:
            image_test=image.copy()
            header_test=header.copy()
    if not fileerr:
        err=np.sqrt(np.abs(image))
    else:
        with fits.open(os.path.join(w_dir,fileerr)) as hdu:
            err=hdu[0].data
    s=image.shape
    if not win:
        win = 0

    if prominence[0] == 0:
            prominence[0]=None
    if prominence[1] == 0:
            prominence[1]=None

    xborders=30
    peaks=[0]*nstripes
    rows=np.linspace(30, s[0]-30, nstripes,dtype=int)
    for ind,row in enumerate(rows):
        cur_ylim = ((int(row - win)), int(row + win))
        if win>1:
            cur_vec = np.nanmedian(image[cur_ylim[0]:cur_ylim[1] + 1, :], axis=0)
        else:
            cur_vec = np.nanmean(image[cur_ylim[0]:cur_ylim[1] + 1, :], axis=0)
        cur_errvec = np.sqrt(np.nansum(err[cur_ylim[0]:cur_ylim[1] + 1, :] ** 2, axis=0))/(2*win+1)
        peaks[ind], _ = signal.find_peaks(signal.resample_poly(cur_vec[xborders:-xborders],oversample,1),
                                          height=(signal.resample_poly(cur_errvec[xborders:-xborders] * snr,oversample,1), None),
                                          distance=peaks_mindist*oversample,prominence=prominence)

    npeaks=[len(p) for p in peaks]
    med_npeaks=int(np.median(npeaks))

    if med_npeaks < 3:
        return False
    rows = np.array([r for ind, r in enumerate(rows) if len(peaks[ind]) == med_npeaks],dtype=float)
    if len(rows) < 5:
        return False
    peaks = np.array([p for p in peaks if len(p) == med_npeaks],dtype=float)/oversample+xborders


    if do_plot:
        add_panel = 0
        if dotmask:
            add_panel = 100
        plt.figure(figsize=(10, 10+add_panel/50))
        plt.subplot(211+add_panel)
        x = np.arange(s[1])
        ind_row=int(np.round(len(rows)/2))
        cur_ylim = ((int(rows[ind_row] - win)), int(rows[ind_row] + win))
        peaks_tmp=np.array(np.round(peaks[ind_row,:]),dtype=int)
        if win > 1:
            cur_vec = np.nanmedian(image[cur_ylim[0]:cur_ylim[1] + 1, :], axis=0)
        else:
            cur_vec = np.nanmean(image[cur_ylim[0]:cur_ylim[1] + 1, :], axis=0)
        cur_errvec = np.sqrt(np.nansum(err[cur_ylim[0]:cur_ylim[1] + 1, :] ** 2, axis=0)) / (2 * win + 1)
        plt.plot(x, cur_vec)
        plt.plot(peaks_tmp, cur_vec[peaks_tmp], "x")
        plt.plot(x, cur_errvec)
        plt.xlim((0,s[1]))
        if prominence[1]:
            plt.ylim((0,prominence[1]))
        plt.title("Line curvature on {}: Nstr={}, Win={}, SNR={}, Pdist={}, Prom.={}".format(filein,nstripes,win,snr,
                                                                                        peaks_mindist,prominence))
        plt.subplot(212+add_panel)

    yy_slit=np.arange(s[0])
    xx_slit=np.arange(s[1])
    lines_coeff=np.ndarray((poly_order+1,med_npeaks),dtype=float)
    lines_xnew=np.ndarray(med_npeaks,dtype=float) # adopted X position in new frame for each line
    yref_new=s[0]/2

    for ind in range(med_npeaks):
        if peaks_mindist:
            rec=np.abs(peaks[:,ind]-np.median(peaks[:,ind]))>peaks_mindist
            peaks[rec,ind]=np.nan
        p, cf = goodpoly(rows, peaks[:,ind], order=poly_order, sigma=3.)
        lines_coeff[:,ind]=np.flip(cf)#np.polyfit(rows,peaks[:,ind],poly_order)
        #p=np.poly1d(lines_coeff[:,ind])
        lines_xnew[ind]=p(yref_new)
        if do_plot:
            plt.plot(peaks[:,ind], rows, ".")
            plt.plot(p(yy_slit),yy_slit,'--')

    if do_plot:
        plt.xlim((0,s[1]))
        if not add_panel:
            plt.xlabel("X, pix")
            plt.draw()
            plt.pause(0.001)

        else:
            plt.subplot(313)

    # Propagate the lines polynomial coefficient for the rest of image
    #   and compute mapping grid matrix where F([0,1],Y,X) =  ([DY,DX]) in respect to original frame
    correspondance_matrix = np.zeros((2, s[0], s[1]), dtype=float)
    distortion_coeff=np.ndarray((poly_order+1,s[1]))
    for cf in range(poly_order+1):
        #p=np.poly1d(np.polyfit(lines_xnew,lines_coeff[cf,:],poly_order+1))
        p,_ = goodpoly(lines_xnew, lines_coeff[cf, :], order=poly_order + 1,sigma=3)
        distortion_coeff[cf,:] = p(xx_slit)
    if do_plot:
        distortion_coeff_lines_for_plot=distortion_coeff.copy()
    for x in xx_slit:
        p=np.poly1d(distortion_coeff[:,x])
        correspondance_matrix[1,:,x]+=(p(yy_slit)-x)



    # IF 13dots is found - then apply distortion and calculate correction along Y
    if dotmask:
        with fits.open(os.path.join(w_dir,dotmask)) as hdu:
            image = hdu[0].data#do_distortion(hdu[0].data, correspondance_matrix)
            # if apply_transform:
            #     image=do_distortion(image,correspondance_matrix)#_ini)
        yborders = 15
        nstripes_dots=40
        win=3
        win_around_peak = 15
        cols_all = np.linspace(30, s[1] - 30, nstripes_dots, dtype=int)
        cols=np.array([],dtype=float)
        peaks=[]
        ndots=13
        for ind, col in enumerate(cols_all):
            cur_xlim = ((int(col - win)), int(col + win))
            if win > 1:
                cur_vec = np.nanmedian(image[:,cur_xlim[0]:cur_xlim[1] + 1], axis=1)
            else:
                cur_vec = np.nanmean(image[:, cur_xlim[0]:cur_xlim[1] + 1], axis=1)
            # peaks_tmp, _ = signal.find_peaks(cur_vec[yborders:-yborders],
            #                     distance=20, prominence=(np.max(cur_vec)/10.,None))
            peaks_tmp, _ = signal.find_peaks(signal.resample_poly(cur_vec[yborders:-yborders],oversample,1),
                                             distance=20*oversample, prominence=(np.max(cur_vec) / 10., None))
            if len(peaks_tmp) == ndots:
                cols=np.append(cols,col)
                peaks.append(peaks_tmp)
        if len(cols) < 5:
            if do_plot:
                plt.draw()
                plt.pause(0.001)
                return False
        peaks = np.array(peaks, dtype=float)/oversample+yborders
        win_around_peak =20
        # bin_over_peak = 0
        dots_coeff = np.ndarray((poly_order + 1, ndots), dtype=float)
        dots_ynew = np.ndarray(ndots, dtype=float)  # adopted Y position in new frame for each line
        xref_new = s[1] / 2

        for ind in range(ndots):
            dots_coeff[:, ind] = np.polyfit(cols, peaks[:, ind], poly_order)
            p = np.poly1d(dots_coeff[:, ind])
            dots_ynew[ind] = p(xref_new)
            if do_plot:
                plt.plot(cols,peaks[:, ind], ".")
                plt.plot(xx_slit,p(xx_slit), '--')

        if do_plot:
            plt.xlim((0, s[1]))
            plt.xlabel("X, pix")
            plt.draw()
            plt.pause(0.001)




        # Propagate the lines polynomial coefficient for the rest of image
        #   and compute mapping grid matrix where F([0,1],Y,X) =  ([DY,DX]) in respect to original frame
        distortion_coeff = np.ndarray((poly_order + 1, s[0]))
        for cf in range(poly_order + 1):
            p = np.poly1d(np.polyfit(dots_ynew, dots_coeff[cf, :], poly_order+1))
            distortion_coeff[cf, :] = p(yy_slit)
        for y in yy_slit:
            p = np.poly1d(distortion_coeff[:, y])
            correspondance_matrix[0, y, :] += (p(xx_slit) - y)

    if do_plot:
        if dotmask:
            add_plotspace=1
        else:
            add_plotspace = 0

        plt.figure(figsize=(6+3*add_plotspace, 9))

        for ind in range(poly_order+1):
            plt.subplot(1+(1+add_plotspace)*ind+10*(1+add_plotspace)+(poly_order+1)*100)
            plt.plot(lines_xnew, lines_coeff[ind, :], '.')
            plt.plot(xx_slit, distortion_coeff_lines_for_plot[ind, :], '-', label='order = {}'.format(poly_order - ind))
            plt.legend()
            if ind == (poly_order + 1):
                plt.xlim((0, s[1]))
                plt.xlabel("X, pix")
            if ind == 0:
                plt.title("Polynome coefficients for calibration lamp")

            if dotmask:
                plt.subplot(2 + (1 + add_plotspace) * ind + 10 * (1 + add_plotspace) + (poly_order + 1) * 100)
                plt.plot(dots_ynew,dots_coeff[ind, :],'.')
                plt.plot(yy_slit,distortion_coeff[ind,:],'-', label='order = {}'.format(poly_order-ind))
                plt.legend()
                if ind == (poly_order+1):
                    plt.xlim((0, s[0]))
                    plt.xlabel("Y, pix")
                if ind == 0:
                    plt.title("Polynome coefficients for 13dots mask")

        plt.draw()
        plt.pause(0.001)


    if trace_star is not None:
        ### IF we set objects to analyse for gradient, then we do that and add it to transform matrix
        message(".... Checking for object slope", noheader=True)
        slope_param=get_specslope(trace_star,transform_matrix=correspondance_matrix,do_plot=do_plot,w_dir=w_dir)
        p=np.poly1d(slope_param)
        if len(slope_param) ==2:
            message(".... Slope calculated: {:.2f}px in Y direction per 1000px in X direction.".format(
                slope_param[0] * 1000.), noheader=True)
            for x in xx_slit:
                correspondance_matrix[0, :, x] += p(x)# - y)
        else:
            message(".... Something went wrong with slope calculation.", noheader=True)

    fits.writeto(os.path.join(w_dir, fileout_transform), np.float32(correspondance_matrix), output_verify='fix',
                         overwrite=True)

    if save_test:
        image_test=do_distortion(image_test,correspondance_matrix)
        fits.writeto(os.path.join(w_dir,save_test),image_test,header_test, overwrite=True)

    return True




def calc_norm_flat(filein, fileerr=None, filein_transform=None, fileout = None, fileout_err = None, w_dir=os.getcwd(),
                                                     nosmo=False, do_plot=True, blue_pix_limit=None, save_warp=None):
    win_y = 10
    win_smo = 31
    with fits.open(os.path.join(w_dir, filein)) as hdu:
        image = hdu[0].data
        head = hdu[0].header
        s = image.shape
        central_row = int(round(s[0] / 2))
        if blue_pix_limit is not None or save_warp is not None:
            low=0
            high=s[0]-1
        else:
            low=np.max([0,central_row-win_y*5])
            high = np.min([s[0]-1, central_row + win_y*5])
        image_cut=image[low:high+1,:]

    if filein_transform:
        with fits.open(os.path.join(w_dir, filein_transform)) as hdu_matrix:
            correspondance_matrix = hdu_matrix[0].data
            image_cut = do_distortion(image_cut, correspondance_matrix[:,low:high+1,:])

    if fileerr:
        with fits.open(os.path.join(w_dir, fileerr)) as hdu:
            err = hdu[0].data
            head_err = hdu[0].header


    xx=np.arange(s[1])
    yy = np.arange(s[0])

    norm_flat=np.ndarray(s,dtype=float)
    central_row_cutted=int(round(image_cut.shape[0] / 2))

    central_flat=np.nanmean(image_cut[central_row_cutted-win_y:central_row_cutted+win_y+1,:],axis=0)
    if not nosmo:
        central_flat_smo=signal.savgol_filter(central_flat, win_smo, 3)
    else:
        central_flat_smo = central_flat.copy()
    norm_flat[:,:]=central_flat_smo[None,:]

    if blue_pix_limit is not None or save_warp is not None:
        image_cut = image_cut / norm_flat
        if save_warp is not None:
            err = do_distortion(err, correspondance_matrix)/norm_flat
            fits.writeto(os.path.join(w_dir, save_warp[0]), np.float32(image_cut), header=head, overwrite=True)
            fits.writeto(os.path.join(w_dir, save_warp[1]), np.float32(err), header=head_err, overwrite=True)
            message(".... Only test distortion corrected flat was created")
            return False

    if filein_transform:
        norm_flat = do_distortion(norm_flat, correspondance_matrix,reverse=True)
        image=image/norm_flat
        err = err / norm_flat


    if blue_pix_limit:
        if blue_pix_limit>=(s[1]-2):
            blue_pix_limit=(s[1]-2)
        red_pix_limit = np.min([blue_pix_limit+int(round(0.03*s[1])),s[1]-1])
        for row in yy:
            image_cut[row,0:blue_pix_limit+1]=np.nanmedian(image_cut[row,blue_pix_limit:red_pix_limit+1])
        if filein_transform:
            image=do_distortion(image_cut,correspondance_matrix,reverse=True)
        else:
            image=image_cut
    image.ravel()[np.isnan(image.ravel()) + ~np.isfinite(image.ravel()) +  (image.ravel()==0)]=1.
    err.ravel()[np.isnan(err.ravel()) + ~np.isfinite(err.ravel()) + (err.ravel() == 0)] = np.nanmax(err.ravel()[np.isfinite(err.ravel())])
    fits.writeto(os.path.join(w_dir, fileout), np.float32(image), header=head, overwrite=True)
    if fileerr:
        fits.writeto(os.path.join(w_dir, fileout_err), np.float32(err), header=head_err, overwrite=True)

    return True


def flat_normilize(filein, fileflat, fileerr=None, flaterr=None, fileout=None, fileerr_out=None, w_dir=os.getcwd()):
    with fits.open(os.path.join(w_dir,fileflat)) as hdu:
        flat=hdu[0].data
    if flaterr is not None:
        if flaterr:
            with fits.open(os.path.join(w_dir,flaterr)) as hdu:
                flaterr=hdu[0].data
    else:
        flaterr=flat*0
    for ind,f in enumerate(filein):
        with fits.open(os.path.join(w_dir,f)) as hdu:
            hdu[0].header.add_history("LONGRED: normalized by flat {}".format(fileflat))
            fits.writeto(os.path.join(w_dir,fileout[ind]),np.float32(hdu[0].data/flat),header=hdu[0].header,overwrite=True, output_verify='fix')
            if fileerr is not None:
                with fits.open(os.path.join(w_dir,fileerr[ind])) as hdu_err:
                    hdu_err[0].header.add_history("LONGRED: normalized by flat {}".format(fileflat))
                    err=np.sqrt(hdu_err[0].data**2/flat**2+(hdu[0].data*flaterr/flat**2)**2)
                    fits.writeto(os.path.join(w_dir,fileerr_out[ind]),np.float32(err),header=hdu_err[0].header,overwrite=True, output_verify='fix')

    return True



def comp_shifts(filein, file_transform=None, do_shift=(True, True), shift_out='shifts.fits',ref_index=None,
                file_numbers=None, w_dir=os.getcwd(), max_shift=(1,100), trace_atm_disp=None,
                plot_atm_disp=True, ypos_atm_disp=None, trace_atm_disp_order=5,save_test=None):
    """
    Derive shifts between different exposures relative to one mentioned in ref_index or to the first.
    Store them to output file 'shift_out'
    Then trace atmospheric dispersion if necessary
    """
    nfiles=len(filein)
    shifts_x=[]
    shifts_y=[]
    if ref_index is not None:
        filein=list(filein)
        file_numbers=list(file_numbers)
        filein.insert(0,filein.pop(ref_index))
        file_numbers.insert(0,file_numbers.pop(ref_index))
    wid=50
    wid_x=50
    wid_y=100
    wid_trace_atm=35
    bin_trace_atm=10
    npoints_trace_atm=100
    borders=20
    #borders_trace_atm=100
    with fits.open(os.path.join(w_dir,file_transform)) as hdu:
        tmatrix=hdu[0].data
    for ind,f in enumerate(filein):
        with fits.open(os.path.join(w_dir,f)) as hdu:
            image=hdu[0].data
            if ind ==0:
                header = hdu[0].header
                s=image.shape
                ccol=int(round(s[1]/2))
                crow = int(round(s[0] / 2))
                borders_x = int(round(s[1]/8))
        img_y=do_distortion(image[:,ccol-wid_x:ccol+wid_x+1],tmatrix[:,:,ccol-wid_x:ccol+wid_x+1])
        vec_y=np.nanmedian(img_y[borders:-borders,wid_x-wid:wid_x+wid+1],axis=1)
        if ind == 0:
            vec_y_ref=vec_y.copy()
            # if ypos_atm_disp is not None:
            #     y_pos_max=[ypos_atm_disp-borders]
            # else:
            y_pos_max, _ = signal.find_peaks(vec_y,distance=vec_y.shape[0]/2,
                                             prominence=[(np.nanmax(vec_y) - np.nanmin(vec_y)) * 0.8, None])
            if len(y_pos_max) > 0:
                y_pos_max=y_pos_max[0]+borders
            else:
                y_pos_max = crow
            if trace_atm_disp and (ypos_atm_disp is None):
                ypos_atm_disp = y_pos_max
            dy=0.
            wid_for_dy = np.min([50,s[0]-1-borders-y_pos_max,y_pos_max-borders])
        else:
            if do_shift[1]:
                dy=derive_vecshift(vec_y[y_pos_max-borders-wid_for_dy:y_pos_max-borders+wid_for_dy+1],
                                   vec_y_ref[y_pos_max-borders-wid_for_dy:y_pos_max-borders+wid_for_dy+1],bin=50,max_ampl=max_shift[1])
            else:
                dy=0.

        if do_shift[0]: # or (trace_atm_disp and ((y_pos_max+wid_trace_atm - wid_y)< ypos_atm_disp_curfile < (y_pos_max + wid_y-wid_trace_atm+1))):
            # cut_bottom = np.max([0, y_pos_max - wid_y])
            # cut_top = np.min([s[0] - 1, y_pos_max + wid_y])
            cut_bottom = np.max([0, crow - wid_y])
            cut_top = np.min([s[0] - 1, crow + wid_y])
            img_x = do_distortion(image[cut_bottom:cut_top + 1,:],
                                  tmatrix[:,cut_bottom:cut_top + 1,:])
            # if trace_atm_disp:
            #     ypos_atm_disp_on_subarray=ypos_atm_disp_curfile-y_pos_max + wid_y
            # if do_shift[0]:
            vec_x = np.nanmedian(img_x[y_pos_max-cut_bottom - wid:y_pos_max-cut_bottom + wid + 1, borders_x:-borders_x], axis=0)
            if ind == 0:
                vec_x_ref=vec_x.copy()
                dx=0.
            else:
                if do_shift[0]:
                    dx=derive_vecshift(vec_x,vec_x_ref,bin=20,max_ampl=max_shift[0])
                else:
                    dx=0.
            # else:
            #     dx=0.
        else:
            dx=0.

        shifts_x.append(dx)
        shifts_y.append(dy)

        if trace_atm_disp:
            # Perfrom transformation to calculate median file and to trace further atm_disp
            cur_tmatrix=tmatrix.copy()
            cur_tmatrix[1,:,:]+=dx
            cur_tmatrix[0, :, :] += dy
            if ind == 0:
                cut_bottom_atm=np.max([0,ypos_atm_disp - wid_y])
                cut_top_atm = np.min([s[0] - 1, ypos_atm_disp + wid_y])

                ypos_atm_disp_on_subarray = ypos_atm_disp - cut_bottom_atm
                wid_trace_atm = np.min([wid_trace_atm, ypos_atm_disp_on_subarray, cut_top_atm - ypos_atm_disp])
                img_for_atmdisp=np.ndarray((nfiles,2*wid_trace_atm+1,cur_tmatrix.shape[2]),dtype=float)

            img_for_atmdisp[ind,:,:]=do_distortion(image[cut_bottom_atm:cut_top_atm + 1, :],
                                  cur_tmatrix[:, cut_bottom_atm:cut_top_atm + 1, :])[ypos_atm_disp_on_subarray-wid_trace_atm:
                                                                             ypos_atm_disp_on_subarray+wid_trace_atm+1,:]

    if trace_atm_disp:
        img_atmdisp_combined=np.median(img_for_atmdisp,axis=0)
        img_atmdisp_combined-=np.nanmin(img_atmdisp_combined.ravel())
        fits.writeto(os.path.join(w_dir,'test.fits'),img_atmdisp_combined,overwrite=True)
        xx_trace_atm = np.linspace(borders_x + wid, img_atmdisp_combined.shape[1] - borders_x - wid,
                                   npoints_trace_atm, dtype=int)
        peaks_trace_atm = np.ndarray(npoints_trace_atm, dtype=float)

        fit=fitting.LevMarLSQFitter()
        g_init=models.Gaussian1D()+models.Polynomial1D(1)
        yy_atm_disp = np.linspace(0, img_atmdisp_combined.shape[0], img_atmdisp_combined.shape[0])
        for xind, x in enumerate(xx_trace_atm):
            cur_vec = np.nanmedian(img_atmdisp_combined[:, x - 30:x + 30 + 1], axis=1)
            # cur_vec -=np.nanmedian([cur_vec[:7],cur_vec[-7:]])
            cur_vec=cur_vec/np.nanmax(cur_vec)
            g_init.mean_0.value=np.argmax(cur_vec)
            g_init.amplitude_0.value = 1.
            g_init.stddev_0.value = 5.
            g=fit(g_init,yy_atm_disp,cur_vec)
            peaks_trace_atm[xind]=g.mean_0.value-wid_trace_atm
            #peaks_trace_atm[xind] = line_center(cur_vec, bin=bin_trace_atm, win=5) - wid_trace_atm
        peaks_trace_atm[np.abs(peaks_trace_atm-np.median(peaks_trace_atm))>10] = np.nan
        # x = np.array(np.tile(xx_trace_atm, len(filein)), dtype=float)
        p, cf = goodpoly(xx_trace_atm, peaks_trace_atm, order=trace_atm_disp_order,sigma=1.5)
        cf[0] -= p(s[1] / 2)
        if plot_atm_disp:
            plt.figure(figsize=(6, 3))
            plt.subplot(111)
            plt.title("Tracing of atmosphere dispersion")
            xx = np.linspace(0, s[1], 100)
            plt.scatter(xx_trace_atm, peaks_trace_atm[:] - p(s[1] / 2), s=4.)
            pp = np.poly1d(np.flip(cf))
            plt.plot(xx, pp(xx), '--')
            plt.xlim(0, s[1])
            plt.ylim(-7, 2)
            plt.xlabel("X, pix")
            plt.ylabel(r"$\Delta$Y, pix")
            plt.draw()
            plt.pause(0.001)
        # cf[0] -= p(s[1] / 2)
        names = ["cf_{}".format(i) for i in range(trace_atm_disp_order + 1)]
        ascii.write(cf, os.path.join(w_dir, trace_atm_disp), names=names, overwrite=True)






        # if trace_atm_disp:
        #     # if not ((y_pos_max+wid_trace_atm - wid_y)< ypos_atm_disp_curfile < (y_pos_max + wid_y-wid_trace_atm+1)):
        #     cut_bottom=np.max([0,ypos_atm_disp_curfile - wid_y])
        #     cut_top = np.min([s[0]-1, ypos_atm_disp_curfile + wid_y])
        #     img_x = do_distortion(image[cut_bottom:cut_top+ 1, :],
        #                     tmatrix[:, cut_bottom:cut_top + 1, :])
        #     ypos_atm_disp_on_subarray = ypos_atm_disp_curfile-cut_bottom
        #     wid_trace_atm=np.min([wid_trace_atm,ypos_atm_disp_on_subarray,cut_top-ypos_atm_disp_curfile])
        #     img_x=img_x[ypos_atm_disp_on_subarray-wid_trace_atm:ypos_atm_disp_on_subarray+wid_trace_atm+1,:]
        #     if ind==0:
        #         xx_trace_atm = np.linspace(borders_x+wid, img_x.shape[1] + borders_x-wid, npoints_trace_atm, dtype=int)
        #         peaks_trace_atm=np.ndarray((nfiles,npoints_trace_atm),dtype=float)
        #         if plot_atm_disp:
        #             plt.figure(figsize=(6,3))
        #             plt.subplot(111)
        #             plt.title("Tracing of atmosphere dispersion")
        #
        #     for xind,x in enumerate(xx_trace_atm):
        #         cur_vec=np.nanmedian(img_x[:,x-30:x+30+1],axis=1)
        #         peaks_trace_atm[ind, xind] = line_center(cur_vec, bin=bin_trace_atm, win=15)-wid_trace_atm-dy+round(dy)


    #
    # if trace_atm_disp:
    #     x = np.array(np.tile(xx_trace_atm, len(filein)), dtype=float)
    #     p,cf = goodpoly(x, peaks_trace_atm.ravel(), order=trace_atm_disp_order)
    #     cf[0] -= p(s[1] / 2)
    #     if plot_atm_disp:
    #         xx = np.linspace(0, s[1], 100)
    #         for ind in range(len(filein)):
    #             plt.scatter(xx_trace_atm, peaks_trace_atm[ind, :]-p(s[1]/2),s=4.)
    #         pp=np.poly1d(np.flip(cf))
    #         plt.plot(xx, pp(xx), '--')
    #         plt.xlim(0, s[1])
    #         plt.ylim(-3, 3)
    #         plt.xlabel("X, pix")
    #         plt.ylabel(r"$\Delta$Y, pix")
    #         plt.draw()
    #         plt.pause(0.001)
    #     #cf[0] -= p(s[1] / 2)
    #     names = ["cf_{}".format(i) for i in range(trace_atm_disp_order + 1)]
    #     ascii.write(cf, os.path.join(w_dir, trace_atm_disp), names=names, overwrite=True)
    file_numbers=np.array(file_numbers)
    shifts_x=np.array(shifts_x)
    shifts_y = np.array(shifts_y)
    srt=np.argsort(file_numbers)
    ascii.write([file_numbers[srt], shifts_x[srt], shifts_y[srt]],os.path.join(w_dir,shift_out),
                names=['filenum', 'dx', 'dy'], overwrite=True)
    if save_test:
        xx=np.arange(s[1])
        if trace_atm_disp: tmatrix[0,:,:]+=np.tile(p(xx),s[0]).reshape((s[0],s[1]))
        image=do_distortion(image,transform_matrix=tmatrix)
        fits.writeto(os.path.join(w_dir,save_test),np.float32(image),header,overwrite=True)

    return True





def calc_dispersion(filein, fileerr=None, filein_transform=None, disper_out = None, w_dir=os.getcwd(),
                                                     do_plot=True,disperser=None,fwhm=5.,wid_y=5,
                                                     poly_order=(5,2), mode="TDS", binning = "1x1", blue_lam_limit=None):

    with fits.open(os.path.join(w_dir, filein)) as hdu:
        image=hdu[0].data
        s=image.shape
        head=hdu[0].header

    if mode == "SAO":
        linesfile = os.path.join(os.path.dirname(os.path.realpath(__file__)), "gratings", "SCORPIO",
                                 "calibration_lines.txt")
        ccd=head.get("detector").lower().replace(" ","")
        if ccd == 'eevccd42-40':
            ccdt=""
        else:
            ccdt="_E2V"
        filedis=os.path.join(os.path.dirname(os.path.realpath(__file__)),"gratings", "SCORPIO", "{}{}.txt".format(disperser,ccdt))
        if not os.path.isfile(filedis):
            filedis=os.path.join(os.path.dirname(os.path.realpath(__file__)),"gratings", "SCORPIO", "{}.txt".format(disperser))
            if not os.path.isfile(filedis):
                message("!!!! Cannot find file with calibration lines for grating {} !!!! ".format(disperser),noheader=True)
                return False

        found=re.findall(r"[VPHGR](\d+)", disperser)
        if len(found) == 1:
            ngr=float(found[0])
        else:
            message("!!!! Cannot identify resolution of the disperser {} !!!! ".format(disperser),noheader=True)
            return False
        lev = (4 - int(np.round(ngr / 600.)))
        if lev < 1:
            lev=1

    elif mode=='TDS':
        t=Time(head['Date'])
        add=""
        min_ampl_mul=1.
        if t.datetime.year<=2020 and t.datetime.month<=1:
            add='_addHe'
            min_ampl_mul = 4.
        if disperser=='G':
            min_ampl_mul=min_ampl_mul/1.6
        linesfile = os.path.join(os.path.dirname(os.path.realpath(__file__)), "gratings", "TDS",
                                 "calibration_lines{}.txt".format(add))
        filedis = os.path.join(os.path.dirname(os.path.realpath(__file__)), "gratings", "TDS",
                               "{}.txt".format(disperser))
        if not os.path.isfile(filedis):
            message("!!!! Cannot find file with calibration lines for grating {} !!!! ".format(disperser),
                    noheader=True)
            return False



    blue_lam_limit=5000


    ref_disp=ascii.read(filedis,header_start=1)

    ref_disp['px']=ref_disp['px'].astype(float)/float(binning.split('x')[0])
    ref_disp['A']=ref_disp['A'].astype(float)
    lines_tab = ascii.read(linesfile,header_start=1)


    # preliminary convertion to wavelength
    p_prelim = np.poly1d(np.polyfit(ref_disp["A"], ref_disp["px"], 3))
    c_revers=np.polyfit(ref_disp["px"], ref_disp["A"], 3)
    p_revers = np.poly1d(c_revers)

    tresh = 2.5

    lambda_beg = p_revers(-5)
    lambda_end = p_revers(s[1] + 150)
    index = (lines_tab["Lam"]>lambda_beg) & (lines_tab["Lam"]<lambda_end)
    lines=lines_tab["Lam"][index]
    intensity = lines_tab["Intens"][index]

    #creation full table spectrum in observed coverage
    pos_previus = p_prelim(lines)
    nlines=len(lines)
    message('.... Number of calibration lines in observed spectral coverage: {}'.format(nlines),noheader=True)

    #creation of tabulated comparison spectra
    x = np.arange(s[1])
    neon_tab = np.zeros(s[1],dtype=float)
    for ind,line in enumerate(lines):
        neon_tab+=gaussian(x, [intensity[ind], pos_previus[ind], fwhm/2.35, 0])
    neon_tab = neon_tab / np.max(neon_tab)

    #=== Read observed calibration file and apply corrections
    if filein_transform:
        with fits.open(os.path.join(w_dir, filein_transform)) as hdu_matrix:
            correspondance_matrix = hdu_matrix[0].data
            image = do_distortion(image, correspondance_matrix)
    if fileerr:
        with fits.open(os.path.join(w_dir, fileerr)) as hdu:
            if filein_transform:
                err = do_distortion(hdu[0].data, correspondance_matrix)
            else:
                err=hdu[0].data
    else:
        err = np.sqrt(np.abs(image))


    #==== Compute shift of observed calibration spectrum from tabulated
    M = 300
    win=5
    image.ravel()[np.isnan(image.ravel()) + ~np.isfinite(image.ravel())]=0

    # fits.writeto(os.path.join(w_dir,'test.fits'),image,overwrite=True)
    # import sys
    # sys.exit()
    neon_obs_vec=np.nanmedian(image[int(round(s[0]/2))-win:int(round(s[0]/2))+win+1,:],axis=0)
    dx = derive_vecshift(neon_tab, neon_obs_vec, bin=20, max_ampl=M)



    lambda_beg = p_revers(-5 + dx)
    lambda_end = p_revers(s[1] + dx)


    if mode == 'SAO':
        index = (lines_tab["Lam"] > lambda_beg) & (lines_tab["Lam"] < lambda_end) & (lines_tab["w"] >= lev)
    elif mode == 'TDS':
        if disperser=='R':
            min_ampl=26*min_ampl_mul/2
        else:
            min_ampl=20*min_ampl_mul
        index = (lines_tab["Lam"] > lambda_beg) & (lines_tab["Lam"] < lambda_end) & (lines_tab["Intens"] >= min_ampl)
    else:
        index = (lines_tab["Lam"] > lambda_beg) & (lines_tab["Lam"] < lambda_end)
    lines = lines_tab["Lam"][index]
    nlines=len(lines)

    message('... Shift= {}; expected lambda = {} - {} '.format(-dx, lambda_beg, lambda_end), noheader=True)
    message('.... Number lines used for creating dispersion curve: {}'.format(nlines),noheader=True)

    #pos_lines = np.ndarray((nlines, s[0]),dtype=float)
    ny_fit=int(np.ceil(s[0]/(2*wid_y+1)))
    pos_lines = np.ndarray((nlines, ny_fit), dtype=float)
    pos_previus = p_prelim(lines)
    for ind in range(nlines):
        pos_lines[ind, :]=pos_previus[ind]-dx
    #if keyword_set(fileobj) then pos_sky=poly(make_array(Ny, val=lsky), c_prelim)+dx

    # determine the accuracy of position of lines
    old = pos_lines.copy()
    w=int(round(fwhm))

    fitter=fitting.LevMarLSQFitter()
    g_init=models.Gaussian1D()

    # for row in range(s[0]):
    row_numbers=np.linspace(wid_y, s[0] - wid_y - 1, ny_fit, dtype=int)
    for row_ind, row in enumerate(row_numbers):
        cur_vec=np.nanmedian(image[row-wid_y:row+wid_y+1,:],axis=0)
        cur_vec_err=np.sqrt(np.nanmedian(err[row-wid_y:row+wid_y+1,:]**2,axis=0))
        for ind in range(nlines):
            pos=int(round(pos_lines[ind,row_ind]))
            if lines[ind] > blue_lam_limit:
                wid=w+2
            else:
                wid=w+8
            xborders=[pos-wid,pos+wid]
            if (xborders[0]>=(s[1]-1)) or (xborders[1] <= 0):
                pos_lines[ind,row]=np.nan
                continue
            if xborders[1]>=(s[1]-1):
                xborders[1]=s[1]-1
            if xborders[0]<=0:
                xborders[0]=0

            Nmax = np.argmax(cur_vec[xborders[0]:xborders[1]+1])
            pos = (Nmax+xborders[0])
            if (pos < w ) or pos > (s[1]-1-w):
                pos_lines[ind, row_ind] = np.nan
                continue
            xx=np.linspace(pos-w,pos+w,2*w+1,dtype=float)
            g_init.mean.value=pos
            g_init.stddev.value=fwhm/2.35428
            g_init.amplitude.value=np.nanmax(cur_vec[pos - w:pos + w+1])-np.nanmin(cur_vec[pos - w:pos + w+1])
            g=fitter(g_init,xx,cur_vec[pos - w:pos + w+1]-np.nanmin(cur_vec[pos - w:pos + w+1]),weights=1./cur_vec_err[pos - w:pos + w+1])
            pos_lines[ind, row_ind]=g.mean.value
            # flux = np.nansum(cur_vec[pos - w:pos + w+1])
            # pos_lines[ind, row] = np.nansum(cur_vec[pos - w:pos+w+1] * x[pos - w:pos + w+1]) / flux


    if do_plot:
        plt.figure(figsize=(18,8))
        plt.subplot(111)
        #row_out = int(round(s[0]/2))
        row_out= row_numbers[ny_fit//2]
        # plot line identification
        vector = np.median(image[row_out-wid_y:row_out+wid_y+1,:],axis=0)
        vector=vector/np.nanmax(vector[int(round(0.05 * s[1])):int(round(0.95 * s[1]))])
        vector[vector<0]=0
        vector = np.sqrt(vector) * 1000

        plt.plot(x, vector,'-',color='black',linewidth=0.8)

        for ind in range(nlines):
            pos = pos_lines[ind, ny_fit//2]
            pos_old = old[ind, ny_fit//2]
            if np.isnan(pos) or not np.isfinite(pos):
                continue
            if ind == 0:
                lab1='Found peak'
                lab2='Expected peak'
            else:
                lab1=''
                lab2=''
            plt.plot([pos]*2, [vector[int(round(pos))], vector[int(round(pos))] + 60],'-', color='red',linewidth=0.5, label=lab1)
            plt.plot([pos_old-fwhm]*2, [vector[int(round(pos))], vector[int(round(pos))] + 60],'-', color='blue',linewidth=0.5, label=lab2)
            plt.plot([pos_old+fwhm]*2, [vector[int(round(pos))], vector[int(round(pos))] + 60],'-', color='blue',linewidth=0.5)
            plt.text(pos, vector[int(round(pos))] + 90, "{:.2f}".format(lines[ind]), rotation=90)
        plt.title("Lines used for dispersion curve based on {} calibration image and {} transform matrix".format(filein, filein_transform))
        plt.subplots_adjust(right=0.99, left=0.08, top=0.96, bottom=0.07)
        plt.legend()
        plt.xlim(0,s[1])
        plt.ylim(0, 1300)
        plt.xlabel("X, pixels")
        plt.ylabel("Y, rel. units")
        plt.draw()
        plt.pause(0.001)

    pos_lines_fit=np.zeros((nlines,s[0]),dtype=float)#pos_lines.copy()
    y_pos=np.arange(s[0])
    if poly_order[1] > 0:
        for ind in range(nlines):
            finite_mask=np.isfinite(pos_lines[ind, :])*(~np.isnan(pos_lines[ind, :]))
            if np.sum(finite_mask) >= ny_fit/10:
                fit_func,_=goodpoly(row_numbers[finite_mask],pos_lines[ind,finite_mask],yerr=None,order=poly_order[1])
                pos_lines_fit[ind,:]=fit_func(y_pos)
            else:
                pos_lines_fit[ind, :]=np.nan
    ##### А у Моисеева здесь идет сглаживание позиции по Y.
    disper_par=np.ndarray((poly_order[0]+2,s[0]),dtype=float)
    lines_calc = np.zeros((nlines,s[0]),dtype=float)
    lines_measured = np.zeros_like(pos_lines)
    for row in range(s[0]):
        finite_mask = np.isfinite(pos_lines_fit[:, row])*(~np.isnan(pos_lines_fit[:, row]))
        fit_func,fit_params = goodpoly(pos_lines_fit[finite_mask, row], lines[finite_mask], yerr=None, order=poly_order[0])
        disper_par[0: poly_order[0]+1, row]=fit_params
        lines_calc[:, row] = fit_func(pos_lines_fit[:, row])
        if row in row_numbers:
            lines_measured[:,list(row_numbers).index(row)]=fit_func(pos_lines[:,list(row_numbers).index(row)])

        disper_par[poly_order[0] + 1, row] = np.nanstd(lines_calc[:, row] - lines)


    fits.writeto(os.path.join(w_dir, disper_out), np.float32(disper_par),overwrite=True)

    if do_plot:
        rms_mean = np.sum(disper_par[poly_order[0] + 1, :]) / s[0]

        currow=int(round(s[0] / 2))
        xx = np.arange(s[1])
        dlam=0.*xx

        for i in range(1,poly_order[0]):
            dlam+=disper_par[i, currow]*i*xx**(i-1)

        plt.figure(figsize=(8, 5))
        plt.subplot(121)
        plt.plot(xx, dlam, '-')
        plt.xlabel("X, pix")
        plt.ylabel(r"d$\lambda$/dX, $\AA$/px")
        plt.subplot(122)
        plt.hist(disper_par[poly_order[0] + 1, :], bins=20, range=[0,2.5])
        plt.xlabel("Error, px")
        plt.ylabel("Nrows")
        plt.tight_layout()


        plt.figure(figsize=(12, 22))
        plt.subplot(111)
        dy = 3.
        for ind in range(nlines):
            rms = np.nanstd(-lines_measured[ind, :] + lines[ind])
            mean = np.nanmean(-lines_measured[ind, :] + lines[ind])
            p = np.poly1d([c_revers[1], c_revers[2] * 2])
            d_lambda = p(lines[ind])
            plt.plot(row_numbers, (-lines_measured[ind, :] + lines[ind] - mean) / d_lambda + dy * ind, '.', color='black',
                     markersize=6)
            plt.plot([0, s[0]], np.array([dy, dy]) * ind, '-', linewidth=1)
            plt.text(s[0] * 1.02, dy * ind, "{:.2f}    {:.2f}±{:.2f}".format(lines[ind], mean, rms))
        plt.text(s[0] * 1.03, dy * (nlines + 0.5), r"$\lambda, \AA$         $\Delta\lambda, \AA$")
        plt.xlim(0, s[0])
        plt.ylim(-dy, nlines * dy)
        plt.subplots_adjust(right=0.85, left=0.08, top=0.96, bottom=0.05)
        plt.xlabel("Y, pixels")
        plt.ylabel(r"$\Delta$X, pixels")
        plt.title(
            'Dispersion curve for disperser {} at binning {}. Mean r.m.s={:.3f}px'.format(disperser, binning, rms_mean))
        plt.draw()
        plt.pause(0.001)

    return True



def linearization(filein_transform=None, fileout_transform= None, dispcurve = None, file_shifts=None,file_atmdisp=None,
                  w_dir=os.getcwd(), lin_parameters={'start':'auto', 'fin': 'auto', 'dlam': 'auto'}, file_number=None):
    with fits.open(os.path.join(w_dir,filein_transform)) as hdu:
            correspondance_matrix=hdu[0].data
            s=correspondance_matrix.shape
    with fits.open(os.path.join(w_dir,dispcurve)) as hdu:
        dispers_param=hdu[0].data[:-1,:]

    if file_shifts is not None:
        shift_tab = ascii.read(os.path.join(w_dir, file_shifts))
        fnum=shift_tab['filenum']
        shifts=np.ndarray((2,len(fnum)))
        shifts[0,:]=shift_tab['dy']
        shifts[1, :] = shift_tab['dx']
    else:
        shifts=np.zeros((2,1),dtype=float)
    if file_atmdisp is not None:
        atmdisp_tab = ascii.read(os.path.join(w_dir, file_atmdisp))
        atmdisp_cf=np.array([atmdisp_tab[key][0] for key in atmdisp_tab.keys()])
    else:
        atmdisp_cf=None

    nfiles = shifts.shape[1]
    row_central=int(round(s[1]/2))
    col_central = int(round(s[2] / 2))
    # dispers_param.ravel()[np.isnan(dispers_param.ravel()) + ~np.isfinite(dispers_param.ravel())] = 0
    dispers_param[:,0]=dispers_param[:,1]
    dispers_param[:, -1] = dispers_param[:, -2]
    p = np.poly1d(np.flip(dispers_param[:,row_central]))

    if lin_parameters['dlam'] == 'auto':
        dlam = round(dispers_param[1,row_central],2)
    else:
        dlam = round(lin_parameters['dlam'],2)
    if lin_parameters['start'] == 'auto':
        lmin=round(p(0),2)
    else:
        lmin =  round(lin_parameters['start'],2)
    if lin_parameters['fin'] == 'auto':
        lfin=p(s[2]-1)
        nxout=int(round((lfin-lmin)/dlam))
    else:
        nxout=int(round(lin_parameters['fin'] -lmin)/dlam)

    xx_out = np.arange(nxout)
    lam_out = xx_out*dlam+lmin
    xx_original = np.arange(s[2])
    correspondance_matrix_new=np.ndarray((nfiles,2,s[1],nxout),dtype=float)
    if atmdisp_cf is not None:
        ap = np.poly1d(np.flip(atmdisp_cf))
        atmdisp_shift=ap(xx_original)
    else:
        atmdisp_shift=np.zeros_like(xx_original)


    for row in range(s[1]):
        p = np.poly1d(np.flip(dispers_param[:, row]))
        lam_original = p(xx_original)
        p_inverse = np.poly1d(np.polyfit(lam_original, xx_original, len(dispers_param[:, row]) - 1))
        xx_corespond = p_inverse(lam_out)
        for file_ind in range(nfiles):
            correspondance_matrix_new[file_ind, 0,row,:]=np.interp(xx_corespond,xx_original,correspondance_matrix[0,row,:]+atmdisp_shift[:])+shifts[0,file_ind]
            correspondance_matrix_new[file_ind, 1, row, :] = (xx_corespond-xx_out)+np.interp(xx_corespond, xx_original, correspondance_matrix[1, row, :]+shifts[1,file_ind])
    header = fits.Header()
    header['NAXIS1'] = nxout
    header['CTYPE1']='AWAV'
    header['CRPIX1'] = 1
    header['CDELT1'] = round(dlam*100)/100.
    header['CRVAL1'] = lam_out[0]

    fits.writeto(os.path.join(w_dir,fileout_transform),np.float32(correspondance_matrix_new),header=header,overwrite=True)
    return True



def correct_geometry(filein, file_transform_matrix, file_err=None, fileout=None, fileout_err=None, w_dir=os.getcwd(), file_number=None):
    with fits.open(os.path.join(w_dir, file_transform_matrix)) as hdu:
        correspondance_matrix=hdu[0].data
        head_matrix=hdu[0].header
        if len(correspondance_matrix.shape) == 3:
            single_matrix=True
            add_message = ""
            correspondance_matrix=correspondance_matrix.reshape((1,correspondance_matrix.shape[0],correspondance_matrix.shape[1],correspondance_matrix.shape[2]))
        else:
            add_message="and linearization "
            single_matrix = False
    for ind,f in enumerate(filein):
        if not single_matrix:
            matrix_ind = ind
        else:
            matrix_ind = 0
        with fits.open(os.path.join(w_dir, f)) as hdu:
            img_process=hdu[0].data
            image = do_distortion(img_process, correspondance_matrix[matrix_ind,:,:,:])
            header = hdu[0].header
            if check_param(head_matrix,["ctype1"],'awav'):
                header['CTYPE1']=head_matrix["CTYPE1"]
                header['CRVAL1'] = head_matrix["CRVAL1"]
                header['CDELT1'] = head_matrix["CDELT1"]
                header['CRPIX1'] = head_matrix["CRPIX1"]
            header.add_history("LONGRED: Geometry correction {}based on transform. {}".format(add_message,file_transform_matrix))
            fits.writeto(os.path.join(w_dir,fileout[ind]),np.float32(image),header=header,output_verify='fix',overwrite=True)
        if file_err is not None:
            if file_err[ind]:
                with fits.open(os.path.join(w_dir, file_err[ind])) as hdu:
                    image = do_distortion(hdu[0].data, correspondance_matrix[matrix_ind,:,:,:])
                    header = hdu[0].header
                    if check_param(head_matrix, ["ctype1"], 'awav'):
                        header['CTYPE1'] = head_matrix["CTYPE1"]
                        header['CRVAL1'] = head_matrix["CRVAL1"]
                        header['CDELT1'] = head_matrix["CDELT1"]
                        header['CRPIX1'] = head_matrix["CRPIX1"]
                    header.add_history("LONGRED: Geometry correction {}based on transform. {}".format(add_message, file_transform_matrix))
                    fits.writeto(os.path.join(w_dir, fileout_err[ind]), np.float32(image), header=header, output_verify='fix',
                                 overwrite=True)
    return True


def remove_sky(filein, file_err=None, fileout=None, fileout_err=None, w_dir=os.getcwd(), file_number=None,
               cr_clean=True,minorder=None,medorder=None,
               order=3, regions=(), file_transform=None, file_shifts=None, do_fft=True, run_test=None,
               file_sunsky=None):
    if filein is None or len(filein) < 1:
        return False
    if file_err is None or len(file_err) == 0: #run_test is not None or
        file_err=None
        fileout_err=None
    else:
        if fileout_err is None or len(fileout_err) == 0:
            file_err = None
            fileout_err = None
    bin = 20
    add_bias=0.
    if do_fft:
        add_bias=1000.
    if run_test:
        filein = [filein[0]]
    if order is None:
        fit_method = "template"
    else:
        fit_method = " polynomial with order={}".format(order)
    if do_fft:
        fit_method = "{} in fft space".format(fit_method)
        if order is None and file_sunsky is not None:
            fit_method = "{} with sunsky".format(fit_method)
    if file_transform:
        with fits.open(os.path.join(w_dir, file_transform)) as hdu:
            correspondance_matrix = hdu[0].data
    if file_shifts:
        shift_tab = ascii.read(os.path.join(w_dir, file_shifts))
        fnum = shift_tab['filenum']
        shifts = np.array(shift_tab['dy'])
    else:
        shifts = np.zeros(len(filein), dtype=float)

    if do_fft and file_sunsky:
        with fits.open(os.path.join(w_dir, file_sunsky)) as hdu:
            sunsky = hdu[0].data
            if file_transform:
                sunsky = do_distortion(sunsky, correspondance_matrix)
            sunsky+=add_bias
    for ind, f in enumerate(filein):
        with fits.open(os.path.join(w_dir, f)) as hdu:
            image = hdu[0].data
            header = hdu[0].header
        if file_err is not None:
            with fits.open(os.path.join(w_dir, file_err[ind])) as hdu:
                err_image_ini = hdu[0].data
                err_header=hdu[0].header
                err_image=err_image_ini.copy()
                rec = np.isnan(err_image.ravel())
                err_image.ravel()[rec]=np.nanmax(err_image)
                if file_transform:
                    err_image = do_distortion(err_image, correspondance_matrix)
        if cr_clean and file_err is not None:
            image_analyse,_ = lacosmic(image, 1.5, 5.,3.,error=err_image_ini)
        else:
            image_analyse = image.copy()
        if file_transform:
            image_analyse = do_distortion(image_analyse, correspondance_matrix)

        s = image_analyse.shape
        xx = np.arange(s[1])
        yy = np.arange(s[0])
        reg = np.zeros(s[0], dtype=bool)
        reg_to_header = ""
        for r in regions:
            r[0]=r[0]-shifts[ind]
            r[1] = r[1]-shifts[ind]
            if (r[0] < (s[0] - 1)) & (r[1] > 0):
                r0 = np.max([int(round(r[0])), 0])
                r1 = np.min([int(round(r[1])), s[0] - 1])
                reg[r0:r1] = True
                reg_to_header = "{}{:d}-{:d}; ".format(reg_to_header, r0, r1)

        spec_analyse = image_analyse.copy()+add_bias
        sky = np.zeros_like(spec_analyse)
        if fileout_err is not None:
            sigma_sky = np.zeros_like(spec_analyse)


        if order is not None:
            err_image.ravel()[np.isnan(err_image.ravel()) + (err_image.ravel() == 0)] = np.max(err_image)
            if do_fft:
                spec_analyse = np.fft.fft(spec_analyse, axis=1)
                sky_fft=np.zeros_like(spec_analyse)
                sigma_sky_fft=np.ndarray((np.sum(reg),s[1]),dtype=np.complex)
            npoints_medorder = 0
            npoints_minorder = 0
            npoints_maxorder = image_analyse.shape[1]
            if file_err is not None:
                rec=np.isnan(err_image[reg,:].ravel())+~np.isfinite(err_image[reg,:].ravel())+(err_image[reg,:].ravel() == 0)
                if np.sum(rec) > 0:
                    err_image[reg,:].ravel()[rec]=np.nanmax(err_image[reg,:].ravel()[rec])
                image_analyse[reg,:].ravel()[np.isnan(image_analyse[reg,:].ravel()) + ~np.isfinite(image_analyse[reg,:].ravel())] = 0
                snrat=np.nanmedian(abs(image_analyse[reg, :]) / err_image[reg, :],axis=0)
                curorder=np.ones_like(snrat,dtype=int)
                # order_med=order-1
                # if order_med<np.min([2,order]):
                #     order_med=np.min([2,order])
                # order_min=np.min([2,order])
                curorder[:]=order
                if medorder is not None:
                    npoints_medorder=np.sum(snrat<medorder[1])
                    npoints_maxorder-=npoints_medorder
                    curorder[snrat<medorder[1]] = medorder[0]
                if minorder is not None:
                    npoints_minorder = np.sum(snrat < minorder[1])
                    if medorder is not None:
                        npoints_medorder -= npoints_minorder
                    else:
                        npoints_maxorder -= npoints_minorder
                    curorder[snrat < minorder[1]] = minorder[0]
                if do_fft:
                    if npoints_medorder>0:
                        sky_fft_med=np.zeros_like(sky_fft)
                        sigma_sky_fft_med = np.zeros_like(sigma_sky_fft)
                    if npoints_minorder > 0:
                        sky_fft_min = np.zeros_like(sky_fft)
                        sigma_sky_fft_min = np.zeros_like(sigma_sky_fft)
                cur_error=err_image[reg,:]
            else:
                cur_error=None
                curorder=np.array(xx.shape[0],dtype=int)
                curorder[:]=order
            for x in xx:
                if do_fft:
                    p_re, _, data_re = simplepoly(yy[reg], spec_analyse[reg, x].real,yerr=None, order=order,get_masked_data=True)
                    p_im, _, data_im = simplepoly(yy[reg], spec_analyse[reg, x].imag,yerr=None, order=order,get_masked_data=True)
                    sky_fft[:, x] = 1j * p_im(yy)
                    sky_fft[:, x] += p_re(yy)
                    if file_err is not None:
                        if npoints_medorder > 0:
                            p_re_med, _, data_re_med = simplepoly(yy[reg], spec_analyse[reg, x].real, yerr=None, order=medorder[0],
                                                              get_masked_data=True)
                            p_im_med, _, data_im_med = simplepoly(yy[reg], spec_analyse[reg, x].imag, yerr=None, order=medorder[0],
                                                              get_masked_data=True)
                            sky_fft_med[:, x] = 1j * p_im_med(yy)
                            sky_fft_med[:, x] += p_re_med(yy)
                        if npoints_minorder > 0:
                            p_re_min, _, data_re_min = simplepoly(yy[reg], spec_analyse[reg, x].real, yerr=None, order=minorder[0],
                                                              get_masked_data=True)
                            p_im_min, _, data_im_min = simplepoly(yy[reg], spec_analyse[reg, x].imag, yerr=None, order=minorder[0],
                                                              get_masked_data=True)

                            sky_fft_min[:, x] = 1j * p_im_min(yy)
                            sky_fft_min[:, x] += p_re_min(yy)


                    if fileout_err is not None:
                        sigma_sky_fft[:,x] = 1j*np.ma.filled(data_im-np.ma.masked_array(p_im(yy[reg]),mask=data_im.mask),0)+\
                                             np.ma.filled(data_re-np.ma.masked_array(p_re(yy[reg]),mask=data_re.mask),0)
                        if npoints_medorder > 0:
                            sigma_sky_fft_med[:, x] = 1j * np.ma.filled(
                            data_im_med - np.ma.masked_array(p_im_med(yy[reg]), mask=data_im_med.mask), 0) + \
                                              np.ma.filled(
                                                  data_re_med - np.ma.masked_array(p_re_med(yy[reg]), mask=data_re_med.mask), 0)
                        if npoints_minorder > 0:
                            sigma_sky_fft_min[:, x] = 1j * np.ma.filled(
                                data_im_min - np.ma.masked_array(p_im_min(yy[reg]), mask=data_im_min.mask), 0) + \
                                                  np.ma.filled(
                                                      data_re_min - np.ma.masked_array(p_re_min(yy[reg]), mask=data_re_min.mask), 0)
                else:
                    p, _, data = simplepoly(yy[reg], spec_analyse[reg, x],yerr=cur_error[:,x], order=curorder[x],get_masked_data=True)
                    sky[:, x] = p(yy)
                    if fileout_err is not None:
                        sigma_sky[:,x]=np.ma.std(data-np.ma.masked_array(p(yy[reg]),mask=data.mask))


            if do_fft:
                sky = np.abs(np.fft.ifft(sky_fft, axis=1))
                if file_err is not None:
                    if npoints_medorder>0:
                        rec=(curorder == medorder[0])
                        sky[:,rec]=np.abs(np.fft.ifft(sky_fft_med, axis=1))[:, rec]
                    if npoints_minorder > 0:
                        rec = (curorder == minorder[0])
                        sky[:, rec] = np.abs(np.fft.ifft(sky_fft_min, axis=1))[:, rec]
                if fileout_err is not None:
                    sigma_sky[:, :] = np.nanstd(np.abs(np.fft.ifft(sigma_sky_fft,axis=1)),axis=0)[None, :]
                    if npoints_medorder>0:
                        rec=(curorder == medorder[0])
                        sigma_sky[:, rec] = np.nanstd(np.abs(np.fft.ifft(sigma_sky_fft_med, axis=1)),axis=0)[None,rec]
                    if npoints_minorder>0:
                        rec = (curorder == minorder[0])
                        sigma_sky[:, rec] = np.nanstd(np.abs(np.fft.ifft(sigma_sky_fft_min, axis=1)),axis=0)[None,rec]
            if minorder is not None or medorder is not None:
                message("...Number of pixels along X fitted by minorder: {}; medorder: {}; maxorder: {}".format(npoints_minorder,
                                                                    npoints_medorder, npoints_maxorder),noheader=True)
        else:
            xx_bin = np.linspace(0, s[1], s[1] * bin, dtype=float) + 0.5
            xx = np.linspace(0, s[1], s[1], dtype=float)
            sky_cut = signal.resample_poly(signal.savgol_filter(spec_analyse[reg, :],7,2,axis=1), bin, 1, axis=1)
            if do_fft:
                sky_over = np.zeros((s[0], s[1] * bin))
                if fileout_err is not None:
                    sigma_over = np.zeros((s[0], s[1] * bin),dtype=float)
                sky_cut_fft=np.fft.fft(sky_cut,axis=1)
                if file_sunsky:
                    sunsky_cut_fft = np.fft.fft(signal.resample_poly(signal.savgol_filter(sunsky[reg, :],7,2,axis=1),bin,1,axis=1),axis=1)
                    sky_rat_i=np.median(sky_cut_fft.imag/sunsky_cut_fft.imag,axis=0)
                    rec=np.argwhere(np.isnan(sky_rat_i) + ~np.isfinite(sky_rat_i))
                    sky_rat_i[rec]=0
                    sky_rat_r = np.median(sky_cut_fft.real / sunsky_cut_fft.real, axis=0)
                    rec = np.argwhere(np.isnan(sky_rat_r) + ~np.isfinite(sky_rat_r))
                    sky_rat_r[rec] = 0
                    sunsky_smo=signal.resample_poly(signal.savgol_filter(sunsky,7,2,axis=1),bin,1,axis=1)
                    sunsky_fft=np.fft.fft(sunsky_smo,axis=1)
                    sky_fft_over_i=np.zeros((s[0],s[1]*bin),dtype=float)
                    sky_fft_over_r = np.ndarray((s[0], s[1] * bin), dtype=float)
                    sky_fft_over_i[:,:]=sunsky_fft.imag[:,:]*sky_rat_i[None,:]
                    sky_fft_over_r[:,:]=sunsky_fft.real[:,:]*sky_rat_r[None,:]
                    rec = np.argwhere(np.isnan(sky_fft_over_r.ravel()) + ~np.isfinite(sky_fft_over_r.ravel()) +
                                      np.abs(sky_fft_over_r.ravel())<1e-10)
                    sky_fft_over_r.ravel()[rec]=0
                    rec = np.argwhere(np.isnan(sky_fft_over_i.ravel()) + ~np.isfinite(sky_fft_over_i.ravel()) +
                                      np.abs(sky_fft_over_i.ravel()) < 1e-10)
                    sky_fft_over_i.ravel()[rec] = 0
                    sky_over=np.real(np.fft.ifft(sky_fft_over_r+1j*sky_fft_over_i,axis=1))

                    if fileout_err is not None:
                        sky_rat_err_fft = np.nanstd(sky_cut_fft.real/sunsky_cut_fft.real, axis=0) + 1j * np.nanstd(sky_cut_fft.imag/sunsky_cut_fft.imag,
                                                                                           axis=0)
                        rec = np.argwhere(np.isnan(sky_rat_err_fft.real) + ~np.isfinite(sky_rat_err_fft.real))
                        sky_rat_err_fft.real[rec] = 0
                        rec = np.argwhere(np.isnan(sky_rat_err_fft.imag) + ~np.isfinite(sky_rat_err_fft.imag))
                        sky_rat_err_fft.imag[rec] = 0
                        sky_sigma_fft_over_i = np.zeros((s[0], s[1] * bin), dtype=float)
                        sky_sigma_fft_over_r = np.ndarray((s[0], s[1] * bin), dtype=float)
                        sky_sigma_fft_over_i[:, :] = sunsky_fft.imag[:, :] * sky_rat_err_fft.imag[None, :]
                        sky_sigma_fft_over_r[:, :] = sunsky_fft.real[:, :] * sky_rat_err_fft.real[None, :]
                        rec = np.argwhere(np.isnan(sky_sigma_fft_over_r.ravel()) + ~np.isfinite(sky_sigma_fft_over_r.ravel()) +
                                          np.abs(sky_sigma_fft_over_r.ravel()) < 1e-10)
                        sky_sigma_fft_over_r.ravel()[rec] = 0
                        rec = np.argwhere(np.isnan(sky_sigma_fft_over_i.ravel()) + ~np.isfinite(sky_sigma_fft_over_i.ravel()) +
                                          np.abs(sky_sigma_fft_over_i.ravel()) < 1e-10)
                        sky_sigma_fft_over_i.ravel()[rec] = 0
                        sigma_over = np.abs(np.fft.ifft(sky_sigma_fft_over_r + 1j * sky_sigma_fft_over_i, axis=1))
                    for i in range(s[0]):
                        t, c, k = interpolate.splrep(xx_bin, sky_over[i,:])
                        spline = interpolate.BSpline(t, c, k, extrapolate=False)
                        sky[i,:] = spline(xx)
                        if fileout_err is not None:
                            t, c, k = interpolate.splrep(xx_bin, sigma_over[i,:])
                            spline = interpolate.BSpline(t, c, k, extrapolate=False)
                            sigma_sky[i, :] = spline(xx)
                else:
                    sky_model_vec_fft=np.median(sky_cut_fft.real,axis=0)+1j*np.median(sky_cut_fft.imag,axis=0)
                    rec = np.argwhere(np.isnan(sky_model_vec_fft.real) + ~np.isfinite(sky_model_vec_fft.real) +
                                      np.abs(sky_model_vec_fft.real) < 1e-10)
                    sky_model_vec_fft.real[rec]=0
                    rec = np.argwhere(np.isnan(sky_model_vec_fft.imag) + ~np.isfinite(sky_model_vec_fft.imag) +
                                      np.abs(sky_model_vec_fft.imag) < 1e-10)
                    sky_model_vec_fft.imag[rec] = 0
                    sky_over[:,:]=np.abs(np.fft.ifft(sky_model_vec_fft))[None,:]
                    t, c, k = interpolate.splrep(xx_bin, sky_over[0, :])
                    spline = interpolate.BSpline(t, c, k, extrapolate=False)
                    sky[:, :] = spline(xx)[None,:]
                    if fileout_err is not None:
                        t, c, k = interpolate.splrep(xx_bin, np.nanstd(sky_cut, axis=0))
                        spline = interpolate.BSpline(t, c, k, extrapolate=False)
                        sigma_sky[:, :] = spline(xx)[None, :]

            else:
                t, c, k = interpolate.splrep(xx_bin, np.median(sky_cut,axis=0))
                spline = interpolate.BSpline(t, c, k, extrapolate=False)
                sky[:, :] = spline(xx)[None,:]
                if fileout_err is not None:
                    t, c, k = interpolate.splrep(xx_bin, np.nanstd(sky_cut, axis=0))
                    spline = interpolate.BSpline(t, c, k, extrapolate=False)
                    sigma_sky[:, :] = spline(xx)[None,:]
        sky-=add_bias
        if run_test:
            fits.writeto(os.path.join(w_dir, run_test), np.float32(image_analyse - sky), header=header,
                         overwrite=True)
            return True
        if file_transform:
            sky = do_distortion(sky, correspondance_matrix,reverse=True)
        image = image - sky
        header.add_history("LONGRED: Sky fitted by {}".format(fit_method))
        header.add_history("... Regions for sky rem.: {}".format(reg_to_header))
        fits.writeto(os.path.join(w_dir, fileout[ind]), np.float32(image), header=header, overwrite=True)
        if file_err is not None and fileout_err is not None:
            if file_transform:
                sigma_sky = do_distortion(sigma_sky, correspondance_matrix,reverse=True)
            err_header.add_history("LONGRED: Sky fitted by {}".format(fit_method))
            err_header.add_history("... Regions for sky rem.: {}".format(reg_to_header))
            fits.writeto(os.path.join(w_dir, fileout_err[ind]), np.float32(np.sqrt(err_image_ini**2+sigma_sky**2)),
                            header=err_header, overwrite=True)
    return True


def exp_combine(filein, file_err=None, fileout=None, fileout_err=None, w_dir=os.getcwd(),
                file_number=None,  sigma=5, mode="Sigma"):
    texp=[]
    for ind,f in enumerate(filein):
        with fits.open(os.path.join(w_dir,f)) as hdu:
            image=hdu[0].data
            header=hdu[0].header
        if ind ==0:
            s=image.shape
            images=np.ndarray((len(filein),s[0],s[1]),dtype=float)
            if file_err is not None:
                err_sqaured=np.zeros((s[0],s[1]),dtype=float)
        if file_err is not None:
            with fits.open(os.path.join(w_dir,file_err[ind])) as hdu:
                header_err=hdu[0].header
                err_sqaured+=hdu[0].data**2
        images[ind,:,:]=image
        texp.append(float(header.get("exptime")))
    ### Normilize for texp differences
    texp=np.array(texp,dtype=float)
    texp_med=np.median(texp)
    images=images*(texp_med/texp[:,None,None])
    texp=np.sum(texp)
    if mode=="Sigma":
        tot_image=combine_sigma(images,cr_treshold=sigma)
        header.add_history("LONGRED: Combined {} exposures using sigma-clip with tr={}".format(len(filein),sigma))
    elif mode=="Median":
        tot_image=np.nanmedian(images,axis=0)*len(filein)
        header.add_history("LONGRED: {} exposures combined using median filter".format(len(filein)))
    elif mode=='Pairs':
        tot_image = combine_pairs(images, cr_treshold=sigma,ncycles=2)
        header.add_history("LONGRED: Combined {} exposures using pairs compar. with tr={}".format(len(filein), sigma))
    header['EXPTIME']=texp
    fits.writeto(os.path.join(w_dir, fileout), np.float32(tot_image), header=header, overwrite=True)
    if file_err is not None:
        header_err.add_history("LONGRED: Combined {} exposures using sigma-clip with tr={}".format(len(filein), sigma))
        header_err['EXPTIME'] = texp
        fits.writeto(os.path.join(w_dir, fileout_err), np.float32(np.sqrt(err_sqaured)), header=header_err, overwrite=True)
    return True

def calc_dqe(filein, file_err=None, fileout='senscuve.fits',standards_dir="/Users/mors/Science/standards/data/",z_star=None,
                dqe_figname='dqe.pdf',mode="TDS",gain=None,extract_window=(None,None),w_dir=os.getcwd(), do_plot=True,star_name=None, smo_win=None):
    with fits.open(os.path.join(w_dir, filein)) as hdu:
        image = hdu[0].data
        header = hdu[0].header
        if extract_window[0] is None:
           y0=0+int(round(image.shape[0]/20))
        else:
            y0=np.max([0,extract_window[0]])
        if extract_window[1] is None:
            y1 = image.shape[0]-1-int(round(image.shape[0]/20))
        else:
            y1 = np.min([image.shape[0]-1, extract_window[1]])
        obs_flux=np.nansum(image[y0:y1+1,:],axis=0)
    if star_name is None:
        star_name=header.get('object')
        if not star_name:
            message("Cannot identify name of the star!")
            return False
    filetab=os.path.join(standards_dir,"m{}.dat".format(re.sub("[\._\-+\s]+", "", star_name.casefold())))
    if not os.path.isfile(filetab):
        message("Cannot find SED for star {}".format(star_name))
        print(filetab)
        return False

    if gain is None:
        gain=header.get('gain')
        if not gain:
            gain=1
    texp=header.get('exptime')
    dlam=header.get('cdelt1')
    if z_star is None:
        z_star=header.get('z')
    if not z_star:
        z_star=20.
        message("... Cannot derive Z(star). Use Z=20 degree.", noheader=True)

    sedtab=ascii.read(filetab)
    wlscale=dlam*(np.arange(image.shape[1])-header.get('crpix1')+1)+header.get('crval1')
    sed_wl=np.array(sedtab['col1'])
    sed_mag = np.array(sedtab['col2'])

    if sed_wl[-1] < wlscale[-1]:
        R = (sed_wl > 8000)
        x = sed_wl[R]
        y = sed_mag[R]#np.log10(sed_flux[R])
        if len(x)>=2:
            f,cf = goodpoly(x, y, order=1)
            N_new = int(np.round((15000. - sed_wl[-1]) / 10.))
            x_new = np.arange(N_new) * 10 + sed_wl[-1] + 10
            y_new = f(x_new)
            sed_wl = np.append(sed_wl, x_new)
            sed_mag = np.append(sed_mag, y_new)

    if mode=="SAO":
        S = 2.51e5 #square of 6-m mirror
    elif mode=="TDS":
        S=42706#49087*0.83
    N_quanta_sed = 10**(-0.4*(sed_mag +calc_ext(sed_wl,z_star)))*948 * S/gain*(5500/sed_wl)**2
    N_quanta_obs = obs_flux/ texp / dlam

    index = ((sed_wl >= (wlscale[0] - 100)) * (sed_wl <= (wlscale[-1] + 100)))
    sed_wl = sed_wl[index]
    N_quanta_sed_original_scale = N_quanta_sed[index]
    f = interpolate.interp1d(sed_wl, N_quanta_sed_original_scale)
    N_quanta_sed = f(wlscale)

    dx=derive_vecshift(signal.savgol_filter(N_quanta_sed-signal.savgol_filter(N_quanta_sed,51,3),5,3),
                       signal.savgol_filter(N_quanta_obs-signal.savgol_filter(N_quanta_obs,51,3),5,3),
                       max_ampl=sed_mag.shape[0]/10.)
    if dx>30:
        message("... Large shift between obs. star and SED. Check this!", noheader=True)
    N_quanta_obs=ndimage.shift(N_quanta_obs,dx,order=1,mode='nearest')


    if mode != "TDS":
        bad_lambda = [3930, 3960, 4130, 4340, 4680, 4861, 5420, 5740, 5770, 6300, 6562,
                  6860, 6900, 6950,  7650, 7600,#7210,
                  7650, 7700, 7750]
    else:
        bad_lambda = [3960,4130, 4340, 4861, 5420, 6562,
                      6860, 6900, 6950, 7210, 7650, 7600,
                      7650, 7700, 7750]
    mask = np.ones(wlscale.shape[0], dtype=bool)
    index = ~(np.isnan(N_quanta_sed) | np.isnan(N_quanta_obs) | (N_quanta_sed == 0) | ~np.isfinite(N_quanta_sed) | ~np.isfinite(N_quanta_obs))
    mask[index]=False
    dqe_all=np.ma.array(N_quanta_obs/N_quanta_sed,mask=mask.copy())
    if mode=="SAO" or mode=="TDS":
        if mode == "SAO": dwin=80
        else: dwin=50
        for blam in bad_lambda:
            index = (abs(wlscale - blam) < dwin)
            mask[index]=True
    dqe=np.ma.array(N_quanta_obs/N_quanta_sed,mask=mask)

    def get_rid_of_nans(y):
        x=np.arange(y.shape[0])
        p=np.poly1d(np.ma.polyfit(x,y,6))
        y=np.ma.filled(y,fill_value=np.nan)
        y[np.isnan(y)]=p(x[np.isnan(y)])
        return y


    if smo_win is None:
        smo_win=int(np.round(obs_flux.shape[0] / 40)) * 2 + 1
    #dqe_smo=signal.savgol_filter(get_rid_of_nans(dqe),smo_win,2)
    # dqe_smo = signal.savgol_filter(signal.savgol_filter(signal.savgol_filter(get_rid_of_nans(dqe), 31, 2),51,2),smo_win,2)
    from statsmodels.nonparametric.smoothers_lowess import lowess
    dqe_smo = lowess(np.ma.filled(dqe,fill_value=np.nan),wlscale,frac=smo_win/obs_flux.shape[0],it=1, is_sorted=True)
    x_smooth = dqe_smo[:, 0]
    y_smooth = dqe_smo[:, 1]
    t, c, k = interpolate.splrep(x_smooth, y_smooth)
    spline = interpolate.BSpline(t, c, k, extrapolate=False)
    dqe_smo = spline(wlscale)
    senscurve= 3.39e-11 * gain / (9.48 * dqe_smo * S)
    #senscurve=np.zeros(wlscale.shape[0],dtype=float)
    # senscurve[0,:]=1./(dqe_smo*2.8e11* S) # ergs
    # senscurve[1,:]=senscurve[0,:]**wlscale**2*3.335e7 # mJy
    header['Z']=z_star
    fits.writeto(os.path.join(w_dir,fileout), senscurve,header=header, overwrite=True)

    plt.figure(figsize=(6,8))
    plt.subplot(211)
    plt.plot(wlscale, N_quanta_obs/np.nanmax(N_quanta_obs),label='Observed')
    plt.plot(sed_wl,N_quanta_sed_original_scale/np.nanmax(N_quanta_sed_original_scale[(sed_wl>wlscale[0])*(sed_wl<wlscale[-1])]),label='Reference')
    plt.legend()
    plt.title("Observed and reference spectra of star {}".format(star_name))
    plt.ylabel(r"Relative flux")
    plt.xlabel(r"$\lambda, \AA$")
    plt.xlim(wlscale[0], wlscale[-1])
    plt.ylim(0, 1.1)
    plt.subplot(212)
    plt.plot(wlscale,dqe_smo*100.)
    plt.plot(wlscale[mask==False],np.ma.filled(dqe_all*100,fill_value=np.nan)[mask==False],'.',color='k',markersize=1.0)
    plt.plot(wlscale[mask == True], np.ma.filled(dqe_all * 100, fill_value=np.nan)[mask == True], '.', color='gray',
             markersize=0.7)
    plt.xlim(wlscale[0],wlscale[-1])
    plt.ylim(0,np.nanmax(dqe_smo*100.)*1.1)
    plt.title("DQE derived with star {}".format(star_name))
    plt.ylabel(r"DQE, %")
    plt.xlabel(r"$\lambda, \AA$")
    plt.tight_layout()
    plt.savefig(os.path.join(w_dir,dqe_figname),dpi=150)
    if do_plot:
        plt.draw()
        plt.pause(0.0001)

    return True

def flux_calib(filein, file_err=None, file_sens=None, fileout=None, z_obj=None, cut_wl=(None,None),
               ypos=None, ydelt=None, fileout_err=None, w_dir=os.getcwd(),files_for_cut=None,
               files_for_cut_out=None, ycut=[0,None]):
    if not os.path.isfile(os.path.join(w_dir,filein)):
        message("... Cannot find input file!",noheader=True)
        return False
    if file_sens is None or not os.path.isfile(os.path.join(w_dir,file_sens)):
        message("... Cannot find sensitivity curve!",noheader=True)
        return False

    with fits.open(os.path.join(w_dir,filein)) as hdu:
        image=hdu[0].data[:,:]#-2]
        header=hdu[0].header
    with fits.open(os.path.join(w_dir,file_sens)) as hdu:
        sens=hdu[0].data
        header_sens=hdu[0].header

    star_name=header_sens.get('object')
    texp_obj = header.get('exptime')
    if z_obj is not None:
        z_obj = header.get('z')
    if not z_obj:
        z_obj=20.
        message("... Cannot derive Z(obj). Use Z=20 degree.", noheader=True)

    dlam=header.get('cdelt1')
    wlscale = dlam * (np.arange(image.shape[1]) - header.get('crpix1') + 1) + header.get('crval1')
    ext_obj = 10 ** (0.4 * calc_ext(wlscale, z_obj))
    z_star = header_sens.get('z')
    if not z_star:
        z_star = 20.
        message("... Cannot derive Z(star). Use Z=20 degree.", noheader=True)
    sens = sens / (texp_obj*dlam)
    # order = int(round(min(np.log10(sens))))
    # if order < 0:
    #     order=order-1
    # output_units = 10.**order #output units in erg/cm^2/s/A
    image[:,:]=image[:,:]*sens[None,:]*ext_obj[None,:]
    header.add_history("LONGRED: flux calibration based on {} standard".format(star_name))
    header["BUNIT"]="erg/cm^2/sec/A"
    if ypos is None:
        ypos=int(np.ceil(image.shape[0]/2))-ycut[0]
    else:
        ypos = ypos - ycut[0]
    if ydelt is None:
        ydelt=1
    header["CRPIX2"]=ypos
    header["CDELT2"] = ydelt
    header["CRVAL2"] = 0
    wlind0 = 0
    wlind1 = len(wlscale) - 1
    if cut_wl[0] is not None:
        wlind0=np.argmin(abs(wlscale-cut_wl[0]))
    if cut_wl[1] is not None:
        wlind1=np.argmin(abs(wlscale-cut_wl[1]))
    image=image[ycut[0]:ycut[1],wlind0:wlind1+1]
    header["CRPIX1"]=1
    header["CRVAL1"]=wlscale[wlind0]

    fits.writeto(os.path.join(w_dir,fileout),np.float32(image),header=header,overwrite=True)
    if file_err is not None and fileout_err is not None:
        with fits.open(os.path.join(w_dir,file_err)) as hdu:
            hdu[0].data = hdu[0].data[:, :] * sens[None, :] * ext_obj[None, :]
            hdu[0].header.add_history("LONGRED: flux calibration based on {} standard".format(star_name))
            hdu[0].header["BUNIT"] = "erg/cm^2/sec/A"
            hdu[0].header["CRPIX2"] = ypos
            hdu[0].header["CDELT2"] = ydelt
            hdu[0].header["CRVAL2"] = 0
            hdu[0].data = np.float32(hdu[0].data[ycut[0]:ycut[1], wlind0: wlind1 + 1])
            hdu[0].header["CRPIX1"] = 1
            hdu[0].header["CRVAL1"] = wlscale[wlind0]
            hdu.writeto(os.path.join(w_dir, fileout_err), overwrite=True)

    if files_for_cut is not None and files_for_cut_out is not None:
        for ind,f in enumerate(files_for_cut):
            if f:
                with fits.open(os.path.join(w_dir, f)) as hdu:
                    hdu[0].header["CRPIX2"] = ypos
                    hdu[0].header["CDELT2"] = ydelt
                    hdu[0].header["CRVAL2"] = 0
                    hdu[0].data = np.float32(hdu[0].data[ycut[0]:ycut[1], wlind0: wlind1 + 1])
                    hdu[0].header["CRPIX1"] = 1
                    hdu[0].header["CRVAL1"] = wlscale[wlind0]
                    hdu.writeto(os.path.join(w_dir, files_for_cut_out[ind]), overwrite=True)



    return True

if __name__ == "__main__":
    print(
        """
        All the main steps of the longslit data reduction pipeline are here
        """
    )
