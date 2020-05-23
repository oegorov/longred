
"""
Different auxiliary programms used for data reduction and files manipulation are located here
==============================================================
"""

import numpy as np
from scipy import signal,ndimage,interpolate
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
import os, errno
from astropy.io import fits, ascii
from astropy.table import Table,vstack
from misc import message
import re
from matplotlib import pyplot as plt
import warnings


def line_center(arr, bin=1, win=None):
    x = np.linspace(0, arr.shape[0], arr.shape[0]*bin)
    if bin>1:
        arr=signal.resample_poly(np.pad(arr,arr.shape[0]//3,mode='reflect'),bin,1)[arr.shape[0]//3*bin:-arr.shape[0]//3*bin]
    if win is None:
        win=int(np.floor((arr.shape[0]-1)/4))
        # center=int(round(arr.shape/2))

    pos = np.argmax(arr)
    if (pos < win) or pos > (arr.shape[0] - 1 - win):
        return np.nan

    flux = np.nansum(arr[pos - win:pos + win + 1])
    return np.nansum(arr[pos - win:pos + win + 1] * x[pos - win:pos + win + 1]) / flux


def simplepoly(x,y, yerr=None,order=3, get_masked_data=False):
    if yerr is not None:
        index = ~(np.isnan(x) | np.isnan(y) | ~np.isfinite(y) | ~np.isfinite(x) | np.isnan(yerr) | ~np.isfinite(yerr)) & (y!=0) & (yerr!=0)
        weights = 1. / yerr[index]
    else:
        index = ~(np.isnan(x) | np.isnan(y) | ~np.isfinite(y) | ~np.isfinite(x)) & (y!=0)
        weights=np.ones_like(y[index])
    if np.sum(index)<5:
        cf=np.zeros(order+1,dtype=float)
        cf[-1]=np.nanmedian(y)
    else:
        cf = np.polyfit(x[index], y[index], order,w=weights)
    fitted_line = np.poly1d(cf)
    if get_masked_data:
        filtered_data = np.ma.masked_array(y, mask=~index)
        return (fitted_line, cf, filtered_data)
    else:
        return (fitted_line, cf)


def goodpoly(x,y,yerr=None, order=3, sigma=3, niter=3, get_masked_data=False):
    """
    Perform polynomial fitting using 3 iteration of sigma-clipping
    Return fit function; its parameters; (also masked data if get_masked_data=True)
    """
    fit = fitting.LinearLSQFitter()
    line_init=models.Polynomial1D(order)
    good_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip, niter=niter, sigma=sigma)
    if yerr is not None:
        weights=1./yerr
    else:
        weights=np.ones_like(y)
    if yerr is not None:
        index = ~(np.isnan(x) | np.isnan(y) | np.isnan(yerr))
    else:
        index = ~(np.isnan(x) | np.isnan(y))
    # fitted_line = fit(line_init, x, y)#,weights=weights)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fitted_line, mask = good_fit(line_init, x[index], y[index], weights=weights[index])

    if get_masked_data:
        filtered_data = np.ma.masked_array(y, mask=False)#mask)
        return (fitted_line,fitted_line.parameters.copy(),filtered_data)
    else:
        return (fitted_line, fitted_line.parameters.copy())



def do_distortion(image, transform_matrix, reverse=False, xchunk_size=500, strip_size=200, pad_border=100):
    """
    Perform geometric transform of the NYxNX image based on the transform_matrix,
    which is (2,NY,NX) array with shifts F([0,1],Y,X) =  ([DY,DX]) of original frame with repect to target
    Simply apply transformation if it is a direct transform.
    Much more complex if backward transform - we need to spline transform matrix and also
    to do it within several chunks because of memory issues
    So, if reverse = True, then one may also set:
    xchunk_size - maximus size of individual block along X-axis
    strip_size - number of overlapping pixels between blocks to exclude edge effects
    pad_border - number of pixels to skip from left and right when filling up final image with those from individual blocks
    transform_matrix should be the same as for direct transform (not inversed!!!)
    """

    def get_direct_trans_coords(xy, tmatrix):
        return xy[0] + tmatrix[0, xy[0], xy[1]], xy[1] + tmatrix[1, xy[0], xy[1]]
    def get_backward_trans_coords(xy,sbs):
        return sbs[0].ev(xy[0], xy[1]), sbs[1].ev(xy[0], xy[1])
    # def get_direct_trans_coords_spline(xy,sbs):
    #     return sbs[0].ev(xy[0], xy[1]), sbs[1].ev(xy[0], xy[1])

    if not reverse:
        xchunk_size=transform_matrix.shape[2]
        strip_size=0
        pad_border=0
    niter = int(np.ceil(float(transform_matrix.shape[2] - strip_size) / (xchunk_size - strip_size)))
    image_out = np.zeros_like(image, dtype=float)
    for curiter in range(niter):
        xstart = curiter * (xchunk_size - strip_size)
        xfinish = np.min([xstart + xchunk_size, transform_matrix.shape[2]])
        cur_tmatrix = transform_matrix[:, :, xstart:xfinish]
        cur_image = image[:, xstart:xfinish]
        xx, yy = np.meshgrid(np.arange(cur_tmatrix.shape[2]), np.arange(cur_tmatrix.shape[1]))
        coords_in_input = get_direct_trans_coords((yy, xx),cur_tmatrix)
        # sbs_y = interpolate.SmoothBivariateSpline(yy.ravel(),xx.ravel(),coords_in_input[0].ravel())
        # sbs_x = interpolate.SmoothBivariateSpline(yy.ravel(),xx.ravel(),coords_in_input[1].ravel())
        # coords_in_input = get_direct_trans_coords_spline((yy, xx), (sbs_y, sbs_x))
        if reverse:
            sbs_y = interpolate.SmoothBivariateSpline(coords_in_input[0].ravel(), coords_in_input[1].ravel(),
                                                      yy.ravel())
            sbs_x = interpolate.SmoothBivariateSpline(coords_in_input[0].ravel(), coords_in_input[1].ravel(),
                                                      xx.ravel())
            coords_in_input = get_backward_trans_coords((yy, xx),(sbs_y,sbs_x))

        # def transform(output_coords):
        #     return get_direct_trans_coords((output_coords[0], output_coords[1]),cur_tmatrix)

        cur_image_out = ndimage.map_coordinates(cur_image, coords_in_input, mode='constant',order=1)#,cval=np.nan)
        # cur_image = ndimage.geometric_transform(cur_image, transform, mode='constant', order=1)  # ,cval=np.nan)
        if niter == 1:
            image_out = cur_image_out
        elif curiter == 0:
            image_out[:, :xfinish - pad_border] = cur_image_out[:, :-pad_border]
        elif curiter == (niter - 1):
            image_out[:, xstart + pad_border:] = cur_image_out[:, pad_border:]
        else:
            image_out[:, xstart + pad_border:xfinish - pad_border] = cur_image_out[:, pad_border:-pad_border]

    return image_out

def calc_ext(wlscale, z):
    """
    Derive atmospheric extinction value (in mag) at each wavelength based on zenith distance
    """
    a = 0.012#0.008
    c = 0.12#0.115
    extin = (a * 1. / ((wlscale/10000.)**4.) + c) * 2.5*np.log10(np.e)
    extin = (extin / np.cos(z*np.pi/180.)*(1-0.0012*(1/np.cos(z*np.pi/180)**2-1)))
    return extin


def normalize(data):
    return (data - np.nanmean(data)) / (np.nanstd(data))

def cos_apod(nsample, perc=10.):
    y=np.ones(nsample)
    nperc=int(np.round(nsample*perc/100))
    x=np.sin(np.pi/2/nperc*np.arange(nperc))
    y[:nperc]=x
    y[-nperc:]=np.flip(x)
    return y


def derive_vecshift(vec, vec_ref, bin=50, max_ampl= None, plot=False):
    """
    Derive shift of 1D-array vec from vec_ref;
    Oversampling is applied using bin
    if max_ampl is set than maximum shift is max_ampl
    """
    nsamples = min([len(vec), len(vec_ref)])
    vec=signal.resample_poly(cos_apod(nsamples)*normalize(vec[:nsamples]), bin, 1)
    vec_ref = signal.resample_poly(cos_apod(nsamples)*normalize(vec_ref[:nsamples]),  bin, 1)
    xcorr = signal.correlate(vec, vec_ref)
    if max_ampl:
        max_ampl=min([(nsamples*bin-1),int(np.floor(max_ampl*bin))])
        xcorr=xcorr[nsamples*bin-(max_ampl+1) : nsamples*bin+max_ampl]
    else:
        max_ampl = nsamples*bin-1
    dt = np.arange( - max_ampl, max_ampl+1)
    shift = dt[xcorr.argmax()]/bin

    if plot:
        plt.figure(figsize=(7,7))
        plt.plot(np.arange(0, nsamples, 1. / bin), vec_ref, label="Reference signal")
        # plt.plot(np.arange(0,nsamples,1./bin),vec,label="Shifted signal")
        # plt.plot(np.arange(0, nsamples, 1. / bin)-shift, vec, label="Corrected signal")
        plt.title("Shift = {:.2f} bins".format(shift))
        plt.legend()
        plt.draw()
        plt.pause(0.001)
    return(shift)



def derive_shift(image_frames, ref_frame=0, bin=50, do_shift = (True,True), max_ampl=(None, None)):
    """
    Определяем относительное смещение изображений в image_frames относительно опорного кадра, заданного ref_frames
    :param image_frames: 3D array with image frames
    :param ref_frame: reference frame number in 3D input array
    :param bin: Do oversampling with corresponding binning
    :param do_shift: whether to do shift on (X; Y)?
    :param max_ampl: maximum amplitude for shift on (X; Y)
    :param
    :return:
    """
    s=image_frames.shape
    if len(s) < 3:
        message("Cannot derive shift => image_frames must be 3D-array!")
        return (None,None)
    if not (do_shift[0] or do_shift[1]):
        return (None,None)
    dx=np.array([],dtype=float)
    dy=np.array([],dtype=float)

    xlim=(int(np.ceil(s[2]*0.3)),int(np.floor(s[2]*0.7)))
    ylim = (int(np.ceil(s[1] * 0.3)), int(np.floor(s[1] * 0.7)))

    if do_shift[0]:
        xvec_ref=np.nanmean(image_frames[ref_frame,ylim[0]:ylim[1]+1,xlim[0]:xlim[1]+1],axis=0)
    if do_shift[1]:
        yvec_ref = np.nanmean(image_frames[ref_frame, ylim[0]:ylim[1] + 1, xlim[0]:xlim[1] + 1], axis=1)
    for ind in range(s[0]):
        if ind == ref_frame:
            dx=np.append(dx,0)
            dy=np.append(dy,0)
            continue
        if do_shift[0]:
            vec = np.nanmean(image_frames[ind, ylim[0]:ylim[1] + 1, xlim[0]:xlim[1] + 1], axis=0)
            dx=np.append(dx,derive_vecshift(vec,xvec_ref,bin=bin,max_ampl=max_ampl[0]))
        else:
            dx=np.append(dx,0)
        if do_shift[1]:
            vec = np.nanmean(image_frames[ind, ylim[0]:ylim[1] + 1, xlim[0]:xlim[1] + 1], axis=1)
            dy=np.append(dy,derive_vecshift(vec,yvec_ref,bin=bin,max_ampl=max_ampl[1]))
        else:
            dy=np.append(dy,0)
    return (dx,dy)

def shift_images(image_frames, dx,dy):
    """
    Смещаем изображения, содержащиеся в image_frames, на значения по x и y, указанные в dx,dy
    """
    s=image_frames.shape
    for ind in range(s[0]):
        if dx[ind] !=0 or dy[ind] !=0:
            image_frames[ind,:,:] = ndimage.shift(image_frames[ind,:,:],[dy[ind],dx[ind]],order=1)
    return image_frames

def ch_clean_pair(image1,image2,t=5,lim=10,gain=1,nonorm=False,do_print=True):
    """
    Чистка космических частиц и артефактов на паре изображений image1 и image2
    """
    norm = 1.
    if not nonorm:
        norm=float(np.median(image1)) / float(np.median(image2))

    if do_print and not nonorm:
        print('Norm factor: {}'.format(norm))
    dif = (image1 - image2 * norm)
    dif=dif.flatten()
    r1 = (dif >= 0)
    counts1=np.sum(r1)
    s = image1.shape
    med = image1.flatten()
    image2=image2.flatten()
    if counts1 > 0:
        med[r1]=image2[r1] * norm
    for i,m in enumerate(med):
        med[i]=abs(m)
        if med[i] < 1:
            med[i]=1
    dif = dif/np.sqrt(gain * med)
    image2.reshape(s)
    dif.reshape(s)

    sig = np.sum(abs(signal.medfilt(dif,kernel_size=3))) / (s[0] * s[1])
    dif=dif.flatten()
    rec1 = (dif > t*sig)
    count1=np.sum(rec1)
    rec2 = (dif <  - t * sig)
    count2=np.sum(rec2)

    mask = image1*0
    mask=mask.flatten()
    if count1 > 0:
        mask[rec1]=100.
    if count2 > 0:
        mask[rec2]=-100.
    mask.reshape(s)
    mask = ndimage.gaussian_filter(mask, 1, mode='nearest')
    mask.flatten()
    if count1 > 0:
        rec1=mask > lim
        count1=np.sum(rec1)
    if count2 > 0:
        rec2=mask < -lim
        count2=np.sum(rec2)
    image1=image1.flatten()
    image2=image2.flatten()
    if do_print:
        print("{} points were replaced".format(count1 + count2))
    if count1 > 0:
        image1[rec1]=image2[rec1] * norm
    if count2 > 0:
        image2[rec2]=image1[rec2] / norm
    return(image1.reshape(s),image2.reshape(s))



def combine_pairs(images, cr_treshold = 10, ncycles = 2):
    """
    Производится попарное сравнение кадров, содержащихся в кубе данных images,
    вычищаются области с отличиями > cr_treshold, после чего возвращается сумма
    Разница в фоне не учитывается!!!
    """
    nfiles=len(images)
    npair = int(nfiles / 2)
    for iter in range(ncycles):
        for i in range(npair):
            im1, im2 = ch_clean_pair(images[i * 2, :, :], images[i * 2 + 1, :, :], t=cr_treshold, nonorm=True)
            images[i * 2, :, :] = im1
            images[i * 2 + 1, :, :] = im2
            del im1, im2
        if nfiles % 2 == 1:
            im1, im2 = ch_clean_pair(images[nfiles - 2, :, :], images[nfiles - 1, :, :], t=cr_treshold,
                                     nonorm=True)
            images[nfiles - 2, :, :] = im1
            images[nfiles - 1, :, :] = im2
        np.roll(images, 1, axis=0)
    return(np.nansum(images, axis=0))



def combine_sigma(images, cr_treshold = 5, return_masked_array=False):
    """
    Производится сложение изображений с использованием sigma-clipping
    Разница в фоне не учитывается, нормировка не выполняется!!!!!
    """
    nfiles=len(images)
    s=images.shape
    images = images.reshape(s[0],-1)
    clip_mask = np.zeros((s[0],s[1]*s[2]),dtype=bool)
    clip_mask[~np.isfinite(images)] = True

    images = np.ma.array(images, mask=clip_mask, fill_value=np.nan)
    images=sigma_clip(images, sigma=cr_treshold, axis=0)
    if not return_masked_array:
        return np.ma.filled(np.ma.mean(images,axis=0)*nfiles,np.nan).reshape(s[1],s[2])
    else:
        return images.reshape(s[0],s[1],s[2])


def remove_overscan(image, mode = "TDS", cut=None):
    """
    Обрезаем область оверскана (или еще чего)
    """
    if cut is not None:
        return image[cut[2]:cut[3],cut[0]:cut[1]]
    else:
        if mode in ["TDS",'TDS_r','TDS_b']:
            return(image)
        # elif mode in ['TDS_r']:
        #     return(image[:497,:2047])
        # elif mode in ['TDS_b']:
        #     return(image[:,140:])
        elif mode == "SAO":
            s=image.shape
            return image[21:,21:]
        else:
            return (image)


def check_overscan(image, mode = "TDS", force_overscan=None):
    """
    Считаем уровень остатка после вычитания bias и dark по оверскану
    Вычитаем этот уровень и возвращаем полученное изображение
    Для TDS оверскана как такового и нет... Хм...
    force_overscan - задаем пиксели для оверскана вместо тех, что тут по умолчанию
    """
    if force_overscan is not None:
        arr=np.array([],dtype=float)
        if force_overscan[0]:
            arr=np.concatenate(image[:,:force_overscan[0]].ravel())
        if force_overscan[1]:
            arr = np.concatenate(image[:, force_overscan[1]].ravel())
        if force_overscan[2]:
            arr = np.concatenate(image[:force_overscan[2],:].ravel())
        if force_overscan[3]:
            arr = np.concatenate(image[force_overscan[3]:,:].ravel())
        med_val=np.median(arr)
        return(image-med_val)
    else:
        if mode in ["TDS","TDS_b","TDS_r"]:
            return(image)
        # elif mode == 'TDS_r':
            # med_val = np.median(image[508:511, :].ravel())
            # return (image - med_val)
        elif mode == "SAO":
            med_val=np.median(np.concatenate((image[0:12,:].ravel(),image[:,0:12].ravel())))
            return (image-med_val)
        else:
            return (image)
def check_param(head, keys=[], val = None):
    """
    Проверяем наличие значения val в одном из параметров keys в шапке (по порядку). Если есть - возвращаем True.
    """
    if len(keys) == 0:
        return False
    for key in keys:
        param_val=head.get(key)
        if (type(param_val) == str) and (type(val) == str):
            param_val=param_val.casefold()
            val = val.casefold()
        if param_val == val:
            return True
    else:
        return False

def get_param(head, keys=[]):
    """
    Ищем первое присутствующее ключевое слово в keys и возвращаем его значение
    """
    if len(keys) == 0:
        return None
    for key in keys:
        try:
            param_val = head[key]
            return param_val
        except (ValueError, TypeError, NameError, KeyError):
            pass
    else:
        return None

def check_if_star(object, standards = "/Users/mors/Science/standards/standards_list.txt"):
    """
    Проверяем, не относится ли текущий объект к списку стандартов. Если да - то возвращаем True
    """
    if not os.path.isfile(standards):
        return False
    stars = ascii.read(standards,format="no_header")['col1']
    if re.sub("[\._\-+\s]+","",object.casefold()) in map(lambda x:re.sub("[\._\-+\s]+","",x.casefold()), stars):
        return True
    else:
        return False

def angle_range(pa):
    """
    Convert input PA to the range (0:360) degrees
    """
    if pa >= 360:
        pa-=360
        pa=angle_range(pa)
    elif pa < 0:
        pa+=360
        pa = angle_range(pa)
    return pa


def parse_redfits(file, mode="TDS"):
    currow = Table(names=('data','filenum','type', 'subtype', 'state', 'disperser', "binning", "subdirectory","filename", "slit", "process"),
                       dtype=('S','S','S', 'S', 'S', 'S', 'S', 'S64','S64', 'S', "?"))
    currow.add_row()
    currow["process"]=True
    #currow["universal"] = False
    all_prefixes={'obj':"obj",'star':"star",'sunsky':"calibration",'neon':"calibration",'flat':"calibration",
                  '13dots':"calibration",'meanbias':"calibration",'dark':"calibration",'transform':"auxiliary",
                  'senscurve':'auxiliary', 'dispcurve':'auxiliary', 'shifts':'auxiliary', 'atmdisp':'auxiliary'}
    currow['filenum'] = None
    for curtype in all_prefixes.keys():
        if file.casefold().startswith(curtype):
            currow['type']=all_prefixes[curtype]
            if curtype not in ['dark', 'meanbias']:
                if curtype in ['transform','shifts', 'atmdisp']:
                    add='_?\w*'
                else:
                    add=''

                prefix_find = re.findall(r"({}{})_s\d\.?\d?_".format(curtype, add), file)
                if len(prefix_find)>0:
                    currow['data'] = prefix_find[0]
                else:
                    currow['data'] = 'unknown'
                if curtype in ['obj', 'star']:
                    ind_find = re.findall(r"{}_s\d\.?\d?_[\w@]+_\dx\d_?(\d?\d?)".format(curtype), file)
                    if len(ind_find) > 0 and ind_find[0] != '':
                        currow['data'] = "{}_{}".format(currow['data'][0],ind_find[0])
                        currow['filenum'] = "{}".format(ind_find[0])
                slit_find = re.findall(r"{}{}_s(\d\.?\d?)_".format(curtype,add), file)
                if len(slit_find) > 0:
                    currow['slit'] = slit_find[0]
                else:
                    currow['slit'] = None
                disp_find = re.findall(r"{}{}_s\d\.?\d?_([\w@]+)_\dx\d".format(curtype,add), file)
                if len(disp_find) > 0:
                    currow['disperser'] = disp_find[0]
                else:
                    currow['disperser'] = None
                bin_find = re.findall(r"{}{}_s\d\.?\d?_[\w@]+_(\dx\d)".format(curtype,add), file)
                if len(bin_find) > 0:
                    currow['binning'] = bin_find[0]
                else:
                    currow['binning'] = None
            else:
                if mode=='TDS':
                    cur_format_bin=r"{}_[\w@]+_(\dx\d)".format(curtype)
                    disp_find = re.findall(r"{}_([\w@]+)_\dx\d".format(curtype), file)
                    if len(disp_find) > 0:
                        if disp_find[0] == 'r':
                            currow['data'] = "{}{}".format(curtype,'_r')
                            currow['disperser'] = 'R'
                        else:
                            currow['data'] = "{}{}".format(curtype, '_b')
                            currow['disperser'] = 'B,G'
                    else:
                        currow['data'] = curtype
                        currow['disperser'] = None
                else:
                    currow['data'] = curtype
                    cur_format_bin = r"{}_(\dx\d)".format(curtype)
                    currow['disperser'] = None
                bin_find = re.findall(cur_format_bin, file)
                if len(bin_find) > 0:
                    currow['binning'] = bin_find[0]
                else:
                    currow['binning'] = None

            break
    else:
        return (None, False)
    currow['filename'] = os.path.basename(file)

    if file.casefold().endswith('_err.fits'):
        currow['subtype'] = 'unc'
    else:
        currow['subtype'] = 'data'

    if file.casefold().endswith('_crmask.fits') or file.casefold().endswith('_test.fits') or \
            file.casefold().endswith('_test_err.fits'):
        currow['state']='unknown'
    elif file.casefold().endswith('_lin.fits') or file.casefold().endswith('_lin_err.fits'):
        currow['state'] = 'linearized'
    elif file.casefold().endswith('_n.fits') or file.casefold().endswith('_n_err.fits'):
        currow['state'] = 'normalized'
    elif file.casefold().endswith('_skyrem.fits') or file.casefold().endswith('_skyrem_err.fits'):
        currow['state'] = 'sky-subtracted'
    elif file.casefold().endswith('_tot.fits') or file.casefold().endswith('_tot_err.fits'):
        currow['state'] = 'combined'
    elif file.casefold().endswith('_abs.fits') or file.casefold().endswith('_abs_err.fits'):
        currow['state'] = 'calibrated'
    elif (currow['type'] == 'calibration' or currow['type'] == 'auxiliary') and \
            (file.casefold().endswith("_{}.fits".format(currow['binning'][0])) or
             file.casefold().endswith("_{}.dat".format(currow['binning'][0])) or
             file.casefold().endswith("_{}_err.fits".format(currow['binning'][0]))):
        currow['state'] = 'ini'
    elif re.match(r'\S*_\d\d.fits', file.casefold()) or re.match(r'\S*_\d\d_err.fits', file.casefold()):
        currow['state'] = 'ini'
    else:
        currow['state'] = 'unknown'

    return (currow,True)

def parse_fits(file, mode = "TDS", standards = "/Users/mors/Science/standards/standards_list.txt"):
    status = False
    currow = Table(names=('type', 'object', 'disperser', "exposure", "binning", "date", "subdirectory","filename", "slit", "pa", "process"),
                   dtype=('S', 'S', 'S', 'f8', 'S', 'S', 'S','S', 'S', 'S', "?"))
    currow.add_row()
    currow["process"]=True
    with fits.open(file) as hdu:
        head = hdu[0].header
        if mode == "TDS":
            #==========================================================
            #===== Derive parameters in each fits for TDS spectra =====
            # ==========================================================
            if not check_param(head,keys = ["instrume"], val = "tds"):
                return (None,False)

            currow['binning'] = "{}x{}".format(head['hbin'], head['vbin'])
            currow['date'] = head['date']
            currow['filename'] = os.path.basename(file)
            currow['exposure']=get_param(head,keys=['exposure','exptime'])


            #===== Check type of disperser
            grism_set = get_param(head,keys=['bdisp','disp','disper','dispers','disperse'])
            channel=get_param(head, keys=["head"])[-1]
            if channel == 'U':
                if grism_set == "G":
                    disperser = "G"
                else:
                    disperser = 'B'
            elif channel == 'V':
                disperser = 'R'
            else:
                message("Unrecognized channel...")
                disperser = None

            # ===== Check slit type (and search for indication of pa)
            slit = 1 #default value
            pa_in_slit = None
            slit_header = get_param(head, keys=["slit"])
            if type(slit_header) == str:
                try:
                    slit=float(slit_header)
                except ValueError:
                    slit_find = re.findall(r"slit[_\s]?([\w])", slit_header)
                    if len(slit_find) > 0:
                        try:
                            slit=float(re.sub("_",".",slit_find[0]))
                        except ValueError:
                            pass

                    pa_find=re.findall(r"pa=?([\w\._]+)", slit_header)
                    if len(pa_find) > 0:
                        try:
                            pa_in_slit=float(re.sub("_",".",pa_find[0]))
                        except ValueError:
                            pa_in_slit=None
            elif (type(slit_header) == int) or (type(slit_header) == float):
                slit=slit_header
            else:
                pass
            if slit is not None:
                slit="{0:.1g}".format(slit).rstrip('0').rstrip('.')
            if check_param(head,keys = ["object","imagetyp"], val = "calibration") or \
                    check_param(head,keys = ["object","imagetyp"], val = "arc"):
                currow['type']='neon'
                currow['object'] = None
                currow['disperser'] = disperser
                currow['slit'] = slit
                currow['pa'] = None
            elif check_param(head,keys = ["object","imagetyp"], val = "bias") or \
                    (currow['exposure'] == 0):
                currow['type'] = 'bias'
                currow['object'] = None
                currow['disperser'] = disperser
                currow['slit'] = None
                currow['pa'] = None
            elif check_param(head,keys = ["object","imagetyp"], val = "dark") and (currow['exposure'] > 0):
                currow['type'] = 'dark'
                currow['object'] = None
                currow['disperser'] = disperser
                currow['slit'] = None
                currow['pa'] = None
            elif check_param(head,keys = ["object","imagetyp"], val = "sunsky"):
                currow['type'] = 'sunsky'
                currow['object'] = None
                currow['disperser'] = disperser
                currow['slit'] = slit
                currow['pa'] = None
            elif check_param(head,keys = ["object","imagetyp"], val = "flat"):
                currow['type'] = 'flat'
                currow['object'] = None
                currow['disperser'] = disperser
                currow['slit'] = slit
                currow['pa'] = None
            else:
                currow['disperser'] = disperser
                currow['slit'] = slit
                curobj =get_param(head,['object'])

                pa_in_header = get_param(head, ['pa'])
                pa_in_obj = re.findall(r"pa[=\s_]?([\w\._]+)", curobj.casefold())
                if len(pa_in_obj) == 1:
                    curobj_splitted = re.split(r'[_\s]+pa[=\s]?[\w]+', curobj, flags=re.IGNORECASE)
                    if len(curobj_splitted) > 1:
                        curobj = curobj_splitted[0]
                    pa_in_obj = float(re.sub("_", ".", pa_in_obj[0]))
                else:
                    pa_in_obj=None
                if pa_in_header is not None:
                    currow['pa']=float(re.sub("_",".",pa_in_header))
                elif pa_in_obj is not None:
                    currow['pa']=pa_in_obj
                else:
                    currow['pa'] = pa_in_slit

                currow['object'] = curobj
                if check_if_star(curobj,standards=standards):
                    currow['type'] = 'star'
                else:
                    currow['type'] = 'obj'

            status = True

        elif mode == "SAO":
            # ==========================================================
            # ===== Derive parameters in each fits for SCORPIO spectra =====
            # ==========================================================
            if not ((check_param(head, keys=["instrume"], val="scorpio") or check_param(head, keys=["instrume"], val="scorpio-1") or
                    check_param(head, keys=["instrume"], val="scorpio-2")) and (check_param(head,keys=["mode"],  val="spectra") or
                    (check_param(head,keys=["mode"],  val="image") and check_param(head,keys=["imagetyp"],  val="bias")))):
                return (None, False)

            currow['binning'] = head['binning']
            currow['date'] = head['date']
            currow['filename'] = os.path.basename(file)
            currow['exposure'] = get_param(head, keys=['exptime','exposure'])

            disperser = get_param(head, keys=['disperse'])
            slit=get_param(head,keys=['slitwid'])
            if slit:
                slit=np.round(slit,1)
            if check_param(head, keys=["instrume"], val="scorpio-1"):
                const = 132.5 + 89.
            else:
                obs_year=int(head['date'][:4])
                if obs_year < 2006:
                    const=131.
                elif obs_year >= 2008:
                    const=132.
                else:
                    const=132.5
                if check_param(head, keys=["instrume"], val="scorpio-2"):
                    const=const+180
            pa = "{0:.02f}".format(angle_range(get_param(head, keys=['parangle']) - get_param(head, keys=['rotangle']) + const)).rstrip('0').rstrip('.')


            if check_param(head,keys = ["imagetyp","object"], val = "neon"):
                currow['type']='neon'
                currow['object'] = None
                currow['disperser'] = disperser
                currow['slit'] = "{0:.1g}".format(slit).rstrip('0').rstrip('.')
                currow['pa'] = None
            elif check_param(head,keys = ["imagetyp","object"], val = "bias") or \
                    (currow['exposure'] == 0):
                currow['type'] = 'bias'
                currow['object'] = None
                currow['disperser'] = None
                currow['slit'] = None
                currow['pa'] = None
            elif check_param(head,keys = ["imagetyp","object"], val = "dark") and (currow['exposure'] > 0):
                currow['type'] = 'dark'
                currow['object'] = None
                currow['disperser'] = None
                currow['slit'] = None
                currow['pa'] = None
            elif check_param(head,keys = ["object","imagetyp"], val = "sunsky"):
                currow['type'] = 'sunsky'
                currow['object'] = None
                currow['disperser'] = disperser
                currow['slit'] = "{0:.1f}".format(slit).rstrip('0').rstrip('.')
                currow['pa'] = None
            elif check_param(head,keys = ["imagetyp","object"], val = "flat"):
                if check_param(head,keys=['slitmask'], val='13 dots'):
                    currow['type'] = '13dots'
                else:
                    currow['type'] = 'flat'
                currow['object'] = None
                currow['disperser'] = disperser
                currow['slit'] = "{0:.1f}".format(slit).rstrip('0').rstrip('.')
                currow['pa'] = None
            else:
                currow['disperser'] = disperser
                currow['slit'] = "{0:.1f}".format(slit).rstrip('0').rstrip('.')
                curobj = get_param(head,['object'])
                #pa_in_header = get_param(head, ['pa'])
                pa_in_obj = re.findall(r"pa[=\s_]?([\w\._]+)", curobj.casefold())
                if len(pa_in_obj) == 1:
                    curobj_splitted = re.split(r'[_\s]+pa[=\s]?[\w]+', curobj, flags=re.IGNORECASE)
                    if len(curobj_splitted) > 1:
                        curobj = curobj_splitted[0]
                    #pa_in_obj = float(re.sub("_", ".", pa_in_obj[0]))
                #else:
                #    pa_in_obj = None

                currow['pa'] = pa
                currow['object'] = curobj
                if check_if_star(curobj, standards=standards):
                    currow['type'] = 'star'
                else:
                    currow['type'] = 'obj'

            status = True

    return(currow, status)


def scan_dir(raw_dir, mode="TDS", standards = "/Users/mors/Science/standards/standards_list.txt", dir_mode="Raw", exist_table=None):
    """
    This routine scan selected directory for object, flat, neon (or calibration), dark, sunsky, bias and other cal. frames.
    It returns a dictionary with lists of filenames for each type
    """
    if dir_mode=='Raw':
        files = {"bias": [], "dark": [], "flat": [], "obj": [], "neon": [], "sunsky": [], "other": [], "star": []}
    else:
        files=None
    all_fits = []
    fnames=[]
    for x in os.walk(raw_dir, followlinks=True):
        if dir_mode == 'Raw':
            all_fits.extend([os.path.join(x[0], f) for f in x[2] if (f.lower().endswith(".fits") or f.lower().endswith(".fts"))])
        else:
            flist=[f for f in x[2] if (f.lower().endswith(".fits")) or (f.lower().endswith(".dat"))]
            fnames.extend(flist)
            all_fits.extend([os.path.join(x[0],f) for f in flist])
    ngood_files=0
    tab=None
    for ind,f in enumerate(all_fits):
        if dir_mode == 'Raw':
            currow, status = parse_fits(f, mode=mode, standards = standards)
        else:
            currow, status = parse_redfits(fnames[ind], mode=mode)
        if status:
            currow['subdirectory'] = f.split(os.path.join(raw_dir, ''))[1].split(currow['filename'][0])[0]
            if exist_table is not None:
                try:
                    ind = list(exist_table['filename']).index(currow['filename'][0])
                    if currow['subdirectory'][0] == exist_table['subdirectory'][ind]:
                        currow['process'][0] = exist_table['process'][ind]
                except ValueError:
                    pass
            if ngood_files == 0:
                tab=currow
            else:
                tab = vstack([tab, currow])
            ngood_files+=1
            if dir_mode=="Raw":
                if currow["type"][0] in files.keys():
                    files[currow["type"][0]].append(f)
                else:
                    files["other"].append(f)

    return(files, tab)



# def scan_w_dir(w_dir, mode="TDS"):
#     """
#     This routine scan selected directory for data in reduction process
#     It returns a dictionary with lists of filenames for each type
#     """
#     files = {"bias": [], "dark": [], "flat": [], "obj": [], "neon": [], "sunsky": [], "other": [], "star": []}
#
#     all_fits = []
#     for x in os.walk(w_dir):
#         all_fits.extend([os.path.join(x[0], f) for f in x[2] if (f.lower().endswith(".fits") and ("test" not in f.lower()))])
#     ngood_files=0
#     tab=None
#     for f in all_fits:
#         currow, status = parse_fits(f, mode=mode)
#         if status:
#             currow['subdirectory'] = f.split(os.path.join(w_dir, ''))[1].split(currow['filename'][0])[0]
#             if ngood_files == 0:
#                 tab=currow
#             else:
#                 tab = vstack([tab, currow])
#             ngood_files+=1
#             if currow["type"][0] in files.keys():
#                 files[currow["type"][0]].append(f)
#             else:
#                 files["other"].append(f)
#
#     return(files, tab)