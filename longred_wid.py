#!/Applications/anaconda3/bin/python
# -*- coding: utf-8 -*-
import tkinter as tk
import re
import tkinter.font as tkFont
import tkinter.filedialog as tkFileDialog
import tkinter.ttk as ttk
from auxiliary import scan_dir,parse_redfits
import os,errno
from astropy.table.pprint import conf
import reduction
from misc import message
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np
from astropy.table import vstack
from operator import itemgetter
matplotlib.use("TkAgg")
from astropy.visualization import PercentileInterval, ZScaleInterval, LinearStretch, \
    AsinhStretch, SqrtStretch, LogStretch, ImageNormalize


#########################
# Geometry parameters
#########################

win_height=700
button_width=15
table_width=820
table_height=380


###################################
# Define some initial parameters
###################################
default_dir="/Users/mors/Science/TDS-red"
os.chdir(default_dir)

w_dir=os.getcwd()
#data_dir=None
data_dir="/Volumes/iData/Science/KGO_data/LS/Raw/TDS20200223/"#"/Volumes/iData/Science/KGO_data/LS/Raw/TDS20200224"

version ='dev-0.1'
standard_dir = "/Users/mors/Science/standards/" # Directory containing information aboud SED standards
standard_data=os.path.join(standard_dir,'data/')

conf.max_lines=-1
conf.max_width=-1






def do_Display():
    if show_summary_mode.get() == 'Red':
        cur_list = red_list
        cur_dir=w_dir
        file_ind = red_list_header['Filename']
        dir_ind = red_list_header['Directory']
    else:
        cur_dir=data_dir
        file_ind=-2
        dir_ind = -3
        if show_datatype.get() == "obj":
            cur_list = obj_list
        elif show_datatype.get() == "flat":
            cur_list = flat_list
        elif show_datatype.get() == "bias":
            cur_list = bias_list
        elif show_datatype.get() == "star":
            cur_list = standard_list
        elif show_datatype.get() == "neon":
            cur_list = neon_list
        elif show_datatype.get() == "dark":
            cur_list = dark_list
        elif show_datatype.get() == "sunsky":
            cur_list = sunsky_list
        elif show_datatype.get() == "other":
            cur_list = other_list
    selection = cur_list.tree.selection()
    if len(selection) > 0:
        selected_file = cur_list.data[cur_list.tree.index(selection[0])][file_ind]
        selected_dir = cur_list.data[cur_list.tree.index(selection[0])][dir_ind]
        display.load_image(filename=os.path.join(cur_dir,selected_dir,selected_file), label=selected_file)



    global mode
    if show_summary_mode.get() == 'Red':
        cur_list=red_list
    else:
        if show_datatype.get() == "obj":
            cur_list = obj_list
        elif show_datatype.get() == "flat":
            cur_list = flat_list
        elif show_datatype.get() == "bias":
            cur_list = bias_list
        elif show_datatype.get() == "stand":
            cur_list = standard_list
    selection=cur_list.tree.item(cur_list.tree.focus())['values']
    if len(selection)==0:
        message("Nothing to display")
        return
    if show_summary_mode.get() == 'Red':
        subdir = ''
        curdir=w_dir
    else:
        curdir=data_dir
        subdir = selection[cur_list.header.index("Directory")]
    filename=os.path.join(curdir,subdir,selection[cur_list.header.index("Filename")])


def get_unique(seq, idfun=None):
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
       if marker in seen: continue
       seen[marker] = 1
       result.append(marker)
   return result


def Load_R_Dir():
    global data_dir, raw_files, summary_tab
    if data_dir:
        ini_dir=data_dir
    else:
        ini_dir=w_dir
    loaded = tkFileDialog.Directory(root, initialdir=ini_dir).show()
    if loaded:
        data_dir=loaded
        field_load_raw_dir.delete(0,'end')
        field_load_raw_dir.insert(0,data_dir)
        raw_files, summary_tab = do_ScanDir(data_dir)
        if summary_tab is not None:
            show_summary_mode.set("Raw")
            ChangeSummaryType()

def Set_R_Dir(ev):
    global data_dir, raw_files, summary_tab
    loaded=field_load_raw_dir.get()
    if loaded:
        if os.path.isdir(loaded):
            global data_dir,raw_files,summary_tab,standard_dir
            data_dir = loaded
            field_load_raw_dir.delete(0, 'end')
            field_load_raw_dir.insert(0, data_dir)
            raw_files, summary_tab = do_ScanDir(data_dir)
            if summary_tab is not None:
                show_summary_mode.set("Raw")
                ChangeSummaryType()

def Load_Log():
    pass

def Save_Log():
    pass


def Load_W_Dir():
    global w_dir,summary_tab_red
    loaded = tkFileDialog.Directory(root, initialdir=w_dir).show()
    if loaded:
        # if loaded == w_dir and summary_tab_red is not None:
        #     field_load_w_dir.delete(0, 'end')
        #     field_load_w_dir.insert(0, w_dir)
        #     return
        w_dir=loaded
        field_load_w_dir.delete(0, 'end')
        field_load_w_dir.insert(0,w_dir)
        if os.path.isdir(w_dir):
            summary_tab_red = do_ScanRedDir(w_dir)
            if summary_tab_red is not None:
                show_summary_mode.set("Red")
                ChangeSummaryType()

def Set_W_Dir(ev):
    global w_dir,summary_tab_red
    if field_load_w_dir.get() == '':
        field_load_w_dir.insert(0, w_dir)
        return
    else:
        selected_dir=field_load_w_dir.get()
        # if selected_dir == w_dir and summary_tab_red is not None:
        #     return
        w_dir=selected_dir
        if os.path.isdir(w_dir):
            summary_tab_red = do_ScanRedDir(w_dir)
            if summary_tab_red is not None:
                show_summary_mode.set("Red")
                ChangeSummaryType()


def do_EditData():
    global summary_tab_red, summary_tab,red_list

    if show_summary_mode.get() == 'Red':
        values_to_edit={'slit':edit_params_slit.get(),'binning':edit_params_binning.get(),
                        'disperser':edit_params_disperser.get()}
        cur_list = red_list
        cur_summary_tab=summary_tab_red
    else:
        values_to_edit = {'slit': edit_params_slit.get(), 'binning': edit_params_binning.get(),'type':edit_params_type.get(),
                          'disperser': edit_params_disperser.get(),'pa': edit_params_pa.get(),'name':edit_params_name.get()}
        cur_summary_tab=summary_tab
        if show_datatype.get() == "obj":
            cur_list = obj_list
        elif show_datatype.get() == "flat":
            cur_list = flat_list
        elif show_datatype.get() == "bias":
            cur_list = bias_list
        elif show_datatype.get() == "star":
            cur_list = standard_list
        elif show_datatype.get() == "neon":
            cur_list = neon_list
        elif show_datatype.get() == "dark":
            cur_list = dark_list
        elif show_datatype.get() == "sunsky":
            cur_list = sunsky_list
        elif show_datatype.get() == "other":
            cur_list = other_list
    for field in [edit_params_pa,edit_params_slit,edit_params_name,edit_params_binning,edit_params_disperser,edit_params_type]:
        field.delete(0,'end')
    if any(values_to_edit.values()):
        selection = cur_list.tree.selection()
        if len(selection) > 0:
            items=[]
            id_in_summary=[]
            for s in selection:
                if show_summary_mode.get()=="Red":
                    selected_file = cur_list.data[cur_list.tree.index(s)][red_list_header['Filename']]
                    selected_dir = cur_list.data[cur_list.tree.index(s)][red_list_header['Directory']]
                else:
                    selected_file = cur_list.data[cur_list.tree.index(s)][-2]
                    selected_dir = cur_list.data[cur_list.tree.index(s)][-3]
                items.append(cur_list.tree.index(s))
                id_in_summary.append((cur_summary_tab['filename'] == selected_file) & (cur_summary_tab['subdirectory'] == selected_dir))
            for id in id_in_summary:
                if show_summary_mode.get()=="Red":
                    old_file=str(cur_summary_tab['filename'][id][0])
                    old_file_prefix=re.findall(r"(.+)_s\d\.?\d?",old_file)[0]
                    old_file_slit = cur_summary_tab['slit'][id][0]
                    old_file_disperser = cur_summary_tab['disperser'][id][0]
                    old_file_binning = cur_summary_tab['binning'][id][0]
                    old_file_suffix = re.findall(r".+_{}(.+)".format(old_file_binning), old_file)[0]
                for k in values_to_edit:
                    if values_to_edit[k]:
                        cur_summary_tab[k][id] = values_to_edit[k]
                if show_summary_mode.get()=="Red":
                    new_file=old_file_prefix
                    if cur_summary_tab['slit'][id]:
                        new_file="{}_s{}".format(new_file,cur_summary_tab['slit'][id][0])
                    if cur_summary_tab['disperser'][id]:
                        new_file="{}_{}".format(new_file,cur_summary_tab['disperser'][id][0])
                    if cur_summary_tab['binning'][id]:
                        new_file="{}_{}".format(new_file,cur_summary_tab['binning'][id][0])
                    new_file="{}{}".format(new_file,old_file_suffix)
                    cur_summary_tab['filename'][id]=new_file
                    if old_file != new_file:
                        os.replace(os.path.join(w_dir,cur_summary_tab['subdirectory'][id][0],old_file),
                                   os.path.join(w_dir,cur_summary_tab['subdirectory'][id][0],new_file))
            if show_summary_mode.get()=="Red":
                do_RedSummary(summary_tab_red)
            else:
                do_Summary(summary_tab)
            child_id = cur_list.tree.get_children()
            for ind, id in enumerate(items):
                if ind == 0:
                    cur_list.tree.selection_set(child_id[id])
                else:
                    cur_list.tree.selection_add(child_id[id])
                cur_list.tree.focus(child_id[id])



def ChangeSummaryType():
    global show_summary_mode,show_datatype,summary_tab_red,red_list
    if show_summary_mode.get() == "Red":
        do_RedSummary(summary_tab_red)
        frame_select_red_datatype.tkraise()
        frame_table_summary['red'].tkraise()
        # do_process_but.grid_remove()
        if (show_red_datatype.get() in ["calibration","auxiliary"]) and (red_list.redsubtype.get() == "data"):
            starobj_selector_button.config(text="Use as Defaul")
            starobj_selector_button.grid()
        else:
            starobj_selector_button.grid_remove()
    else:
        frame_select_datatype.tkraise()
        frame_table_summary[show_datatype.get()].tkraise()
        do_process_but.grid()
        if show_datatype.get() == 'obj':
            starobj_selector_button.config(text="Use as Standard")
            starobj_selector_button.grid()
        elif show_datatype.get() == 'star':
            starobj_selector_button.config(text="Use as Object ")
            starobj_selector_button.grid()
        else:
            starobj_selector_button.grid_remove()
        #do_process_but.config(state='normal')


# def SaveFile(ev):
#     fn = tkFileDialog.SaveAs(root, filetypes=[('*.txt files', '.txt')]).show()
#     if fn == '':
#         return
#     if not fn.endswith(".txt"):
#         fn += ".txt"
#     open(fn, 'wt').write(textbox.get('1.0', 'end'))


def get_reduced_data(summary_tab_red, types=['obj'],slit=None, disperser=None, binning=None,qual='.fits',
                 state_in='ini', state_out=None, process_unc=True, w_dir=os.getcwd(), process_numbers=True, get_index=False):
    """
    Select all data from summary_tab_red that correspond to mentioned criteria.
    Return the lists of:
    filein, fileout;
    filein_err, fileout_err (if process_unc is true, else skip);
    filenum (if process_numbers is True, else skip)
    """


    all_states={'ini':'','normalized':'_n','linearized':'_lin',
                'sky-subtracted':'_skyrem','combined':'_tot','calibrated':'_abs'}

    skip_binning, skip_slit, skip_disperser=(False,False,False)
    if binning is None:
        skip_binning=True
    if slit is None:
        skip_slit=True
    if disperser is None:
        skip_disperser=True
    if type(types) != str:
        types = tuple(types)
    mask = [idx for idx, val in enumerate(summary_tab_red['data']) if
            val.startswith(types) and (summary_tab_red['process'][idx] == True) and
            any([(summary_tab_red['binning'][idx] == binning), skip_binning]) and
            any([(summary_tab_red['slit'][idx] == slit), skip_slit]) and
            any([(summary_tab_red['disperser'][idx] == disperser), skip_disperser]) and
            (summary_tab_red['subtype'][idx] == 'data') and (summary_tab_red['state'][idx] == state_in)]
    filein = []
    fileout = []
    if process_unc:
        fileerr = []
        fileerr_out = []
    if process_numbers:
        filenum = []
    if get_index:
        index_out=[]
    if len(mask) > 0:
        if process_numbers:
            mask = [x for _, x in sorted(zip(summary_tab_red['filenum'][mask], mask))]
        for m in mask:
            if os.path.isfile(os.path.join(w_dir,summary_tab_red['subdirectory'][m],
                                            summary_tab_red['filename'][m])):
                filein.append(os.path.join(summary_tab_red['subdirectory'][m],
                                                summary_tab_red['filename'][m]))
                if get_index:
                    index_out.append(m)
                if process_numbers:
                    filenum.append(summary_tab_red['filenum'][m])

                found = re.findall(r"(.+){}\.fi?ts".format(all_states[state_in]), summary_tab_red['filename'][m])
                if state_out is not None:
                    fileout.append(os.path.join(summary_tab_red['subdirectory'][m],"{}{}.fits".format(found[0],
                                                                    all_states[state_out])))
                else:
                    fileout.append(None)
                if process_unc:
                    f="{}{}_err.fits".format(found[0],all_states[state_in])
                    try:
                        ind=list(summary_tab_red['filename']).index(f)
                        if summary_tab_red['process'][ind] == True and os.path.isfile(os.path.join(w_dir,
                                                                                summary_tab_red['subdirectory'][ind], f)):
                            fileerr.append(os.path.join(summary_tab_red['subdirectory'][ind], f))
                            if state_out is not None:
                                fileerr_out.append(os.path.join(summary_tab_red['subdirectory'][ind],
                                                    "{}{}_err.fits".format(found[0],all_states[state_out])))
                            else:
                                fileerr_out.append(None)
                        else:
                            fileerr.append(None)
                            fileerr_out.append(None)
                    except ValueError:
                        fileerr.append(None)
                        fileerr_out.append(None)

    output=[filein,fileout]
    if process_unc:
        output.append(fileerr)
        output.append(fileerr_out)
    if process_numbers:
        output.append(filenum)
    if get_index:
        output.append(index_out)
    return output




def get_all_spec_params(summary_tab, selectors=(), selector_types=(), file_types=()):
    # Return all possible values of parameters set by selectors or presented in summary_tab.
    # Order to check: binning, disperser, slit
    # All other parameters should be checked by yourself!!!
    output = []
    for id_sel, selector in enumerate(selectors):
        if selector.get() != "All":
            output.append(selector.get())
        else:
            mask = [idx for idx, val in enumerate(summary_tab['type']) if (val in file_types) and (summary_tab['process'][idx] == True)]
            output.append(get_unique(summary_tab[selector_types[id_sel]][mask]))
    return output



def get_all_redspec_params(summary_tab_red, file_types=()):
    # Return all possible values of parameters for files marked to process.
    # Order to check: binning, disperser, slit
    # All other parameters should be checked by yourself!!!
    output = []
    for selector_type in ['binning','disperser','slit']:
        if type(file_types) !=str:
            file_types=tuple(file_types)
        mask = [idx for idx, val in enumerate(summary_tab_red['data']) if (val.startswith(file_types)
                        and (summary_tab_red['process'][idx] == True))]

        output.append(get_unique(summary_tab_red[selector_type][mask]))
    return output


def get_all_reduced_params(summary_tab, selectors=(), selector_types=(), file_types=(), w_dir=os.getcwd()):
    all_fits = []
    for x in os.walk(w_dir, followlinks=True):
        all_fits.extend([f for f in x[2] if (f.lower().endswith(".fits") or f.lower().endswith(".fts"))])

    binnings=[]
    slits=[]
    dispersers=[]
    for f in all_fits:
        for ftype in file_types:
            if not f.lower().startswith(ftype):
                continue
            slit_find = re.findall(r"{}_s(\d\.?\d?)_".format(ftype), f)
            if len(slit_find) > 0:
                slits.append(slit_find[0])
            disp_find = re.findall(r"{}_s\d\.?\d?_([\w@]+)_\dx\d".format(ftype), f)
            if len(disp_find) > 0:
                dispersers.append(disp_find[0])
            bin_find = re.findall(r"{}_s\d\.?\d?_[\w@]+_(\dx\d)".format(ftype), f)
            if len(bin_find) > 0:
                binnings.append(bin_find[0])
    return (get_unique(binnings),get_unique(dispersers),get_unique(slits))

def search_files_in_dir(w_dir, rootname=None, exclude=["None"], suffix='', suffix_out=''):
    if not rootname:
        return (None,None,None,None,None)
    if exclude[0] == "All":
        exclude=["_abs","_lin","_n","_skyrem", "_rskyrem","_tot","_tot-sky"] # all types except initial frames
    all_fits = []
    filein=[]
    fileerr=[]
    fileout = []
    fileerr_out = []
    filenum= []
    for x in os.walk(w_dir, followlinks=True):
        all_fits.extend([f for f in x[2] if (f.startswith(rootname) and
                                             (f.lower().endswith("{}.fits".format(suffix)) or
                                              f.lower().endswith("{}.fts".format(suffix))))])
    for f in all_fits:
        found=re.findall(r"(.+){}\.fi?ts".format(suffix), f)
        if len(found)>0:
            if ("crmask" in found[0]) or found[0].endswith(tuple(exclude)) or \
                    ("test" in found[0]) or ('err' in found[0]):
                continue
            found_num=re.findall(r"_(\d\d)",found[0])
            if len(found_num)>0:
                filenum.append(found_num[-1])
            filein.append(f)
            fileout.append("{}{}.fits".format(found[0],suffix_out))
            fileerr.append(None)
            fileerr_out.append(None)
            ferr=["{}{}_err{}".format(found[0],suffix,qual) for qual in ['.fits','.fts']]
            for fe in ferr:
                if os.path.isfile(os.path.join(w_dir,fe)):
                    fileerr[-1] = fe
                    fileerr_out[-1] = "{}{}_err.fits".format(found[0],suffix_out)
    sort_ind=np.argsort(np.array(filenum))
    return (np.array(filein)[sort_ind],np.array(fileout)[sort_ind],np.array(fileerr)[sort_ind],
            np.array(fileerr_out)[sort_ind], np.array(filenum)[sort_ind])





def do_Reduction():
    global data_dir, raw_files, summary_tab, mode, w_dir,summary_tab_red
    select_state=list(steps_buttons.state())
    for index,option in enumerate(steps.items()):
        option[1][1]=select_state[index]



    if any(select_state):
        if summary_tab_red is None or w_dir != field_load_w_dir.get():
            Set_W_Dir(None)
            ### Check whether w_dir exists. Create if not
            try:
                os.makedirs(w_dir)
                message("Working directory '{:s}' is created".format(w_dir))
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise

    if steps["cre_bias"][1] and summary_tab is not None:
        #======== Creates meanbias ===========
        binnings, = get_all_spec_params(summary_tab, selectors=(selected_binning,), selector_types=("binning",),
                                        file_types=['bias'])
        for cur_binning in binnings:
            if mode.get() == "TDS":
                filein_r = []
                filein_b = []
                mask = [idx for idx, val in enumerate(summary_tab['type']) if (val == 'bias') and (summary_tab['process'][idx] == True) and
                        (summary_tab['binning'][idx] == cur_binning) and (summary_tab['disperser'][idx] == "R")]
                for m in mask:
                    filein_r.append(os.path.join(summary_tab['subdirectory'][m],summary_tab['filename'][m]))
                mask = [idx for idx, val in enumerate(summary_tab['type']) if (val == 'bias') and (summary_tab['process'][idx] == True) and
                          (summary_tab['binning'][idx] == cur_binning) and ((summary_tab['disperser'][idx] == "G") or
                                                                            (summary_tab['disperser'][idx] == "B"))]
                for m in mask:
                    filein_b.append(os.path.join(summary_tab['subdirectory'][m],summary_tab['filename'][m]))
                filein = (filein_r,filein_b)
                fileout = ("meanbias_r_{}.fits".format(cur_binning), "meanbias_b_{}.fits".format(cur_binning))
                #filenoise = ("readnoise_r_{}.fits".format(cur_binning), "readnoise_b_{}.fits".format(cur_binning))
            else:
                mask = [idx for idx, val in enumerate(summary_tab['type']) if ((val == 'bias') and (summary_tab['process'][idx] == True) and
                        (summary_tab['binning'][idx] == cur_binning))]

                filein=[]
                fileout="meanbias_{}.fits".format(cur_binning)
                #filenoise = "readnoise_{}.fits".format(cur_binning)
                for m in mask:
                    filein.append(os.path.join(summary_tab['subdirectory'][m], summary_tab['filename'][m]))
            status=reduction.cre_bias_or_dark(filein,fileout=fileout,w_dir=w_dir,raw_dir=data_dir, mode=mode.get())#, filenoise=filenoise)
            if status:
                if mode.get() != "TDS":
                    fileout=[fileout]
                do_update_RedSummary(add=fileout)
                message("Bias frames were prepared for binning {}".format(cur_binning))
            else:
                message("Something goes wrong with preparing of bias frame for binning {}".format(cur_binning))

    if steps['cre_dark'][1] and summary_tab is not None:
        # ======== Creates dark frames ===========
        binnings, = get_all_spec_params(summary_tab, selectors=(selected_binning,),selector_types=("binning",),file_types=['dark'])
        cur_overscan = None
        fields_check=[field_cutframe_x0.get(), field_cutframe_x1.get(), field_cutframe_y0.get(), field_cutframe_y1.get()]
        if cut_as_overcan_val.get() and any(fields_check):
            cur_overscan=[None,None,None,None]
            for ind, field in enumerate(fields_check):
                if field != "":
                    cur_overscan[ind] = int(field)

        for cur_binning in binnings:
            if mode.get() == "TDS":
                filein_r = []
                filein_b = []
                mask = [idx for idx, val in enumerate(summary_tab['type']) if (val == 'dark') and (summary_tab['process'][idx] == True) and
                        (summary_tab['binning'][idx] == cur_binning) and (summary_tab['disperser'][idx] == "R")]
                for m in mask:
                    filein_r.append(os.path.join(summary_tab['subdirectory'][m], summary_tab['filename'][m]))
                mask = [idx for idx, val in enumerate(summary_tab['type']) if (val == 'dark') and (summary_tab['process'][idx] == True) and
                        (summary_tab['binning'][idx] == cur_binning) and ((summary_tab['disperser'][idx] == "G") or
                                                                          (summary_tab['disperser'][idx] == "B"))]
                for m in mask:
                    filein_b.append(os.path.join(summary_tab['subdirectory'][m], summary_tab['filename'][m]))
                filein = (filein_r, filein_b)
                fileout = ("dark_r_{}.fits".format(cur_binning),
                           "dark_b_{}.fits".format(cur_binning))
                #filebias = ("meanbias_r_{}.fits".format(cur_binning), "meanbias_b_{}.fits".format(cur_binning))
                mask = [idx for idx, val in enumerate(summary_tab_red['data']) if
                        (val == 'meanbias_r') and (summary_tab_red['process'][idx] == True) and
                        (summary_tab_red['binning'][idx] == cur_binning)]
                if len(mask)==1:
                    filebias_r=os.path.join(summary_tab_red['subdirectory'][mask[0]],
                                            summary_tab_red['filename'][mask[0]])
                else: filebias_r=None
                mask = [idx for idx, val in enumerate(summary_tab_red['data']) if
                        (val == 'meanbias_b') and (summary_tab_red['process'][idx] == True) and
                        (summary_tab_red['binning'][idx] == cur_binning)]

                if len(mask) == 1:
                    filebias_b = os.path.join(summary_tab_red['subdirectory'][mask[0]],
                                              summary_tab_red['filename'][mask[0]])
                else: filebias_b=None
                filebias=(filebias_r,filebias_b)
            else:
                mask = [idx for idx, val in enumerate(summary_tab['type']) if ((val == 'dark') and (summary_tab['process'][idx] == True) and
                                                                               (summary_tab['binning'][
                                                                                    idx] == cur_binning))]
                filein = []
                fileout = "dark_{}.fits".format(cur_binning)
                for m in mask:
                    filein.append(os.path.join(summary_tab['subdirectory'][m], summary_tab['filename'][m]))
                # filebias = "meanbias_{}.fits".format(cur_binning)
                mask = [idx for idx, val in enumerate(summary_tab_red['data']) if
                        (val == 'meanbias') and (summary_tab_red['process'][idx] == True) and
                        (summary_tab_red['binning'][idx] == cur_binning)]

                if len(mask) == 1:
                    filebias = os.path.join(summary_tab_red['subdirectory'][mask[0]],
                                              summary_tab_red['filename'][mask[0]])
                else: filebias=None

            status=reduction.cre_bias_or_dark(filein,fileout=fileout,w_dir=w_dir,raw_dir=data_dir, mode=mode.get(), filebias=filebias, overscan=cur_overscan)
            if status:
                if mode.get() != "TDS":
                    fileout=[fileout]
                do_update_RedSummary(add=fileout)
                message("Dark frames were prepared for binning {}".format(cur_binning))
            else:
                message("Something goes wrong with preparing of dark frame for binning {}".format(cur_binning))

    if steps['cre_calibration'][1] and summary_tab is not None:
        seek_types=['flat', 'sunsky', 'neon', '13dots']
        binnings, dispersers, slits = get_all_spec_params(summary_tab, selectors=(
            selected_binning, selected_disperser, selected_slit),
                                                          selector_types=("binning", "disperser", "slit"),
                                                          file_types=seek_types)
        if auto_gain_val.get() or not float(field_gain.get()):
            gain=None
        else:
            gain=float(field_gain.get())
        if auto_readnoise_val.get() or not float(field_readnoise.get()):
            readnoise=None
        else:
            readnoise=float(field_readnoise.get())
        cut_borders = [None, None, None, None]
        for ind, field in enumerate([field_cutframe_x0, field_cutframe_x1, field_cutframe_y0, field_cutframe_y1]):
            if field.get() != "":
                cut_borders[ind] = int(field.get())
        cut_is_overscan = cut_as_overcan_val.get()

        for cur_binning in binnings:
            for cur_disperser in dispersers:
                for cur_slit in slits:
                    for cur_type in seek_types:
                        filein = [os.path.join(summary_tab['subdirectory'][idx],summary_tab['filename'][idx]) for idx, val in
                                  enumerate(summary_tab['type']) if (val == cur_type) and (summary_tab['process'][idx] == True) and
                            summary_tab['disperser'][idx] == cur_disperser and summary_tab['binning'][idx] == cur_binning and
                                summary_tab['slit'][idx] == cur_slit]

                        fileout = "{}_s{}_{}_{}.fits".format(cur_type,cur_slit,cur_disperser,cur_binning)
                        file_error_out = "{}_s{}_{}_{}_err.fits".format(cur_type,cur_slit, cur_disperser, cur_binning)
                        if mode.get() == "TDS":
                            if cur_disperser == "G" or cur_disperser == "B":
                                add_mode="_b"
                            else:
                                add_mode="_r"
                        else:
                            add_mode=""

                        mask = [idx for idx, val in enumerate(summary_tab_red['data']) if
                                ((val == 'meanbias' and mode.get()!="TDS") or
                                 (val == 'meanbias_r' and mode.get()=="TDS" and cur_disperser=="R") or
                                 (val == 'meanbias_b' and mode.get()=="TDS" and (cur_disperser=="B" or cur_disperser=="G"))) and (summary_tab_red['process'][idx] == True) and
                                (summary_tab_red['binning'][idx] == cur_binning)]
                        if len(mask) == 1:
                            filebias = os.path.join(summary_tab_red['subdirectory'][mask[0]],
                                                    summary_tab_red['filename'][mask[0]])
                        else:
                            filebias = None

                        mask = [idx for idx, val in enumerate(summary_tab_red['data']) if
                                ((val == 'dark' and mode.get() != "TDS") or
                                 (val == 'dark_r' and mode.get() == "TDS" and cur_disperser == "R") or
                                 (val == 'dark_b' and mode.get() == "TDS" and (
                                             cur_disperser == "B" or cur_disperser == "G"))) and (
                                            summary_tab_red['process'][idx] == True) and
                                (summary_tab_red['binning'][idx] == cur_binning)]
                        if len(mask) == 1:
                            filedark = os.path.join(summary_tab_red['subdirectory'][mask[0]],
                                                    summary_tab_red['filename'][mask[0]])
                        else:
                            filedark = None
                        status = reduction.cre_calibration(filein,fileout=fileout,w_dir=w_dir,raw_dir=data_dir,mode="{}{}".format(mode.get(),add_mode),
                                                      filebias=filebias,filedark=filedark, file_error_out=file_error_out,gain=gain,
                                                      imgtype = cur_type, readnoise=readnoise, cut=cut_borders, overscan=cut_is_overscan)
                        if status:
                            do_update_RedSummary(add=[fileout,file_error_out])
                            message("{} frames were prepared for slit {}, disperser {} and binning {}".format(cur_type, cur_slit,cur_disperser,cur_binning))
                        else:
                            message("Something goes wrong with preparing {} for slit {}, disperser {} and binning {}".format(cur_type, cur_slit,cur_disperser,cur_binning))

    if steps['cre_ini'][1] and summary_tab is not None:
        seek_types = ['obj', 'star']
        if auto_gain_val.get() or not float(field_gain.get()):
            gain=None
        else:
            gain=float(field_gain.get())
        if auto_readnoise_val.get() or not float(field_readnoise.get()):
            readnoise=None
        else:
            readnoise=float(field_readnoise.get())

        cut_borders=[None,None,None,None]
        for ind,field in enumerate([field_cutframe_x0,field_cutframe_x1,field_cutframe_y0,field_cutframe_y1]):
            if field.get()!="":
                cut_borders[ind]=int(field.get())
        cut_is_overscan=cut_as_overcan_val.get()

        binnings, dispersers, slits = get_all_spec_params(summary_tab, selectors=(
            selected_binning, selected_disperser, selected_slit),
                                                          selector_types=("binning", "disperser", "slit"),
                                                          file_types=seek_types)
        for cur_binning in binnings:
            for cur_disperser in dispersers:
                for cur_slit in slits:
                    for cur_type in seek_types:
                        filein = [os.path.join(summary_tab['subdirectory'][idx],summary_tab['filename'][idx]) for idx, val in
                                  enumerate(summary_tab['type']) if (val == cur_type) and (summary_tab['process'][idx] == True) and
                            summary_tab['disperser'][idx] == cur_disperser and summary_tab['binning'][idx] == cur_binning and
                                summary_tab['slit'][idx] == cur_slit]
                        if len(filein)>0:
                            message("######==== Start processing of {} data.....".format(cur_type))
                        fileout_prefix = "{}_s{}_{}_{}".format(cur_type,cur_slit,cur_disperser,cur_binning)
                        fileout_suffix = ""
                        errorout_suffix = "_err"
                        if mode.get() == "TDS":
                            if cur_disperser == "G" or cur_disperser == "B":
                                add_mode = "_b"
                            else:
                                add_mode = "_r"
                        else:
                            add_mode = ""

                        mask = [idx for idx, val in enumerate(summary_tab_red['data']) if
                                ((val == 'meanbias' and mode.get() != "TDS") or
                                 (val == 'meanbias_r' and mode.get() == "TDS" and cur_disperser == "R") or
                                 (val == 'meanbias_b' and mode.get() == "TDS" and (
                                             cur_disperser == "B" or cur_disperser == "G"))) and (
                                            summary_tab_red['process'][idx] == True) and
                                (summary_tab_red['binning'][idx] == cur_binning)]
                        if len(mask) == 1:
                            filebias = os.path.join(summary_tab_red['subdirectory'][mask[0]],
                                                    summary_tab_red['filename'][mask[0]])
                        else:
                            filebias = None

                        mask = [idx for idx, val in enumerate(summary_tab_red['data']) if
                                ((val == 'dark' and mode.get() != "TDS") or
                                 (val == 'dark_r' and mode.get() == "TDS" and cur_disperser == "R") or
                                 (val == 'dark_b' and mode.get() == "TDS" and (
                                         cur_disperser == "B" or cur_disperser == "G"))) and (
                                        summary_tab_red['process'][idx] == True) and
                                (summary_tab_red['binning'][idx] == cur_binning)]
                        if len(mask) == 1:
                            filedark = os.path.join(summary_tab_red['subdirectory'][mask[0]],
                                                    summary_tab_red['filename'][mask[0]])
                        else:
                            filedark = None

                        file_to_del, _, fileerr_to_del,_,ind_to_del = get_reduced_data(summary_tab_red, types=[cur_type],
                                                                 slit=cur_slit, disperser=cur_disperser,
                                                                 binning=cur_binning, state_in='ini', w_dir=w_dir,
                                                                 get_index=True,process_unc=True,process_numbers=False)
                        if len(ind_to_del) > 0:
                            for f in file_to_del+fileerr_to_del:
                                if f is not None:
                                    os.remove(os.path.join(w_dir,f))
                            summary_tab_red.remove_rows(ind_to_del)




                        status,all_prepared_files = reduction.cre_ini(filein,fileout_prefix=fileout_prefix, fileout_suffix=fileout_suffix,w_dir=w_dir,
                                        raw_dir=data_dir,mode="{}{}".format(mode.get(),add_mode),filebias=filebias,
                                        cr_clean=use_lacosmic_val.get(),filedark=filedark,cut=cut_borders,
                                        overscan=cut_is_overscan,
                                        errorout_suffix=errorout_suffix,gain=gain, readnoise=readnoise)
                        if not status:
                            message("No any raw frames was found => cannot contunue.")
                        else:
                            do_update_RedSummary(add=all_prepared_files)
                            message("Initial frames were prepared")

    if summary_tab_red is None and any([steps['calc_geometry'][1],steps['norm_flat'][1],
                                       steps['calc_shifts'][1], steps['calc_dispersion'][1],
                                       steps['lin_and_align'][1],steps['do_transform'][1],
                                       steps['sky_subtract'][1], steps['combine'][1],
                                       steps['calc_dqe'][1],steps['flux_cal'][1]]):
        message("Cannot proceed selected steps.")
        return

    if steps['calc_geometry'][1]:
        seek_types=['neon','13dots']
        binnings, dispersers, slits = get_all_redspec_params(summary_tab_red,file_types=seek_types)
        for cur_binning in binnings:
            for cur_disperser in dispersers:
                for cur_slit in slits:
                    if geometry_trace_star_val.get():
                        #=== Test slope from object and star frames
                        filein_star, _ = get_reduced_data(summary_tab_red, types=['star'],
                        slit=cur_slit,disperser=cur_disperser, binning=cur_binning, state_in='ini',
                                        w_dir=w_dir, process_numbers=False, process_unc=False)
                        if len(filein_star) == 0:
                            filein_star=None
                    else:
                        filein_star=None

                    cur_type='neon'
                    filein,_,fileerr,_=get_reduced_data(summary_tab_red,types=[cur_type],
                                        slit=cur_slit,disperser=cur_disperser,
                                        binning=cur_binning, state_in='ini',w_dir=w_dir,
                                                  process_numbers=False)
                    if len(filein) == 0:
                        continue
                    if mode.get()=='SAO':
                        dotmask, _ = get_reduced_data(summary_tab_red, types=['13dots'],
                                                            slit=cur_slit, disperser=cur_disperser,
                                                            binning=cur_binning, state_in='ini', w_dir=w_dir,
                                                            process_numbers=False, process_unc=False)
                        if len(dotmask) == 0:
                                dotmask = [None]
                    else:
                        dotmask=[None]
                    fileout_transform = "transform_s{}_{}_{}.fits".format(cur_slit, cur_disperser, cur_binning)
                    if geometry_do_test_val.get():
                        save_test = "neon_s{}_{}_{}_warp_test.fits".format(cur_slit, cur_disperser, cur_binning)
                    else:
                        save_test = None
                    if geometry_prominance_1.get() == 0:
                        prom=None
                    else:
                        prom = float(geometry_prominance_1.get())
                    status = reduction.calc_geometry(filein[0], fileerr=fileerr[0],
                                            fileout_transform=fileout_transform, w_dir=w_dir,
                                            do_plot=geometry_plot_results_val.get(), dotmask=dotmask[0],
                                            nstripes=int(geometry_nstripes.get()), win=int(geometry_win.get()),
                                            snr=float(geometry_snr.get()), peaks_mindist=float(geometry_peakdist.get()),
                                            poly_order=int(geometry_poly_lines.get()),
                                            prominence=[float(geometry_prominance_0.get()), prom],
                                            oversample=50,trace_star=filein_star, save_test=save_test)
                    if status:
                        if not dotmask:
                            add_mess = 'file'
                        else:
                            add_mess="and 13dots files"
                        message("Geometry model was constructed using neon {} with binning {}, disperser {} and slit {}".format(
                            add_mess,cur_binning,cur_disperser,cur_slit))
                        do_update_RedSummary(add=[fileout_transform,save_test])

                    else:
                        message("Something went wrong with calculation of geometry for binning {}, disperser {} and slit {}".format(
                            cur_binning,cur_disperser,cur_slit))


    if steps['cre_norm_flat'][1]:
        message("Start creating normalized flat")
        seek_types = ['flat','sunsky']
        binnings, dispersers, slits = get_all_redspec_params(summary_tab_red,file_types=seek_types)
        for cur_binning in binnings:
            for cur_disperser in dispersers:
                if os.path.isfile(os.path.join(w_dir,default_transform.get())):
                    cur_default_disperser=re.findall(r"\S+_s\d?\.?\d?_([\w@]+)_\d+x\d+",default_transform.get())
                    cur_default_bin = re.findall(r"\S+_s\d?\.?\d?_[\w@]+_(\d+x\d+)", default_transform.get())
                    if len(cur_default_disperser)==0:
                        cur_default_disperser=None
                    if len(cur_default_bin) == 0:
                        cur_default_bin = None
                    if cur_binning == cur_default_bin and cur_disperser == cur_default_disperser:
                        cur_default_transform=default_transform.get()
                    else:
                        cur_default_transform = None
                cur_default_transform=None
                for cur_slit in slits:
                    if default_transform_force.get() and cur_default_transform is not None:
                        filein_transform=[cur_default_transform]
                    else:
                        filein_transform, _, = get_reduced_data(summary_tab_red, types=['transform'],
                                                             slit=cur_slit, disperser=cur_disperser,
                                                             binning=cur_binning, state_in='ini', w_dir=w_dir,
                                                             process_numbers=False, process_unc=False)
                    # filein_transform = "transform_s{}_{}_{}.fits".format(cur_slit, cur_disperser, cur_binning)
                    # if not os.path.isfile(os.path.join(w_dir, filein_transform)):
                    if len(filein_transform) == 0 and cur_default_transform is None:
                        message("!!!! Can't find transformation matrix computed at previous stages. Check this!",
                                noheader=True)
                        filein_transform = [None]
                    elif len(filein_transform) == 0:
                        filein_transform=[cur_default_transform]
                    nosmo=False
                    for cur_type in seek_types:
                        filein, _,fileerr,_ = get_reduced_data(summary_tab_red, types=[cur_type],
                                                                slit=cur_slit, disperser=cur_disperser,
                                                                binning=cur_binning, state_in='ini', w_dir=w_dir,
                                                                process_numbers=False)
                        file_normflat = "flat_s{}_{}_{}_n.fits".format(cur_slit, cur_disperser, cur_binning)
                        file_normflat_err = "flat_s{}_{}_{}_n_err.fits".format(cur_slit, cur_disperser, cur_binning)
                        if len(filein)==0:
                            message(".... Cannot find {} for disperser {}, binning {} and slit {}".format(cur_type,
                                                            cur_disperser,cur_binning,cur_slit),noheader=True)
                        else:
                            if cur_type == 'sunsky':
                                message(".... Use sunsky as flat for disperser {}, binning {} and slit {}".format(cur_disperser,
                                                                                    cur_binning,cur_slit), noheader=True)
                                nosmo=True
                            break
                    else:
                        message("!!!! Cannot proceed with normalized flat creation for disperser {}, binning {} and slit {}".format(cur_disperser,
                                                                                                              cur_binning,
                                                                                                              cur_slit), noheader=True)
                        continue
                    #if cur_binning != '2x4':
                    status=reduction.calc_norm_flat(filein[0], fileerr=fileerr[0], fileout=file_normflat, fileout_err=file_normflat_err, nosmo=nosmo,
                                                   filein_transform=filein_transform[0], w_dir=w_dir, do_plot=True, blue_pix_limit=None)
                    # else:
                    #     status=True
                    if not status:
                        message(
                            "Something went wrong with construction of normalized flat for binning {}, disperser {} and slit {}".format(
                                cur_binning, cur_disperser, cur_slit))
                    else:
                        do_update_RedSummary(add=[file_normflat,file_normflat_err])
        message("Creating of normalized flat finished.")

    if steps['norm_flat'][1]:
        seek_types = ["obj",'neon','star','sunsky']
        binnings, dispersers, slits = get_all_redspec_params(summary_tab_red, file_types=seek_types)
        for cur_binning in binnings:
            for cur_disperser in dispersers:
                if os.path.isfile(os.path.join(w_dir,default_flat.get())):
                    cur_default_disperser=re.findall(r"\S+_s\d?\.?\d?_([\w@]+)_\d+x\d+",default_flat.get())
                    cur_default_bin = re.findall(r"\S+_s\d?\.?\d?_[\w@]+_(\d+x\d+)", default_flat.get())
                    if len(cur_default_disperser)==0:
                        cur_default_disperser=None
                    if len(cur_default_bin) == 0:
                        cur_default_bin = None
                    if cur_binning == cur_default_bin[0] and cur_disperser == cur_default_disperser[0]:
                        cur_default_flat=default_flat.get()
                        cur_default_flat_err = "{}_err.fits".format(re.findall(r"(\S+).fits", default_flat.get())[0])
                        if not os.path.isfile(os.path.join(w_dir,cur_default_flat_err)):
                            cur_default_flat_err=None
                    else:
                        cur_default_flat = None
                        cur_default_flat_err = None
                else:
                    cur_default_flat=None
                    cur_default_flat_err = None
                for cur_slit in slits:
                    if cur_default_flat is not None and default_flat_force.get():
                        file_normflat=[cur_default_flat]
                        file_normflat_err=[cur_default_flat_err]
                    else:
                        file_normflat,_,file_normflat_err,_=get_reduced_data(summary_tab_red, types=['flat'],
                                                    slit=cur_slit, disperser=cur_disperser,process_numbers=False,
                                                    binning=cur_binning, state_in='normalized',
                                                                         w_dir=w_dir)
                        if len(file_normflat)==0 and cur_default_flat is not None:
                            file_normflat = [cur_default_flat]
                            file_normflat_err = [cur_default_flat_err]
                        elif len(file_normflat)==0:
                            message("...Cannot normalize files: norm.flat is not available",noheader=True)
                            continue
                    for cur_type in seek_types:
                        message("Normalize {} frame to flat".format(cur_type))
                        if cur_type not in ['sunsky','neon']:
                            filein, fileout, fileerr, fileerr_out,_ = get_reduced_data(summary_tab_red, types=[cur_type],
                                                                         slit=cur_slit, disperser=cur_disperser,
                                                                         binning=cur_binning, state_in='ini',state_out='normalized',
                                                                         w_dir=w_dir)

                        else:
                            filein, fileout, fileerr, fileerr_out = get_reduced_data(summary_tab_red,
                                        types=[cur_type],slit=cur_slit,disperser=cur_disperser,binning=cur_binning,
                                        state_in='ini',state_out='normalized',w_dir=w_dir,process_numbers=False)
                        if len(filein) > 0:
                            file_to_del, _, fileerr_to_del, _, ind_to_del = get_reduced_data(summary_tab_red,
                                                                                                 types=[cur_type],
                                                                                                 slit=cur_slit,
                                                                                                 disperser=cur_disperser,
                                                                                                 binning=cur_binning,
                                                                                                 state_in='normalized',
                                                                                                 w_dir=w_dir,
                                                                                                 get_index=True,
                                                                                                 process_unc=True,
                                                                                                 process_numbers=False)
                            if len(ind_to_del) > 0:
                                for f in file_to_del + fileerr_to_del:
                                    if f is not None:
                                        os.remove(os.path.join(w_dir, f))
                                summary_tab_red.remove_rows(ind_to_del)
                            message("Normalize {} frames to flat from {}".format(cur_type,file_normflat[0]))
                            status = reduction.flat_normilize(filein, file_normflat[0], fileerr=fileerr, flaterr=file_normflat_err[0],
                                                                fileout=fileout, fileerr_out=fileerr_out, w_dir=w_dir)
                            if status:
                                do_update_RedSummary(add=fileout+fileerr_out)
                            else:
                                message("Something went wrong with dividing to flat")
        message("Flat-fielding finished.")

    if steps['calc_shifts'][1]:
        seek_types = ['obj', 'star']
        message("Begin to calculate relative shifts and atmospheric dispersion")
        max_shifts=(float(atmdisp_maxshift_x.get()), float(atmdisp_maxshift_y.get()))
        binnings, dispersers, slits = get_all_redspec_params(summary_tab_red,file_types=seek_types)
        for cur_binning in binnings:
            for cur_disperser in dispersers:
                if os.path.isfile(os.path.join(w_dir,default_transform.get())):
                    cur_default_disperser=re.findall(r"\S+_s\d?\.?\d?_([\w@]+)_\d+x\d+",default_transform.get())
                    cur_default_bin = re.findall(r"\S+_s\d?\.?\d?_[\w@]+_(\d+x\d+)", default_transform.get())
                    if len(cur_default_disperser)==0:
                        cur_default_disperser=None
                    if len(cur_default_bin) == 0:
                        cur_default_bin = None
                    if cur_binning == cur_default_bin[0] and cur_disperser == cur_default_disperser[0]:
                        cur_default_transform=default_transform.get()
                    else:
                        cur_default_transform = None
                else:
                    cur_default_transform=None
                for cur_slit in slits:
                    for cur_type in seek_types:
                        do_shift_x = False
                        do_shift_y=False
                        if default_transform_force.get() and cur_default_transform is not None:
                            filein_transform = [cur_default_transform]
                        else:
                            file_transform, _ = get_reduced_data(summary_tab_red, types=['transform'],
                                                                        slit=cur_slit,disperser=cur_disperser,
                                                                    binning=cur_binning, state_in='ini',
                                                        w_dir=w_dir, process_unc=False,process_numbers=False)
                        if len(file_transform)==0 and cur_default_transform is None:
                            file_transform=[None]
                        elif len(file_transform)==0:
                            file_transform=[cur_default_transform]
                        if cur_type != 'sunsky':
                            if (cur_type == 'obj' and atmdisp_findshift_x_obj_val.get()) or \
                                    (cur_type == 'star' and atmdisp_findshift_x_star_val.get()):
                                do_shift_x=True
                            if (cur_type == 'obj' and atmdisp_findshift_y_obj_val.get()) or \
                                    (cur_type == 'star' and atmdisp_findshift_y_star_val.get()):
                                do_shift_y=True
                            if atmdisp_do_test_val.get():
                                save_test = "{}_s{}_{}_{}_atmdisp_test.fits".format(cur_type, cur_slit, cur_disperser,
                                                                                    cur_binning)
                            else:
                                save_test = None
                            if (cur_type == 'obj'or cur_type == 'star'):
                                atm_disp_out = "atmdisp_{}_s{}_{}_{}.dat".format(cur_type, cur_slit, cur_disperser,
                                                                                 cur_binning)
                                if cur_type=='obj' and int(atmdisp_obj_ypos.get()) !=0:
                                    atm_disp_ypos = int(atmdisp_obj_ypos.get())
                                elif cur_type=='star' and int(atmdisp_star_ypos.get()) !=0:
                                    atm_disp_ypos = int(atmdisp_star_ypos.get())
                                else:
                                    atm_disp_ypos = None
                            else:
                                atm_disp_out = None
                                atm_disp_ypos = None
                        else:
                            atm_disp_out=None
                            save_test = None
                            atm_disp_ypos=None

                        shift_out = "shifts_{}_s{}_{}_{}.dat".format(cur_type, cur_slit, cur_disperser,
                                                                          cur_binning)

                        filein, _, filenum = get_reduced_data(summary_tab_red, types=[cur_type],
                                                             slit=cur_slit, disperser=cur_disperser,
                                                             binning=cur_binning, state_in='normalized',
                                                             w_dir=w_dir, process_unc=False, process_numbers=True)

                        if len(filein) > 0 and file_transform[0] is not None:
                            if atm_disp_out is not None:
                                if os.path.isfile(os.path.join(w_dir,atm_disp_out)):
                                    os.remove(os.path.join(w_dir,atm_disp_out))
                            if not ((cur_type == 'obj' and atmdisp_obj_apply_val.get()) or
                                    (cur_type == 'star' and atmdisp_star_apply_val.get())):
                                atm_disp_out=None
                            message("Calculate relative shifts between different exposures of {} for slit {}, binning {} and disperser {}".format(
                                cur_type,cur_slit,cur_binning,cur_disperser))

                            if atmdisp_ref_index.get() != '':
                                ref_index=int(atmdisp_ref_index.get())
                            else:
                                ref_index=None
                            status=reduction.comp_shifts(filein, file_transform=file_transform[0],ypos_atm_disp=atm_disp_ypos,
                                                  do_shift=(do_shift_x, do_shift_y),max_shift=max_shifts,
                                                  shift_out=shift_out, w_dir=w_dir, file_numbers=filenum, trace_atm_disp=atm_disp_out,
                                                  plot_atm_disp=atmdisp_plot_val.get(),ref_index=ref_index,
                                                  trace_atm_disp_order=int(atmdisp_order.get()), save_test=save_test)
                            if status:
                                do_update_RedSummary(add=[shift_out, atm_disp_out, save_test])
                            else:
                                message("Something went wrong with calculation of relative shifts and atmospheric dispersion")
                        else:
                            message("Cannot proceed atm.disp and shifts correction of {} for slit {}, binning {} and disperser {}".format(
                                cur_type,cur_slit,cur_binning,cur_disperser))
                        # elif len(filein) == 1:
                        #     ascii.write([filenum[0], 0, 0], os.path.join(w_dir,shift_out),
                        #                 names=['filenum', 'dx', 'dy'], overwrite=True)
        message("Calculation of relative shifts and atmospheric dispersion is finished")



    if steps['calc_dispersion'][1]:
        cur_type = 'neon'
        message("Start to compute the dispersion curve")
        binnings, dispersers, slits = get_all_redspec_params(summary_tab_red,file_types=[cur_type])

        for cur_binning in binnings:
            for cur_disperser in dispersers:
                if os.path.isfile(os.path.join(w_dir,default_transform.get())):
                    cur_default_disperser=re.findall(r"\S+_s\d?\.?\d?_([\w@]+)_\d+x\d+",default_transform.get())
                    cur_default_bin = re.findall(r"\S+_s\d?\.?\d?_[\w@]+_(\d+x\d+)", default_transform.get())
                    if len(cur_default_disperser)==0:
                        cur_default_disperser=None
                    if len(cur_default_bin) == 0:
                        cur_default_bin = None
                    if cur_binning == cur_default_bin[0] and cur_disperser == cur_default_disperser[0]:
                        cur_default_transform=default_transform.get()
                    else:
                        cur_default_transform = None
                else:
                    cur_default_transform=None
                for cur_slit in slits:
                    if default_transform_force.get() and cur_default_transform is not None:
                        filein_transform=[cur_default_transform]
                    else:
                        filein_transform, _ = get_reduced_data(summary_tab_red, types=['transform'],
                                                          slit=cur_slit, disperser=cur_disperser,
                                                          binning=cur_binning, state_in='ini',
                                                          w_dir=w_dir, process_unc=False, process_numbers=False)
                    if len(filein_transform) == 0 and cur_default_transform is None:
                        message("!!!! Can't find transformation matrix computed at previous stages. Check this!",noheader=True)
                        filein_transform=[None]
                    elif len(filein_transform) == 0:
                        filein_transform=[cur_default_transform]
                    dispcurve_out= "dispcurve_s{}_{}_{}.fits".format(cur_slit, cur_disperser, cur_binning)
                    filein, _,fileerr,_ = get_reduced_data(summary_tab_red, types=[cur_type],
                                                           slit=cur_slit, disperser=cur_disperser,
                                                           binning=cur_binning, state_in='normalized',
                                                           w_dir=w_dir, process_unc=True, process_numbers=False)
                    if len(filein)>0:
                        status = reduction.calc_dispersion(filein[0], fileerr=fileerr[0], filein_transform=filein_transform[0],
                                                         disper_out=dispcurve_out, w_dir=w_dir,
                                                         do_plot=True,mode=mode.get(),wid_y=int(dispcurve_window_y.get()),
                                                         poly_order=(int(dispcurve_order_lam.get()),int(dispcurve_order_y.get())),
                                                           disperser=cur_disperser, binning = cur_binning)
                        if status:
                            message("Disperson curve calculated for binning {}, disperser {} and slit {}".format(
                                    cur_binning, cur_disperser, cur_slit))
                            do_update_RedSummary(add=[dispcurve_out])
                        else:
                            message("Something went wrong with construction of disperson curve for binning {}, disperser {} and slit {}".format(
                                cur_binning,cur_disperser,cur_slit))
        message("Calculation of dispersion curve is finished")

    if steps['lin_and_align'][1]:
        message("Start to prepare final transformation to linearize spectra")
        seek_types = ['obj', 'star','sunsky']
        binnings, dispersers, slits = get_all_redspec_params(summary_tab_red,file_types=seek_types)

        if linearizarion_disp_auto_val.get():
            dl='auto'
        else:
            dl=float(linearizarion_disp.get())
        if linearizarion_lam0_auto_val.get():
            l0='auto'
        else:
            l0=float(linearizarion_lam0.get())
        if linearizarion_lam1_auto_val.get():
            l1 = 'auto'
        else:
            l1 = float(linearizarion_lam1.get())
        lin_parameters = {"start": l0, 'dlam': dl, 'fin': l1}
        for cur_binning in binnings:
            for cur_disperser in dispersers:
                if os.path.isfile(os.path.join(w_dir,default_transform.get())):
                    cur_default_disperser=re.findall(r"\S+_s\d?\.?\d?_([\w@]+)_\d+x\d+",default_transform.get())
                    cur_default_bin = re.findall(r"\S+_s\d?\.?\d?_[\w@]+_(\d+x\d+)", default_transform.get())
                    if len(cur_default_disperser)==0:
                        cur_default_disperser=None
                    if len(cur_default_bin) == 0:
                        cur_default_bin = None
                    if cur_binning == cur_default_bin[0] and cur_disperser == cur_default_disperser[0]:
                        cur_default_transform=default_transform.get()
                    else:
                        cur_default_transform = None
                else:
                    cur_default_transform=None
                if os.path.isfile(os.path.join(w_dir,default_disp_curve.get())):
                    cur_default_disperser=re.findall(r"\S+_s\d?\.?\d?_([\w@]+)_\d+x\d+",default_disp_curve.get())
                    cur_default_bin = re.findall(r"\S+_s\d?\.?\d?_[\w@]+_(\d+x\d+)", default_disp_curve.get())
                    if len(cur_default_disperser)==0:
                        cur_default_disperser=None
                    if len(cur_default_bin) == 0:
                        cur_default_bin = None
                    if cur_binning == cur_default_bin[0] and cur_disperser == cur_default_disperser[0]:
                        cur_default_disp_curve=default_disp_curve.get()
                    else:
                        cur_default_disp_curve = None
                else:
                    cur_default_disp_curve=None
                for cur_slit in slits:
                    if default_transform_force.get() and cur_default_transform is not None:
                        filein_transform=[cur_default_transform]
                    else:
                        filein_transform, _ = get_reduced_data(summary_tab_red, types=['transform'],
                                                             slit=cur_slit, disperser=cur_disperser,
                                                             binning=cur_binning, state_in='ini',
                                                             w_dir=w_dir, process_unc=False, process_numbers=False)
                    if default_disp_curve_force.get() and cur_default_disp_curve is not None:
                        dispcurve = [cur_default_disp_curve]
                    else:
                        dispcurve, _ = get_reduced_data(summary_tab_red, types=['dispcurve'],
                                                           slit=cur_slit, disperser=cur_disperser,
                                                           binning=cur_binning, state_in='ini',
                                                           w_dir=w_dir, process_unc=False, process_numbers=False)
                    if len(filein_transform) == 0 and cur_default_transform is not None:
                        filein_transform=[cur_default_transform]
                    if len(dispcurve) == 0 and cur_default_disp_curve is not None:
                        dispcurve=[cur_default_disp_curve]
                    if len(dispcurve)>0:
                        if len(filein_transform)==0:
                            message("!!!! Cannot find transformation matrix computed at previous stages for slit {}, binning {}, disperser {}  => cannot continue!".format(
                                cur_slit,cur_binning,cur_disperser),noheader=True)

                        else:
                            for cur_type in seek_types:
                                filein, _,filenum = get_reduced_data(summary_tab_red, types=[cur_type],
                                                                       slit=cur_slit, disperser=cur_disperser,
                                                                       binning=cur_binning, state_in='normalized',
                                                                       w_dir=w_dir, process_unc=False,
                                                                       process_numbers=True)
                                if len(filein)>0 or cur_type=="sunsky":
                                    if linearization_do_test_val.get() and (cur_type != "sunsky"):
                                        run_test = "{}_s{}_{}_{}_lin_test.fits".format(cur_type, cur_slit,
                                                                                           cur_disperser,
                                                                                           cur_binning)
                                    else:
                                        run_test = None
                                    if cur_type != "sunsky":
                                        file_atmdisp, _ = get_reduced_data(summary_tab_red,
                                                                          types=['atmdisp_{}'.format(cur_type)],
                                                                          slit=cur_slit, disperser=cur_disperser,
                                                                          binning=cur_binning,
                                                                          state_in='ini',
                                                                          w_dir=w_dir, process_unc=False,
                                                                          process_numbers=False, qual='.dat')
                                        if len(file_atmdisp)==0:
                                            file_atmdisp=[None]
                                    else:
                                        file_atmdisp=[None]
                                    file_shifts, _ = get_reduced_data(summary_tab_red,
                                                        types=['shifts_{}'.format(cur_type)],
                                                        slit=cur_slit, disperser=cur_disperser,binning=cur_binning,
                                                        state_in='ini',w_dir=w_dir, process_unc=False,
                                                        process_numbers=False, qual='.dat')
                                    if len(file_shifts)==0:
                                        message(
                                                "!!!! Can't find file with relative shifts computed at previous stages. Check this!",
                                                noheader=True)
                                        file_shifts = [None]

                                    fileout_transform = "transform_{}_s{}_{}_{}_lin.fits".format(cur_type, cur_slit, cur_disperser,
                                                                                            cur_binning)
                                    message("Prepare final geometric correction for {} files of {} type for slit {}, binning {}, disperser {}".format(
                                        len(filein), cur_type, cur_slit,cur_binning,cur_disperser))
                                    status=reduction.linearization(filein_transform=filein_transform[0],fileout_transform=fileout_transform,dispcurve=dispcurve[0],
                                                            file_shifts=file_shifts[0], lin_parameters=lin_parameters, w_dir=w_dir, file_number=filenum,
                                                            file_atmdisp=file_atmdisp[0])
                                    if status:
                                        message(".... Done",noheader=True)
                                        do_update_RedSummary(add=[fileout_transform])
                                        if run_test:
                                            reduction.correct_geometry([filein[0]], fileout_transform, file_err=None,
                                                        fileout=[run_test],w_dir=w_dir)
                                            do_update_RedSummary(add=[run_test])
                    else:
                        message("!!!! Cannot find disperison curve for disperser {}, binning {}, slit {} => cannot continue!".format(
                            cur_disperser,cur_binning,cur_slit))

    if steps['do_transform'][1]:
        message("Apply transformation to all data")
        seek_types = ['obj', 'star','sunsky']
        binnings, dispersers, slits = get_all_redspec_params(summary_tab_red,file_types=seek_types)
        for cur_binning in binnings:
            for cur_disperser in dispersers:
                cur_default_transform = None
                if os.path.isfile(os.path.join(w_dir,default_transform.get())):
                    #cur_default_state = re.findall(r"_\d+x\d+_?(\S+).fits", default_transform.get())
                    #if len(cur_default_state) > 0 and cur_default_state[0] == 'lin':
                    #cur_default_type=re.findall(r"\S+_(\w+)_s\d?\.?\d?_[\w@]+_\d+x\d+",default_transform.get())
                    #if len(cur_default_type)>0 and cur_default_type==cur_type:
                    cur_default_disperser=re.findall(r"\S+_s\d?\.?\d?_([\w@]+)_\d+x\d+",default_transform.get())
                    cur_default_bin = re.findall(r"\S+_s\d?\.?\d?_[\w@]+_(\d+x\d+)", default_transform.get())
                    if len(cur_default_disperser)>0 and cur_disperser == cur_default_disperser[0] and \
                        len(cur_default_bin) > 0 and cur_binning == cur_default_bin[0]:
                            cur_default_transform=default_transform.get()
                for cur_type in seek_types:
                    for cur_slit in slits:
                        if default_transform_force.get() and cur_default_transform is not None:
                            file_transform = [cur_default_transform]
                        else:
                            file_transform, _ = get_reduced_data(summary_tab_red, types=['transform_{}'.format(cur_type)],
                                                              slit=cur_slit, disperser=cur_disperser,
                                                              binning=cur_binning, state_in='linearized',
                                                              w_dir=w_dir, process_unc=False,
                                                              process_numbers=False)
                            if len(file_transform)==0 and cur_default_transform is not None:
                                file_transform = [cur_default_transform]
                            elif len(file_transform)==0:
                                message("Cannot find transformation matrix computed at previous stages for type {}, slit {}, binning {}, disperser {}  => cannot continue with final transformation!".format(
                                    cur_type,cur_slit,cur_binning,cur_disperser))
                                continue
                        if cur_type is not "sunsky":
                            filein, fileout, fileerr, fileerrout, filenum = get_reduced_data(summary_tab_red,
                                                                 types=[cur_type],
                                                                 slit=cur_slit, disperser=cur_disperser,
                                                                 binning=cur_binning, state_in='normalized',
                                                                 state_out='linearized',w_dir=w_dir, process_unc=True,
                                                                 process_numbers=True)

                        else:
                            filein, fileout, fileerr, fileerrout = get_reduced_data(summary_tab_red,
                                                                                             types=[cur_type],
                                                                                             slit=cur_slit,
                                                                                             disperser=cur_disperser,
                                                                                             binning=cur_binning,
                                                                                             state_in='normalized',
                                                                                             state_out='linearized',
                                                                                             w_dir=w_dir,
                                                                                             process_unc=True,
                                                                                             process_numbers=False)
                            filenum=None
                        if len(filein) > 0:
                            file_to_del, _, fileerr_to_del, _, ind_to_del = get_reduced_data(summary_tab_red,
                                                                                             types=[cur_type],
                                                                                             slit=cur_slit,
                                                                                             disperser=cur_disperser,
                                                                                             binning=cur_binning,
                                                                                             state_in='linearized',
                                                                                             w_dir=w_dir,
                                                                                             get_index=True,
                                                                                             process_unc=True,
                                                                                             process_numbers=False)
                            if len(ind_to_del) > 0:
                                for f in file_to_del + fileerr_to_del:
                                    if f is not None:
                                        os.remove(os.path.join(w_dir, f))
                                summary_tab_red.remove_rows(ind_to_del)
                            message("Perform geometric transformation for {} {} frames with binning {}, disperser {} and slit {}".format(
                                    len(filein), cur_type, cur_binning, cur_disperser, cur_slit))
                            status=reduction.correct_geometry(filein,file_transform[0],file_err=fileerr, fileout=fileout,
                                                        fileout_err=fileerrout, w_dir=w_dir, file_number=filenum)
                            if status:
                                do_update_RedSummary(add=fileout+fileerrout)
                            else:
                                message("Something went wrong with linearization")
        message("Finished with linearization")

    if steps['sky_subtract'][1]:
        seek_types=[]
        if skyrem_obj_process_val.get():
            seek_types.append('obj')
        if skyrem_star_process_val.get():
            seek_types.append('star')
        if len(seek_types) == 0:
            message("Nothing to do with sky subtraction")
        else:
            message("Begin to model and subtract night-sky lines")
            skyrem_apply_to_raw = False
            if not skyrem_apply_to_raw:
                suffix_in='_lin'
                cur_state='linearized'
                suffix_out = '_skyrem'
            else:
                suffix_in = '_n'
                cur_state = 'normalized'
                suffix_out = '_rskyrem'
            binnings, dispersers, slits = get_all_redspec_params(summary_tab_red,file_types=seek_types)

            if skyrem_fft_val.get():
                do_fft=True
            else:
                do_fft=False
            for cur_binning in binnings:
                for cur_disperser in dispersers:
                    for cur_slit in slits:
                        file_sunsky,_,fileerr_sunsky,_ = get_reduced_data(summary_tab_red,
                                                            types=['sunsky'],slit=cur_slit,
                                                        disperser=cur_disperser,
                                                        binning=cur_binning,
                                                        state_in=cur_state,
                                                        w_dir=w_dir, process_unc=True,
                                                        process_numbers=False)
                        if not do_fft or skyrem_mode.get() != "Template" or not skyrem_sunsky_val.get() \
                                or (len(file_sunsky) == 0):
                            file_sunsky=None
                        else:
                            file_sunsky=file_sunsky[0]
                        for cur_type in seek_types:
                            order = None
                            cr_clean = True
                            medorder=None
                            minorder=None
                            if cur_type=="obj":
                                if skyrem_mode.get() != "Template":
                                    order = int(skyrem_poly_obj.get())
                                    if skyrem_minorder_obj.get()!='' and skyrem_minsn_obj.get() != '':
                                            minorder=(int(skyrem_minorder_obj.get()),float(skyrem_minsn_obj.get()))
                                    if skyrem_medorder_obj.get()!='' and skyrem_medsn_obj.get() != '':
                                            medorder=(int(skyrem_medorder_obj.get()),float(skyrem_medsn_obj.get()))
                                if window_skyreg_object is not None:
                                    sky_regions=window_skyreg_object.regions
                                else:
                                    sky_regions=[]
                                if len(sky_regions)==0:
                                    sky_regions = [[20, 200], [750,990]]

                                try:
                                    id=list(summary_tab_red['filename']).index('obj_s{}_{}_{}_crmask.fits'.format(cur_slit,cur_disperser,cur_binning))
                                    if summary_tab_red['process'][id] == False:
                                        cr_clean=True
                                    else:
                                        cr_clean=False
                                except ValueError:
                                    cr_clean = True
                            else:
                                if window_skyreg_star is not None:
                                    sky_regions = window_skyreg_star.regions
                                else:
                                    sky_regions = []
                                if len(sky_regions) == 0:
                                    sky_regions = [[20, 200], [750,990]]

                                try:
                                    id = list(summary_tab_red['filename']).index(
                                        'star_s{}_{}_{}_crmask.fits'.format(cur_slit, cur_disperser, cur_binning))
                                    if summary_tab_red['process'][id] == False:
                                        cr_clean = True
                                    else:
                                        cr_clean=False
                                except ValueError:
                                    cr_clean = True

                                if skyrem_mode.get() != "Template":
                                    order = int(skyrem_poly_star.get())
                                    if skyrem_minorder_star.get()!='' and skyrem_minsn_star.get() != '':
                                            minorder=(int(skyrem_minorder_star.get()),float(skyrem_minsn_star.get()))
                                    if skyrem_medorder_star.get()!='' and skyrem_medsn_star.get() != '':
                                            medorder=(int(skyrem_minorder_star.get()),float(skyrem_medsn_star.get()))
                            filein, fileout, fileerr, fileerrout, filenum = get_reduced_data(summary_tab_red,
                                                types=[cur_type], slit=cur_slit,disperser=cur_disperser,
                                                binning=cur_binning,state_in=cur_state,state_out='sky-subtracted',
                                                w_dir=w_dir, process_unc=True,process_numbers=True)
                            if len(filein) > 0:
                                if skyrem_do_test_val.get():
                                    run_test="{}_s{}_{}_{}_skyrem_test.fits".format(cur_type, cur_slit, cur_disperser,
                                                                                     cur_binning)
                                else:
                                    run_test=None
                                message("Perform sky subtraction for {} {} frames with binning {}, disperser {} and slit {}".format(
                                        len(filein), cur_type, cur_binning, cur_disperser, cur_slit))
                                if skyrem_apply_to_raw:
                                    file_transform, _ = get_reduced_data(summary_tab_red,types=['transform'],
                                                slit=cur_slit,disperser=cur_disperser,binning=cur_binning,state_in='ini',
                                                w_dir=w_dir,process_unc=False,process_numbers=False)
                                    file_shifts, _ = get_reduced_data(summary_tab_red, types=['shifts_{}'.format(cur_type)],
                                                                         slit=cur_slit, disperser=cur_disperser,
                                                                         binning=cur_binning, state_in='ini',
                                                                         w_dir=w_dir, process_unc=False,
                                                                         process_numbers=False,qual='.dat')
                                    if len(file_transform) == 0:
                                        file_transform=[None]
                                    if len(file_shifts)==0:
                                        file_shifts = [None]
                                else:
                                    file_transform=[None]
                                    file_shifts=[None]

                                if not run_test:
                                    file_to_del, _, fileerr_to_del, _, ind_to_del = get_reduced_data(summary_tab_red,
                                                                                                     types=[cur_type],
                                                                                                     slit=cur_slit,
                                                                                                     disperser=cur_disperser,
                                                                                                     binning=cur_binning,
                                                                                                     state_in='sky-subtracted',
                                                                                                     w_dir=w_dir,
                                                                                                     get_index=True,
                                                                                                     process_unc=True,
                                                                                                     process_numbers=False)
                                    if len(ind_to_del) > 0:
                                        for f in file_to_del + fileerr_to_del:
                                            if f is not None:
                                                os.remove(os.path.join(w_dir, f))
                                        summary_tab_red.remove_rows(ind_to_del)
                                status=reduction.remove_sky(filein, file_err=fileerr, fileout=fileout,
                                                     file_transform=file_transform[0], do_fft=do_fft,
                                                    fileout_err=fileerrout, w_dir=w_dir, file_number=filenum,
                                                    order=order, regions=sky_regions, run_test=run_test,
                                                    file_shifts=file_shifts[0],cr_clean=cr_clean,
                                                    minorder=minorder,medorder=medorder, file_sunsky=file_sunsky)
                                if status:
                                    if run_test:
                                        do_update_RedSummary(add=[run_test])
                                    else:
                                        do_update_RedSummary(add=fileout+fileerrout)
                                    message("..... Done",noheader=True)
        message("Finish with sky subtraction")
    if steps['combine'][1]:
        message("Start combining different exposures")
        seek_types = ['obj', 'star']
        binnings, dispersers, slits = get_all_redspec_params(summary_tab_red,file_types=seek_types)
        for cur_binning in binnings:
            for cur_disperser in dispersers:
                for cur_slit in slits:
                    for cur_type in seek_types:
                        filein, _, fileerr, _, filenum = get_reduced_data(summary_tab_red, types=[cur_type],
                                                             slit=cur_slit, disperser=cur_disperser,
                                                             binning=cur_binning, state_in='sky-subtracted',
                                                             w_dir=w_dir, process_unc=True, process_numbers=True)
                        if len(filein) > 0:
                            message("Summing individual exposures for {} {} frames with binning {}, disperser {} and slit {}".format(
                                    len(filein), cur_type, cur_binning, cur_disperser, cur_slit))
                            #filein=[filein[0]]
                            fileout="{}_s{}_{}_{}_tot.fits".format(cur_type,cur_slit,cur_disperser,cur_binning)
                            fileout_err = "{}_s{}_{}_{}_tot_err.fits".format(cur_type,cur_slit,cur_disperser,cur_binning)
                            status=reduction.exp_combine(filein, file_err=fileerr, fileout=fileout,mode=combine_mode.get(),
                                                  fileout_err=fileout_err, w_dir=w_dir, file_number=filenum,
                                                  sigma=float(combine_treshold.get()))
                            if status:
                                do_update_RedSummary(add=[fileout, fileout_err])
                            else:
                                message("Something went wrong with combination of different exposures")
        message("Finish with combination of different exposures")

    if steps['calc_dqe'][1]:
        message("Start to compute spectral sensitivity curve")
        seek_types = ['star']
        binnings, dispersers, slits = get_all_redspec_params(summary_tab_red,file_types=seek_types)

        do_plot = dqe_plot_val.get()
        star_name = dqe_starname.get()
        if star_name is '':
            star_name = None
        if dqe_smooth.get() is '':
            smo_win = None
        else:
            smo_win = int(dqe_smooth.get()) // 2 * 2 + 1
        if dqe_starz.get() is '':
            z_star = None
        else:
            z_star = float(dqe_starz.get())
        if dqe_extract_bottom.get() is '':
            ew0 = None
        else:
            ew0 = int(dqe_extract_bottom.get())
        if dqe_extract_top.get() is '':
            ew1 = None
        else:
            ew1 = int(dqe_extract_top.get())
        extract_window=(ew0,ew1)
        for cur_binning in binnings:
            for cur_disperser in dispersers:
                for cur_slit in slits:
                    filein, _, filein_err, _ = get_reduced_data(summary_tab_red, types=['star'],
                                                                      slit=cur_slit, disperser=cur_disperser,
                                                                      binning=cur_binning, state_in='combined',
                                                                      w_dir=w_dir, process_unc=True,
                                                                      process_numbers=False)
                    if len(filein) ==0:
                        message("Can't find star for binning {}, slit {} and disperser {}".format(cur_binning,cur_slit,cur_disperser),noheader=True)
                        continue
                    file_dqe="senscurve_s{}_{}_{}.fits".format(cur_slit, cur_disperser, cur_binning)
                    dqe_figname="dqe_s{}_{}_{}.pdf".format(cur_slit, cur_disperser, cur_binning)
                    status=reduction.calc_dqe(filein[0], file_err=filein_err[0], fileout=file_dqe,mode=mode.get(),
                                       standards_dir=standard_data,w_dir=w_dir, do_plot=do_plot,
                                       star_name=star_name, smo_win=smo_win, z_star=z_star,
                                       extract_window=extract_window,dqe_figname=dqe_figname)
                    if status:
                        do_update_RedSummary(add=[file_dqe])
                    else:
                        message("Something went wrong with calculation of sensitivity curve")
        message("Spectral sensitivity curve is computed")

    if steps['flux_cal'][1]:
        seek_types = ['obj']
        message("Start flux calibration")
        binnings, dispersers, slits = get_all_redspec_params(summary_tab_red,file_types=seek_types)

        if fluxcal_objz.get() is "":
            z_obj=None
        else:
            z_obj=float(fluxcal_objz.get())
        if fluxcal_ypos.get() is "":
            ypos = None
        else:
            ypos = float(fluxcal_ypos.get())
        if fluxcal_ydelt.get() is "":
            ydelt = None
        else:
            ydelt = float(fluxcal_ydelt.get())
        if fluxcal_lam0.get() is '':
            lam0 = None
        else:
            lam0 = float(fluxcal_lam0.get())
        if fluxcal_lam1.get() is '':
            lam1 = None
        else:
            lam1 = float(fluxcal_lam1.get())
        cut_wl=(lam0,lam1)
        ycut=[0,None]
        if fluxcal_ycut0.get():
            ycut[0]=int(fluxcal_ycut0.get())
        if fluxcal_ycut1.get():
            ycut[1]=int(fluxcal_ycut1.get())
        for cur_binning in binnings:
            for cur_disperser in dispersers:
                cur_default_senscurve = None
                if os.path.isfile(os.path.join(w_dir, default_sens_curve.get())):
                    cur_default_disperser = re.findall(r"\S+_s\d?\.?\d?_([\w@]+)_\d+x\d+", default_sens_curve.get())
                    cur_default_bin = re.findall(r"\S+_s\d?\.?\d?_[\w@]+_(\d+x\d+)", default_sens_curve.get())
                    if len(cur_default_disperser) > 0 and cur_disperser == cur_default_disperser[0] and \
                            len(cur_default_bin) > 0 and cur_binning == cur_default_bin[0]:
                        cur_default_senscurve = default_sens_curve.get()
                for cur_slit in slits:
                    if default_sens_curve_force.get() and cur_default_senscurve is not None:
                        file_sens = [cur_default_senscurve]
                    else:
                        file_sens, _ = get_reduced_data(summary_tab_red, types=['senscurve'],
                                                    slit=cur_slit, disperser=cur_disperser, binning=cur_binning,
                                                    state_in='ini',
                                                    w_dir=w_dir, process_unc=False, process_numbers=False)
                    if len(file_sens) == 0 and cur_default_senscurve is None:
                        message("Sensitivity curve is not available for binning {}, disperser {} and slit {}".format(
                            cur_binning,
                            cur_disperser,
                            cur_slit))
                        continue
                    elif len(file_sens) == 0:
                        file_sens=[cur_default_senscurve]
                    filein, fileout, filein_err, fileout_err = get_reduced_data(summary_tab_red, types=['obj'],
                                    slit=cur_slit, disperser=cur_disperser,binning=cur_binning, state_in='combined',
                                    state_out='calibrated',w_dir=w_dir, process_unc=True,process_numbers=False)
                    if len(filein)>0:
                        file_sunsky, file_sunsky_out, file_sunsky_err, file_sunsky_err_out = get_reduced_data(summary_tab_red, types=['sunsky'],
                                                                                    slit=cur_slit,
                                                                                    disperser=cur_disperser,
                                                                                    binning=cur_binning,
                                                                                    state_in='linearized',
                                                                                    state_out='calibrated', w_dir=w_dir,
                                                                                    process_unc=True,
                                                                                    process_numbers=False)
                        if len(file_sunsky)==0:
                            file_sunsky=[None]
                            file_sunsky_out = [None]
                            file_sunsky_err = [None]
                            file_sunsky_err_out = [None]
                        message("Perform flux calibration for binning {}, disperser {} and slit {}".format(cur_binning,
                                                                                                           cur_disperser,
                                                                                                           cur_slit))
                        status=reduction.flux_calib(filein[0], file_err=filein_err[0], file_sens=file_sens[0], fileout=fileout[0],
                                    cut_wl=cut_wl,fileout_err=fileout_err[0], w_dir=w_dir, z_obj=z_obj, ypos=ypos, ydelt=ydelt,
                                                    files_for_cut=[file_sunsky[0],file_sunsky_err[0]],ycut=ycut,
                                                    files_for_cut_out=[file_sunsky_out[0],file_sunsky_err_out[0]])
                        if status:
                            message(".... Done", noheader=True)
                            do_update_RedSummary(add=[fileout[0],fileout_err[0],file_sunsky_out[0],file_sunsky_err_out[0]])
                        else:
                            message(".... Something went wrong", noheader=True)
        message("Finish with flux calibration")

def update_all_droplists(summary_tab=None, objects=[],dispersers=[], pas=[], slits=[], binnings=[]):
    global flat_list, bias_list, obj_list, standard_list, sunsky_list, neon_list, other_list, dark_list
    if summary_tab:
        mask = [idx for idx, val in enumerate(summary_tab['type']) if val in ["obj", "star"]]
        objects = get_unique(summary_tab['object'][mask])
        mask = [idx for idx, val in enumerate(summary_tab['type']) if val in ["obj", "neon", 'star', 'flat', 'sunsky']]
        dispersers = get_unique(summary_tab['disperser'][mask])
        mask = [idx for idx, val in enumerate(summary_tab['type']) if
                val in ["obj", "neon", 'star', 'flat', 'sunsky', 'bias', 'dark']]
        binnings = get_unique(summary_tab['binning'][mask])
        mask = [idx for idx, val in enumerate(summary_tab['type']) if
                val in ["obj", "neon", 'star', 'flat', 'sunsky']]
        slits = get_unique(summary_tab['slit'][mask])
        mask = [idx for idx, val in enumerate(summary_tab['type']) if val in ["obj", 'star']]
        pas = get_unique(summary_tab['pa'][mask])

    bias_list.update_droplists(objects=objects,dispersers=dispersers, pas=pas, slits=slits, binnings=binnings)
    flat_list.update_droplists(objects=objects,dispersers=dispersers, pas=pas, slits=slits, binnings=binnings)
    obj_list.update_droplists(objects=objects,dispersers=dispersers, pas=pas, slits=slits, binnings=binnings)
    sunsky_list.update_droplists(objects=objects,dispersers=dispersers, pas=pas, slits=slits, binnings=binnings)
    dark_list.update_droplists(objects=objects,dispersers=dispersers, pas=pas, slits=slits, binnings=binnings)
    standard_list.update_droplists(objects=objects,dispersers=dispersers, pas=pas, slits=slits, binnings=binnings)
    other_list.update_droplists(objects=objects,dispersers=dispersers, pas=pas, slits=slits, binnings=binnings)
    neon_list.update_droplists(objects=objects,dispersers=dispersers, pas=pas, slits=slits, binnings=binnings)




def do_update_RedSummary(add=[]):
    global mode,w_dir, summary_tab_red
    if summary_tab_red is None:
        summary_tab_red=do_ScanRedDir(w_dir)
    else:
        for f in add:
            if f is None:
                continue
            currow, _ = parse_redfits(f, mode=mode.get())
            currow['subdirectory'] = f.split(currow['filename'][0])[0]
            try:
                ind=list(summary_tab_red['filename']).index(currow['filename'][0])
                if currow['subdirectory'][0] == summary_tab_red['subdirectory'][ind]:
                    summary_tab_red[ind] = currow[0]
                else:
                    summary_tab_red = vstack([summary_tab_red, currow])
            except ValueError:
                summary_tab_red = vstack([summary_tab_red, currow])
            do_RedSummary(summary_tab_red)


def do_ScanRedDir(w_dir):
    global mode,summary_tab,summary_tab_red
    sd = scan_dir(w_dir, mode=mode.get(), dir_mode='Red', exist_table=summary_tab_red)
    if not sd[1]:
        summary_tab_red=None
        if summary_tab is None:
            button_run.config(state="disabled")
            do_process_but.config(state='disabled')
        message("I can't find any appropriate fits file in w_dir")
        do_RedSummary(None)
    else:
        summary_tab_red=sd[1]
        button_run.config(state="normal")
        do_process_but.config(state='normal')
        starobj_selector_button.config(state="normal")
        do_RedSummary(summary_tab_red)
    return summary_tab_red

def do_ScanDir(data_dir):
    global mode, standard_dir,summary_tab_red
    sd = scan_dir(data_dir, mode=mode.get(), standards = os.path.join(standard_dir,'standards_list.txt'))
    if not sd[1]:
        summary_tab=None
        raw_files=None
        if summary_tab_red is None:
            button_run.config(state="disabled")
            do_process_but.config(state='disabled')
        starobj_selector_button.config(state="disabled")
        message("I can't find any appropriate fits file in data_dir")
        do_Summary(None)
        update_all_droplists()
    else:
        raw_files=sd[0]
        summary_tab=sd[1]
        button_run.config(state="normal")
        starobj_selector_button.config(state="normal")
        do_process_but.config(state='normal')
        do_Summary(summary_tab)
        update_all_droplists(summary_tab)
    return (raw_files,summary_tab)


def do_RedSummary(summary):
    global red_list, show_red_datatype
    if summary is None:
        red_list.update([''], None)
        return
    datatype = show_red_datatype.get()
    redstate= red_list.redstate.get()
    frametype=red_list.redsubtype.get()

    data_show = []
    process_show = []


    for row in summary:
        if (datatype in ["All", row['type']]) and (redstate in ["All", row['state']]) and (frametype in ["All", row['subtype']]):
            data_show.append((row['data'], row['slit'],row['disperser'],  row['binning'],
                                   row['type'], row['state'], row['subtype'],
                                  row['subdirectory'],row['filename'],row['process']))
            process_show.append(row['process'])

    red_list.update(data_show, process_show)


def do_Summary(summary, object="All", disperser="All", binning = "All", slit = "All", pa="All"):
    global flat_list, bias_list, obj_list, standard_list, sunsky_list, neon_list, other_list, dark_list
    if not summary:
        flat_list.update([''],None)
        bias_list.update([''],None)
        obj_list.update([''],None)
        neon_list.update([''], None)
        dark_list.update([''], None)
        sunsky_list.update([''], None)
        other_list.update([''], None)
        standard_list.update([''], None)
        return
    data_flat=[]
    data_bias = []
    data_obj=[]
    data_neon = []
    data_sunsky = []
    data_dark = []
    data_other = []
    data_standard = []
    process_flat = []
    process_bias = []
    process_obj = []
    process_standard = []
    process_neon = []
    process_sunsky = []
    process_dark = []
    process_other = []

    for row in summary:
        if row['type']=='bias':
            if (binning in ["All", row['binning']]):
                data_bias.append((row['binning'],row['date'],row['subdirectory'],row['filename'],row['process']))
            process_bias.append(row['process'])
        elif row['type']=='obj':
            if (binning in ["All", row['binning']]) and (object in ["All", row['object']]) and \
                    (pa in ["All", row['pa']]) and (slit in ["All", row['slit']]) and (disperser in ["All",row['disperser']]):
                data_obj.append((row['object'],row['disperser'],row['pa'],row['slit'],row['exposure'],row['binning'],
                                  row['date'], row['subdirectory'], row['filename'],row['process']))
            process_obj.append(row['process'])
        elif row['type'] =='star':
            if (binning in ["All", row['binning']]) and (object in ["All", row['object']]) and \
                    (pa in ["All", row['pa']]) and (slit in ["All", row['slit']]) and (disperser in ["All", row['disperser']]):
                data_standard.append((row['object'], row['disperser'], row['pa'], row['slit'], row['exposure'], row['binning'],
                             row['date'], row['subdirectory'], row['filename'],row['process']))
            process_standard.append(row['process'])
        elif row['type'] =='flat':
            if (binning in ["All", row['binning']]) and (slit in ["All", row['slit']]) and (disperser in ["All", row['disperser']]):
                data_flat.append((row['disperser'], row['slit'], row['exposure'], row['binning'],
                             row['date'], row['subdirectory'], row['filename'],row['process']))
            process_flat.append(row['process'])
        elif row['type'] =='neon':
            if (binning in ["All", row['binning']]) and (slit in ["All", row['slit']]) and (disperser in ["All", row['disperser']]):
                data_neon.append((row['disperser'], row['slit'], row['exposure'], row['binning'],
                             row['date'], row['subdirectory'], row['filename'],row['process']))
            process_neon.append(row['process'])
        elif row['type'] =='dark':
            if (binning in ["All", row['binning']]):
                data_dark.append((row['exposure'], row['binning'],
                             row['date'], row['subdirectory'], row['filename'],row['process']))
            process_dark.append(row['process'])
        elif row['type'] =='sunsky':
            if (binning in ["All", row['binning']]) and (slit in ["All", row['slit']]) and (disperser in ["All", row['disperser']]):
                data_sunsky.append((row['disperser'], row['slit'], row['exposure'], row['binning'],
                             row['date'], row['subdirectory'], row['filename'],row['process']))
            process_sunsky.append(row['process'])
        else:
            if (binning in ["All", row['binning']]) and (slit in ["All", row['slit']]) and (disperser in ["All", row['disperser']]):
                data_other.append((row['type'],row['disperser'], row['slit'], row['exposure'], row['binning'],
                                row['date'], row['subdirectory'], row['filename'],row['process']))
            process_other.append(row['process'])

    bias_list.update(data_bias,process_bias)
    flat_list.update(data_flat,process_flat)
    obj_list.update(data_obj,process_obj)
    standard_list.update(data_standard, process_standard)
    neon_list.update(data_neon, process_neon)
    dark_list.update(data_dark, process_dark)
    sunsky_list.update(data_sunsky, process_sunsky)
    other_list.update(data_other,process_other)


def ChangeProcess():
    if show_summary_mode.get() == 'Red':
        cur_list = red_list
        cur_summary_tab=summary_tab_red
    else:
        cur_summary_tab=summary_tab
        if show_datatype.get() == "obj":
            cur_list = obj_list
        elif show_datatype.get() == "flat":
            cur_list = flat_list
        elif show_datatype.get() == "bias":
            cur_list = bias_list
        elif show_datatype.get() == "star":
            cur_list = standard_list
        elif show_datatype.get() == "neon":
            cur_list = neon_list
        elif show_datatype.get() == "dark":
            cur_list = dark_list
        elif show_datatype.get() == "sunsky":
            cur_list = sunsky_list
        elif show_datatype.get() == "other":
            cur_list = other_list
    selection = cur_list.tree.selection()
    if len(selection) > 0:
        items=[]
        for s in selection:
            if show_summary_mode.get()=="Red":
                selected_file = cur_list.data[cur_list.tree.index(s)][red_list_header['Filename']]
                selected_dir = cur_list.data[cur_list.tree.index(s)][red_list_header['Directory']]
            else:
                selected_file = cur_list.data[cur_list.tree.index(s)][-2]
                selected_dir = cur_list.data[cur_list.tree.index(s)][-3]
            items.append(cur_list.tree.index(s))
            id = ((cur_summary_tab['filename'] == selected_file) & (cur_summary_tab['subdirectory'] == selected_dir))
            cur_summary_tab['process'][id] = bool(do_process_val.get())

        if show_summary_mode.get()=="Red":
            do_RedSummary(summary_tab_red)
        else:
            do_Summary(summary_tab)
        child_id = cur_list.tree.get_children()
        for ind, id in enumerate(items):
            if ind == 0:
                cur_list.tree.selection_set(child_id[id])
            else:
                cur_list.tree.selection_add(child_id[id])
            cur_list.tree.focus(child_id[id])
        # if cur_list.do_process:
            #     cur_list.do_process[cur_list.tree.index(s)]=do_process_val.get()

class RegTable(tk.Frame):
    def __init__(self, parent, nrows):
        tk.Frame.__init__(self, parent)
        self.cells={}
        self.buttons={}
        self.parent=parent
        tk.Label(self,text="Pix Start:").grid(row=0,column=0)
        tk.Label(self, text="Pix End:").grid(row=0, column=1)
        self.cells_frame=tk.Frame(self)
        self.cells_frame.grid(row=1,column=0,columnspan=2,rowspan=12)
        for i in range(nrows):  # Rows
            for j in range(2):  # Columns
                b = tk.Entry(self.cells_frame, text="", validate='key')
                b["validatecommand"] = (b.register(IntFieldTestVal), '%P', '%d')
                b.grid(row=i, column=j, sticky='ns')
                self.cells[(i, j)] = b
            self.buttons[i]=tk.Button(self.cells_frame, text="X", command=lambda c=i: self.delrow(c))
            self.buttons[i].grid(row=i, column=j + 1)


    def regrid(self):
        nrows = int(len(self.cells) / 2)
        for i in range(nrows):
            for j in range(2):
                self.cells[(i,j)].grid_forget()
                self.cells[(i,j)].grid(row=i, column=j, sticky='ns')
            self.buttons[i].grid_forget()
            self.buttons[i].grid(row=i, column=j + 1)


    def delrow(self, ind):
        if len(self.cells)/2 > 1:
            cells_new={}
            but_new={}
            for i in range(int(len(self.cells)/2)):
                if i<ind:
                    cells_new[(i,0)]=self.cells[(i,0)]
                    cells_new[(i, 1)] = self.cells[(i, 1)]
                    but_new[i]=self.buttons[i]
                elif i>ind:
                    cells_new[(i-1, 0)] = self.cells[(i, 0)]
                    cells_new[(i-1, 1)] = self.cells[(i, 1)]
                    but_new[i-1] = self.buttons[i]
                    but_new[i - 1].config(command=lambda c=(i - 1): self.delrow(c))
            self.cells[(ind, 0)].grid_forget()
            self.cells[(ind, 1)].grid_forget()
            self.buttons[ind].grid_forget()
            del self.cells[(ind, 0)], self.cells[(ind, 1)], self.buttons[ind]
            self.cells=cells_new.copy()
            self.buttons=but_new.copy()
            del cells_new,but_new
            self.regrid()


    def addrow(self):
        i = (int(len(self.cells) / 2))
        for j in range(2):  # Columns
            b = tk.Entry(self.cells_frame, text="", validate="key")
            b["validatecommand"] = (b.register(IntFieldTestVal), '%P', '%d')
            b.grid(row=i, column=j, sticky='ns')
            self.cells[(i, j)] = b
        self.buttons[i] = tk.Button(self.cells_frame, text="X", command=lambda c=i: self.delrow(c))
        self.buttons[i].grid(row=i, column=j + 1, sticky='ns')

class SkyRegWindow(tk.Toplevel):
    def __init__(self, parent, title="Sky regions", sky_regions=[], image=None, summary_data=None, mode=None):
        tk.Toplevel.__init__(self,parent)
        self.parent=parent
        self.image=image
        self.mode=mode
        self.regions=sky_regions
        self.summary_data=summary_data
        self.title(title)
        self.protocol('WM_DELETE_WINDOW', self.hidewin)
        self.regtab=RegTable(self,nrows=summary_data["nregions"].get()+1)
        self.regtab.grid(row=0, column=0, columnspan=2, rowspan=4,sticky=tk.NSEW)
        tab_button = tk.Frame(self)
        tab_button.grid(row=4, column=0,columnspan=2)
        tk.Button(tab_button, text="Load reference image", command=self.get_image).grid(row=0, column=0, padx=5)
        tk.Button(tab_button, text="Add region", command=self.regtab.addrow).grid(row=0, column=1, padx=5)
        tk.Button(tab_button, text="Save and close", command=self.savereg).grid(row=0, column=2, padx=5)
        tk.Button(tab_button, text="Close", command=self.hidewin).grid(row=0, column=3, padx=5)

        fig = Figure(figsize=(12, 7))
        self.ax = fig.add_subplot(111)
        tab_image = tk.Frame(self)
        tab_image.grid(row=0, column=2,columnspan=4,rowspan=5,sticky=tk.NSEW)
        self.canvas = FigureCanvasTkAgg(fig, master=tab_image)  # A tk.DrawingArea.
        self.canvas.get_tk_widget().grid(row=0, column=0, columnspan=4, rowspan=5)
        toolbarFrame = tk.Frame(tab_image)
        toolbarFrame.grid(row=7, column=0,sticky=tk.NSEW,columnspan=4)
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbarFrame)
        self.toolbar.update()
        fig.canvas.callbacks.connect('button_press_event', self.on_click)

    def on_click(self, event):
        if self.image is not None:
            if event.inaxes is not None:
                row=int(len(self.regtab.cells)/2)
                if event.button == 1:
                    self.regtab.cells[(row - 1, 0)].delete(0, 'end')
                    self.regtab.cells[(row-1,0)].insert(0,(int(round(event.ydata))))
                elif event.button == 3:
                    self.regtab.cells[(row - 1, 1)].delete(0, 'end')
                    self.regtab.cells[(row-1,1)].insert(0,(int(round(event.ydata))))


    def savereg(self):
        nrows=int(len(self.regtab.cells)/2)
        self.regions=[]
        npix=0
        add_lowborder_sign=""
        for row in reversed(range(nrows)):
            try:
                v0=int(self.regtab.cells[(row,0)].get())
            except ValueError:
                v0=0
            try:
                v1 = int(self.regtab.cells[(row, 1)].get())
            except ValueError:
                v1=0
            if ((v0 == 0) and (v1 == 0)) or ((v1 < v0) and v1 != 0):
                if len(self.regtab.cells)/2 > 1:
                    self.regtab.delrow(row)
                else:
                    self.regtab.cells[(0,0)].delete(0,'end')
                    self.regtab.cells[(0, 1)].delete(0, 'end')
            else:
                if v1 == 0:
                    v1=10000
                if v1!=10000:
                    npix+=(v1-v0+1)
                else:
                    add_lowborder_sign=">"
                self.regions.append([v0,v1])
        self.summary_data["nregions"].set(len(self.regions))
        self.summary_data["npix"].set(npix)
        self.summary_data["label"].config(text="N regs:{:2d}; N pix:{}{:4d}".format(self.summary_data["nregions"].get(),
                                                            add_lowborder_sign,self.summary_data["npix"].get()))
        self.hidewin()


    def hidewin(self):
        self.withdraw()
        self.parent.withdraw()
        self.parent.deiconify()

    def get_image(self):
        global w_dir
        loaded = tk.filedialog.askopenfilename(master=self,initialdir=w_dir,
                                               filetypes=(("Fits files", ["*.fits",".fts"]),("All files", "*.*")))
        if os.path.isfile(loaded):
            with fits.open(loaded) as hdu:
                self.image=hdu[0].data
            self.show_image()

    def show_image(self):
        if self.image is not None:
            interval = PercentileInterval(97.)
            f_stretch = LinearStretch()
            transform_r = f_stretch + interval
            im = transform_r(self.image)
            self.ax.imshow(im, origin='lower', cmap=plt.cm.Greys, extent=(0,im.shape[1],0,im.shape[0]), aspect='auto')
            self.ax.set_xlabel("X, pix")
            self.ax.set_ylabel("Y, pix")
            self.canvas.draw()
            self.toolbar.update()



def SetSkyRegions(mode):
    global window_skyreg_object,window_skyreg_star,sky_regions_obj,sky_regions_star,\
        skyrem_nregions_label_star,skyrem_nregions_label_obj, \
            skyrem_nregions_star,skyrem_nregions_obj,skyrem_npix_star,skyrem_npix_obj

    if mode=='obj':
        window=window_skyreg_object
        sky_regions=sky_regions_obj
        summary_data={"nregions": skyrem_nregions_obj, "npix": skyrem_npix_obj, "label": skyrem_nregions_label_obj}
        title="Regions for sky model (for object frames)"
    else:
        window = window_skyreg_star
        sky_regions = sky_regions_star
        summary_data = {"nregions": skyrem_nregions_star, "npix": skyrem_npix_star, "label": skyrem_nregions_label_star}
        title="Regions for sky model (for star frames)"
    if window is not None: window.deiconify()
    try:
        window.deiconify()
        window.focus()
    except (NameError, AttributeError):
        window = SkyRegWindow(parent=root,title=title, sky_regions=sky_regions,
                             summary_data=summary_data, mode=mode)
        if mode == 'obj':
            window_skyreg_object=window
        else:
            window_skyreg_star=window


def ChangeTypeOfFile():
    if show_summary_mode.get() == "Red":
        selection = red_list.tree.selection()
        if len(selection) == 1:
            selected_file = red_list.data[red_list.tree.index(selection[0])][red_list_header['Filename']]
            selected_dir = red_list.data[red_list.tree.index(selection[0])][red_list_header['Directory']]
            selected_type = red_list.data[red_list.tree.index(selection[0])][red_list_header['Type']]
            selected_data = red_list.data[red_list.tree.index(selection[0])][red_list_header['Data']]
            selected_state = red_list.data[red_list.tree.index(selection[0])][red_list_header['State']]
            if selected_type == 'calibration' and selected_state == 'normalized' and selected_data == 'flat':
                default_flat.delete(0,'end')
                default_flat.insert(0,os.path.join(selected_dir,selected_file))
            elif selected_type == 'auxiliary' and selected_data.casefold().startswith('transform'):
                default_transform.delete(0,'end')
                default_transform.insert(0,os.path.join(selected_dir,selected_file))
            elif selected_type == 'auxiliary' and selected_data == 'dispcurve':
                default_disp_curve.delete(0,'end')
                default_disp_curve.insert(0,os.path.join(selected_dir,selected_file))
            elif selected_type == 'auxiliary' and selected_data == 'senscurve':
                default_sens_curve.delete(0,'end')
                default_sens_curve.insert(0,os.path.join(selected_dir,selected_file))
    else:
        if show_datatype.get() == "obj":
            cur_list = obj_list
            new_type='star'
        elif show_datatype.get() == "star":
            cur_list = standard_list
            new_type='obj'
        selection = cur_list.tree.selection()
        if len(selection) > 0:
            for s in selection:
                selected_file = cur_list.data[cur_list.tree.index(s)][-2]
                selected_dir = cur_list.data[cur_list.tree.index(s)][-3]
                id=((summary_tab['filename'] == selected_file) & (summary_tab['subdirectory'] == selected_dir))
                summary_tab['type'][id] = new_type
            do_Summary(summary_tab)

def ChangeParamSensibility():
    if skyrem_fft_val.get() and (skyrem_mode.get() == "Template"):
        skyrem_sunsky.config(state='normal')
    else:
        skyrem_sunsky.config(state='disabled')
    if not skyrem_obj_process_val.get():
        skyrem_poly_obj.config(state='disabled')
        skyrem_setregion_obj.config(state='disabled')
    else:
        skyrem_poly_obj.config(state='normal')
        skyrem_setregion_obj.config(state='normal')
    if not skyrem_star_process_val.get():
        skyrem_poly_star.config(state='disabled')
        skyrem_setregion_star.config(state='disabled')
    else:
        skyrem_poly_star.config(state='normal')
        skyrem_setregion_star.config(state='normal')
    if not auto_gain_val.get():
        field_gain.config(state='normal')
    else:
        field_gain.config(state='disabled')
    if not auto_readnoise_val.get():
        field_readnoise.config(state='normal')
    else:
        field_readnoise.config(state='disabled')
    if not atmdisp_obj_apply_val.get():
        atmdisp_obj_ypos.config(state='disabled')
    else:
        atmdisp_obj_ypos.config(state='normal')
    if not atmdisp_star_apply_val.get():
        atmdisp_star_ypos.config(state='disabled')
    else:
        atmdisp_star_ypos.config(state='normal')
    if not linearizarion_disp_auto_val.get():
        linearizarion_disp.config(state='normal')
    else:
        linearizarion_disp.config(state='disabled')
    if not linearizarion_lam0_auto_val.get():
        linearizarion_lam0.config(state='normal')
    else:
        linearizarion_lam0.config(state='disabled')
    if not linearizarion_lam1_auto_val.get():
        linearizarion_lam1.config(state='normal')
    else:
        linearizarion_lam1.config(state='disabled')

def ChangeObsDevice():
    if mode.get() != "SAO":
        cut_as_overcan.config(state="normal")
        cut_as_overcan_val.set(True)
    else:
        cut_as_overcan_val.set(False)
        cut_as_overcan.config(state="disabled")




#===============  Begining of root frame description ==================================




###################################
# Define some initial parameters
###################################

raw_files={}
sky_regions_obj=[]
sky_regions_star=[]
window_skyreg_object=None
window_skyreg_star=None

summary_tab=None
summary_tab_red=None

##########################
# Setup widgets hierarchy
##########################

root = tk.Tk()
root.title('Long-slit Data Reduction (v{})'.format(version))
root['bg'] = 'white'
wid_head_font = ("Helvetica",16)
wid_head_color = 'black'
wid_subhead_font = ("Arial",14)
wid_subhead_color = 'black'

frame_control = tk.Frame(root,bg='white',highlightbackground="black",highlightthickness=1)#,height=win_height,width=frames_width[0])
frame_data = tk.Frame(root,bg='white',highlightbackground="black",highlightthickness=1)#,height=win_height,width=frames_width[1])
frame_parameters = tk.Frame(root,bg='white',highlightbackground="black",highlightthickness=1)
frame_display = tk.Frame(root,bg='white',highlightbackground="black",highlightthickness=1)#,height=win_height,width=frames_width[2])
#frame_summary.grid_propagate(False)
frame_control.grid(row=0,column=0,columnspan=1,rowspan=2,sticky="nsew",padx=3)
frame_data.grid(row=0,column=1,columnspan=3,rowspan=2,sticky="nsew",padx=3)
frame_parameters.grid(row=0,column=4,columnspan=1,rowspan=4,sticky="nsew",padx=3)
frame_display.grid(row=2,column=0,columnspan=4,rowspan=2,sticky="nsew",padx=3)



########################
# CONTROL SECTION
#######################

mode_section=tk.Frame(frame_control)
mode_section.grid(row=0,column=0,sticky='n',pady=5)


all_modes = [("TDS data", "TDS"),
        ("SCORPIO data", "SAO"),
        ("Other data", "OTHER")]

mode = tk.StringVar()
mode.set("TDS") # initialize
i=0
for text, curmode in all_modes:
    b = tk.Radiobutton(mode_section, text=text,
                        variable=mode, value=curmode, command = ChangeObsDevice)
    b.grid(column=i,row=0,padx=10)
    i+=1


load_section=tk.Frame(frame_control)
load_section.grid(row=1,column=0,pady=5)#pack(side="top")
button_load_raw_dir=tk.Button(load_section, text="Load Raw Dir.",command=Load_R_Dir)
field_load_raw_dir=tk.Entry(load_section)
button_load_w_dir=tk.Button(load_section, text="Load Work Dir.",command=Load_W_Dir)
field_load_w_dir=tk.Entry(load_section)
field_load_w_dir.insert(0,w_dir)
button_load_log=tk.Button(load_section, text="Load Log file",command=Load_Log)
button_save_log=tk.Button(load_section, text="Save Log file",command=Save_Log)
if data_dir:
    field_load_raw_dir.insert(0, data_dir)

button_load_log.grid(row=0,column=0,padx=1,pady=2)
button_save_log.grid(row=0,column=1,pady=2)
button_load_raw_dir.grid(row=1,column=0,padx=1,pady=2)
field_load_raw_dir.grid(row=1,column=1,pady=2)
button_load_w_dir.grid(row=2,column=0,padx=1,pady=2)
field_load_w_dir.grid(row=2,column=1,pady=2)


class Checkcolumn(tk.Frame):
    def __init__(self, parent=None, picks={}, side=tk.TOP, anchor=tk.W):
        tk.Frame.__init__(self, parent)
        self.vars = []
        for pick in picks.items():
            var = tk.BooleanVar()
            chk = tk.Checkbutton(self, text=pick[1][0], variable=var)
            if pick[1][1]:
                chk.select()
            chk.pack(anchor=anchor, expand=tk.YES)
            self.vars.append(var)

    def state(self):
        return map((lambda var: var.get()), self.vars)


frame_steps=tk.Frame(frame_control)
frame_steps.grid(row=2,column=0, sticky='n',padx=5,pady=10)#pack(side='top')



steps={"cre_bias":["Create meanbias (if any)",False],
       "cre_dark":["Create dark (if any)",False],
       "cre_calibration": ["Create calibrations",False],
       "cre_ini":["Prepare obj and star frames",False],
       "calc_geometry":["Trace geometry",False],
       "cre_norm_flat": ["Create normalized flat",False],
       "norm_flat": ["Correct for normalized flat",False],
       "calc_shifts":["Calculate shifts and atm. dispersion",False],
       "calc_dispersion":["Create dispersion curve",False],
       "lin_and_align": ["Prepare final transformation",False],
       "do_transform":["Apply all transformations",False],
       "sky_subtract":["Model and subtract night sky",False],
       "combine":["Combine exposures",False],
       "calc_dqe":["Derive DQE",True],
       "flux_cal":["Flux Calibration",False]}

steps_buttons=Checkcolumn(frame_steps,steps)
steps_buttons.grid(row=0,column=0, sticky='nw')
steps_buttons.config(relief=tk.GROOVE, bd=2)



frame_buttons=tk.Frame(frame_control,bd=2)
frame_buttons.grid(row=3,column=0,pady=10)
button_run=tk.Button(frame_buttons,text="Run reduction",command=do_Reduction,state="disabled",width=button_width)
button_run.grid(row=0,column=0,padx=5)#.pack(side='top')


button_quit=tk.Button(frame_buttons, text="Quit",command=lambda:root.destroy(),width=button_width)
button_quit.grid(row=0,column=1,padx=5)#.pack(side='bottom')





#####################
# ===== Data section
#####################


frame_summary_type=tk.Frame(frame_data)
frame_summary_type.grid(row=0,column=0,pady=3)

show_summary_mode = tk.StringVar()
show_summary_mode.set("Raw") # initialize
all_modes = [("Raw data", "Raw"),
        ("Reduced data", "Red")]
i=0
for text, curmode in all_modes:
    b = tk.Radiobutton(frame_summary_type, text=text,
                        variable=show_summary_mode, value=curmode, command=ChangeSummaryType)
    b.grid(column=i,row=0,padx=10)#b.pack(side='left')
    i+=1


# frame_select_datatype_blank=tk.Frame(frame_data)
# frame_select_datatype_blank.grid(row=1,column=0,sticky='news')
frame_select_red_datatype=tk.Frame(frame_data)
frame_select_red_datatype.grid(row=1,column=0,sticky='news')


frame_select_datatype=tk.Frame(frame_data)
frame_select_datatype.grid(row=1,column=0,sticky='news')


show_datatype = tk.StringVar()
show_datatype.set("obj") # initialize
all_modes = [("Obj data", "obj"),
        ("Standard star", "star"),
        ("Sunsky", "sunsky"),
        ("Cal. lamp", "neon"),
        ("Flat data", "flat"),
        ("Bias data", "bias"),
        ("Dark data", "dark"),
        ("Other calib.", "other")]
i=0
for text, curmode in all_modes:
    b = tk.Radiobutton(frame_select_datatype, text=text,
                        variable=show_datatype, value=curmode, command=ChangeSummaryType)
    b.grid(column=i,row=0,padx=4,in_=frame_select_datatype)
    i+=1

show_red_datatype = tk.StringVar()
show_red_datatype.set("obj") # initialize
all_modes = [("Obj data", "obj"),
        ("Standard star", "star"),
        ("Calibration", "calibration"),
        ("Auxiliary", "auxiliary")]
i=0
for text, curmode in all_modes:
    b = tk.Radiobutton(frame_select_red_datatype, text=text,
                        variable=show_red_datatype, value=curmode, command=ChangeSummaryType)
    b.grid(column=i,row=0,padx=4,in_=frame_select_red_datatype)
    i+=1


class McListBox(object):
    """use a ttk.TreeView as a multicolumn ListBox"""

    def __init__(self, label='',header=[],data=[],rootframe=None,table_width=table_width,
                 table_height=table_height,dispers_selector=None,dispers_values=["All"],
                 binning_selector=None, binning_values=["All"],
                 pa_selector=None, pa_values=["All"],
                 slit_selector=None, slit_values=["All"],
                 obj_selector=None,obj_values=["All"], reddata=False):
        self.label=label
        self.frame_width=table_width
        self.frame_height = table_height
        self.header=header
        self.data=data
        self.sort_col=None
        self.sort_descending=0
        self.root=rootframe
        self.tree = None
        self.obj_selector=obj_selector
        self.obj_values=obj_values
        self.obj_list = None
        self.dispers_selector = dispers_selector
        self.dispers_values = dispers_values
        self.dispers_list=None
        self.binning_selector = binning_selector
        self.binning_values = binning_values
        self.binning_list = None
        self.slit_selector = slit_selector
        self.slit_values = slit_values
        self.slit_list = None
        self.pa_selector = pa_selector
        self.pa_values = pa_values
        self.pa_list = None
        self.do_process=None
        self._setup_widgets(reddata)
        self._build_tree()


    def set_droplist(self, cur_var, value):
        global summary_tab
        cur_var.set(value)
        def get_selector(selector):
            if selector:
                return selector.get()
            else:
                return "All"
        do_Summary(summary_tab, disperser=get_selector(self.dispers_selector),
                   object=get_selector(self.obj_selector), slit=get_selector(self.slit_selector), pa=get_selector(self.pa_selector),
                   binning=get_selector(self.binning_selector))
    def update_droplists(self, dispersers=[], objects=[], binnings=[], slits=[], pas=[]):
        if self.obj_selector:
            menu = self.obj_list["menu"]
            menu.delete(0, "end")
            for string in ["All"]+objects:
                menu.add_command(label=string,
                                 command=lambda value=string: self.set_droplist(self.obj_selector, value))
        if self.dispers_selector:
            menu = self.dispers_list["menu"]
            menu.delete(0, "end")
            for string in ["All"]+dispersers:
                menu.add_command(label=string,
                                 command=lambda value=string: self.set_droplist(self.dispers_selector, value))
        if self.binning_selector:
            menu = self.binning_list["menu"]
            menu.delete(0, "end")
            for string in ["All"] + binnings:
                menu.add_command(label=string,
                                 command=lambda value=string: self.set_droplist(self.binning_selector, value))
        if self.slit_selector:
            menu = self.slit_list["menu"]
            menu.delete(0, "end")
            for string in ["All"] + slits:
                menu.add_command(label=string,
                                 command=lambda value=string: self.set_droplist(self.slit_selector, value))
        if self.pa_selector:
            menu = self.pa_list["menu"]
            menu.delete(0, "end")
            for string in ["All"] + pas:
                menu.add_command(label=string,
                                 command=lambda value=string: self.set_droplist(self.pa_selector, value))

    def update(self, data=[],do_process=None):
        for el in self.tree.get_children():
            self.tree.delete(el)
        self.data = data
        self.do_process=do_process
        self._build_tree()
        if self.sort_col is not None and self.data:
            self.sortby(col=self.sort_col,descending=self.sort_descending)
    def _setup_widgets(self,reddata=False):
        cover_frame=tk.Frame(self.root,width=self.frame_width,height=self.frame_height)
        cover_frame.pack_propagate(False)
        cover_frame.pack(fill='x',expand=True)

        selector_frame = tk.Frame(cover_frame, width=self.frame_width, height=30)
        selector_frame.pack(fill='x')
        selector_frame.pack_propagate(False)

        selector_du_frame = tk.Frame(cover_frame, width=self.frame_width, height=30)
        selector_du_frame.pack(fill='x')
        selector_du_frame.pack_propagate(False)
        if self.obj_selector:
            msg = tk.Label(selector_frame, text="Object:")
            msg.pack(side='left')
            self.obj_list = tk.OptionMenu(selector_frame, self.obj_selector, *self.obj_values)
            self.obj_list.pack(side='left', padx=(10, 50))
        if self.dispers_selector:
            msg = tk.Label(selector_frame, text="Disperser:")
            msg.pack(side='left')
            self.dispers_list = tk.OptionMenu(selector_frame, self.dispers_selector, *self.dispers_values)
            self.dispers_list.pack(side='left', padx=(10, 50))
        if self.binning_selector:
            msg = tk.Label(selector_frame, text="Binning:")
            msg.pack(side='left')
            self.binning_list = tk.OptionMenu(selector_frame, self.binning_selector, *self.binning_values)
            self.binning_list.pack(side='left', padx=(10, 50))
        if self.pa_selector:
            msg = tk.Label(selector_frame, text="PA:")
            msg.pack(side='left')
            self.pa_list = tk.OptionMenu(selector_frame, self.pa_selector, *self.pa_values)
            self.pa_list.pack(side='left', padx=(10, 50))
        if self.slit_selector:
            msg = tk.Label(selector_frame, text="Slit:")
            msg.pack(side='left')
            self.slit_list = tk.OptionMenu(selector_frame, self.slit_selector, *self.slit_values)
            self.slit_list.pack(side='left', padx=(10, 50))
        if reddata:
            self.redstate = tk.StringVar()
            self.redstate.set("ini")  # initialize
            all_modes = [("Initial", "ini"),
                         ("Normalized", "normalized"),
                         ("Linearized", "linearized"),
                         ("Sky Subtracted", "sky-subtracted"),
                         ("Combined", "combined"),
                         ("Calibrated", "calibrated"),
                         ("Other","unknown")]
            msg = tk.Label(selector_frame, text="Process state:")
            msg.pack(side='left', padx=(5, 10))
            for text, curmode in all_modes:
                b = tk.Radiobutton(selector_frame, text=text,
                                   variable=self.redstate, value=curmode, command=ChangeSummaryType)
                b.pack(side='left', padx=(0, 5))

            self.redsubtype = tk.StringVar()
            self.redsubtype.set("data")  # initialize
            all_modes = [("Data", "data"),
                         ("Uncertainties", "unc")]
            msg = tk.Label(selector_du_frame, text="File type:")
            msg.pack(side='left', padx=(5, 10))
            for text, curmode in all_modes:
                b = tk.Radiobutton(selector_du_frame, text=text,
                                   variable=self.redsubtype, value=curmode, command=ChangeSummaryType)
                b.pack(side='left', padx=(0, 5))

        msg = ttk.Label(cover_frame,wraplength="4i", justify="center", anchor="center",
            padding=(2, 2, 2, 2), text=self.label)
        msg.pack(side='top',fill='x',expand=True)
        container = ttk.Frame(cover_frame)
        container.pack(fill='x', expand=True, side='left')
        # create a treeview with dual scrollbars
        style = ttk.Style()
        # style.configure(".", font=('Helvetica', 8), foreground="white")
        # style.configure("Treeview", foreground='red')
        style.configure("Treeview.Heading", foreground='blue')

        self.tree = ttk.Treeview(container,columns=self.header, show="headings",height=table_height)
        self.tree.pack(fill='x')
        vsb = ttk.Scrollbar(orient="vertical",
            command=self.tree.yview)
        hsb = ttk.Scrollbar(orient="horizontal",
            command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set,
            xscrollcommand=hsb.set)
        self.tree.grid(column=0, row=0, sticky='nsew', in_=container)
        vsb.grid(column=1, row=0, sticky='ns', in_=container)
        hsb.grid(column=0, row=1, sticky='ew', in_=container)
        container.grid_columnconfigure(0, weight=1)
        container.grid_rowconfigure(0, weight=1)
    def _build_tree(self):
        for col in self.header:
            self.tree.heading(col, text=col.title(),
                command=lambda c=col: self.sortby(c, 0))
            # adjust the column's width to the header string
            self.tree.column(col,
                width=tkFont.Font().measure(col.title()))
        for item in self.data:
            self.tree.insert('', 'end', values=item)

            # adjust column's width if necessary to fit each value
            for ix, val in enumerate(item):
                col_w = tkFont.Font().measure(val)+20
                if tkFont.Font().measure(self.header[ix].title()) > col_w:
                    col_w=tkFont.Font().measure(self.header[ix].title())
                #if self.tree.column(self.header[ix],width=None)<col_w:
                self.tree.column(self.header[ix], width=col_w)
        self.tree.bind('<ButtonRelease>',self.update_selection)
        self.tree.bind('<Down>', self.update_selection)
        self.tree.bind('<Up>', self.update_selection)

    def update_selection(self,event):
        global do_process_val
        if self.do_process:
            item = self.tree.focus()
            do_process_val.set(str(self.do_process[self.tree.index(item)]))
#
    def sortby(self, col, descending):
        """sort tree contents when a column header is clicked on"""
        self.sort_col=col
        # grab values to sort
        data = [(self.tree.set(child, col), child) \
            for child in self.tree.get_children('')]
        # if the data to be sorted is numeric change to float
        #data =  change_numeric(data)
        # now sort the data in place
        srtarr=[srtval for srtval in sorted(enumerate(data),key=lambda i:i[1],reverse=descending)]
        #data.sort(reverse=descending)
        sort_order=[]
        for ix, item in enumerate(srtarr):
            sort_order.append(item[0])
            self.tree.move(item[1][1], '', ix)

        if self.data:
            self.data = list(itemgetter(*sort_order)(self.data))
        if self.do_process:
            self.do_process = list(itemgetter(*sort_order)(self.do_process))
        # switch the heading so it will sort in the opposite direction
        self.sort_descending = descending
        new_desc=int(not self.sort_descending)
        self.tree.heading(col, command=lambda col=col: self.sortby(col, new_desc))



frame_table_summary={"obj":None,"flat":None,"bias":None,"dark":None,"neon":None,"star":None,"red":None,"sunsky":None,
                     "other":None}

for key in frame_table_summary.keys():
    frame_table_summary[key]=tk.Frame(frame_data)
    frame_table_summary[key].grid(row=2,column=0)


selected_object = tk.StringVar()
selected_object.set("All") # default value
selected_disperser = tk.StringVar()
selected_disperser.set("All")
selected_binning = tk.StringVar()
selected_binning.set("All")
selected_slit = tk.StringVar()
selected_slit.set("All")
selected_pa = tk.StringVar()
selected_pa.set("All")

obj_list=McListBox(rootframe=frame_table_summary["obj"],
                   label='Object frames',header=["Object","Disperser","Pos.angle","Slit","Exposure","Binning",'Date',"Directory","Filename","Process"],
                   data=[""], obj_selector=selected_object,dispers_selector=selected_disperser,
                        binning_selector=selected_binning, pa_selector=selected_pa, slit_selector = selected_slit)


flat_list=McListBox(rootframe=frame_table_summary["flat"],label='Flat frames',
                    header=["Disperser","Slit","Exposure","Binning",'Date',"Directory","Filename","Process"],
                    data=[""], obj_selector=selected_object,dispers_selector=selected_disperser,
                    binning_selector=selected_binning,slit_selector=selected_slit)
bias_list = McListBox(rootframe=frame_table_summary["bias"],label='Bias frames',binning_selector=selected_binning,
                      header=["Binning",'Date',"Directory","Filename","Process"],data=[""])
dark_list = McListBox(rootframe=frame_table_summary["dark"],label='Dark frames',binning_selector=selected_binning,
                      header=["Exposure","Binning",'Date',"Directory","Filename","Process"],data=[""],)
neon_list = McListBox(rootframe=frame_table_summary["neon"],label='Cal.lamp frames',
                      header=["Disperser","Slit","Exposure","Binning",'Date',"Directory","Filename","Process"],data=[""],
                      dispers_selector=selected_disperser,binning_selector=selected_binning,slit_selector=selected_slit)
sunsky_list = McListBox(rootframe=frame_table_summary["sunsky"],label='Sunsky frames',
                        header=["Disperser","Slit","Exposure","Binning",'Date',"Directory","Filename","Process"],data=[""],
                        obj_selector=selected_object,dispers_selector=selected_disperser,
                        binning_selector=selected_binning, slit_selector = selected_slit)
other_list = McListBox(rootframe=frame_table_summary["other"],label='Other calibrations',
                      header=["Type","Disperser","Slit","Exposure","Binning",'Date',"Directory","Filename","Process"],data=[""],
                       obj_selector=selected_object, dispers_selector=selected_disperser,
                       binning_selector=selected_binning, slit_selector=selected_slit)
standard_list=McListBox(rootframe=frame_table_summary["star"],label='Standard frames',
                        header=["Object","Disperser","Pos.angle","Slit", "Exposure","Binning",'Date',"Directory","Filename","Process"],
                        data=[""], obj_selector=selected_object,pa_selector=selected_pa,dispers_selector=selected_disperser,
                        binning_selector=selected_binning, slit_selector = selected_slit)

red_list_header={"Data":0,"Slit":1,"Disperser":2,"Binning":3,"Type":4,"State":5,"Subtype":6,'Directory':7,"Filename":8,"Process":9}
red_list=McListBox(rootframe=frame_table_summary["red"],label='In-processing data',
                   header=sorted(red_list_header, key=red_list_header.get),
                   data=[""],reddata=True)




view_and_mask_panel=tk.Frame(frame_data)
view_and_mask_panel.grid(row=3,column=0)

view_buttion=tk.Button(view_and_mask_panel,text="Display Image",command=do_Display,width=button_width)
view_buttion.grid(row=0,column=0,columnspan=2,padx=5)



starobj_selector_button=tk.Button(view_and_mask_panel,text="Use as Standard",command=ChangeTypeOfFile,
                                state="disabled",width=button_width)
starobj_selector_button.grid(row=0,column=3,columnspan=2,padx=5)

do_process_val = tk.BooleanVar()
do_process_val.set(True)
do_process_but = tk.Checkbutton(view_and_mask_panel, text="Process", variable=do_process_val,
                                   command=ChangeProcess)
do_process_but.grid(row=0,column=5,columnspan=2,padx=3)
do_process_but.configure(state='disabled')

ChangeSummaryType()




def IntFieldTestVal(inStr,acttyp):
    if acttyp == '1': #insert
        if not inStr.isdigit():
            return False
    return True

def FloatFieldTestVal(inStr,acttyp):
    if acttyp == '1': #insert
        try:
            float(inStr)
            return True
        except ValueError:
            return False
    return True

edit_params_frame=tk.Frame(frame_data)
edit_params_frame.grid(row=4,column=0)
tk.Label(edit_params_frame,text="Edit parameters of selected files (if necessary):",font=wid_subhead_font).grid(column=0,row=0,columnspan=13)
r=1
c=0
tk.Label(edit_params_frame,text="Disperser: ").grid(column=c,row=r)
c+=1
edit_params_disperser=tk.Entry(edit_params_frame,width=8)
edit_params_disperser.grid(column=c,row=r)
c+=1
tk.Label(edit_params_frame,text="Slit: ").grid(column=c,row=r)
c+=1
edit_params_slit=tk.Entry(edit_params_frame,width=4,validate='key')
edit_params_slit.grid(column=c,row=r)
c+=1
tk.Label(edit_params_frame,text="Binning: ").grid(column=c,row=r)
c+=1
edit_params_binning=tk.Entry(edit_params_frame,width=4)
edit_params_binning.grid(column=c,row=r)
c+=1
tk.Label(edit_params_frame,text="Object: ").grid(column=c,row=r)
c+=1
edit_params_name=tk.Entry(edit_params_frame,width=8)
edit_params_name.grid(column=c,row=r)
c+=1
tk.Label(edit_params_frame,text="PA: ").grid(column=c,row=r)
c+=1
edit_params_pa=tk.Entry(edit_params_frame,width=5,validate='key')
edit_params_pa.grid(column=c,row=r)
c+=1
tk.Label(edit_params_frame,text="DataType: ").grid(column=c,row=r)
c+=1
edit_params_type=tk.Entry(edit_params_frame,width=5)
edit_params_type.grid(column=c,row=r)
c+=1
edit_params_button=tk.Button(edit_params_frame,text="Edit parameters",command=do_EditData)
edit_params_button.grid(column=c,row=r)

edit_params_slit["validatecommand"]=(edit_params_slit.register(FloatFieldTestVal), '%P', '%d')
edit_params_pa["validatecommand"]=(edit_params_pa.register(FloatFieldTestVal), '%P', '%d')
#
# #############################
# # Parameters of data reduction section
# ############################


#==== Initial preparation (overscan, gain, lacosmic) =====
ccd_and_iniproc_section=tk.Frame(frame_parameters, highlightbackground="silver",highlightthickness=1)
tk.Label(ccd_and_iniproc_section,text="Initial CCD and reduction parameters:",
         font=wid_head_font,fg=wid_head_color).grid(column=0,row=0)
overscan_section=tk.Frame(ccd_and_iniproc_section)
overscan_section.grid(row=1, column=0)
tk.Label(overscan_section,text="Cut frames:", font=wid_subhead_font,fg=wid_subhead_color).grid(column=0,row=0,columnspan=10)
tk.Label(overscan_section,text="From left:").grid(column=0,row=1)
tk.Label(overscan_section,text="right:").grid(column=2,row=1)
tk.Label(overscan_section,text="bottom:").grid(column=4,row=1)
tk.Label(overscan_section,text="top:").grid(column=6,row=1)
field_cutframe_x0=tk.Entry(overscan_section,width=5, validate="key")
field_cutframe_x0.insert(0,0)
field_cutframe_x1=tk.Entry(overscan_section,width=5, validate="key")
field_cutframe_x1.insert(0,0)
field_cutframe_y0=tk.Entry(overscan_section,width=5, validate="key")
field_cutframe_y0.insert(0,0)
field_cutframe_y1=tk.Entry(overscan_section,width=5, validate="key")
field_cutframe_y1.insert(0,0)
for field in [field_cutframe_x0,field_cutframe_x1,field_cutframe_y1,field_cutframe_y1]:
    field["validatecommand"]=(field.register(IntFieldTestVal), '%P', '%d')

cut_as_overcan_val = tk.BooleanVar()
cut_as_overcan_val.set(False)
cut_as_overcan = tk.Checkbutton(overscan_section, text="Use as overscan", variable=cut_as_overcan_val)
cut_as_overcan.configure(state='disabled')

for wid, widget in enumerate([field_cutframe_x0,field_cutframe_x1,field_cutframe_y0,
                              field_cutframe_y1]):
    widget.grid(column=2*wid+1, row=1, padx=3)
cut_as_overcan.grid(column=8,row=1, columnspan=2)

iniprocess_section=tk.Frame(ccd_and_iniproc_section)
iniprocess_section.grid(row=2, column=0,sticky='w')
tk.Label(iniprocess_section,text="Initial preparation:", font=wid_subhead_font,
                    fg=wid_subhead_color).grid(column=0,row=0,columnspan=8)
tk.Label(iniprocess_section,text="Gain:").grid(column=0,row=1)
tk.Label(iniprocess_section,text="RNe:").grid(column=3,row=1)
field_gain=tk.Entry(iniprocess_section,width=3, validate="key")
field_gain.grid(row=1,column=1)
field_gain.insert(0,1)
field_gain["validatecommand"]=(field_gain.register(FloatFieldTestVal), '%P', '%d')
field_readnoise=tk.Entry(iniprocess_section,width=3, validate="key")
field_readnoise.grid(row=1,column=4)
field_readnoise.insert(0,2.5)
field_readnoise["validatecommand"]=(field_gain.register(FloatFieldTestVal), '%P', '%d')

auto_gain_val = tk.BooleanVar()
auto_gain_val.set(True)
auto_gain = tk.Checkbutton(iniprocess_section, text="Get gain from header", variable=auto_gain_val, command=ChangeParamSensibility)
auto_gain.grid(row=2,column=0,columnspan=2)
field_gain.configure(state='disabled')

auto_readnoise_val = tk.BooleanVar()
auto_readnoise_val.set(True)
auto_readnoise = tk.Checkbutton(iniprocess_section, text="Get RNe from bias", variable=auto_readnoise_val, command=ChangeParamSensibility)
auto_readnoise.grid(row=2,column=3,columnspan=2)
field_readnoise.configure(state='disabled')

use_lacosmic_val = tk.BooleanVar()
use_lacosmic_val.set(True)
use_lacosmic = tk.Checkbutton(iniprocess_section, text="Use LACOSMIC for CR clean", variable=use_lacosmic_val)
use_lacosmic.grid(row=1,column=6,columnspan=2)


ccd_and_iniproc_section.grid(column=0,row=0,pady=5)

#==== Geometry =====
geometry_section=tk.Frame(frame_parameters,bg='white',highlightbackground="silver",highlightthickness=1)
lab_geometry=tk.Label(geometry_section,text="Geometry control parameters:", font=wid_head_font,fg=wid_head_color)
lab_geometry.grid(column=0,row=0,columnspan=8)

r=1
c=0
tk.Label(geometry_section, text="Find peaks of calibration lines:",
      font=wid_subhead_font,fg=wid_subhead_color).grid(column=0,row=r,columnspan=6)
r+=1
tk.Label(geometry_section, text="S/N rat:").grid(column=c,row=r)
c+=1
geometry_snr=tk.Entry(geometry_section,width=3, validate="key")
geometry_snr.insert(0,20)
geometry_snr.grid(column=c,row=r)
c+=1
tk.Label(geometry_section, text="Win. for sum:").grid(column=c,row=r)
c+=1
geometry_win=tk.Entry(geometry_section,width=3, validate="key")
geometry_win.insert(0,3)
geometry_win.grid(column=c,row=r)
c+=1
tk.Label(geometry_section, text="Num. of rows:").grid(column=c,row=r)
c+=1
geometry_nstripes=tk.Entry(geometry_section,width=3, validate="key")
geometry_nstripes.insert(0,20)
geometry_nstripes.grid(column=c,row=r)

c=0
r+=1
tk.Label(geometry_section, text="Peaks prominence from:").grid(column=c,row=r)#,columnspan=3)
c+=1
geometry_prominance_0=tk.Entry(geometry_section,width=5, validate="key")
geometry_prominance_0.insert(0,1)
geometry_prominance_0.grid(column=c,row=r)
c+=1
tk.Label(geometry_section, text="to (0 for Any):").grid(column=c,row=r)
c+=1
geometry_prominance_1=tk.Entry(geometry_section,width=5, validate="key")
geometry_prominance_1.insert(0,0)
geometry_prominance_1.grid(column=c,row=r)
c+=1
tk.Label(geometry_section, text="Min. distance:").grid(column=c,row=r)
c+=1
geometry_peakdist=tk.Entry(geometry_section,width=3, validate="key")
geometry_peakdist.insert(0,30)
geometry_peakdist.grid(column=c,row=r)
c=0
r+=1
tk.Label(geometry_section, text="Extend geometry for whole spectrum:",
        font=wid_subhead_font,fg=wid_subhead_color).grid(column=0,row=r,columnspan=6)
r+=1
c=0
tk.Label(geometry_section, text="Lines polinomial order:").grid(column=c,row=r)
c+=1
geometry_poly_lines=tk.Entry(geometry_section,width=3, validate="key")
geometry_poly_lines.insert(0,2)
geometry_poly_lines.grid(column=c,row=r)
c+=1
geometry_trace_star_val = tk.BooleanVar()
geometry_trace_star_val.set(True)
geometry_trace_star = tk.Checkbutton(geometry_section, text="Calc gradient using star", variable=geometry_trace_star_val)
geometry_trace_star.grid(column=c,row=r)

c+=2
geometry_plot_results_val = tk.BooleanVar()
geometry_plot_results_val.set(True)
geometry_plot_results = tk.Checkbutton(geometry_section, text="Plot results", variable=geometry_plot_results_val)
geometry_plot_results.grid(column=c,row=r)
c+=1
geometry_do_test_val = tk.BooleanVar()
geometry_do_test_val.set(True)
geometry_do_test = tk.Checkbutton(geometry_section, text="Save test", variable=geometry_do_test_val)
geometry_do_test.grid(column=c,row=r)

for field in [geometry_poly_lines, geometry_nstripes,geometry_win]:
    field["validatecommand"]=(field.register(IntFieldTestVal), '%P', '%d')
for field in [geometry_snr,geometry_prominance_0,geometry_prominance_1]:
    field["validatecommand"]=(field.register(FloatFieldTestVal), '%P', '%d')

geometry_section.grid(column=0,row=1,pady=5)



#==== Relative shifts and atmospheric dispersion =====
atmdisp_section=tk.Frame(frame_parameters,bg='white',highlightbackground="silver",highlightthickness=1)
lab_atmdisp=tk.Label(atmdisp_section,text="Relative shifts and atmospheric dispersion:", font=wid_head_font,fg=wid_head_color)
lab_atmdisp.grid(column=0,row=0,columnspan=10)
#
r=1
c=0

tk.Label(atmdisp_section, text="Max shift along X axis:").grid(column=c,row=r,columnspan=2)
c+=2
atmdisp_maxshift_x=tk.Entry(atmdisp_section,width=3, validate="key")
atmdisp_maxshift_x.insert(0,2)
atmdisp_maxshift_x.grid(column=c,row=r)
c+=1
tk.Label(atmdisp_section, text="Y axis:").grid(column=c,row=r, columnspan=1)
c+=1
atmdisp_maxshift_y=tk.Entry(atmdisp_section,width=3, validate="key")
atmdisp_maxshift_y.insert(0,5)
atmdisp_maxshift_y.grid(column=c,row=r)
c+=1
tk.Label(atmdisp_section,text='Index of reference file:').grid(column=c,row=r,columnspan=2)
c+=2
atmdisp_ref_index=tk.Entry(atmdisp_section,width=3, validate="key")
atmdisp_ref_index.grid(column=c,row=r)
r+=1
c=0
tk.Label(atmdisp_section, text="Derive X-shift for:").grid(column=c,row=r,columnspan=2)
c+=2
atmdisp_findshift_x_obj_val = tk.BooleanVar()
atmdisp_findshift_x_obj_val.set("True")
atmdisp_findshift_x_star_val = tk.BooleanVar()
atmdisp_findshift_x_star_val.set("True")
atmdisp_findshift_x_obj=tk.Checkbutton(atmdisp_section, text="object", variable=atmdisp_findshift_x_obj_val)
atmdisp_findshift_x_obj.grid(column=c,row=r)
c+=1
atmdisp_findshift_x_star=tk.Checkbutton(atmdisp_section, text="star", variable=atmdisp_findshift_x_star_val)
atmdisp_findshift_x_star.grid(column=c,row=r)
c+=1
tk.Label(atmdisp_section, text="Y-shift for:").grid(column=c,row=r,columnspan=1)
c+=1
atmdisp_findshift_y_obj_val = tk.BooleanVar()
atmdisp_findshift_y_obj_val.set("True")
atmdisp_findshift_y_star_val = tk.BooleanVar()
atmdisp_findshift_y_star_val.set("True")
atmdisp_findshift_y_obj=tk.Checkbutton(atmdisp_section, text="object", variable=atmdisp_findshift_y_obj_val)
atmdisp_findshift_y_obj.grid(column=c,row=r)
c+=1
atmdisp_findshift_y_star=tk.Checkbutton(atmdisp_section, text="star", variable=atmdisp_findshift_y_star_val)
atmdisp_findshift_y_star.grid(column=c,row=r)

r+=1
atmdisp_obj_apply_val = tk.BooleanVar()
atmdisp_obj_apply_val.set("True")
atmdisp_star_apply_val = tk.BooleanVar()
atmdisp_star_apply_val.set("True")
atmdisp_plot_val = tk.BooleanVar()
atmdisp_plot_val.set("True")
atmdisp_do_test_val=tk.BooleanVar()
atmdisp_do_test_val.set("True")
c=0
tk.Label(atmdisp_section,text="Trace atm. disp. with for", font=wid_subhead_font,fg=wid_subhead_color).grid(column=c,row=r,columnspan=2)
c+=2
atmdisp_obj_apply=tk.Checkbutton(atmdisp_section, text="object", variable=atmdisp_obj_apply_val, command=ChangeParamSensibility)
atmdisp_obj_apply.grid(column=c,row=r)
c+=1
atmdisp_star_apply=tk.Checkbutton(atmdisp_section, text="star", variable=atmdisp_star_apply_val, command=ChangeParamSensibility)
atmdisp_star_apply.grid(column=c,row=r, columnspan=1)
c+=2
atmdisp_plot=tk.Checkbutton(atmdisp_section, text="Plot results", variable=atmdisp_plot_val)
atmdisp_plot.grid(column=c,row=r,columnspan=1)
c+=1
atmdisp_do_test=tk.Checkbutton(atmdisp_section, text="Save test", variable=atmdisp_do_test_val)
atmdisp_do_test.grid(column=c,row=r,columnspan=1)

r+=1
c=0
tk.Label(atmdisp_section, text="Pos of object (0 for auto) on OBJ frame:").grid(column=c,row=r, columnspan=3)
c+=3
atmdisp_obj_ypos=tk.Entry(atmdisp_section,width=4, validate="key")
atmdisp_obj_ypos.insert(0,0)
atmdisp_obj_ypos.grid(column=c,row=r)

c+=1
tk.Label(atmdisp_section, text="on STAR frame:").grid(column=c,row=r)
c+=1
atmdisp_star_ypos=tk.Entry(atmdisp_section,width=4, validate="key")
atmdisp_star_ypos.insert(0,0)
atmdisp_star_ypos.grid(column=c,row=r, columnspan=1)

c+=1
tk.Label(atmdisp_section, text="Poly order:").grid(column=c,row=r)
c+=1
atmdisp_order=tk.Entry(atmdisp_section,width=2, validate="key")
atmdisp_order.insert(0,3)
atmdisp_order.grid(column=c,row=r)



for field in [atmdisp_star_ypos, atmdisp_obj_ypos, atmdisp_order,atmdisp_ref_index]:
    field["validatecommand"]=(field.register(IntFieldTestVal), '%P', '%d')
for field in [atmdisp_maxshift_x, atmdisp_maxshift_y]:
    field["validatecommand"]=(field.register(FloatFieldTestVal), '%P', '%d')
atmdisp_section.grid(column=0,row=2,pady=5)



#==== Dispersion curve ===
dispcurve_section=tk.Frame(frame_parameters,bg='white',highlightbackground="silver",highlightthickness=1)
lab_dc=tk.Label(dispcurve_section,text="Parameters of dispersion curve:", font=wid_head_font,fg=wid_head_color)
lab_dc.grid(column=0,row=0,columnspan=8)

r=1
c=0
tk.Label(dispcurve_section, text="Polynomial order along dispersion:").grid(column=c,row=r, columnspan=3)
c+=3
dispcurve_order_lam=tk.Entry(dispcurve_section,width=3, validate="key")
dispcurve_order_lam.insert(0,5)
dispcurve_order_lam.grid(column=c,row=r)
c+=1
tk.Label(dispcurve_section, text="along slit:").grid(column=c,row=r, columnspan=2)
c+=2
dispcurve_order_y=tk.Entry(dispcurve_section,width=3, validate="key")
dispcurve_order_y.insert(0,2)
dispcurve_order_y.grid(column=c,row=r)
c+=1
tk.Label(dispcurve_section, text="Window along slit:").grid(column=c,row=r, columnspan=2)
c+=2
dispcurve_window_y=tk.Entry(dispcurve_section,width=3, validate="key")
dispcurve_window_y.insert(0,5)
dispcurve_window_y.grid(column=c,row=r)
# c+=1
# Label(geometry_section, text="Num. of rows:").grid(column=c,row=r)
# c+=1

for field in [dispcurve_order_y,dispcurve_order_lam, dispcurve_window_y]:
    field["validatecommand"]=(field.register(IntFieldTestVal), '%P', '%d')

dispcurve_section.grid(column=0,row=3,pady=5)


#==== Linearization =====
linearizarion_section=tk.Frame(frame_parameters,bg='white',highlightbackground="silver",highlightthickness=1)
lab_lin=tk.Label(linearizarion_section,text="Linearization parameters:", font=wid_head_font,fg=wid_head_color)
lab_lin.grid(column=0,row=0,columnspan=9)
#
r=1
c=0
tk.Label(linearizarion_section, text="Dispersion A/px:").grid(column=c,row=r, columnspan=2)
c+=2
linearizarion_disp=tk.Entry(linearizarion_section,width=5, validate="key")
linearizarion_disp.insert(0,0.8)
linearizarion_disp.grid(column=c,row=r)
c+=1
tk.Label(linearizarion_section, text="Wavelength range:").grid(column=c,row=r, columnspan=2)
c+=2
linearizarion_lam0=tk.Entry(linearizarion_section,width=6, validate="key")
linearizarion_lam0.insert(0,3600)
linearizarion_lam0.grid(column=c,row=r)
c+=1
tk.Label(linearizarion_section, text=" - ").grid(column=c,row=r)
c+=1
linearizarion_lam1=tk.Entry(linearizarion_section,width=6, validate="key")
linearizarion_lam1.insert(0,7200)
linearizarion_lam1.grid(column=c,row=r)

c+=1
linearization_do_test_val = tk.BooleanVar()
linearization_do_test_val.set(True)
linearization_do_test = tk.Checkbutton(linearizarion_section, text="Apply for test", variable=linearization_do_test_val)
linearization_do_test.grid(column=c,row=r)



linearizarion_disp_auto_val = tk.BooleanVar()
linearizarion_lam0_auto_val = tk.BooleanVar()
linearizarion_lam1_auto_val = tk.BooleanVar()
linearizarion_disp_auto_val.set(True)
linearizarion_lam0_auto_val.set(True)
linearizarion_lam1_auto_val.set(True)
r+=1
c=2
linearizarion_disp_auto = tk.Checkbutton(linearizarion_section, text="Auto", variable=linearizarion_disp_auto_val, command=ChangeParamSensibility)
linearizarion_disp_auto.grid(column=c,row=r)
c=5
linearizarion_lam0_auto = tk.Checkbutton(linearizarion_section, text="Auto", variable=linearizarion_lam0_auto_val, command=ChangeParamSensibility)
linearizarion_lam0_auto.grid(column=c, row=r)
c=7
linearizarion_lam1_auto = tk.Checkbutton(linearizarion_section, text="Auto", variable=linearizarion_lam1_auto_val,command=ChangeParamSensibility)
linearizarion_lam1_auto.grid(column=c,row=r)




for field in [linearizarion_disp,linearizarion_lam1,linearizarion_lam0]:
    field["validatecommand"]=(field.register(FloatFieldTestVal), '%P', '%d')



linearizarion_section.grid(column=0,row=4,pady=5)

linearizarion_disp.config(state='disabled')
linearizarion_lam0.config(state='disabled')
linearizarion_lam1.config(state='disabled')




#==== Sky subtraction =====
skyrem_section=tk.Frame(frame_parameters,bg='white',highlightbackground="silver",highlightthickness=1)
lab_sky=tk.Label(skyrem_section,text="Night sky model parameters:", font=wid_head_font,fg=wid_head_color)
lab_sky.grid(column=0,row=0,columnspan=8)
#
r=1
skyrem_mode = tk.StringVar()
skyrem_mode.set("Poly") # initialize
all_skyrem_modes = [("Template", "Template"),
        ("Polynom.", "Poly")]
tk.Label(skyrem_section,text="Mode:", font=wid_subhead_font,fg=wid_subhead_color).grid(column=0,row=r)
c=1
for text, curmode in all_skyrem_modes:
    b = tk.Radiobutton(skyrem_section, text=text,
                        variable=skyrem_mode, value=curmode, command=ChangeParamSensibility)
    b.grid(column=c,row=r,padx=10,columnspan=1)
    c+=1

skyrem_fft_val=tk.BooleanVar()
skyrem_fft_val.set(False)
skyrem_fft=tk.Checkbutton(skyrem_section, text="in Fourier space", variable=skyrem_fft_val, command=ChangeParamSensibility)
skyrem_fft.grid(column=c,row=r,columnspan=1)

c+=1

skyrem_sunsky_val=tk.BooleanVar()
skyrem_sunsky_val.set(False)
skyrem_sunsky=tk.Checkbutton(skyrem_section, text="Use sunsky", variable=skyrem_sunsky_val)
skyrem_sunsky.grid(column=c,row=r,padx=10)
skyrem_sunsky.config(state='disabled')
c+=1
skyrem_do_test_val = tk.BooleanVar()
skyrem_do_test_val.set(True)
skyrem_do_test = tk.Checkbutton(skyrem_section, text="Run test only", variable=skyrem_do_test_val)
skyrem_do_test.grid(column=c,row=r,padx=10)

r+=1
c=0
tk.Label(skyrem_section, text="Polynomial order (Obj):").grid(column=c,row=r, columnspan=2)
c+=2
skyrem_poly_obj=tk.Entry(skyrem_section,width=3, validate="key")
skyrem_poly_obj.insert(0,3)
skyrem_poly_obj.grid(column=c,row=r)
c+=1
skyrem_setregion_obj=tk.Button(skyrem_section,text="Set Sky Reg. for Obj",command=lambda c="obj": SetSkyRegions(c),state="normal")#,width=button_width)
skyrem_setregion_obj.grid(row=r,column=c,columnspan=1,padx=4)
c+=1
skyrem_nregions_obj=tk.IntVar(0)
skyrem_npix_obj=tk.IntVar(0)
skyrem_nregions_label_obj=tk.Label(skyrem_section, text="N regs:{:2d}; N pix:{:4d}".format(skyrem_nregions_obj.get(),skyrem_npix_obj.get()))
skyrem_nregions_label_obj.grid(column=c,row=r, columnspan=1)
c+=1
skyrem_obj_process_val = tk.BooleanVar()
skyrem_obj_process_val.set(True)
skyrem_obj_process = tk.Checkbutton(skyrem_section, text="Process", variable=skyrem_obj_process_val, command=ChangeParamSensibility)
skyrem_obj_process.grid(column=c,row=r)

r+=1
c=0
tk.Label(skyrem_section, text="Polynomial order (Star):").grid(column=c,row=r, columnspan=2)
c+=2
skyrem_poly_star=tk.Entry(skyrem_section,width=3, validate="key")
skyrem_poly_star.insert(0,3)
skyrem_poly_star.grid(column=c,row=r)
c+=1
skyrem_setregion_star=tk.Button(skyrem_section,text="Set Sky Reg. for Star",command=lambda c="star": SetSkyRegions(c),state="normal")#,width=button_width)
skyrem_setregion_star.grid(row=r,column=c,columnspan=1, padx=4)
c+=1
skyrem_nregions_star=tk.IntVar(0)
skyrem_npix_star=tk.IntVar(0)
skyrem_nregions_label_star=tk.Label(skyrem_section, text="N regs:{:2d}; N pix:{:4d}".format(skyrem_nregions_star.get(),skyrem_npix_star.get()))
skyrem_nregions_label_star.grid(column=c,row=r, columnspan=1)
c+=1
skyrem_star_process_val = tk.BooleanVar()
skyrem_star_process_val.set(True)
skyrem_star_process = tk.Checkbutton(skyrem_section, text="Process", variable=skyrem_star_process_val, command=ChangeParamSensibility)
skyrem_star_process.grid(column=c,row=r)

r+=1
c=0
skyrem_min_med_poly_frame=tk.Frame(skyrem_section)
skyrem_min_med_poly_frame.grid(column=c,row=r,columnspan=12)
r=0
c=0
tk.Label(skyrem_min_med_poly_frame, text="Minimal order and limiting S/N (Obj):").grid(column=c,row=r, columnspan=3)
c+=3
skyrem_minorder_obj=tk.Entry(skyrem_min_med_poly_frame,width=3, validate="key")
skyrem_minorder_obj.grid(column=c,row=r)
c+=1
skyrem_minsn_obj=tk.Entry(skyrem_min_med_poly_frame,width=3, validate="key")
skyrem_minsn_obj.grid(column=c,row=r)
c+=1
tk.Label(skyrem_min_med_poly_frame, text="Medium order and limiting SN (Obj):").grid(column=c,row=r, columnspan=3)
c+=3
skyrem_medorder_obj=tk.Entry(skyrem_min_med_poly_frame,width=3, validate="key")
skyrem_medorder_obj.grid(column=c,row=r)
c+=1
skyrem_medsn_obj=tk.Entry(skyrem_min_med_poly_frame,width=3, validate="key")
skyrem_medsn_obj.grid(column=c,row=r)
c=0
r+=1
tk.Label(skyrem_min_med_poly_frame, text="Minimal order and limiting S/N (Star):").grid(column=c,row=r, columnspan=3)
c+=3
skyrem_minorder_star=tk.Entry(skyrem_min_med_poly_frame,width=3, validate="key")
skyrem_minorder_star.grid(column=c,row=r)
c+=1
skyrem_minsn_star=tk.Entry(skyrem_min_med_poly_frame,width=3, validate="key")
skyrem_minsn_star.grid(column=c,row=r)
c+=1
tk.Label(skyrem_min_med_poly_frame, text="Medium order and limiting S/N (Star):").grid(column=c,row=r, columnspan=3)
c+=3
skyrem_medorder_star=tk.Entry(skyrem_min_med_poly_frame,width=3, validate="key")
skyrem_medorder_star.grid(column=c,row=r)
c+=1
skyrem_medsn_star=tk.Entry(skyrem_min_med_poly_frame,width=3, validate="key")
skyrem_medsn_star.grid(column=c,row=r)




for field in [skyrem_poly_obj,skyrem_poly_star,skyrem_medorder_star,skyrem_minorder_star,
              skyrem_medorder_obj,skyrem_minorder_obj]:
    field["validatecommand"]=(field.register(IntFieldTestVal), '%P', '%d')
for field in [skyrem_medsn_star,skyrem_minsn_star,skyrem_medsn_obj,skyrem_minsn_obj]:
    field["validatecommand"]=(field.register(FloatFieldTestVal), '%P', '%d')

skyrem_section.grid(column=0,row=5,pady=5)


#=== Combine exposures calculation ====

combine_section=tk.Frame(frame_parameters,bg='white',highlightbackground="silver",highlightthickness=1)
tk.Label(combine_section,text="Combine exposures parameters: ", font=wid_head_font,
         fg=wid_head_color).grid(column=0,row=0,columnspan=12)

r=1
c=0
tk.Label(combine_section,text="Method:",font=wid_subhead_font,fg=wid_subhead_color).grid(column=c,row=r)
c+=1
combine_mode = tk.StringVar()
combine_mode.set("Sigma") # initialize
all_skyrem_modes = [("Sigma Clip.", "Sigma"),("Median", "Median"),("Comp. pairs", "Pairs")]
c=1
for text, curmode in all_skyrem_modes:
    b = tk.Radiobutton(combine_section, text=text,
                        variable=combine_mode, value=curmode)
    b.grid(column=c,row=r,padx=10,columnspan=1)
    c+=1
tk.Label(combine_section,text="Treshold:").grid(column=c,row=r)
c+=1
combine_treshold=tk.Entry(combine_section,width=4,validate='key')
combine_treshold.grid(row=r,column=c)
combine_treshold.insert(0,2.1)

for field in [combine_treshold]:
    field["validatecommand"]=(field.register(FloatFieldTestVal), '%P', '%d')


combine_section.grid(column=0,row=6,pady=5)


#=== DQE calculation ====

dqe_section=tk.Frame(frame_parameters,bg='white',highlightbackground="silver",highlightthickness=1)
tk.Label(dqe_section,text="DQE calculation parameters (leave empty for auto-deriving): ", font=wid_head_font,
         fg=wid_head_color).grid(column=0,row=0,columnspan=12)

r=1
c=0
tk.Label(dqe_section,text="Star name:").grid(column=c,row=r)
c+=1
dqe_starname=tk.Entry(dqe_section,width=10)
dqe_starname.grid(row=r,column=c)
c+=1
tk.Label(dqe_section,text="Ypos:").grid(column=c,row=r)
c+=1
dqe_extract_bottom=tk.Entry(dqe_section,width=4,validate='key')
dqe_extract_bottom.grid(row=r,column=c)
c+=1
tk.Label(dqe_section,text="-").grid(column=c,row=r)
c+=1
dqe_extract_top=tk.Entry(dqe_section,width=4,validate='key')
dqe_extract_top.grid(row=r,column=c)
c+=1
tk.Label(dqe_section,text="Z Star:").grid(column=c,row=r)
c+=1
dqe_starz=tk.Entry(dqe_section,width=3,validate='key')
dqe_starz.grid(row=r,column=c)
c+=1
tk.Label(dqe_section,text="Smooth Win:").grid(column=c,row=r)
c+=1
dqe_smooth=tk.Entry(dqe_section,width=4,validate='key')
dqe_smooth.grid(row=r,column=c)
c+=2

dqe_plot_val=tk.BooleanVar()
dqe_plot_val.set(True)
dqe_plot=tk.Checkbutton(dqe_section, text="Show DQE", variable=dqe_plot_val)
dqe_plot.grid(column=c,row=r,columnspan=1)

for field in [dqe_smooth,dqe_extract_bottom,dqe_extract_top]:
    field["validatecommand"]=(field.register(IntFieldTestVal), '%P', '%d')
dqe_starz["validatecommand"]=(dqe_starz.register(FloatFieldTestVal), '%P', '%d')

dqe_section.grid(column=0,row=7,pady=5)

#=== Flux calibration ====

fluxcal_section=tk.Frame(frame_parameters,bg='white',highlightbackground="silver",highlightthickness=1)
tk.Label(fluxcal_section,text="Flux calibration parameters (leave empty for auto-deriving):", font=wid_head_font,
         fg=wid_head_color).grid(column=0,row=0,columnspan=10)

r=1
c=0
tk.Label(fluxcal_section,text="Z Obj:").grid(column=c,row=r)
c+=1
fluxcal_objz=tk.Entry(fluxcal_section,width=3,validate='key')
fluxcal_objz.grid(row=r,column=c)
c+=1
tk.Label(fluxcal_section,text="Ypos of zero:").grid(row=r,column=c)
c+=1
fluxcal_ypos=tk.Entry(fluxcal_section,width=4,validate='key')
fluxcal_ypos.grid(row=r,column=c)
c+=1
tk.Label(fluxcal_section,text="Pixsize:").grid(row=r,column=c)
c+=1
fluxcal_ydelt=tk.Entry(fluxcal_section,width=4,validate='key')
fluxcal_ydelt.grid(row=r,column=c)
r+=1
c=1
tk.Label(fluxcal_section,text="Cut Wavelength:").grid(column=c,row=r)
c+=1
fluxcal_lam0=tk.Entry(fluxcal_section,width=5,validate='key')
fluxcal_lam0.grid(row=r,column=c)
c+=1
tk.Label(fluxcal_section,text=" - ").grid(column=c,row=r)
c+=1
fluxcal_lam1=tk.Entry(fluxcal_section,width=5,validate='key')
fluxcal_lam1.grid(row=r,column=c)
c+=1
tk.Label(fluxcal_section,text="Cut along slit:").grid(column=c,row=r)
c+=1
fluxcal_ycut0=tk.Entry(fluxcal_section,width=5,validate='key')
c+=1
fluxcal_ycut0.grid(row=r,column=c)
c+=1
tk.Label(fluxcal_section,text=" - ").grid(column=c,row=r)
c+=1
fluxcal_ycut1=tk.Entry(fluxcal_section,width=5,validate='key')
c+=1
fluxcal_ycut1.grid(row=r,column=c)
for field in [fluxcal_ydelt,fluxcal_ypos, fluxcal_lam0, fluxcal_lam1, fluxcal_objz]:
    field["validatecommand"]=(field.register(FloatFieldTestVal), '%P', '%d')
for field in [fluxcal_ycut0,fluxcal_ycut1]:
    field["validatecommand"]=(field.register(IntFieldTestVal), '%P', '%d')

fluxcal_section.grid(column=0,row=8,pady=5)


#=== Set calibration to be used as default for all slits ====

default_calibrations_section=tk.Frame(frame_parameters,bg='white',highlightbackground="silver",highlightthickness=3)
tk.Label(default_calibrations_section,text="Default calibration (to be used for any slit if necessary):", font=wid_head_font,
         fg=wid_head_color).grid(column=0,row=0,columnspan=10)

r=1
c=0
tk.Label(default_calibrations_section,text="Transform:").grid(column=c,row=r)
c+=1
default_transform=tk.Entry(default_calibrations_section,width=10)
default_transform.grid(row=r,column=c)
c+=1
tk.Label(default_calibrations_section,text="Flat:").grid(column=c,row=r)
c+=1
default_flat=tk.Entry(default_calibrations_section,width=10)
default_flat.grid(row=r,column=c)
c+=1
tk.Label(default_calibrations_section,text="Disp.curve:").grid(column=c,row=r)
c+=1
default_disp_curve=tk.Entry(default_calibrations_section,width=10)
default_disp_curve.grid(row=r,column=c)
c+=1
tk.Label(default_calibrations_section,text="Sens.curve:").grid(column=c,row=r)
c+=1
default_sens_curve=tk.Entry(default_calibrations_section,width=10)
default_sens_curve.grid(row=r,column=c)

default_sens_curve_force=tk.BooleanVar()
default_sens_curve_force.set(True)
default_disp_curve_force=tk.BooleanVar()
default_disp_curve_force.set(False)
default_transform_force=tk.BooleanVar()
default_transform_force.set(False)
default_flat_force=tk.BooleanVar()
default_flat_force.set(False)

r+=1
c=0
tk.Checkbutton(default_calibrations_section,text="Force use transform",
               variable=default_transform_force).grid(row=r,column=c,columnspan=2)
c+=2
tk.Checkbutton(default_calibrations_section,text="Force use flat",
               variable=default_flat_force).grid(row=r,column=c,columnspan=2)
c+=2
tk.Checkbutton(default_calibrations_section,text="Force use disp.curve",
               variable=default_disp_curve_force).grid(row=r,column=c,columnspan=2)
c+=2
tk.Checkbutton(default_calibrations_section,text="Force use sens.curve",
               variable=default_sens_curve_force).grid(row=r,column=c,columnspan=2)

default_calibrations_section.grid(column=0,row=9,pady=5)

#############################
# Display section
############################

class DisplayPanel(tk.Frame):
    def redraw(self):
        if self.image is None:
            return
        if self.im_scale.get() == "Linear":
            f_stretch = LinearStretch()
        elif self.im_scale.get() == "Log":
            f_stretch = LogStretch()
        elif self.im_scale.get() == "Asinh":
            f_stretch = AsinhStretch()
        elif self.im_scale.get() == "Sqrt":
            f_stretch = SqrtStretch()

        vmin = float(self.im_minval.get())
        vmax = float(self.im_maxval.get())
        norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=f_stretch)
        self.ax.cla()
        self.ax.imshow(self.image, origin='lower', cmap=plt.cm.Greys, norm=norm,
                       extent=(0, self.image.shape[1], 0, self.image.shape[0]),
                       aspect='auto')
        self.ax.set_xlabel("X, pix")
        self.ax.set_ylabel("Y, pix")
        self.ax.set_title(self.label)
        self.canvas.draw()
        self.snap_cursor = self.SnaptoCursor(self.ax)
        self.ax.figure.canvas.callbacks.connect('motion_notify_event', self.snap_cursor.mouse_move)
        self.toolbar.update()

    def load_image(self, filename, label=""):
        if not os.path.isfile(filename):
            return
        with fits.open(filename) as hdu:
            objname=hdu[0].header.get('object')
            if len(hdu[0].data.shape) == 2:
                self.image = np.array(hdu[0].data,dtype=float)
                if self.cube is not None:
                    self.cube_button_incr.grid_remove()
                    self.cube_button_decr.grid_remove()
                    self.cube_channel_label.grid_remove()
                self.cube= None
                self.cube_channel = 0

            elif len(hdu[0].data.shape) == 3:
                if self.cube is None:
                    self.cube_button_incr.grid()
                    self.cube_button_decr.grid()
                    self.cube_channel_label.grid()
                self.cube=np.array(hdu[0].data,dtype=float)
                self.image=self.cube[0,:,:]
                self.cube_channel=0
                self.cube_channel_label.config(text="Channel: {:02d}".format(self.cube_channel))
            elif len(hdu[0].data.shape) == 1:
                if self.cube is not None:
                    self.cube_button_incr.grid_remove()
                    self.cube_button_decr.grid_remove()
                    self.cube_channel_label.grid_remove()
                self.cube=None
                self.image=np.array(hdu[0].data,dtype=float).reshape((1,hdu[0].data.shape[0]))
                self.cube_channel = 0
            else:
                message("Too many dimensions in this fits to show.")
                return
        self.label=label
        if objname: self.label="{} ({})".format(self.label,objname)
        if filename.casefold().endswith('_crmask.fits'):
            minval=0
            maxval=1
        else:
            minval, maxval = ZScaleInterval().get_limits(self.image)
        self.im_minval.delete(0,'end')
        self.im_minval.insert(0,minval)
        self.im_maxval.delete(0,'end')
        self.im_maxval.insert(0,maxval)
        self.redraw()

    def nextChannel(self, increment):
        if self.cube is not None:
            if 0 <= self.cube_channel+increment <= (self.cube.shape[0]-1):
                self.cube_channel+= increment
                self.image=self.cube[self.cube_channel,:,:]
                self.cube_channel_label.config(text="Channel: {:02d}".format(self.cube_channel))
                self.redraw()


    class SnaptoCursor(object):
        def __init__(self, ax):
            # Have to draw the canvas once beforehand to cache the renderer
            ax.figure.canvas.draw()
            self.bg = ax.figure.canvas.copy_from_bbox(ax.bbox)
            self.ax = ax
            self.lx = ax.axhline(color='green')  # the horiz line
            self.ly = ax.axvline(color='green')  # the vert line
            self.x = 0
            self.y = 0

        def mouse_move(self, event):
            if not event.inaxes:
                return
            x, y = event.xdata, event.ydata
            # update the line positions
            self.ax.figure.canvas.restore_region(self.bg)
            self.lx.set_ydata(y)
            self.ly.set_xdata(x)
            self.ax.draw_artist(self.lx)
            self.ax.draw_artist(self.ly)
            self.ax.figure.canvas.blit(self.ax.bbox)

    def __init__(self, parent, label="", image=None, cube=None):
        global wid_head_font
        tk.Frame.__init__(self,parent)
        self.parent=parent
        self.image=image
        self.cube = cube
        self.cube_channel = 0
        self.label=label
        self.snap_cursor=None

        display_panel = tk.Frame(self)
        display_panel.grid(row=0, column=0, columnspan=16,rowspan=8)
        control_panel = tk.Frame(self)
        control_panel.grid(row=8, column=0, columnspan=16, rowspan=1, padx=20, pady=10)

        tk.Label(control_panel, text="Brightness Min:").grid(row=1, column=0)
        tk.Label(control_panel, text="Max:").grid(row=1, column=2)
        self.im_minval = tk.Entry(control_panel, width=7, validate='key')
        self.im_maxval = tk.Entry(control_panel, width=7, validate='key')
        for field in [self.im_minval, self.im_maxval]:
            field["validatecommand"] = (field.register(FloatFieldTestVal), '%P', '%d')
            field.bind("<Return>", lambda x=None: self.redraw())
        self.im_minval.grid(row=1, column=1)
        self.im_maxval.grid(row=1, column=3)
        tk.Label(control_panel, text="      ").grid(row=1, column=4, columnspan=1)
        tk.Label(control_panel, text="Scale:", font=wid_head_font).grid(row=1, column=5, columnspan=1)
        self.im_scale = tk.StringVar()
        self.im_scale.set("Linear")  # initialize
        all_modes = ["Linear", "Log", "Asinh", "Sqrt"]
        i = 6
        for curmode in all_modes:
            b = tk.Radiobutton(control_panel, text=curmode,
                               variable=self.im_scale, value=curmode, command=self.redraw)
            b.grid(row=1, column=i, columnspan=1, sticky='nsew')
            i += 1


        fig = Figure(figsize=(11.3, 4))
        self.ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.05, top=0.93, right=0.99, bottom=0.12)

        #fig.tight_layout(pad=2)
        self.canvas = FigureCanvasTkAgg(fig, master=display_panel)  # A tk.DrawingArea.
        self.canvas.get_tk_widget().grid(row=0, column=1, columnspan=15, rowspan=5)

        toolbarFrame = tk.Frame(control_panel)
        toolbarFrame.grid(row=0, column=0, sticky=tk.NSEW, columnspan=10)
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbarFrame)
        self.toolbar.update()

        self.cube_panel=tk.Frame(display_panel)
        self.cube_panel.grid(row=0,column=0,sticky='ew')
        self.cube_channel_label=tk.Label(self.cube_panel, text="Channel: {:02d}".format(0))
        self.cube_channel_label.grid(row=0, column=0,sticky=tk.NSEW)
        self.cube_button_incr=tk.Button(self.cube_panel,text="Z+",command=lambda c=1: self.nextChannel(c),state="normal",width=5)
        self.cube_button_incr.grid(row=1,column=0)
        self.cube_button_decr=tk.Button(self.cube_panel,text="Z-",command=lambda c=-1: self.nextChannel(c),state="normal",width=5)
        self.cube_button_decr.grid(row=2,column=0)
        self.cube_button_incr.grid_remove()
        self.cube_button_decr.grid_remove()
        self.cube_channel_label.grid_remove()

display=DisplayPanel(frame_display)
display.grid(row=0,column=0)



#==========================
# Make all to be resizable
#==========================
root.resizable(height = True, width = True)
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)
frame_control.grid_rowconfigure(0, weight=2)
frame_control.grid_columnconfigure(0, weight=2)
#frame_data.grid_rowconfigure(0, weight=1)
frame_data.grid_columnconfigure(0, weight=1)
frame_display.grid_rowconfigure(0, weight=1)
frame_display.grid_columnconfigure(0, weight=1)
frame_parameters.grid_rowconfigure(0, weight=1)
frame_parameters.grid_columnconfigure(0, weight=1)
frame_buttons.grid_rowconfigure(0, weight=1)
frame_buttons.grid_columnconfigure(0, weight=1)
frame_steps.grid_rowconfigure(0, weight=1)
frame_steps.grid_columnconfigure(0, weight=1)
#frame_summary_type.grid_columnconfigure(0,weight=1)
#frame_summary_type.grid_rowconfigure(0,weight=1)
load_section.grid_columnconfigure(0,weight=1)
load_section.grid_rowconfigure(0,weight=1)
mode_section.grid_columnconfigure(0,weight=1)
mode_section.grid_columnconfigure(1,weight=1)
mode_section.grid_columnconfigure(2,weight=1)
mode_section.grid_rowconfigure(0,weight=1)

# #############################
# #  Define widgets behaviour
# #############################


field_load_raw_dir.bind("<Return>",Set_R_Dir)
field_load_w_dir.bind("<Return>",Set_W_Dir)
# field_load_calib_dir.bind("<Return>",Set_Cal_Dir)

root.mainloop()

