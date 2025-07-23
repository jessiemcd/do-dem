#Path to top-level do-dem directory - edit for your system.
path_to_dodem = '/Users/jmdunca2/do-dem/'
from sys import path as sys_path
sys_path.append(path_to_dodem+'/dodem/')

import nustar_utilities as nuutil

import pathlib
import zipfile

import matplotlib.pyplot as plt
import numpy as np
import xspec

import os

import response as rres

#Path to top-level do-dem directory - edit for your system.
path_to_dodem = '/Users/jmdunca2/do-dem/'
from sys import path as sys_path
sys_path.append(path_to_dodem+'/pyxspec_extension/')


from pyxspec_extension.interface import XSPECInterface
from pyxspec_extension.plotter import ModelPlotter
from pyxspec_extension.config import DEFAULT_PARAMETER_FILE, DEFAULT_PILEUP_PARAMETER_FILE

abundance_file = '/Users/jmdunca2/do-dem/reference_files/feld92a_coronal0.txt'


def path_obsid(thepath):

    files = thepath.glob('*A06_0_4_p_sr.pha')
    fff = [f for f in files]
    obsid = str(fff[0])[-27:-16]

    return obsid

def path_timestring(thepath):

    return str(thepath).split('/')[-1]


def all_srms(obsidstr):
    """
    As may be obvious, needs to be run in the directory with the pha and response files.
    
    """

    fpms=['A', 'B']
    grade_exps = ['0_4', '21_24']
    for fpm in fpms:
        for ge in grade_exps:
            rres.make_srm_file(
                'nu'+obsidstr+fpm+'06_'+ge+'_p_sr.srm',
                'nu'+obsidstr+fpm+'06_'+ge+'_p_sr.rmf',
                arf_file = 'nu'+obsidstr+fpm+'06_'+ge+'_p_sr.arf',
                data_file = 'nu'+obsidstr+fpm+'06_'+ge+'_p_sr.pha'
            )


def check_enough_10bins(hdul):

    csum=np.cumsum(hdul[1])
    thresh=10
    firstind=0
    indices=[]
    for i in range(0, len(csum)):
        if csum[i] > thresh:
            #print(csum[i])
            indices.append(np.s_[firstind:i])
            firstind=i+1
            thresh+=10
    
    finebins=0
    for ii in indices:
        #print(ii)
        if not np.any(np.array(hdul[0][ii]) < 2.5):
            finebins+=1
    
    #print(finebins)
    return (finebins > 1), finebins





def do_xspectral_fit(pha_path, expression='const*vapec',
                           model_name = 'isothermal_with_pileup',
                           fit_slope=True, slopes = {},
                           plot=False, dopileup=True):

    out_path = pha_path / 'xspec_out/'

    #Save the current path to return later.
    current_path = os.getcwd()
    print(current_path)

    #Go to the data directory.
    os.chdir(pha_path)
    print(os.getcwd())

    obsid = path_obsid(pha_path)



    #====================================================================
    #XSPEC time!
    #====================================================================
    
    xspec.Xset.abund = f'file {abundance_file}'
    interface = XSPECInterface()
    
    fpmA_pha = 'nu'+obsid+'A06_0_4_p_sr_grp.pha'
    fpmA_pileup = 'nu'+obsid+'A06_21_24_p_sr_grp.pha'
    fpmB_pha = 'nu'+obsid+'B06_0_4_p_sr_grp.pha'
    fpmB_pileup = 'nu'+obsid+'B06_21_24_p_sr_grp.pha'

    #Check if there is enough of a pile-up component to use in both FPM.
    hdul = nuutil.read_pha(fpmA_pileup, return_dat_hdr=False)
    FPMAcheck = check_enough_10bins(hdul)
    hdul = nuutil.read_pha(fpmB_pileup, return_dat_hdr=False)
    FPMBcheck = check_enough_10bins(hdul)
    
    if not np.all([FPMAcheck[0], FPMBcheck[0]]):
        print("INSUFFICIENT PILE-UP COUNTS TO DO PILE UP FITTING - NOT INCLUDING THOSE MODELS!")
        print("FPMA groups: ", FPMAcheck[1], '; FPMB groups: ', FPMBcheck[1])
        model_name = model_name+'_nm_re_the_pileup'
        dopileup=False
    
    if dopileup:
        interface.add_instrument(
            name = 'FPM A',
            signal_file = fpmA_pha,
            pileup_file = fpmA_pileup)
        interface.add_instrument(
            name = 'FPM B',
            signal_file = fpmB_pha,
            pileup_file = fpmB_pileup)
    else:
        interface.add_instrument(
            name = 'FPM A',
            signal_file = fpmA_pha)
        interface.add_instrument(
            name = 'FPM B',
            signal_file = fpmB_pha)        
    
    interface.read_data(pha_path)
    #xspec.AllData.ignore('**-2.5 10.0-**') # Specify which energy channels to ignore.
    xspec.AllData.ignore('bad') 
    xspec.AllData.ignore('**-2.5 15.0-**') 

    if dopileup:
        pileup_models = interface.set_pileup_model('expmodgauss')
        p_lim_file = DEFAULT_PILEUP_PARAMETER_FILE
    else:
        p_lim_file = DEFAULT_PARAMETER_FILE
        
    model = interface.add_component(
        model_name = model_name,
        expression = expression,
        parameter_limits_file = p_lim_file,
        out_dir = out_path)

    if expression=='const*(vapec+bknpower)':
        print('')
        print('LOOK HERE')
        comp2 = model.bknpower
        par1 = comp2.PhoIndx1
        par1.values = 2
        par1.frozen = True

        #value, fit delta, min, bot, top, max
        par2 = comp2.BreakE
        par2.values = [6, 0.1, 4, 4, 20, 20]

        par3 = comp2.PhoIndx2
        par3.values = [2, 0.1, 0, 0, 15, 15]

    if dopileup:
        # Modify the pileup norm to accurately represent the grade 0-4 pileup (1.25x).
        for instrument, pileup_model in zip(interface.instruments, pileup_models):
                pileup_component = pileup_model.componentNames[0]
                signal_model = interface.instruments[instrument].get_signal_model(
                    model_name)
                pileup_model_name = interface.instruments[instrument].pileup_model_name
                signal_model.__dict__[pileup_component].norm.link = '1.25 * ' \
                    f'{pileup_model_name}:p4'
                signal_model.__dict__[pileup_component].norm.link = '1.25 * ' \
                    f'{pileup_model_name}:p4'


    if not slopes:
        slopes = {
            'FPM A': (1, 0.01, 0.94, 0.94, 1.0, 1.0), # Follows XSPEC's convention: (initial, delta, hard min, soft min, soft max, hard max)
            'FPM B': (1, 0.01, 0.94, 0.94, 1.0, 1.0)
        }
    
    offsets = {'A': 0, 'B': 0}
    interface.set_gain(
        slopes = slopes, offsets = offsets,
        fit_slope = fit_slope, fit_offset = False)
    
    interface.fit(
        num_iterations = 1000,
        critical_delta = 0.01,
        fit_statistic = 'cstat',
        fit_error = True # XSPEC will determine the error on all unfrozen parameters
    )
    
    interface.archive_previous()
    interface.archive.save(out_path / f'{model_name}_archive.pkl')
    archive = interface.archive

    if plot:
        fig, axs = plt.subplots(
            2, 1, figsize=(5,7),
            sharex=True, layout='constrained',
            gridspec_kw={'height_ratios': [3, 1], 'hspace': 0})
        
        plotter = ModelPlotter(archive)
        ax0, ax1 = plotter.make_xspec_plot(model_name, axs=axs)
        ax1.axhline(0, color='gray', ls=':')
        ax0.legend()
        ax1.legend()
        plt.savefig(out_path / f'{model_name}.png')


    #print(archive.instruments['FPM A'].MODEL_NAME.response_parameters)
    #print(archive.instruments['FPM B'].MODEL_NAME.response_parameters)

    os.chdir(current_path)
    print(os.getcwd())

    return archive, model_name


