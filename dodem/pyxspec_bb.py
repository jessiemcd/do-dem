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
path_to_spec = '/Users/jmdunca2/do-dem/spectroscopy/'
from sys import path as sys_path
sys_path.append(path_to_spec+'/pyxspec_extension/')


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

    print(model.vapec)

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
        par3.values = [2, 0.1, 1, 1, 20, 20]

    #Default is 30 MK; if you want to change it change it here.
    upper_lim_mk=30

    ulkt=upper_lim_mk/11.6
    normlo = 1e41*3.5557e-42
    normhi = 1e49*3.5557e-42

    if expression=='const*vapec' or expression=='const*(vapec+vapec)':
        par1 = model.vapec.kT
        par2 = model.vapec.norm
        par1.values = par1.values[0:4]+[ulkt, ulkt]
        par2.values = par2.values[0:2]+[normlo, normlo, normhi, normhi]
    
    if expression=='const*(vapec+vapec)':
        if dopileup:
            par1 = model.vapec_4.kT
            par2 = model.vapec_4.norm
        else:
            par1 = model.vapec_3.kT
            par2 = model.vapec_3.norm
        par1.values = par1.values[0:4]+[ulkt, ulkt]
        par2.values = par2.values[0:2]+[normlo, normlo, normhi, normhi]

    

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



def convert_component_params(component_parameters, labstring):

    """
    Takes XSPEC output parameters and makes a dictionary of parameters (label + value),
    coordinate conversions are completed for the thermal parameters. 

    Inputs
    ---------
    labstring            - label to be added to the parameter labels (e.g. indicating the name of the model/fit).
    component_parameters - free_parameters attribute of a specific instrument/model in an archive object from 
                            pyxspec_extension, as in:
                                    component_parameters = archive.instruments[{instrument}][{model_name}].free_parameters


    Outputs
    -------
    converted_dict - dictionary of parameters, with converted thermal parameter values
    
    """
    from astropy import units as u
    
    
    PARAMETER_CONVERSIONS = dict(
    nlapec=dict(
        kT=11.6045 * u.MK / u.keV,
        norm=u.cm**-3 / 3.5557e-42
    ),
    vapec=dict(
        kT=11.6045 * u.MK / u.keV,
        norm=u.cm**-3 / 3.5557e-42
    )
    )

    converted_dict = {}


    for component, parameters in component_parameters.items():
        for parameter in parameters:
            value = parameter.quantity
            errors = np.array([*parameter.error[0:2]]) * parameter.unit
            errors -= value
    
            if component.base_name in PARAMETER_CONVERSIONS:
                conversions = PARAMETER_CONVERSIONS[component.base_name]
                if parameter.name in conversions:
                    value *= conversions[parameter.name]
                    errors *= conversions[parameter.name]


            #print(parameter.name)
            #if parameter.name == 'kT':
            #    if value > 15*u.MK:
            #        print('High temp!: ', value, errors)
    

            if parameter.name+labstring in converted_dict:
                converted_dict[parameter.name+'_2'+labstring] = [value, errors]
            else:
                converted_dict[parameter.name+labstring] = [value, errors]

    return converted_dict


def print_model_params(out_path, model_name, modelstring, pileup=True):

    """
    For a given path containing a saved pyxspec_extension archive for a given model

    Inputs
    -------
    out_path - path to archive (.pkl)
    model_name - string name of fit model (found in archive .pkl filename, and used internally 
                    in the archive structure). Must be the model name from the original fit.
    modelstring - shorter identifying model string for organizing files/parameters. Can be 
                    made up (not pre-existing in archive structure). 
    pileup - set True if the model in question included a pileup component.


    Outputs
    -------
    all_model_dict - dictionary containing parameter values for both instruments (FPMA, B),
                        pile-up (if relevant) and regular.
                    
    
    """

    interface = XSPECInterface()

    try:
        archive = interface.archive.load(out_path / f'{model_name}_archive.pkl')
    except FileNotFoundError:
        #print('No file found for ', model_name, ' at ', out_path)
        return

    all_model_dict = {}

    component_parameters = archive.instruments['FPM A'][model_name].free_parameters
    cd = convert_component_params(component_parameters, '_m'+modelstring)
    all_model_dict.update(cd)
    #print(all_model_dict.keys())

    if pileup:
        component_parameters = archive.instruments['FPM A']['pileupFPMA_'+model_name].free_parameters
        cd = convert_component_params(component_parameters, '_puA'+modelstring)
        all_model_dict.update(cd)
        #print(all_model_dict.keys())

    #print('')
    
    component_parameters = archive.instruments['FPM B'][model_name].free_parameters
    cd = convert_component_params(component_parameters, '_m'+modelstring)
    all_model_dict.update(cd)
    #print(all_model_dict.keys())

    if pileup:
        #print('')
        component_parameters = archive.instruments['FPM B']['pileupFPMB_'+model_name].free_parameters
        cd = convert_component_params(component_parameters, '_puB'+modelstring)
        all_model_dict.update(cd)
        #print(all_model_dict.keys())



    all_model_dict.update({'CSTAT_FPMA'+modelstring: archive.instruments['FPM A'][model_name].statistic,
                           'CSTAT_FPMB'+modelstring: archive.instruments['FPM B'][model_name].statistic
                          })

    #print('')
    #print('')

    return all_model_dict
    


def check_rez_pyxspec(filelist, shush=False):
    
    all_finishes=[]
    nonefiles=[]
    incompletes=[]
    fullsets=[]
    nopileup=[]
    for f in filelist:
        if not shush:
            print(f)
        finishes=np.zeros(3)
        npfinishes=np.zeros(3)
        thepath = pathlib.Path('/'.join(f.split('/')[0:-1])+'/')
        timestring = path_timestring(thepath)
        isotfile = str(thepath)+'/xspec_out/'+'isothermal_with_pileup_nogain.png'
        isotfile2 = str(thepath)+'/xspec_out/'+'isothermal_with_pileup_nogain_nm_re_the_pileup.png'
        if os.path.exists(isotfile):
            if not shush:
                print('yes isot finish')
            finishes[0]=1
        elif os.path.exists(isotfile2):
            if not shush:
                print('yes isot finish')
            npfinishes[0]=1
            
        #twotfile = str(thepath)+'/xspec_out/'+'two_thermal_with_pileup.png'
        #twotfile1 = str(thepath)+'/xspec_out/'+'two_thermal_with_pileup_nm_re_the_pileup.png'
        twotfile2 = str(thepath)+'/xspec_out/'+'two_thermal_with_pileup_nogain.png'
        twotfile3 = str(thepath)+'/xspec_out/'+'two_thermal_with_pileup_nogain_nm_re_the_pileup.png'
        # if os.path.exists(twotfile):
        #     print('yes twot finish')
        #     finishes[1]=1
            
        # elif os.path.exists(twotfile1):
        #         print('yes twot finish')
        #         finishes[1]=1
        # elif os.path.exists(twotfile2):
        #         print('yes twot finish')
        #         finishes[1]=1
        # elif os.path.exists(twotfile3):
        #         print('yes twot finish')
        #         finishes[1]=1
        if os.path.exists(twotfile2):
            if not shush:
                print('yes twot finish')
            finishes[1]=1
        elif os.path.exists(twotfile3):
            if not shush:
                print('yes twot finish')
            npfinishes[1]=1
            
        #tntfile = str(thepath)+'/xspec_out/'+'thermal_nonthermal_with_pileup.png'
        #tntfile1 = str(thepath)+'/xspec_out/'+'thermal_nonthermal_with_pileup_nm_re_the_pileup.png'
        tntfile2 = str(thepath)+'/xspec_out/'+'thermal_nonthermal_with_pileup_nogain.png'
        tntfile3 = str(thepath)+'/xspec_out/'+'thermal_nonthermal_with_pileup_nogain_nm_re_the_pileup.png'
        # if os.path.exists(tntfile):
        #     print('yes tnt finish')  
        #     finishes[2]=1
            
        # elif os.path.exists(tntfile1):
        #         print('yes tnt finish')  
        #         finishes[2]=1
        # elif os.path.exists(tntfile2):
        #         print('yes tnt finish')  
        #         finishes[2]=1

        # elif os.path.exists(tntfile3):
        #         print('yes tnt finish')  
        #         finishes[2]=1

        if os.path.exists(tntfile2):
            if not shush:
                print('yes tnt finish')  
            finishes[2]=1
        elif os.path.exists(tntfile3):
            if not shush:
                print('yes tnt finish')  
            npfinishes[2]=1
        
        if np.sum(finishes)==0:
            nonefiles.append(f)
            
        if np.sum(finishes) < 3 and np.sum(finishes) > 0:
            incompletes.append(f)
            
        if np.sum(finishes) == 3:
            fullsets.append(f)

        if np.sum(npfinishes) == 3:
            nopileup.append(f)

        all_finishes.append(finishes)
        if not shush:
            print('')

    return np.array(all_finishes), nonefiles, incompletes, fullsets, nopileup




def save_SPEX_params(filelist_, extrastring, pileup=True):

    """
    Takes a list of files, and an extra string which may be used to label them. Acquires dictionaries of fit parameters
    for three spectral models, and saves them in a dictionary which is then added to the result file dictionary as
    'SPEX_dict'.
    
    Inputs
    -------
    filelist_ - list of result files in expected locations (same directory as an xspec results dictionary
    extrastring - additional model label (e.g. "no pileup") 
    pileup - set True if pileup is included in the models desired. 
    
    
    Outputs
    -------
    None - edits result file.
    
    """

    import pickle
    
    #empty lists for isothermal
    kt_its = []
    em_its = []
    
    #empty lists for two-thermal
    kt1_tts = []
    em1_tts = []
    kt2_tts = []
    em2_tts = []
    
    for i in range(0, len(filelist_)):
    #for i in range(0,1):
        f_=filelist_[i]
        #print(f_)
        #result file xspec directory path (location of archive files)
        thepath = pathlib.Path('/'.join(f_.split('/')[0:-1])+'/xspec_out/')
    
        #empty dictionary for all fit parameters (all models)
        spec_dict = {}
    
        #adding isothermal model parameters from archive
        model_name = 'isothermal_with_pileup'+extrastring
        amd = print_model_params(thepath, model_name, '_it', pileup=pileup)
        if amd:
            spec_dict.update(amd)
    
        #adding double thermal model parameters from archive
        model_name = 'two_thermal_with_pileup'+extrastring
        amd = print_model_params(thepath, model_name, '_tt', pileup=pileup)
        #print('tt amd: ', amd)
        if amd:
            spec_dict.update(amd)
    
        #adding thermal + nonthermal parameters from archive
        model_name = 'thermal_nonthermal_with_pileup'+extrastring
        amd = print_model_params(thepath, model_name, '_nt', pileup=pileup)
        if amd:
            spec_dict.update(amd)
    
        #print(len(spec_dict.keys()))
        #print(spec_dict.keys())
    
        try:
            kt_its.append(spec_dict['kT_m_it'][0])
            em_its.append(spec_dict['norm_m_it'][0])
            ('It worked: ', f_)
        except KeyError:
            pass
            #print('No archive found for any models for this result with this string in title:', extrastring)
            #print('')
    
    
        try:
            if spec_dict['kT_m_tt'][0] > spec_dict['kT_2_m_tt'][0]:
                kt1_tts.append(spec_dict['kT_m_tt'][0])
                em1_tts.append(spec_dict['norm_m_tt'][0])
                kt2_tts.append(spec_dict['kT_2_m_tt'][0])
                em2_tts.append(spec_dict['norm_2_m_tt'][0])
            else:
                kt1_tts.append(spec_dict['kT_2_m_tt'][0])
                em1_tts.append(spec_dict['norm_2_m_tt'][0])
                kt2_tts.append(spec_dict['kT_m_tt'][0])
                em2_tts.append(spec_dict['norm_m_tt'][0])
                
        except KeyError:
            pass
    
        if spec_dict:
    
            with open(f_, 'rb') as f:
                data = pickle.load(f)

            if 'SPEX_dict' in data.keys():
                thespexdict = data['SPEX_dict']
                thespexdict.update(spec_dict)
                data.update({'SPEX_dict': thespexdict})
            else:
                data['SPEX_dict'] = spec_dict
        
            with open(f_, 'wb') as f:
                pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)  
    

def prep_for_spec(files, out_dir_name='xspec_out_/'):

    #Making SRM files for all of the quiet region/times. Also, making a little xspec results dictionary. 
    import subprocess
    
    for f in files:
        the_path = pathlib.Path('/'.join(f.split('/')[0:-1])+'/')
    
        #Make an xspec out directory in the data directory if one does not already exist.
        if not the_path.exists():
            print('Path wrong!!!!!!!', f)
            continue
    
        #Make an xspec out directory in the data directory if one does not already exist.
        out_path = the_path / out_dir_name
        if not out_path.exists():
            out_path.mkdir()
    
        obsidstr = path_obsid(the_path)
    
    
        #Save the current path to return later.
        current_path = os.getcwd()
        print(current_path)
    
        #Go to the data directory.
        os.chdir(the_path)
        print(os.getcwd())


        #Make grouped PHA files (only if they don't already exist
        fpms=['A', 'B']
        grade_exps = ['0_4', '21_24']
        for fpm in fpms:
            for ge in grade_exps:
                pha_file = 'nu'+obsidstr+fpm+'06_'+ge+'_p_sr.pha'
                grp_file = 'nu'+obsidstr+fpm+'06_'+ge+'_p_sr_grp.pha'
                if not os.path.exists(str(the_path)+'/'+grp_file):
                    logname = fpm+'_'+ge+'_grppha.log'
                    print(pha_file, grp_file, logname)
    
                    status = subprocess.call('grppha infile='+pha_file+' outfile='+grp_file+' chatter=1 comm="group min 10 & exit" >& '+logname, 
                                 shell=True) 
    
        #Make SRM response files
        all_srms(obsidstr)
        
        
        os.chdir(current_path)
        print(os.getcwd())
        
    
        print('')


def make_pyxspec_scripts(filelist, plot=True):

    """
    This is the old version where have one big python script for each case (that runs all three spectral fits).

    See make_pyxspec_subprocess_scripts() for the version that makes a script for each spectral fit individually 
    and sets a time limit on how long they are allowed to run before we quit them.
    
    """

    
    templatefile = './template_pyxspec_run.py'
    pystrings = []

    for i in range(0, len(filelist)):
        f=filelist[i]
        thepath = pathlib.Path('/'.join(f.split('/')[0:-1])+'/')
        timestring = path_timestring(thepath)
        pyfile = str(thepath)+'/xspec_out/run_pyxspec.py'
        pystring = 'python '+pyfile+' > '+' '+str(thepath)+'/xspec_out/pyxspec_out.txt &'
        #print(pystring)
        pystrings.append(pystring)
        if i % 5 == 0 and i > 0:
            pystrings.extend(['', 'wait', '', 'echo "finished up to '+str(i)+'"', ''])
    
        with open(templatefile, 'r') as f:
            lines = f.read()
            llist = lines.split('\n')
            llist[8] = 'pathstring = "'+str(thepath)+'"'
            if plot:
                llist[9] = 'plot = True'
            else:
                llist[9] = 'plot = False'
            newlist = '\n'.join(llist)
            #print(newlist)
        
            with open(pyfile, 'w') as file_out:
                file_out.seek(0)
                file_out.write(newlist)
                file_out.truncate()
    

    with open('run_quiet_spectra_template.sh', 'r') as f:
        lines = f.read()
        llist = lines.split('\n')
        llist.extend(['','',])
        llist.extend(pystrings)
        llist.extend(['wait', '', 'echo "all scripts finished"'])
        pylist = '\n'.join(llist)
        print(pylist)

        with open('run_quiet_spectra.sh', 'w') as file_out:
            file_out.seek(0)
            file_out.write(pylist)
            file_out.truncate()




def make_pyxspec_subprocess_scripts(filelist, xspec_dir_str='xspec_out',
                                    plot=True):
    """
    This one makes a specific script for each spectral fit to be run as a subprocess so that we 
    can set a time limit.
    
    """
    
    templatefile = './template_pyxspec_subrun.py'
    templateisot = './template_pyxspec_isot.py'
    templatetwot = './template_pyxspec_twot.py'
    templatetnt = './template_pyxspec_tnt.py'
    pystrings = []

    for i in range(0, len(filelist)):
        f=filelist[i]
        thepath = pathlib.Path('/'.join(f.split('/')[0:-1])+'/')
        timestring = path_timestring(thepath)
        pyfile = str(thepath)+'/'+xspec_dir_str+'/subrun_pyxspec.py'
        pystring = 'python '+pyfile+' > '+' '+str(thepath)+'/'+xspec_dir_str+'/pyxspec_out.txt &'
        #print(pystring)
        pystrings.append(pystring)
        if i % 5 == 0 and i > 0:
            pystrings.extend(['', 'wait', '', 'echo "finished up to '+str(i)+'"', ''])
    
        with open(templatefile, 'r') as f:
            lines = f.read()
            llist = lines.split('\n')
            llist[3] = 'pathstring = "'+str(thepath)+'/'+xspec_dir_str+'"'
            newlist = '\n'.join(llist)
            #print(newlist)
        
            with open(pyfile, 'w') as file_out:
                file_out.seek(0)
                file_out.write(newlist)
                file_out.truncate()

        pyfile = str(thepath)+'/'+xspec_dir_str+'/isot_pyxspec.py'

        with open(templateisot, 'r') as f:
            lines = f.read()
            llist = lines.split('\n')
            llist[8] = 'pathstring = "'+str(thepath)+'"'
            if plot:
                llist[9] = 'plot = True'
            else:
                llist[9] = 'plot = False'
            newlist = '\n'.join(llist)
            #print(newlist)
        
            with open(pyfile, 'w') as file_out:
                file_out.seek(0)
                file_out.write(newlist)
                file_out.truncate()

        pyfile = str(thepath)+'/'+xspec_dir_str+'/twot_pyxspec.py'

        with open(templatetwot, 'r') as f:
            lines = f.read()
            llist = lines.split('\n')
            llist[8] = 'pathstring = "'+str(thepath)+'"'
            if plot:
                llist[9] = 'plot = True'
            else:
                llist[9] = 'plot = False'
            newlist = '\n'.join(llist)
            #print(newlist)
        
            with open(pyfile, 'w') as file_out:
                file_out.seek(0)
                file_out.write(newlist)
                file_out.truncate()

        pyfile = str(thepath)+'/'+xspec_dir_str+'/tnt_pyxspec.py'

        with open(templatetnt, 'r') as f:
            lines = f.read()
            llist = lines.split('\n')
            llist[8] = 'pathstring = "'+str(thepath)+'"'
            if plot:
                llist[9] = 'plot = True'
            else:
                llist[9] = 'plot = False'
            newlist = '\n'.join(llist)
            #print(newlist)
        
            with open(pyfile, 'w') as file_out:
                file_out.seek(0)
                file_out.write(newlist)
                file_out.truncate()
    
    
    with open('run_quiet_spectra_template.sh', 'r') as f:
        lines = f.read()
        llist = lines.split('\n')
        llist.extend(['','',])
        llist.extend(pystrings)
        llist.extend(['', 'wait', '', 'echo "all scripts finished"'])
        pylist = '\n'.join(llist)
        print(pylist)

        with open('run_quiet_spectra.sh', 'w') as file_out:
            file_out.seek(0)
            file_out.write(pylist)
            file_out.truncate()







