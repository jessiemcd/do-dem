import numpy as np
from astropy import units as u
import pickle
import pathlib

path_to_dodem = '/Users/jmdunca2/do-dem/'
from sys import path as sys_path
sys_path.append(path_to_dodem+'/dodem/')

import time_interval_selection as tis
import gauss2D as g2d





def goes_string(flux):

    #goes_conversion: [1e-8, 1e-7, 1e-6, 1e-5, 1e-4] = ['A', 'B', 'C', 'M', 'X']
    goes_levels = np.array([-12,-11,-10, -9, -8, -7, -6, -5, -4])
    goes_labels = np.array(['A0.0001', 'A0.001', 'A0.01','A0.1', 'A', 'B', 'C', 'M', 'X'])

    order_ = omag(flux.value)
    ind = np.where(goes_levels == order_)[0]
    if order_ < -8:
        round_ = -1*(order_+7)
        numberstr = str(round(flux.value/(10**(-8)), round_)) 
        gstring = 'A'+numberstr
    else:
        gstring = goes_labels[ind][0]+str(round(flux.value/(10**(order_)), 1)) 

    return gstring

def goes_minmax(vals):
    """
    Takes in a list of GOES XRSB values and returns the minimum and maximum values, 
    both as numbers and nicely formatted strings.
    
    """
    gmin, gmax = np.nanmin(vals), np.nanmax(vals)

    if gmin < 0:
        print('Negative?', gmin)
    if gmax < 0:
        print('Negative?', gmax)

    
    goesrange_vals = [gmin.value, gmax.value]*gmin.unit
    goesrange_strings = [goes_string(gmin),
                         goes_string(gmax)]
    
    print('GOES Class range: ', goesrange_strings[0], \
          ' to ', goesrange_strings[1])

    return goesrange_vals, goesrange_strings


def try_get_goes(datapath, satellite=13):
    """
    Get GOES data + return minimum and maximum GOES class during a NuSTAR orbit
    """

    #GOES16 data seems to start in february 2017 
    
    #Path to top-level do-dem directory - edit for your system.
    path_to_dodem = '/Users/jmdunca2/do-dem/'
    from sys import path as sys_path
    sys_path.append(path_to_dodem+'/dodem/')
    
    import initial_analysis as ia
    import nustar_utilities as nuutil
    import lightcurves as lc
    import astropy.time

    import datetime
    
    fpm='A'
    #Get first and last times associated with an event in the cleaned event file.
    evt_data, hdr = ia.return_submap(datapath=datapath, fpm=fpm, return_evt_hdr=True)
    time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]

    time0_ = time0.tt.datetime
    time1_ = time1.tt.datetime

    #print(time0, time1)

    timerange = [time0_,time1_]

    #from datetime import timezone
    #timerange = [t.replace(tzinfo=timezone.utc) for t in timerange]

    lc.get_goes(timerange, satellite=satellite, peek=False)
    instrument='GOES'
    data = lc.load_lightcurves(instrument)
    
    ylabel = data['GOES flux label']
    goestimes = data['GOES Times']
    xrsbcounts = data['XRSB counts']
    xrsblabel = data['XRSB label']
    gts = np.array([t.datetime for t in goestimes])

    int_inds = np.where(np.logical_and(gts > time0_, gts < time1_))[0]
    int_counts = xrsbcounts[int_inds]

    return goes_minmax(int_counts)


def all_obs_goes(datapaths, satellite=13):
    """
    Compile GOES info for all orbits in a given AR observation.
    """

    all_goes_vals=[]
    goes_per_orbit=[]
    goes_per_orbit_strings=[]
    for id in datapaths:
        goesrange_vals, goesrange_strings = try_get_goes(id, satellite=satellite)
        all_goes_vals.extend(goesrange_vals)
        goes_per_orbit.append(goesrange_vals)
        goes_per_orbit_strings.append(goesrange_strings)


    print('')
    print('all-day:')
    all_goes_vals_ = [vv.value for vv in all_goes_vals]*all_goes_vals[0].unit
    allminmax, allstrings = goes_minmax(all_goes_vals_)

    return allminmax, allstrings, goes_per_orbit, goes_per_orbit_strings

    
    

def omag(x):
    return int(np.floor(np.log10(np.abs(x))))


def get_durations(datapaths, fpm='A'):


    #Path to top-level do-dem directory - edit for your system.
    path_to_dodem = '/Users/jmdunca2/do-dem/'
    from sys import path as sys_path
    sys_path.append(path_to_dodem+'/dodem/')
    
    import initial_analysis as ia
    import nustar_utilities as nuutil

    import importlib
    from astropy import units as u
    import glob
    from astropy.io import fits

    import astropy.time
    

    importlib.reload(ia)

    durations = []
    lvttotals = []
    total = 0.*u.min
    for id in datapaths:
        #print(id)
        #Get first and last times associated with an event in the cleaned event file.
        evt_data, hdr = ia.return_submap(datapath=id, fpm=fpm, return_evt_hdr=True)
        time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
        #Duration: difference between them.
        durations.append((time1-time0).to(u.min))

        #Load in housekeeping data
        hk = glob.glob(id+'/hk/*'+fpm+'_fpm.hk')
        hdulist = fits.open(hk[0])
        dat = hdulist[1].data
        hdr = hdulist[1].header
        hdulist.close()

        #Get the livetimes from between the time0,1 found above.
        livetimes = dat['livetime']
        mjd_ref_time=astropy.time.Time(hdr['mjdrefi'],format='mjd')
        livetime_times = astropy.time.Time(mjd_ref_time+dat['time']*u.s,format='mjd')
        durbins = np.where(np.logical_and(livetime_times > time0, livetime_times < time1))[0]
        livetimes_ = np.array(livetimes)[durbins]
        
        #print(id[-12:-1], (np.sum(livetimes_)*u.s)/((time1-time0).to(u.s)))

        lvttotals.append((np.sum(livetimes_)*u.s).to(u.min))

    return durations, sum(durations), sum(lvttotals)


def print_times(files):

    #Path to top-level do-dem directory - edit for your system.
    path_to_dodem = '/Users/jmdunca2/do-dem/'
    from sys import path as sys_path
    sys_path.append(path_to_dodem+'/dodem/')
    
    import initial_analysis as ia
    import nustar_utilities as nuutil
    
    for d in files:
        evt_data, hdr = ia.return_submap(datapath=d, fpm='A', return_evt_hdr=True)
        time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
        print(time0, time1)



def load_nufiles(f):
    """
    """
    hdulist = fits.open(f)
    dat = hdulist[1].data
    hdr = hdulist[1].header
    hdulist.close()
        
    return dat, hdr




def single_gauss_prep(key, plot=True):


    with open('all_targets.pickle', 'rb') as f:
        data = pickle.load(f)
    
    ARDict = data[key]
    
    id_dirs = ARDict['datapaths']
    obsids = ARDict['obsids']
    working_dir = ARDict['working_dir']
    
    #Make a new working directory for prepped data/etc if it doesn't yet exist
    save_path = pathlib.Path(working_dir)
    if not save_path.exists():
        save_path.mkdir()

    gauss_stats=[]
    for i in range(0, len(id_dirs)):
        #guess, fast_min_factor 
        res = g2d.per_orbit_onegauss_params(id_dirs[i], guess=[], plot=plot)
        gauss_stats.append(res)


    ARDict['gauss_stats'] = gauss_stats

    data[key] = ARDict
    
    with open('all_targets.pickle', 'wb') as f:
             # Pickle the 'data' dictionary using the highest protocol available.
             pickle.dump(data, f, pickle.HIGHEST_PROTOCOL) 

    #where: where to find templates + place scripts.
    tis.make_tis_scripts(obsids, key, where='./')   




def do_key_dem(key, missing_last=False, missing_orbit=4):

    """
    Set missing_last=True to trim time interval list to exclude the last interval in an orbit (missing_orbit)
    (useful due to some NCCS AIA data intervals slightly shorter than NuSTAR intervals). 
    
    """

    import orbit_auto as oa

    with open('all_targets.pickle', 'rb') as f:
        data = pickle.load(f)

    ARDict = data[key]
    
    id_dirs = ARDict['datapaths']
    obsids = ARDict['obsids']
    working_dir = ARDict['working_dir']
    prepped_aia_dir = ARDict['prepped_aia']
    
    all_time_intervals, all_time_intervals_list = oa.find_all_intervals(working_dir, shush=True, 
                                                                        missing_last=missing_last, missing_orbit=missing_orbit)

    #What instruments are you using?
    #---------------------------------
    aia=True
    #---------------------------------
    eis=False
    xrt=True
    #This is where I'm putting my XRT level-1 data and grade maps:
    xrt_path=working_dir+'/XRT_for_DEM/'
    xrt=True
    plot_xrt=False
    from astropy import units as u
    exposure_dict={'Be_thin': [1*u.s, 10*u.s],
                    'Be_med': [],
                  'Al_poly': [0.1*u.s, 1*u.s]}
    #---------------------------------
    plot=False
    #---------------------------------
    nustar=True
    #If nustar is being used, here are the chosen energy ranges:
    nuenergies=[[2.5,3.5], [3.5,6.], [6.,10.]]
    nuradius=150
    onegauss=True
    #---------------------------------
    
    #---------------------------------
    #---------------------------------
    #What temperature range would you like to use? (units: log(T))
    minT=5.6
    maxT=7.2
    
    #Would you prefer to plot temperatures in MK, or the default (logT)
    plotMK=False
    #---------------------------------
    #---------------------------------
    
    name=key
    
    #Path to top-level do-dem directory - edit for your system.
    path_to_dodem = '/Users/jmdunca2/do-dem/'
    
    import dodem
    
    for o in range(0, len(obsids)):
    
        datapath=id_dirs[o]
        gtifile=datapath+'event_cl/nu'+obsids[o]+'A06_gti.fits'
        regfile=path_to_dodem+'starter_region.reg'
    
        orbit_aia_dir = prepped_aia_dir+'/orbit_'+obsids[o]+'/'
    
        time_intervals = all_time_intervals[o]
    
        guess = ARDict['gauss_stats'][o][0]
    
        xrt_path=path_to_dodem+'other_idl/'+obsids[o]+'_coobs/XRT_for_DEM/'

        if not pathlib.Path(xrt_path).is_dir():
            xrt=False
    
        for time in time_intervals:
    
            data, bl, tr, region_input = oa.read_interval_dicts(time, where=orbit_aia_dir, bltr=True)

            dodem.dodem(time, bl, tr, xrt=xrt, aia=aia, nustar=nustar, name=name,
                                        plotMK=plotMK, minT=minT, maxT=maxT,
                                        plotresp=False, working_directory=working_dir,
                                        default_err=0.2, path_to_dodem=path_to_dodem,
                
                                        #demreg related
                                        rgt_fact=1.2, max_iter=30,
                                        reg_tweak=1, gloci=1, mc_in=True, mc_rounds=100, 
                                        
                                        #nustar related 
                                        combine_fpm=True, nuenergies=nuenergies, 
                                        datapath=datapath, gtifile=gtifile,
                                        nuradius=nuradius, guess=guess, onegauss=onegauss,
                                        adjacent_grades=True, pile_up_corr=True,
                
                                        #aia related
                                        load_prepped_aia=data, 
    
                                        #xrt related
                                       xrtmethod='Average', real_xrt_err=True, xrt_path=xrt_path,
                                        xrt_exposure_dict=exposure_dict,
                                        input_xrt_region="circle", input_xrt_region_dict=region_input)
        print('')

















    
    