import numpy as np
from astropy import units as u
import pickle
import pathlib

path_to_dodem = '/Users/jmdunca2/do-dem/'
from sys import path as sys_path
sys_path.append(path_to_dodem+'/dodem/')

import time_interval_selection as tis
import visualize_dem_results as viz
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


def get_exposures(target_dict, dogoes=False):

    keys = target_dict.keys()
    
    all_effective_exposure=0*u.min
    all_duration=0*u.min
    all_all_goes=[]
    for key in keys:
        ARdict = target_dict[key]
        durations, total, lvttotal = get_durations(ARdict['datapaths'])
        print(key)
        print('Duration: ', np.round(total, 2), '| Effective Exposure: ', np.round(lvttotal,2), '| Livetime %: ', np.round((lvttotal/total),2))
    
        target_dict[key]['orbit durations'] = durations
        target_dict[key]['total duration'] = total
        target_dict[key]['total livetime'] = lvttotal
    
        all_effective_exposure += lvttotal
        all_duration += total
    
        if dogoes:
            allminmax, allstrings, goes_per_orbit, goes_per_orbit_strings = ana.all_obs_goes(ARdict['datapaths'], satellite=ARdict['goes_satellite'])
            target_dict[key]['AR GOES min, max vals'] = allminmax
            target_dict[key]['AR GOES min, max strings'] = allstrings
            target_dict[key]['orbit GOES min, max vals'] = goes_per_orbit
            target_dict[key]['orbit GOES min, max strings'] = goes_per_orbit_strings
            print('')
            all_all_goes.append(allminmax)
    
    
    print('')
    print('Total Dataset Effective Exposure: ', all_effective_exposure.to(u.h))
    print('Total Dataset Observation Duration: ', all_duration.to(u.h))
    
    if dogoes:
        print('GOES - all-dataset: ')
        all_all_goes_vals_ = [vv.value for vv in all_all_goes]*all_all_goes[0].unit
        allminmax, allstrings = goes_minmax(all_all_goes_vals_)




def single_gauss_prep(key, plot=True, guess=[]):


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
        res = g2d.per_orbit_onegauss_params(id_dirs[i], guess=guess, plot=plot)
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
    plot_xrt=True
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
                                        xrt_exposure_dict=exposure_dict, plot_xrt=plot_xrt,
                                        input_xrt_region="circle", input_xrt_region_dict=region_input)
        print('')




def get_above10s(key='', all=True, plot=False, time_weighted=False, seconds_per=5):
    
    import glob
    
    import pandas as pd
    import astropy.time
    
    df = pd.read_csv('fpmA.csv')
    starts = df['flare_start'].values
    stops = df['flare_end'].values
    
    from astropy import units as u
    early_starts = [(astropy.time.Time(s)-2*u.min) for s in starts]
    late_stops = [(astropy.time.Time(s)+2*u.min) for s in stops]

    if all:
        res_files = glob.glob('./compact_results/*')
    elif key:
        res_files = glob.glob('./compact_results/*'+key+'*.pickle')
    else:
        print('Either set all=True, or select a key!')
    res_files.sort()
    
    all_above10s_flares = []
    all_above10s_non = []
    all_above10s=[]
    times = []
    durations =  []
    
    for f in res_files:
        flare=False
    
        data, timestring, time = viz.load_DEM(f)
        times.append(time)
        dur = (time[1]-time[0]).to(u.s)
        durations.append(dur)
        time_mult = round(dur.value/seconds_per)
    
        b4 = [s < time[0] for s in early_starts]
        ea = [s > time[0] for s in late_stops]
        es = np.where(np.logical_and(b4, ea))
    
        if es[0].size > 0:
            flare=True
    
        b4 = [s > time[0] for s in early_starts]
        ea = [s < time[1] for s in early_starts]
        es = np.where(np.logical_and(b4, ea))
    
        if es[0].size > 0:
            flare=True
        
        res = viz.get_DEM_params(f)
    
        m1, max1, above5_, above7_, above10_, \
            above_peak, below_peak, above_635, below_635, \
               data['chanax'], data['dn_in'], data['edn_in'], \
                    powerlaws, EMT_all, EMT_thresh = res


    
        if flare:
            if time_weighted:
                #print(time_mult)
                #print([above10_[0] for i in range(0,time_mult)])
                all_above10s_flares.extend([above10_[0] for i in range(0,time_mult)])
            else:
                all_above10s_flares.append(above10_[0])
        else:
            if time_weighted:
                #print(time_mult)
                #print([above10_[0] for i in range(0,time_mult)])
                all_above10s_non.extend([above10_[0] for i in range(0,time_mult)])
            else:  
                all_above10s_non.append(above10_[0])
        
        if time_weighted:
            all_above10s.extend([above10_[0] for i in range(0,time_mult)])
        else:
            all_above10s.append(above10_[0])

    # #print(durations)
    # for d in durations:
    #     print(round(d.value/5))

    if plot:

        area_i = 100**2
        area_m = np.pi*150**2
        #print(area_i, area_m)
        factor = area_m/area_i
        factor = 1
        
        from matplotlib import pyplot as plt 
        
        logbins = np.geomspace(np.min(all_above10s), np.max(all_above10s), 50)
        
        fig, ax = plt.subplots(figsize=(15,4), tight_layout = {'pad': 1})
        
        ax.hist(all_above10s, bins=logbins, color='green', edgecolor='black', label='All bins')
        ax.hist(np.array(all_above10s_non)*factor, bins=logbins, color='skyblue', edgecolor='black', label='Non-flare time bins')
        ax.hist(np.array(all_above10s_flares)*factor, bins=logbins, color='purple', edgecolor='black', label="Bins during Reed's flares")
        ax.set_xscale('log')
        ax.axvline(1.8e22, color='Red')
        ax.axvline(1.5e23, color='Red')
        ax.axvspan(1.8e22, 1.5e23, alpha=0.3, color='Red', label='Ishikawa (2017) 95% Interval')
        ax.set_ylabel('Number of intervals')
        ax.set_xlabel('EM Integrated >10 MK')
        if key:
            ax.set_title(key)
        else:
            ax.set_title('All')
        ax.legend()
        if time_weighted:
            plt.savefig(key+'_'+str(seconds_per)+'s_bin_time_weighted_above10splot.png')
        else:
            plt.savefig(key+'_above10splot.png')
        

    return all_above10s, all_above10s_flares, all_above10s_non












    
    