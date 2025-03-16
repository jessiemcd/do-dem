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
import images_and_coalignment as iac
import nustar_dem_prep as nu




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

    import nustar_utilities as nuutil
    import lightcurves as lc
    import astropy.time

    import datetime
    
    fpm='A'
    #Get first and last times associated with an event in the cleaned event file.
    evt_data, hdr = nu.return_submap(datapath=datapath, fpm=fpm, return_evt_hdr=True)
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
    
    import nustar_utilities as nuutil

    from astropy import units as u
    import glob
    from astropy.io import fits

    import astropy.time
    


    durations = []
    lvttotals = []
    total = 0.*u.min
    for id in datapaths:
        #print(id)
        #Get first and last times associated with an event in the cleaned event file.
        evt_data, hdr = nu.return_submap(datapath=id, fpm=fpm, return_evt_hdr=True)
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
    
    import nustar_utilities as nuutil
    
    for d in files:
        evt_data, hdr = nu.return_submap(datapath=d, fpm='A', return_evt_hdr=True)
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
            allminmax, allstrings, \
                goes_per_orbit, goes_per_orbit_strings = all_obs_goes(ARdict['datapaths'], satellite=ARdict['goes_satellite'])
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




def single_gauss_prep(key, plot=True, guess=[], make_scripts=True,
                     plotregion=[], plotgaussregions=False):


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
        res = g2d.per_orbit_onegauss_params(id_dirs[i], guess=guess, plot=plot,
                                           plotregion=plotregion, plotgaussregions=plotgaussregions)
        gauss_stats.append(res)


    ARDict['gauss_stats'] = gauss_stats

    data[key] = ARDict
    
    with open('all_targets.pickle', 'wb') as f:
             # Pickle the 'data' dictionary using the highest protocol available.
             pickle.dump(data, f, pickle.HIGHEST_PROTOCOL) 

    if make_scripts:
        #where: where to find templates + place scripts.
        tis.make_tis_scripts(obsids, key, where='./scripts/')  



def double_gauss_prep(key, plot=True, guess=[], guess2=[], sep_axis='SN', make_scripts=True,
                      plotregion=[], write_input_regions=False,
                      plotgaussregions=False, write_regions=False):


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
        res = g2d.per_orbit_twogauss_params(id_dirs[i], sep_axis=sep_axis, guess=guess, guess2=guess2, plot=plot,
                                           plotregion=plotregion, plotgaussregions=plotgaussregions,
                                            write_input_regions=write_input_regions,
                                            write_regions=write_regions, region_dir=working_dir)
                        
        gauss_stats.append(res)
        print('')

    print(gauss_stats)


    ARDict['gauss_stats'] = gauss_stats

    data[key] = ARDict
    
    with open('all_targets.pickle', 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL) 

    if make_scripts:
        ##where: where to find templates + place scripts.
        tis.make_tis_scripts(obsids, key, where='./scripts/', tworegion=True)   
    



def manual_prep(key, plot=True, guess=[], guess2=[], make_scripts=True,
                      plotregion=[], write_input_regions=True,
                      plotgaussregions=False):

    """
    key - key for all target dictionary (where to get information about the nustar data, region, etc)
    plot - set True for plots to be made in general
    guess, guess2 - for tweaking the double gaussian fit used for context. Not related to final regions 
                    saved, etc. Optional. 

    make_scripts - set True with you've finalized your regions and are ready to write corresponding scripts for TIS
    plotregion - list of region dictionaries, of the form:
                    
                    plotregion = [{'centerx': 950, 'centery': -325, 'radius': 150},
                               {'centerx': 900, 'centery': -50, 'radius': 150}]

                (You can have as many regions as you want. Values in arcseconds from solar center.)
                
    write_input_regions - set True to write a .reg file for every region in plotregion. This is needed to run TIS 
                            (scripts written in make_scripts will call functions that will look for these regions).

    plotgaussregions - set True to plot 150" circles centered at gaussian fit result locations, if useful for visualization.

    
    """


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

    region_stats=[]
    for i in range(0, len(id_dirs)):
        #guess, fast_min_factor 
        res = g2d.per_orbit_manual_params(id_dirs[i], guess=guess, guess2=guess2, plot=plot,
                                           plotregion=plotregion, plotgaussregions=plotgaussregions,
                                            write_input_regions=write_input_regions,
                                            region_dir=working_dir)
                        
        region_stats.append(res)
        print('')


    ARDict['region_stats'] = region_stats

    data[key] = ARDict
    
    with open('all_targets.pickle', 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL) 

    if make_scripts:
        ##where: where to find templates + place scripts.
        tis.make_tis_scripts(obsids, key, where='./scripts/', manualregion=True)  
    



def do_key_dem(key, missing_last=False, missing_orbit=4, plot_xrt=True, method='fit',
              high_temp_analysis=False):

    """
    Set missing_last=True to trim time interval list to exclude the last interval in an orbit (missing_orbit)
    (useful due to some NCCS AIA data intervals slightly shorter than NuSTAR intervals). 
    
    """


    #Path to top-level do-dem directory - edit for your system.
    path_to_dodem = '/Users/jmdunca2/do-dem/'

    with open('all_targets.pickle', 'rb') as f:
        data = pickle.load(f)

    ARDict = data[key]
    
    id_dirs = ARDict['datapaths']
    obsids = ARDict['obsids']
    working_dir = ARDict['working_dir']
    prepped_aia_dir = ARDict['prepped_aia']
    if method=='double':
        gauss_stats = ARDict['gauss_stats']
        sep_axis = gauss_stats[0][0]
    else:
        sep_axis = ''
        

    if method in ['input', 'double']:
        directories = get_region_directories(key, method=method)
        all_all_time_intervals, fixit = tis.region_time_intervals(directories, id_dirs, shush=True)

    if method=='fit':
        onegauss=True
        regfile=path_to_dodem+'starter_region.reg'
        all_time_intervals, all_time_intervals_list = tis.find_all_intervals(working_dir, shush=True, 
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
    
    import dodem
    import glob
    
    for o in range(0, len(obsids)):
    
        datapath=id_dirs[o]
        xrt_path=path_to_dodem+'other_idl/'+obsids[o]+'_coobs/XRT_for_DEM/'
        if not pathlib.Path(xrt_path).is_dir():
            xrt=False
        gtifile=datapath+'event_cl/nu'+obsids[o]+'A06_gti.fits'
        orbit_aia_dir = prepped_aia_dir+'/orbit_'+obsids[o]+'/'
        obsid=obsids[o]

        if method=='fit':
            guess = ARDict['gauss_stats'][o][0]
            time_intervals = all_time_intervals[o]    
            for time in time_intervals:

                print(orbit_aia_dir)
                data, bl, tr, region_input = iac.read_interval_dicts(time, where=orbit_aia_dir, bltr=True)
                #print('data: ', data)
                #print(data['aia_dn_s_px'])

                if high_temp_analysis:
                    dodem.high_temp_analysis(time, bl, tr, xrt=xrt, aia=aia, nustar=nustar, name2=name,
                                                   plotMK=plotMK, highT=7.2, #(minT, maxT are varied)
                                                   working_directory=working_dir, #(plotresp set false in high_temp_analysis)
                                                   default_err=0.2, path_to_dodem=path_to_dodem,
                                                   demmethod='DEMREG',
                            
                                                   #demreg/xrt_iterative related
                                                   rgt_fact=1.2, max_iter=30, 
                                                   reg_tweak=1, #mc_in, mc_rounds hard coded in high_temp_analysis (same as below)
                            
                                                   #nustar=related
                                                   nuenergies = nuenergies, #combine_fpm hard coded in high_temp_analysis (same as below)
                                                   datapath=datapath, gtifile=gtifile, 
                                                   nuradius=nuradius, guess=guess, onegauss=onegauss,
                                                   adjacent_grades=True, pile_up_corr=True,
                                                   force_nustar=True,
                                             
                                                   #aia related
                                                   load_prepped_aia=data, 
                            
                                                   #xrt related
                                                   xrtmethod='Average', real_xrt_err=True, xrt_path=xrt_path,
                                                   xrt_exposure_dict = exposure_dict, plot=plot_xrt, #(this plots xrt)
                                                   input_xrt_region="circle", input_xrt_region_dict=region_input)
                    


                else:
                    dodem.dodem(time, bl, tr, xrt=xrt, aia=aia, nustar=nustar, name=name,
                                            plotMK=plotMK, minT=minT, maxT=maxT,
                                            plotresp=False, working_directory=working_dir,
                                            default_err=0.2, path_to_dodem=path_to_dodem,
                    
                                            #demreg related
                                            rgt_fact=1.2, max_iter=30,
                                            reg_tweak=1, mc_in=True, mc_rounds=100, 
                                            
                                            #nustar related 
                                            combine_fpm=True, nuenergies=nuenergies, 
                                            datapath=datapath, gtifile=gtifile,
                                            nuradius=nuradius, guess=guess, onegauss=onegauss,
                                            adjacent_grades=True, pile_up_corr=True,
                                            force_nustar=True,
                    
                                            #aia related
                                            load_prepped_aia=data, 
        
                                            #xrt related
                                           xrtmethod='Average', real_xrt_err=True, xrt_path=xrt_path,
                                            xrt_exposure_dict=exposure_dict, plot_xrt=plot_xrt,
                                            input_xrt_region="circle", input_xrt_region_dict=region_input)
        print('')

        if method in ['input', 'double']:
            fpm='A'
            if method=='input':
                regfiles = glob.glob(working_dir+'gauss_cen_'+obsid+'_'+fpm+'_user_input*.reg')
            if method=='double':
                regfiles = glob.glob(working_dir+'gauss_cen_'+obsid+'_'+fpm+'_*.reg')
                
            regfiles.sort()
            
            for i in range(0, len(directories)):
                #Time intervals for this region, orbit
                time_intervals = all_all_time_intervals[i][o]

                regfile = regfiles[i]
                print(directories[i], regfiles[i])

                for time in time_intervals:
                    res = iac.read_interval_dicts(time, where=orbit_aia_dir, bltr=True)
                    if res is None:
                        continue
                    datas, bl, tr, xrt_region_inputs = res
                    data = datas['region'+str(i)]
                    region_input = xrt_region_inputs[i]

                    if high_temp_analysis:
                        dodem.high_temp_analysis(time, bl, tr, xrt=xrt, aia=aia, nustar=nustar, name2=name,
                                                plotMK=plotMK, highT=7.2, #(minT, maxT are varied)
                                                working_directory=working_dir, #(plotresp set false in high_temp_analysis)
                                                default_err=0.2, path_to_dodem=path_to_dodem,
                                                demmethod='DEMREG',

                                 
                                                #demreg/xrt_iterative related
                                                rgt_fact=1.2, max_iter=30, 
                                                reg_tweak=1, #mc_in, mc_rounds hard coded in high_temp_analysis (same as below)

                            
                                                #nustar=related
                                                nuenergies = nuenergies, #combine_fpm hard coded in high_temp_analysis (same as below)
                                                datapath=datapath, gtifile=gtifile, 
                                                nuradius=nuradius, edit_regfile=False,
                                                regfile=regfile,
                                                adjacent_grades=True, pile_up_corr=True,
                                                force_nustar=True,

                                                #aia related
                                                load_prepped_aia=data, 

                                                #xrt related
                                                xrtmethod='Average', real_xrt_err=True, xrt_path=xrt_path,
                                                xrt_exposure_dict=exposure_dict, plot=plot_xrt, #(this plots xrt)
                                                input_xrt_region="circle", input_xrt_region_dict=region_input)



                                                 
                    else:
                        dodem.dodem(time, bl, tr, xrt=xrt, aia=aia, nustar=nustar, name=name,
                                    plotMK=plotMK, minT=minT, maxT=maxT,
                                    plotresp=False, working_directory=directories[i],
                                    default_err=0.2, path_to_dodem=path_to_dodem,
            
                                    #demreg related
                                    rgt_fact=1.2, max_iter=30,
                                    reg_tweak=1, mc_in=True, mc_rounds=100, 
                                    
                                    #nustar related 
                                    combine_fpm=True, nuenergies=nuenergies, 
                                    datapath=datapath, gtifile=gtifile,
                                    nuradius=nuradius, edit_regfile=False,
                                    regfile=regfile,
                                    adjacent_grades=True, pile_up_corr=True,
                                    force_nustar=True,
            
                                    #aia related
                                    load_prepped_aia=data, 
        
                                    #xrt related
                                   xrtmethod='Average', real_xrt_err=True, xrt_path=xrt_path,
                                    xrt_exposure_dict=exposure_dict, plot_xrt=plot_xrt,
                                    input_xrt_region="circle", input_xrt_region_dict=region_input)




def get_key_resultfiles(key, withparams=False):
    
    import glob

    if withparams:
        res_files = glob.glob('./compact_results/*'+key+'*withparams.pickle')
    else:
        rfs = glob.glob('./compact_results/*'+key+'*.pickle')
        res_files = [f for f in rfs if 'withparams' not in f]
        #print(res_files)
        
    res_files.sort()
        
    if 'region' in res_files[0]:
        rez = []
        zeros = [f for f in res_files if 'region_0' in f]
        if zeros:
            rez.append(zeros)
        ones = [f for f in res_files if 'region_1' in f]
        if ones:
            rez.append(ones)


        if zeros and ones:
            return rez, True
        else:
            #print(rez[0])
            return rez[0], False


    else:
        return res_files, False
        

def get_above10s(key='', all=True, plot=False, time_weighted=False, seconds_per=5, return_loc=False,
                EMT=False, regions_return=False):
    
    import glob

    if all:
        res_files = glob.glob('./compact_results/*')
    elif key:
        res_files, tworegion = get_key_resultfiles(key, withparams=False)
        if tworegion:
            res=[]
            if res_files[0]:
                res0 = extract_and_plot_above10s(res_files[0], key=key, 
                                             plot=plot, time_weighted=time_weighted, seconds_per=seconds_per, plotlabel='Region 0',
                                                EMT=EMT)
                res.append(res0)
                
            if res_files[1]:
                res1 = extract_and_plot_above10s(res_files[1], key=key, 
                                             plot=plot, time_weighted=time_weighted, seconds_per=seconds_per, plotlabel='Region 1',
                                                EMT=EMT)  
                res.append(res1)

            if regions_return:
                return res 

            res_files = res_files[0] + res_files[1]
       
    else:
        print('Either set all=True, or select a key!')
        return

    res = extract_and_plot_above10s(res_files, key=key, plot=plot, time_weighted=time_weighted, seconds_per=seconds_per)

    if return_loc:
        with open('all_targets.pickle', 'rb') as f:
            data = pickle.load(f)

        ARDict = data[key]
        loc = ARDict['loc']
        
        return res, loc
    else:
        return res

    
def check_for_flare(time, starts, stops):

    flare=False
    #If there is a flare that starts before time0 and stops after time0
    b4 = [s < time[0] for s in starts]
    ea = [s > time[0] for s in stops]
    es = np.where(np.logical_and(b4, ea))

    if es[0].size > 0:
        flare=True

    #If there is a flare that starts between time0 and time1
    b4 = [s > time[0] for s in starts]
    ea = [s < time[1] for s in starts]
    es = np.where(np.logical_and(b4, ea))

    if es[0].size > 0:
        flare=True

    return flare


def get_saved_flares(add_stdv_flares=True):

    import pandas as pd
    import astropy.time as time

    risefactor=1
    fallfactor=2
    
    df = pd.read_csv('fpmA.csv')
    starts = df['flare_start'].values
    stops = df['flare_end'].values
    num=len(starts)
    
    from astropy import units as u
    durs = [(time.Time(stops[i])-time.Time(starts[i])).to(u.s) for i in range(0,num)]
    early_starts = [(time.Time(starts[i])-risefactor*durs[i]).datetime for i in range(0,num)]
    late_stops = [(time.Time(stops[i])+fallfactor*durs[i]).datetime for i in range(0,num)]

    #early_starts = [(time.Time(s)-2*u.min) for s in starts]
    #late_stops = [(time.Time(s)+4*u.min) for s in stops]

    if add_stdv_flares:
        with open('stdv_flares.pickle', 'rb') as f:
            data = pickle.load(f)

        flarewindows = np.array(data['stdv_flares'])
        num=len(flarewindows)
        durs = [time.Time(flarewindows[i,1])-time.Time(flarewindows[i,0]) for i in range(0, num)]
        early_starts.extend([(time.Time(flarewindows[i,0])-risefactor*durs[i]).datetime for i in range(0, num)])
        late_stops.extend([(time.Time(flarewindows[i,1])+fallfactor*durs[i]).datetime for i in range(0, num)])
        #early_starts.extend([(time.Time(s)-2*u.min).datetime for s in flarewindows[:,0]])
        #late_stops.extend([(time.Time(s)+4*u.min).datetime for s in flarewindows[:,1]]) 

    return [early_starts, late_stops]

def extract_and_plot_above10s(res_files, key='', plot=False, time_weighted=False, seconds_per=5, plotlabel='', show=False,
                             EMT=False):


    import pandas as pd
    import astropy.time

    flare_res = get_saved_flares(add_stdv_flares=True)
    early_starts = flare_res[0]
    late_stops = flare_res[1]
    
    all_above10s_flares = []
    all_above10s_non = []
    all_above10s=[]
    all_above10ls_flares = []
    all_above10ls_non = []
    all_above10ls=[]
    all_above10hs_flares = []
    all_above10hs_non = []
    all_above10hs=[]    
    times = []
    durations =  []
    time1 = 0

    EMTall = []
    EMTflare = []
    EMTnon = []
    
    for f in res_files:
        #flare=False
    
        data, timestring, time = viz.load_DEM(f)
        times.append(time)
        dur = (time[1]-time[0]).to(u.s)
        if time1==0:
            time1=time[0]
        durations.append(dur)
        time_mult = round(dur.value/seconds_per)

        flare = check_for_flare(time, early_starts, late_stops)
        res = viz.get_DEM_params(f)

        #print(res)
    
        m1, max1, above5_, above7_, above10_, \
            above_peak, below_peak, above_635, below_635, \
               chanax, dn_in, edn_in, \
                    powerlaws, EMT_all, EMT_thresh = res

        #if above10_[0] > 1e24:
        #    print(f, above10_)

        if len(data['dn_in']) <= 6:
            print('Suspect missing nustar:')
            print(chanax)
            print(f)
            continue


    
        if flare:
            if time_weighted:
                #print(time_mult)
                #print([above10_[0] for i in range(0,time_mult)])
                all_above10s_flares.extend([above10_[0] for i in range(0,time_mult)])
                all_above10ls_flares.extend([above10_[1] for i in range(0,time_mult)])
                all_above10hs_flares.extend([above10_[2] for i in range(0,time_mult)])
                EMTflare.extend([EMT_thresh for i in range(0,time_mult)])
            else:
                all_above10s_flares.append(above10_[0])
                all_above10ls_flares.append(above10_[1])
                all_above10hs_flares.append(above10_[2])
                EMTflare.append(EMT_thresh)
                
        else:
            if time_weighted:
                #print(time_mult)
                #print([above10_[0] for i in range(0,time_mult)])
                all_above10s_non.extend([above10_[0] for i in range(0,time_mult)])
                all_above10ls_non.extend([above10_[1] for i in range(0,time_mult)])
                all_above10hs_non.extend([above10_[2] for i in range(0,time_mult)])
                EMTnon.extend([EMT_thresh for i in range(0,time_mult)])
            else:  
                all_above10s_non.append(above10_[0])
                all_above10ls_non.append(above10_[1])
                all_above10hs_non.append(above10_[2])
                EMTnon.append(EMT_thresh)
        
        if time_weighted:
            all_above10s.extend([above10_[0] for i in range(0,time_mult)])
            all_above10ls.extend([above10_[1] for i in range(0,time_mult)])
            all_above10hs.extend([above10_[2] for i in range(0,time_mult)])
            EMTall.extend([EMT_thresh for i in range(0,time_mult)])
        else:
            all_above10s.append(above10_[0])
            all_above10ls.append(above10_[1])
            all_above10hs.append(above10_[2])
            EMTall.append(EMT_thresh)

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
        
        logbins = np.geomspace(1e17, 1e25, 55)
        
        fig, axes = plt.subplots(2, 1, figsize=(15,6), tight_layout = {'pad': 1})

        ax=axes[0]
        #ax.hist(all_above10s, bins=logbins, color='green', edgecolor='black', label='All bins')
        ax.hist(np.array(all_above10ls_non)*factor, bins=logbins, color='aliceblue', edgecolor='black', 
                label='Uncertainty Lowbound', alpha=0.9, linestyle='dotted')
        ax.hist(np.array(all_above10hs_non)*factor, bins=logbins, color='lightcyan', edgecolor='black', 
                label='Uncertainty Highbound', alpha=0.9, linestyle='dashed')
        ax.hist(np.array(all_above10s_non)*factor, bins=logbins, color='skyblue', edgecolor='black', label='Non-flare time bins', alpha=0.9)
        #ax.hist(np.array(all_above10s_flares)*factor, bins=logbins, color='purple', edgecolor='black', label="Bins during Reed's flares")
        ax.set_xscale('log')
        ax.set_xlim([1e17,1e25])
        ax.axvline(1.8e22, color='Red')
        ax.axvline(1.5e23, color='Red')
        ax.axvspan(1.8e22, 1.5e23, alpha=0.3, color='Red', label='Ishikawa (2017) 95% Interval')
        if time_weighted:
            ax.set_ylabel('Number of '+str(seconds_per)+'s intervals')
        ax.set_xlabel('EM Integrated >10 MK')
        
        if key:
            ax.set_title(key+' '+plotlabel)
        else:
            ax.set_title('All')
        ax.legend()

        ax=axes[1]
        #ax.hist(all_above10s, bins=logbins, color='green', edgecolor='black', label='All bins')
        ax.hist(np.array(all_above10ls_flares)*factor, bins=logbins, color='lavender', edgecolor='black', 
                label='Uncertainty Lowbound', alpha=0.9, linestyle='dotted')
        ax.hist(np.array(all_above10hs_flares)*factor, bins=logbins, color='thistle', edgecolor='black', 
                label='Uncertainty Highbound', alpha=0.9, linestyle='dashed')
        ax.hist(np.array(all_above10s_flares)*factor, bins=logbins, color='purple', edgecolor='black', label="Bins during flares")
        ax.set_xscale('log')
        ax.set_xlim([1e17,1e25])
        ax.axvline(1.8e22, color='Red')
        ax.axvline(1.5e23, color='Red')
        ax.axvspan(1.8e22, 1.5e23, alpha=0.3, color='Red', label='Ishikawa (2017) 95% Interval')
        if time_weighted:
            ax.set_ylabel('Number of '+str(seconds_per)+'s intervals')
        ax.set_xlabel('EM Integrated >10 MK')
        ax.legend()
        
        if time_weighted:
            plt.savefig('figures_etc/'+key+'_'+str(seconds_per)+'s_bin_time_weighted_above10splot_'+plotlabel+'.png')
        else:
            plt.savefig('figures_etc/'+key+'_above10splot_'+plotlabel+'.png')

        if not show:
            plt.close()


    EMT = [EMTall, EMTflare, EMTnon]

    averages = [np.mean(all_above10s), np.mean(all_above10s_flares), np.mean(all_above10s_non)]


    return all_above10s, all_above10s_flares, all_above10s_non, np.array(averages), time1, EMT




def make_orbit_plots(working_dir, key, minT=5.6, maxT=7.2, show=False):

    
    all_time_intervals, all_time_intervals_list = tis.find_all_intervals(working_dir, shush=True, 
                                                                        missing_last=False)
    print(len(all_time_intervals))
    
    
    print('Number of orbits:', len(all_time_intervals))
    for k in range(0, len(all_time_intervals)):
        time_intervals = all_time_intervals[k]
    

        vals = viz.get_DEM_timeseries(time_intervals, working_dir, minT, maxT, key)    

        if len(vals['result_time_intervals'])==0:
            print('No sucessful DEM intervals.')
            continue
        
        
        dn_inss=vals['dn_ins']
        labels= ['NuSTAR 6-10 keV Emission', 'NuSTAR 3.5-6 keV Emission', 'NuSTAR 2.5-3.5 keV Emission']
        backcolors=['turquoise', 'lightcyan']
        color='darkcyan'
        for i in range(1,4):
            ind = i*-1
            nuvals = [din[ind] for din in dn_inss]    
            viz.pretty_orbit_timeseries(vals['result_time_intervals'], nuvals, 'Cts/s', labels[i-1],
                                color, backcolors, working_dir=working_dir, plot_flares=True, show=show)
        
        
        peaks=vals['peaks']
        peaksmk = [10**m1/1e6 for m1 in peaks]    
        
        backcolors=['pink', 'lavenderblush']
        color='Red'
            
        viz.pretty_orbit_timeseries(vals['result_time_intervals'], peaksmk, 'DEM Peak Temperature (MK)', 'DEM Peak Temperature',
                                color, backcolors, working_dir=working_dir, show=show)
        
        
        backcolors=['powderblue', 'aliceblue']
        color='Blue'
        
        above10s=np.array(vals['above10s'])
        above10s_=above10s[:,0]
        
        viz.pretty_orbit_timeseries(vals['result_time_intervals'], above10s_, 'EM (cm^-5)', 'Total EM >10 MK',
                                color, backcolors, error=True, quantity_low=above10s[:,1], quantity_high=above10s[:,2], 
                                ylog=True, comparisonbar=True, comp_band=[1.8e22, 1.5e23, 'Ishikawa (2017) 95%'],
                                    working_dir=working_dir, plot_flares=True, show=show)
        
        backcolors=['powderblue', 'aliceblue']
        color='Green'
        
        above7s=np.array(vals['above7s'])
        above7s_=above7s[:,0]
        
        viz.pretty_orbit_timeseries(vals['result_time_intervals'], above7s_, 'EM (cm^-5)', 'Total EM >7 MK',
                                color, backcolors, error=True, quantity_low=above7s[:,1], quantity_high=above7s[:,2], 
                                ylog=True, working_dir=working_dir, show=show)
        
        backcolors=['powderblue', 'aliceblue']
        color='Purple'
        
        above5s=np.array(vals['above5s'])
        above5s_=above5s[:,0]
        
        viz.pretty_orbit_timeseries(vals['result_time_intervals'], above5s_, 'EM (cm^-5)', 'Total EM >5 MK',
                                color, backcolors, error=True, quantity_low=above5s[:,1], quantity_high=above5s[:,2], 
                                ylog=True, working_dir=working_dir, show=show)
        
        
        backcolors=['khaki', 'lemonchiffon']
        color='Orange'
        
        val=np.array(vals['low_powers'])
        
        viz.pretty_orbit_timeseries(vals['result_time_intervals'], val, 'Index', 'Lower Power Law',
                                color, backcolors, error=False, working_dir=working_dir, show=show)
        
        
        backcolors=['khaki', 'lemonchiffon']
        color='Red'
        
        val=np.array(vals['hi_powers'])*-1
        
        viz.pretty_orbit_timeseries(vals['result_time_intervals'], val, 'Index', 'Upper Power Law',
                                color, backcolors, error=False, working_dir=working_dir, show=show)






def get_region_directories(key, method='input'):

    import glob
    import os

    with open('all_targets.pickle', 'rb') as f:
        data = pickle.load(f)

    ARDict = data[key]
    working_dir = ARDict['working_dir']
    
    if method=='fit':
        directories=[working_dir]


    if method=='input':
        directories = [f+'/' for f in glob.glob(working_dir+'/region_*') if os.path.isdir(f)]
        directories.sort()

    if method=='double':
        gauss_stats = ARDict['gauss_stats']
        sep_axis = gauss_stats[0][0]
        if sep_axis=='EW':
            directions = ['east', 'west']
        elif sep_axis=='SN':
            directions = ['south', 'north']
            
        directories=[]    
        for r in directions:
            directories.append(working_dir+'/'+r+'/')

    return directories


def do_stdv_analysis(key, method='input', show=True):

    import nustar_utilities as nuutil
    import glob

    with open('all_targets.pickle', 'rb') as f:
        data = pickle.load(f)

    ARDict = data[key]
    
    id_dirs = ARDict['datapaths']
    obsids = ARDict['obsids']
    directories = get_region_directories(key, method=method)

    all_newwindows=[]

    for dd in directories:
    
        all_time_intervals, all_time_intervals_list = tis.find_all_intervals(dd, shush=True, 
                                                                            missing_last=False)
    
        for ind in range(0, len(all_time_intervals)):
            datapath=id_dirs[ind]
            obsid=obsids[ind]
            print(datapath)

            evt_data, hdr = nu.return_submap(datapath=datapath, fpm='A', return_evt_hdr=True)
            time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
            timerange = [time0.tt.datetime, time1.tt.datetime]
            from datetime import timezone
            
            #Comment second line if you're not using this same example nustar orbit
            #Edit it to include only the desired time interval (default- all times in file) once you've run this once
            #timerange=[]
            #timerange=[datetime.datetime(2018, 5, 29, 16, 00), datetime.datetime(2018, 5, 29, 16, 55)]
            timerange = [t.replace(tzinfo=timezone.utc) for t in timerange]


            evtA = glob.glob(datapath+'/event_cl/*A06_cl.evt')
            evtB = glob.glob(datapath+'/event_cl/*B06_cl.evt')
            hkA  = glob.glob(datapath+'/hk/*A_fpm.hk')
            hkB  = glob.glob(datapath+'/hk/*B_fpm.hk')

            import lightcurves as lc
            
            lc.prepare_nustar_lightcurves(evtA, evtB, hkA, hkB, timebin=15, erange=[2.,4.], 
                                          livetime_corr=False, save_dir=dd, event_stats=True)
            lc.prepare_nustar_lightcurves(evtA, evtB, hkA, hkB, timebin=15, erange=[4.,6.], 
                                          livetime_corr=False, save_dir=dd, event_stats=True)
            lc.prepare_nustar_lightcurves(evtA, evtB, hkA, hkB, timebin=15, erange=[6.,10.], 
                                          livetime_corr=False, save_dir=dd, event_stats=True)
            
            newwindows = lc.plot_with_stdv(aia_inputs=[], fexviii=False, nustar_inputs=[[2.,4.],[4.,6.],[6.,10.]], 
                   goes_inputs=[], smooth=24, 
                   plot_each=True, plot_logic=True, remove_fexviii_max=False, 
                   analyze_transients=True, transient_number=3,
                   timerange=timerange,
                  excluded_range=[], save_dir=dd, savestring=key+'_'+obsid, show=show)
            

            all_newwindows.extend(newwindows)


    return all_newwindows





def make_summary_lcs(key, method='input', show=True):

    
    from matplotlib import pyplot as plt
    import copy
    import matplotlib.dates as mdates
    import nustar_utilities as nuutil
    import lightcurves as lc
    import glob

    minT=5.6
    maxT=7.2
    nustar_acc_color='xkcd:apple'
    nustar_cts_color=['xkcd:sky blue', 'xkcd:cerulean', 'xkcd:periwinkle', 'xkcd:cyan']

    with open('all_targets.pickle', 'rb') as f:
        data = pickle.load(f)

    ARDict = data[key]
    
    id_dirs = ARDict['datapaths']
    directories = get_region_directories(key, method=method)

    for dd in directories:
    
        all_time_intervals, all_time_intervals_list = tis.find_all_intervals(dd, shush=True, 
                                                                            missing_last=False)
    
        for ind in range(0, len(all_time_intervals)):
            
            datapath=id_dirs[ind]
        
            time_intervals = all_time_intervals[ind]
            vals = viz.get_DEM_timeseries(time_intervals, dd, minT, maxT, key)   
            time_intervals = vals['result_time_intervals']
            
            if not time_intervals:
                print('This key/region/orbit had no sucessful DEMs.')
                print(key)
                print(dd)
                print(ind)
                print('')
                return
            
            evt_data, hdr, obsid = nu.return_submap(datapath=datapath, fpm='A', return_evt_hdr=True, return_obsid=True)
            time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
            timerange = [time0.datetime, time1.datetime]
            print(timerange[0].strftime('%H-%M-%S'), timerange[1].strftime('%H-%M-%S'))
            
            
            from datetime import timezone
            newtimerange = [t.replace(tzinfo=timezone.utc) for t in timerange]
            
            axiscolor='black'
            
            labelfontsize=10
            tickfontsize=17
            
            
            fig, axes = plt.subplots(6, 1, figsize=(15, 15), sharex=True)
            plt.subplots_adjust(hspace=0)
            
            
            
            # AIA PLOTS + NuSTAR DEM INPUT PLOTS + >10 MK EM plots #======================================================================
            
            lw=2
            times = [t[0].datetime for t in time_intervals]
            midtimes=[(t[0]+(t[1]-t[0]).to(u.s)/2).datetime for t in time_intervals]
            
            
            starttime = (time_intervals[0][0]-120*u.s).datetime
            stoptime = (time_intervals[-1][1]+120*u.s).datetime
            
            times_ = copy.deepcopy(times)
            times_.append(time_intervals[-1][1].datetime)
            
            
            dn_ins = vals['dn_ins']
            chanaxs = vals['chanaxs']
            
            allcolors = lc.make_colors(10)
            
            ax=axes[0]
            normin=1
            for i in range(0, 6):  
                aiavals = [din[i] for din in dn_ins] 
                label = chanaxs[0][i]
                color = allcolors[i]
                normvals = np.array(aiavals)/np.max(aiavals)
                if np.min(normvals) < normin:
                    normin=np.min(normvals)
                ax.stairs(normvals, times_, linewidth=lw, color=color, 
                          label=label,
                         baseline=None)
            
            ax.set_ylim([normin*0.95, 1.01])
            
            ax=axes[1]
            normin=1
            for i in range(0, 3): 
                ii = (i+1)*-1
                aiavals = [din[ii] for din in dn_ins] 
                label = chanaxs[0][ii]
                color = allcolors[ii]
                normvals = np.array(aiavals)/np.max(aiavals)
                if np.min(normvals) < normin:
                    normin=np.min(normvals)
                ax.stairs(normvals, times_, linewidth=lw, color=color, 
                          label=label,
                         baseline=None)
            
            ax.set_ylim([normin*0.95, 1.01])
    
    
            ax=axes[5]
            above10s=np.array(vals['above10s'])
            above10s_=above10s[:,0]
            label='Total EM >10 MK'
            ylabel='EM (cm^-5)'
            color='Blue'
            ax.stairs(above10s_, times_, linewidth=lw, color=color, 
                          label=label,
                         baseline=None)
            
            quantity_low=above10s[:,1]
            quantity_high=above10s[:,2]
            
            error=True
            if error:
                lp = np.hstack([quantity_low, quantity_low[-1]])
                hp = np.hstack([quantity_high, quantity_high[-1]])
    
                fill = ax.fill_between(times_, lp, hp, step="post", 
                                     color=color, alpha=0.1) 
            
            ax.set_ylabel(ylabel)
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=10))
            ax.set_yscale('log')
            
            comp_band=[1.8e22, 1.5e23, 'Ishikawa (2017) 95%']
            comparisonbar=True
            if comparisonbar:
                ax.axhspan(comp_band[0], comp_band[1], color='brown', alpha=0.5, label=comp_band[2])   
                ax.legend()
            
            
            
            #==================================================================================================================================
            
            #NUSTAR - ALL-FOV, FINE TIME RESOLUTION ==========================================================================================
            
            tr=newtimerange
                
            evtA = glob.glob(datapath+'/event_cl/*A06_cl.evt')
            evtB = glob.glob(datapath+'/event_cl/*B06_cl.evt')
            hkA  = glob.glob(datapath+'/hk/*A_fpm.hk')
            hkB  = glob.glob(datapath+'/hk/*B_fpm.hk')
            
            #Load in the evt file (has the list of photons)
            evtdataA, hdrA = lc.load_nufiles(evtA[0])
            # Load in the hk file (has the livetime info)
            lvdataA, lvhdrA = lc.load_nufiles(hkA[0])
            evtdataB, hdrB = lc.load_nufiles(evtB[0])
            lvdataB, lvhdrB = lc.load_nufiles(hkB[0])
    
            
            eranges=[[2.5,3.5],[3.5,6.], [6.,10.], [10.,15.]]
            labels = ['NuSTAR 2.5-3.5 keV', 'NuSTAR 3.5-6 keV', 'NuSTAR 6-10 keV', 'NuSTAR 10-15 keV']
            eind=0
            #     acc=0
            for erange in eranges:
                kevA = evtdataA['PI']*0.04+1.6
                erange_evtdataA = evtdataA[np.where(np.logical_and(kevA > erange[0],kevA < erange[1]))]
                kevB = evtdataB['PI']*0.04+1.6
                erange_evtdataB = evtdataB[np.where(np.logical_and(kevB > erange[0],kevB < erange[1]))]
            
            
                res = lc.get_a_nustar_lightcurve(erange_evtdataA, hdrA, lvdataA, lvhdrA, timebin=10, livetime_corr=True, event_stats=True)
                times_converted, countrate, lvt, counts, acc_sample, rej_sample, all_sample = res
              
                keepinds = np.nonzero(np.logical_and(times_converted > tr[0], times_converted < tr[1]))
                times_converted=times_converted[keepinds]
                countrateA=countrate[keepinds]   

                acc_sample=acc_sample[keepinds]
                rej_sample=rej_sample[keepinds]
                evsumA=acc_sample/(acc_sample+rej_sample)*100

                    
                res = lc.get_a_nustar_lightcurve(erange_evtdataB, hdrB, lvdataB, lvhdrB, timebin=10, livetime_corr=True, event_stats=True)
                times_converted, countrate, lvt, counts, acc_sample, rej_sample, all_sample = res
    
            
                keepinds = np.nonzero(np.logical_and(times_converted > tr[0], times_converted < tr[1]))
                times_converted=times_converted[keepinds]
                countrateB=countrate[keepinds]

                acc_sample=acc_sample[keepinds]
                rej_sample=rej_sample[keepinds]
        
                evsumB=acc_sample/(acc_sample+rej_sample)*100


                evsum = (np.array(evsumA)+np.array(evsumB))/2.
                axes[3].plot(times_converted, evsum, label='Accepted Event %', color=nustar_acc_color)
    
                totalrate = countrateA+countrateB
                totalrate_vals = totalrate[np.isfinite(totalrate)]
                maximum = np.nanmax(totalrate_vals) 
                axes[2].plot(times_converted, totalrate/maximum, label=labels[eind], color=nustar_cts_color[eind]) 

            
                eind+=1
                
            axes[2].set_ylim([0,1.01])
            
            #==================================================================================================================================
            
            #Add Goes Curve (both days) #=========================================================================================================
            tr=timerange

            goes_number = ARDict['goes_satellite']
            
            #print(tr)
            lc.get_goes(tr, satellite=goes_number)
            instrument='GOES'
            data = lc.load_lightcurves(instrument)
            
            ylabel = data['GOES flux label']
            goestimes = data['GOES Times']
            xrsbcounts = data['XRSB counts']
            xrsblabel = data['XRSB label']
            gts = [t.datetime for t in goestimes]
            
            axg = axes[4] #.twinx()
            #xrs=[x.value for x in all_xrsb]
            axg.plot(gts, xrsbcounts, label='GOES '+xrsblabel, color='red')
            axg.set_yscale('log')
            #axg.set_ylim([1e-7,1e-6])
            #axg.yaxis.set_ticks([1e-6, 1e-5, 1e-4], labels=['C','M','X'], color='red', fontsize=labelfontsize)
            axg.grid(True, which='major', axis='y', linestyle=':', color='red')
            
            #axg.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            #axg.xaxis.set_minor_locator(mdates.MinuteLocator(interval=10))
            
            #==================================================================================================================================
            
            # #Plot Flares from Reed's list #====================================================================================================
            # plot_flares=True
            # if plot_flares:
            #     import pandas as pd
            #     import astropy.time
            #     df = pd.read_csv('fpmA.csv')
            #     starts = df['flare_start'].values
            #     stops = df['flare_end'].values
                
            
            #     early_starts = [(astropy.time.Time(s)-2*u.min).datetime for s in starts]
            #     late_stops = [(astropy.time.Time(s)+4*u.min).datetime for s in stops]  
            
            #     for j in range(0, len(early_starts)):
            #         for ax in axes:
            #             ax.axvspan(early_starts[j], late_stops[j], alpha=.25, color='purple')
            #             ax.set_xlim(timerange[0], timerange[1])
            

            # #Plot Flares stdv analysis list #====================================================================================================
            # plot_stdv_flares=True
            # if plot_stdv_flares:
            #     #print('got here')

            #     with open('stdv_flares.pickle', 'rb') as f:
            #         data = pickle.load(f)

            #     flarewindows = np.array(data['stdv_flares'])
                
            #     early_starts = [(astropy.time.Time(s)-2*u.min).datetime for s in flarewindows[:,0]]
            #     late_stops = [(astropy.time.Time(s)+4*u.min).datetime for s in flarewindows[:,1]]  
            
            #     for j in range(0, len(early_starts)):
            #         for ax in axes:
            #             ax.axvspan(early_starts[j], late_stops[j], alpha=.25, color='green')
            #             ax.set_xlim(timerange[0], timerange[1])
                        
            
            plot_flares=True
            plot_stdv_flares=True
            
            if plot_flares:
                flare_res = get_saved_flares()
                early_starts = flare_res[0]
                late_stops = flare_res[1]

                for j in range(0, len(early_starts)):
                    for ax in axes:
                        ax.axvspan(early_starts[j], late_stops[j], alpha=.25, color='blue')
                        ax.set_xlim(timerange[0], timerange[1])

                if plot_stdv_flares:
                    flare_res = get_saved_flares(add_stdv_flares=True)
                    early_starts = flare_res[0]
                    late_stops = flare_res[1]
                    for j in range(0, len(early_starts)):
                        for ax in axes:
                            ax.axvspan(early_starts[j], late_stops[j], alpha=.25, color='red')
                            ax.set_xlim(timerange[0], timerange[1])
                    

            
            
            #==================================================================================================================================


            #Plot SAAs and Pointing Shifts #====================================================================================================
            plot_badtimes=True
            if plot_badtimes:
                import use_correlator as uc
                import astropy.time
                t_corr=2
                #minimum length of suborbit
                min_t = 30
                #minimum number of steps needed in a sufficiently-sized suborbit, based on the above. 
                min_step = int(np.ceil(min_t/t_corr))
                res = uc.get_suborbits(id_dirs[ind], t_corr, min_step, plot=False)
                bad_grouptimes = res[3]
                
            
                for j in range(0, len(bad_grouptimes)):
                    for ax in axes:
                        ax.axvspan(bad_grouptimes[j][0], bad_grouptimes[j][1], alpha=.25, color='orange')
                        ax.set_xlim(timerange[0], timerange[1])
            
            #==================================================================================================================================
            
            
            leg0 = axes[0].legend()
            leg1 = axes[1].legend()#(loc='best')
            leg2 = axes[2].legend()#loc='best')
            leg3 = axes[3].legend()
            leg4 = axes[4].legend()
            leg5 = axes[5].legend()
            
            
            plt.savefig(dd+obsid+'_lc_summary.png', transparent=False)

            if not show:
                plt.close()
                
        print('')

    

def plot_with_discriminator(axes, resfiles, discriminator, options, label, paramx, paramy, optionnumbers, caselabels,
                  xindex='', yindex='', flare=False, startstop=[], update_flarecheck=False):

    """

    Takes in plot axes, a set of result files from a given key/region, and selections re what parameters we would like to plot.

    paramx, paramy - to be plotted on the horizontal and vertical axes, respectively. Should be keys in the result file dictionaries. 

    xindex, yindex - for result params that come in lists/arrays (say, [value, lowbound, highbound]), indicates which to use. If set to
                        something other than an integer, will not be used (assume single value). 

    label - uniquely identifies the key/region; used for plot labels.

    Also takes in DISCRIMINATOR – a value of the key/region with which we would like to separate out different cases and plot them
    on different axes. The number of axes should be equal to the number of possible values of discriminator. For recordkeeping, we 
    also need: 
    
    optionnumbers – array (of the same length as the number of possible values of the discriminator), containing the number of each which
                    has been found so far (we'll add to the number associated with the value of discriminator each time this is run). 

    caselabels - list (of the same length as the number of possible values of the discriminator) containting lists of the LABELS of each
                    key/region (we'll add the present label value to the list for the present value of the discriminator.)
    
    If startstop is not set, take all intervals.
    If startstop contains lists of paired start and stop times for flares, then:
        if Flare == True: take only times intersecting with a flare
        if Flare == False: take only times that do not intersect with a flare

    update_flarecheck - set True to reevaluate the resultfiles to check if they are for time intervals intersecting with a flare. This
                            should be done if the flarelist changes. 
                            
    """

    from astropy import units as u

    greenz = ['green', 'darkgreen', 'mediumseagreen', 'lightgreen', 'forestgreen', 'seagreen', 'darkseagreen', 'lime', 'springgreen', 
          'mediumaquamarine', 'aquamarine', 'lawngreen', 'green', 'darkgreen', 'mediumseagreen', 'lightgreen', 'forestgreen', 'seagreen', 
              'darkseagreen', 'lime', 'springgreen', 
          'mediumaquamarine', 'aquamarine', 'lawngreen']

    bluez = ['blue', 'royalblue', 'mediumblue', 'darkblue', 'navy', 'cornflowerblue', 'dodgerblue', 'deepskyblue', 'skyblue',
         'lightblue', 'lightskyblue', 'steelblue', 'lightsteelblue', 'powderblue', 'darkslateblue', 'slateblue', 'lightslategray',
         'slategray']


    reds = ['lightcoral', 'indianred', 'brown', 'red', 'salmon', 'tomato']

    purples = ['plum', 'violet', 'purple', 'darkviolet', 'darkorchid', 'thistle']

    colorz = [greenz, bluez, reds, purples]

    first = 0

    #If we're using a list of flaretimes to decide whether to plot contents of these resultfiles...
    if startstop:
        early_starts = startstop[0]
        late_stops = startstop[1]

    valuex=[]
    valuey=[]   
    #For each result...
    for r in resfiles:
        #load in file
        with open(r, 'rb') as f:
            data = pickle.load(f)

        #If we're checking for a flare...
        if startstop:
            time = data['time_interval']
            #If there's already a note about whether it's a flare in the result file (and we aren't re-checking)
            if 'flare?' in data.keys() and update_flarecheck == False:
                flarecheck = data['flare?'] 
            else:
                #Check if this time interval intersects with a flare
                flarecheck = check_for_flare(time, early_starts, late_stops)
                #If set, update the result dictionary to make a note of the result.
                if update_flarecheck:
                    data['flare?'] = flarecheck
                    with open(r, 'wb') as f:
                         # Pickle the 'data' dictionary using the highest protocol available.
                         pickle.dump(data, f, pickle.HIGHEST_PROTOCOL) 
                        
            #If this time interval does not satisfy our condition re flaring, don't plot.            
            if flarecheck != flare:
                #print(time[0].strftime('%H-%M-%S'))
                continue
                
        #Check if x,yindex set; in either case, fetch the values needed.  
        if type(xindex) == int:
            valuex.append(data[paramx][xindex])
        else:
            valuex.append(data[paramx])
            
        if type(yindex) == int:
            valuey.append(data[paramy][yindex])
        else:
            valuey.append(data[paramy])

    #Check to see if the discriminator is equal to any of the options.
    for o in range(0, len(options)):
        #If so, plot it in the corresponding place.
        if discriminator == options[o]:
            ax=axes[o]
            #If this is the first time we're plotting a result from this specific key/region, add the label to our list.
            if label not in caselabels[o]:
                optionnumbers[o]+=1
                caselabels[o].append(label)
                first = 1
            thecolor = colorz[o][int(optionnumbers[o])]
            #If this is the first time we're plotting a result from this specific key/region, plot with label
            if first:
                ax.scatter(valuex, valuey, color=thecolor, label=label)
                ax.set_title(options[o]+' regions, '+paramy+' vs. '+paramx)
                firstl=0
            else:
                ax.scatter(valuex, valuey, color=thecolor)             


    return optionnumbers, caselabels








    
    