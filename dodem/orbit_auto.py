"""
Code related to automizing (as much as possible) full-orbit DEM analysis.
"""
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import pathlib
import glob
import os

from astropy import units as u 

from operator import itemgetter
from itertools import groupby

from nustar_tools.trackers import CoordinateTracker as ct
from nustar_tools.trackers.TrackerCorrelator import TrackerCorrelator
from nustar_tools.utils import utilities

import time_interval_selection as tis
import nustar_dem_prep as nu
import initial_analysis as ia
import nustar_utilities as nuutil


def get_correlator(id_dir, fpm, t_corr):

    """
    Wrapper for getting a correlator object for a given NuSTAR OBSID + fpm. As demonstrated in
    https://github.com/masek014/nustar_tools/blob/main/nustar_tools/trackers/examples/correlator.ipynb
    """

    time_range = ct.get_observation_time(id_dir, fpm)
    
    att_tracker = ct.AttitudeTracker(id_dir)
    att_tracker.read_data(time_range)

    raw_tracker = ct.CentroidTracker(id_dir, fpm, 3, 'RAW')
    raw_tracker.initialize()

    correlator = TrackerCorrelator(att_tracker, raw_tracker)
    #The timestep set here will be the resulting timestep of the correlator output
    correlator.correlate_trackers(time_step=t_corr)
    correlator.make_overview()

    return correlator
    
    
def groupem(indices, min_length):
    """
    Helper. 
    
    Make groups of consecutive values from a longer monotonically increasing
    (but not necessarily consecutive) list of intergers (indices). Save only
    groups of length greater than min_length. 
    """
    
    groups=[]
    data = indices
    for key, group in groupby(enumerate(data), lambda i: i[0] - i[1]):
        group = list(map(itemgetter(1), group))
        if len(group) >= min_length:
            groups.append(group)
    return groups

def get_group_info(groups, times, printtimes=False):

    """
    Helper.
    
    For a list of groups of indices (e.g. output of groupem() ) and a list, "times" to which the indices refer
    (i.e. the times corresponding to the grouped indices), return a list of all indices in all the groups (groupinds)
    as well as a list of start/stop times for each group (grouptimes). 
    """
    
    #get start/end time for each group; get only the "good" indices that are PART OF GROUPS
    grouptimes, groupinds=[], []
    for i in range(0, len(groups)):
        timez = [times[groups[i][0]], times[groups[i][-1]]]
        if printtimes:
            print(timez[0])
            print(timez[1])
            print('')
        grouptimes.append(timez)
        groupinds.extend(groups[i])

    return groupinds, grouptimes


def get_correlator_groups(correlator, min_step, printtimes=False):

    """
    Wrapper: for a given correlator object (for an OBSID, fpm), identify groups of "good" times, where the pointing is
    predicted to be stable. Exclude times including SAA and pointing shifts (the latter as well as possible). 
    
    We want "good" indices: so amplitude is NOT above threshold (indicates pointing shift), and IS non-zero (zero indicates
    SAA passage). These will be times for which we can run the DEM time interval selection process. 

    """
    corramp = correlator.amplitude 
    goodarr = np.logical_and(corramp > 0, corramp < correlator.threshold)
    good_inds = np.nonzero(goodarr)[0]

    #Let's identify groups of consecutive good indices of length greater than our minimum length.
    groups=groupem(good_inds, min_step)
    #Let's get all indices in groups, and start/end times
    groupinds, grouptimes = get_group_info(groups, correlator.times, printtimes=printtimes)
    #Let's make a boolean array covering the whole time, checking if each index is in a group
    grouparr = np.isin(np.arange(0,len(corramp),1), groupinds)

    #print(correlator.times[0], correlator.times[-1])

    return groups, grouptimes, grouparr, correlator.times


def get_suborbits(id_dir, t_corr, min_step, plot=False):

    """
    
    We want to split the data time interval into "suborbits" that each contain no correlator-identified pointing shifts, 
    or SAAs. The idea is that a pointing correction (vs AIA) for one interval within the suborbit is good for the whole 
    suborbit. There will be some times excluded: SAAs, and also some amount of time around NuSTAR pointing shifts. These 
    non-sub-orbits should also be saved, for possible later use (analysis methods would be needed which are okay with the 
    NuSTAR pointing being in-motion). 

    There will be some minimum suborbit time for which we can make a useful DEM (i.e. sufficient NuSTAR statistics). 
    This will depend on observation livetime, so we want to make sure it can be adjusted (see min_t below).  

    FPMA and FPMB may have slightly different identified pointing shift times. Since we are adding data together in most 
    cases, we need to make a unified list of good suborbits. To be conservative, we will only take times that both FPMA, B 
    don't identify as being part of a pointing shift, etc. 

    PROCESS:
    
    – Make a correlator object for each fpm, find groups of good times.
    – Combine the two fpm: find times where BOTH fpm are "good", then make groups of these.
    – Return:
        both_groupinds – all good indices (both fpm)
        both_grouptimes - list of start/stop times for good groups
        bad_groups – groups of bad indices (list)
        bad_ids – boolean list corresponding to bad groups. True if there are at least some non-SAA 
            times in a given bad group.

    """


    fpms=['A', 'B']
    colors=['red', 'blue']
    if plot:
        fig = plt.figure(figsize=(15,5))
        axs = fig.add_subplot(1,1,1)
        i=0
        ylim=[0,1]
    goodarrs = []
    times = []
    for j in range(0,2):
        evt_file = glob.glob(id_dir+'/event_cl/*'+fpms[j]+'06_cl.evt')[0]
        print(evt_file)
        nu.convert_wrapper(evt_file)
        
        correlator = get_correlator(id_dir, fpms[j], t_corr)
        groups, grouptimes, goodarr, timez = get_correlator_groups(correlator, min_step)  
        goodarrs.append(goodarr)
        times.append(timez)
        if plot:
            for g in grouptimes:
                #print(g[0])
                #print(g[1])
                #print('')
                axs.fill_between(g, *ylim, color=colors[j], alpha=0.2)
                axs.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
                axs.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))


    rez = goodarrs[0]*1
    rez2 = goodarrs[1]*1

    #The following handles one special case where the fpma, b lists are different lengths. Other cases with break it.
    #Perhaps then would be a good idea to figure out why this is happening.
    if len(times[0]) != len(times[1]):
        if len(times[1]) > len(times[0]):
            if times[1][-1] > times[0][-1]:
                times[1] = times[1][0:len(times[0])]
                rez2 = rez2[0:len(times[0])]
    
    good_inds = np.nonzero((rez+rez2)==2)[0]
    groups = groupem(good_inds, min_step)
    #note, from here on we will be using FPMB correlator only – but it is okay because times, SAA times are the same.
    both_groupinds, both_grouptimes = get_group_info(groups, correlator.times, printtimes=True)
    
    #Indices NOT part of groups
    bad_inds = np.delete(np.arange(0,len(correlator.times),1), both_groupinds)
    #Groups of indices not part of groups:
    bad_groups = groupem(bad_inds, 0) 
    bad_groupinds, bad_grouptimes = get_group_info(bad_groups, correlator.times, printtimes=False)
    #For each group, is True if there are at least some non-SAA times in the bad group (False if pure SAA).
    bad_ids = [(max(correlator.amplitude[b]) > 0) for b in bad_groups]


    return both_groupinds, both_grouptimes, bad_groups, bad_grouptimes, bad_ids




def get_suborbit_intervals(both_grouptimes, id_dir, working_dir, erange=[6.,10],
                          lctype='corr54',fast_min_factor=2, countmin=10,
                          minimum_seconds=30, centroid_region=True, 
                        twogauss=False, direction='', guess=[],
                          force_both_fpm_always=False, shush=False,
                          nuradius=150):

    """
    Having found sub-orbits with get_suborbits(), do time interval selection for each.

    For time interval selection keywords, see tis.find_time_intervals_plus() documentation.
    
    """

    bad_suborbits=[]
    all_intervals=[]

    for g in both_grouptimes:

        time = g
        timestring = time[0].strftime('%H-%M-%S')
        stopstring = time[1].strftime('%H-%M-%S')
        timestring=timestring+'_'+stopstring
    
        #Make suborbit directory
        suborbit_dir=working_dir+'/suborbit_'+timestring
    
        #Make a new working directory for prepped data/etc if it doesn't yet exist
        save_path = pathlib.Path(suborbit_dir)
        if not save_path.exists():
            save_path.mkdir()
        
        #uncomment to plot lightcurves for each suborbit
        #tis.real_count_lightcurves(id_dir, g, working_dir, erange)
    
        datapath=id_dir
        timerange=g
        res = tis.find_time_intervals_plus(datapath, timerange, working_dir, erange=erange, 
                                   lctype=lctype, fast_min_factor=fast_min_factor, countmin=countmin,
                                  minimum_seconds=minimum_seconds, centroid_region=centroid_region,
                                           force_both_fpm_always=force_both_fpm_always, shush=shush,
                                           twogauss=twogauss, direction=direction, guess=guess,
                                          nuradius=nuradius)
        if not res[0]:
            bad_suborbits.append(res[1])
        else:
            all_intervals.extend(res[0])


    return bad_suborbits, all_intervals


def find_region_dirs(working_dir):
    
    region_dirs = [f+'/' for f in glob.glob(working_dir+'/region_*') if os.path.isdir(f)]
    region_dirs.sort()
    
    return region_dirs

def find_direction_dirs(working_dir, sep_axis):

    if sep_axis=='EW':
        directions = ['east', 'west']
    elif sep_axis=='SN':
        directions = ['south', 'north']

    region_dirs = []
    for d in directions:
        region_dirs.append(working_dir+d+'/')
    
    return region_dirs

def find_all_intervals(working_dir, shush=False, missing_last=False, missing_orbit=0):

    """
    Just to take a look at all the intervals found for all the suborbits.
    
    """
    
    all_intervals = glob.glob(working_dir+'*intervals.pickle')
    all_intervals.sort()
    timerange='hi' #not using this, as supplying custom file

    all_time_intervals, all_time_intervals_list = [], []
    orbit=0
    for tt in all_intervals:
        time_intervals = tis.get_saved_intervals(timerange, custom_file=tt)
        if missing_last:
            if orbit == missing_orbit:
                time_intervals = time_intervals[0:-1]
        all_time_intervals.append(time_intervals)
        all_time_intervals_list.extend(time_intervals)
        count=0 
        if shush==False:
            for t in time_intervals:
                print(orbit,'-',count, t[0].strftime('%H-%M-%S'), t[1].strftime('%H-%M-%S'))
                count+=1
            
        orbit+=1
        if shush==False:
            print('')
        
    return all_time_intervals, all_time_intervals_list


def check_conseq(working_dir):

    all_intervals = glob.glob(working_dir+'*intervals.pickle')
    all_intervals.sort()
    print(all_intervals)
    timerange='hi' #not using this, as supplying custom file

    all_time_intervals, all_time_intervals_list = [], []
    orbit=0
    for tt in all_intervals:
        time_intervals, full_interval = tis.get_saved_intervals(timerange, custom_file=tt, return_full_range=True)
        print('')
        print('')
        print('')
        print(full_interval[0].strftime('%H-%M-%S'), full_interval[1].strftime('%H-%M-%S'))
        for i in range(0, len(time_intervals)-1):
            now = time_intervals[i]
            next = time_intervals[i+1]            
            if i==0:
                if now[0] == full_interval[0]:
                    print('good:')
                else:
                    print('first bad:')
                    print(full_interval[0].strftime('%H-%M-%S'), full_interval[1].strftime('%H-%M-%S'))
                    print(now[0].strftime('%H-%M-%S'), now[1].strftime('%H-%M-%S'))
                    print('')
                    
                    
            if now[1] == next[0]:
                print('good:')
                #print(now[0].strftime('%H-%M-%S'), now[1].strftime('%H-%M-%S'))
                #print(next[0].strftime('%H-%M-%S'), next[1].strftime('%H-%M-%S'))
                #print('')
            else:
                print('bad:')
                print(now[0].strftime('%H-%M-%S'), now[1].strftime('%H-%M-%S'))
                print(next[0].strftime('%H-%M-%S'), next[1].strftime('%H-%M-%S'))
                print('')
            if i==(len(time_intervals)-2):
                if time_intervals[i+1][1] == full_interval[1]:
                    print('good:')
                else:
                    print('last bad:')
                    print(full_interval[0].strftime('%H-%M-%S'), full_interval[1].strftime('%H-%M-%S'))
                    print(time_intervals[i+1][0].strftime('%H-%M-%S'), time_intervals[i+1][1].strftime('%H-%M-%S'))
                    print('')
                
                
    


def regfile_to_regdict(regfile, time):

    import region_fitting as rf

    offset, rad = rf.read_regfile(regfile, time[0], time[1], regRAunit='hourangle')

    regiondict = {'radius': rad.value,
              'centerx': offset[0],
              'centery': offset[1]}

    return regiondict


def check_region_emission(all_time_intervals, working_dir, grade='0_4', plot=True, efilter=[],
                         keep_problem_intervals=True):

    import region_fitting as rf

    all_percentAs, all_percentBs = [], []
    problem_plot_intervals=[]
    
    for at in all_time_intervals:
        percentAs, percentBs = [],[]
        for time_interval in at:

            time = time_interval
            timestring = time[0].strftime('%H-%M-%S')
            stopstring = time[1].strftime('%H-%M-%S')
            timestring=timestring+'_'+stopstring
        
    
            if grade=='0':
                specific_time_evt = glob.glob(working_dir+timestring+'/'+'*0_p_cl_sunpos.evt')
                
            if grade=='0_4':
                specific_time_evt = glob.glob(working_dir+timestring+'/'+'*0_4_p_cl_sunpos.evt')

            problem_plot = glob.glob(working_dir+timestring+'/'+'*_twogauss_problem_plot.png')
            if problem_plot:
                print('This time interval had an issue with two-gaussian fitting.')
                print(time_interval)
                problem_plot_intervals.append(time_interval)
                if not keep_problem_intervals:
                    continue


            specific_time_evt.sort() 
            #print(specific_time_evt)
            if not specific_time_evt:
                print('not - ', timestring)
                continue
                
            evtA=specific_time_evt[0]
            evtB=specific_time_evt[1]
    
            regionfileA = glob.glob(working_dir+timestring+'/'+'*A*.reg')[0]
            regionfileB = glob.glob(working_dir+timestring+'/'+'*B*.reg')[0]
        
            percentA = rf.check_region(evtA, time[0], time[1], regfile=True, file=regionfileA, shush=True, get_percent=True, efilter=efilter)
            percentB = rf.check_region(evtB, time[0], time[1], regfile=True, file=regionfileB, shush=True, get_percent=True, efilter=efilter)
    
            percentAs.append(percentA)
            percentBs.append(percentB)
    
    
        all_percentAs.append(percentAs)
        all_percentBs.append(percentBs)



    pcts = [all_percentAs, all_percentBs]
    colors = ['pink', 'red']
    labels = ['FPMA', 'FPMB']
    
    if plot:
        fig = plt.figure(figsize=(15,5))
        axs = fig.add_subplot(1,1,1)

    allallpercents=[]
    allalltimes=[]
    for fpm in range(0,2):
        once=0
        allpercent = pcts[fpm]
        for t in range(0, len(allpercent)):
            if plot:
                times=[]
                for tt in all_time_intervals[t]:
                    times.append(tt[0].datetime)
        
                times.append(tt[1].datetime)
                #print(len(times))
                #print(len(allpercent[t]))
                if once==0:
                    axs.stairs(allpercent[t], np.array(times), color=colors[fpm], baseline=None, label=labels[fpm])
                    once+=1
                else:
                    axs.stairs(allpercent[t], np.array(times), color=colors[fpm], baseline=None)
            
            allallpercents.extend(allpercent[t])
            allalltimes.extend(times[0:-1])
    if plot:  
        axs.set_ylim(np.min(allallpercents)*0.75, 1)
        axs.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        axs.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))
        axs.legend()

        axs.axhline(np.mean(allallpercents), linestyle='-.', color='Black')
        axs.axhline(np.mean(allallpercents)+0.05, linestyle='-.', color='Black')
        axs.axhline(np.mean(allallpercents)-0.05, linestyle='-.', color='Black')

    print('E-filter was: ', efilter)
    print('Bins >5% above mean:', len(np.nonzero(np.array(allallpercents) > np.mean(allallpercents)+0.05)[0]), ' out of ', len(allallpercents))
    print('Bins >5% below mean:', len(np.nonzero(np.array(allallpercents) < np.mean(allallpercents)-0.05)[0]), ' out of ', len(allallpercents))
    print('')

    
    allallpercents=100*np.array(allallpercents)
    print('Mean emision included: ', np.mean(allallpercents).round(2), '%')
    print('Minimum emission included: ', np.min(allallpercents).round(2), '%')
    print('Maximum emission included: ', np.max(allallpercents).round(2), '%')
    print('STDV: ', np.std(allallpercents).round(2), '%')    

    #print(len(allallpercents))
    #print(len(allalltimes))
    wheremin=np.where(allallpercents == np.min(allallpercents))[0]
    print('time of minimum: ', allalltimes[wheremin[0]])
    print('')
    


    

    return all_percentAs, all_percentBs



def get_orbit_aiamaps(aia_dir, id_dirs, wave=94):

    import sunpy.map
    from aiapy.calibrate.util import get_correction_table, get_pointing_table
    from aiapy.calibrate import register, update_pointing, degradation, estimate_error

    aiamaps = []
    for id in id_dirs:
        print(id)
        evt_data, hdr = ia.return_submap(datapath=id, fpm='A', return_evt_hdr=True)
        time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
        start = str(time0)
        files = glob.glob(aia_dir+'*'+str(wave)+'A_'+start[0:10]+'T'+start[11:13]+'*')
        print(start)
        print(files)
        amap=sunpy.map.Map(files[0])
        #ptab = get_pointing_table(amap.date - 12 * u.h, amap.date + 12 * u.h)
        #m_temp = update_pointing(amap, pointing_table=ptab)
        try:
            m = register(amap)
        except TypeError:
            amap.meta.pop('crpix1')
            amap.meta.pop('crpix2')
            print('CRPIX issue on ', files)
            m = register(amap)
    
        aiamaps.append(m)

    return aiamaps
    

def nu_aia_coalign(time_interval, working_dir, nushift, regionmethod='fit',
                   obsid='', region_dir='',
                   input_aia=[], 
                   grade='0_4', justCOM=False, save_dict=False,
                  savefigdir=[]):
    """
    nushift in x, y arcseconds
    """
    
    import pickle

    time = time_interval
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring

    if regionmethod=='fit':
        regionfileA = glob.glob(working_dir+timestring+'/'+'*A*sunpos*.reg')
        regionfileB = glob.glob(working_dir+timestring+'/'+'*B*sunpos*.reg')
        specific_time_evt = glob.glob(working_dir+timestring+'/'+'*cl.evt')

    if regionmethod=='input' or regionmethod=='double':
        if not obsid or not region_dir:
            print('This method requires you specify the obsid, and the region directory.')
            print('(set obsid and region_dir)')
            return
            
        fpm='A' #WHEN DOING MANUAL INPUT, REGIONS ARE THE SAME!
        if regionmethod=='input':
            regionfileA = glob.glob(working_dir+'gauss_cen_'+obsid+'_'+fpm+'_user_input*.reg')
            regionfileB = glob.glob(working_dir+'gauss_cen_'+obsid+'_'+fpm+'_user_input*.reg')
        if regionmethod=='double':
            regionfileA = glob.glob(working_dir+'gauss_cen_'+obsid+'_'+fpm+'_*.reg')
            regionfileB = glob.glob(working_dir+'gauss_cen_'+obsid+'_'+fpm+'_*.reg')
            
        specific_time_evt = glob.glob(region_dir+timestring+'/'+'*cl.evt')

    
    specific_time_evt.sort()      
        
        

    regiondictA = [regfile_to_regdict(rA, time) for rA in regionfileA]
    regiondictB = [regfile_to_regdict(rB, time) for rB in regionfileB]

    #print(regiondictA)
    #print(regionfileA)
    #regionsavename=working_dir+'/'+timestring+'/'+timestring
    if not savefigdir:
        if region_dir:
            savefigdir=region_dir+timestring
        else:
            savefigdir=working_dir+timestring


    #print(specific_time_evt)

    #if grade=='0':
    #    evtA=specific_time_evt[1]
    #    evtB=specific_time_evt[4]
        
    if grade=='0_4':
        evtA=specific_time_evt[0]
        evtB=specific_time_evt[2]
        #evtB=specific_time_evt[1]

    if grade=='21_24':
        evtA=specific_time_evt[1]
        evtB=specific_time_evt[3]



    if input_aia:
        m, nu_smap, aiareg, COMxy = ia.nuevtplot(evtA=evtA, evtB=evtB,
              savefigdir=savefigdir, AIA94=True, input_aia=input_aia,
              regiondictA=regiondictA, regiondictB=regiondictB,
             regionsave=False, #regionsavename=regionsavename, 
                             overlimb=True, nushift=nushift) 
    else:
        m, nu_smap, aiareg, COMxy = ia.nuevtplot(evtA=evtA, evtB=evtB,
              savefigdir=savefigdir, AIA94=True,
              regiondictA=regiondictA, regiondictB=regiondictB,
             regionsave=False, #regionsavename=regionsavename, 
                             overlimb=True, nushift=nushift)  


    if save_dict:
        dict = {'aiaregdicts': aiareg,
            'map': m,
            'nuCOM': COMxy,
            'nushift': nushift}
    
        if region_dir:
            file=region_dir+timestring+'/'+timestring+'_aia_region.pickle'
        else:
            file=working_dir+timestring+'/'+timestring+'_aia_region.pickle'
        with open(file, 'wb') as f:
            pickle.dump(dict, f, pickle.HIGHEST_PROTOCOL)

        return dict, file
    


    return aiareg, m, COMxy



def region_time_intervals(region_dirs, id_dirs, shush=True, list_=False):

    """
    Wrapper for find_all_intervals that deals with the case that TIS may have suceeded/failed for different
    sets of orbits for different regions (when using more than one region). 
    """
    #Find all time intervals for each region and make a big nested list
    all_all_time_intervals=[]
    all_all_time_intervals_list=[]
    starts=[]
    fixit=False
    for r in region_dirs:
        all_time_intervals, all_time_intervals_list = find_all_intervals(r, shush=shush)
        if len(all_time_intervals) != len(id_dirs):
            print('TIS failed on at least one orbit. Orbits completed: ', len(all_time_intervals))
            print('Orbits total: ', len(id_dirs))
            print('Region was: ', r)
            fixit=True
            
        starts_reg=[]
        for at in all_time_intervals:
            starts_reg.append(at[0][0])
        starts.append(starts_reg)
        all_all_time_intervals.append(all_time_intervals)
        all_all_time_intervals_list.append(all_time_intervals_list)
    
    #THE FOLLOWING is activated when TIS has failed for at least one orbit for at least one region (no saved 
    #intervals file).
    #It adds an empty string placeholder to the all_all_time_intervals array for that orbit(s) for that region.
    #It seems to work, but if your regions are looking very off the culprit may be incorrect performance here. 

    #print(starts)
    if fixit:
        ls = [len(s) for s in starts]
        longest = starts[np.argmax(ls)]
        #Looping over # orbits in region with the most
        for i in range(0, len(longest)):
            print('')
            #Looping over regions
            for j in range(0, len(starts)):
                #print('longest: ', longest[i])
                try:
                    #print('starts val: ', starts[j][i])
                    test=starts[j][i]
                except IndexError:
                    test=''
                if test != longest[i]:
                    #print('no,', test, longest[i])
                    all_all_time_intervals[j].insert(i, '')
                    starts[j].insert(i, '')

    #print(starts)
    #print(all_all_time_intervals)

    if list_:
        return all_all_time_intervals, fixit, all_all_time_intervals_list
        
    return all_all_time_intervals, fixit



def per_orbit_region_adjustment(working_dir, id_dirs, obsids, orbit_ind, aiamaps, nushift=[20,0],
                               method='input', shush=False, sep_axis=''):

    """
    Takes a given orbit (orbit_ind) for a given AR observation which has completed time interval selection.
    Plots AIA + NuSTAR data together, showing analysis regions. Allows adjustment of the needed shift for
    NuSTAR/AIA coalignment, and then saves an aia region pickle file for the first interval in the orbit. 

    Do this for every orbit, then run make_all_aia_dicts() to make aia region files for all time intervals 
    in all orbits (ready for upload to the NCCS, or wherever AIA data is being prepped based on region inputs. 
    """
    
    import copy

    if method=='input' or method =='double':
        #Find the top level directory for each region
        if method=='input':
            region_dirs = find_region_dirs(working_dir)
        if method=='double':
            region_dirs = find_direction_dirs(working_dir, sep_axis)
            #print(region_dirs)
        all_all_time_intervals, fixit = region_time_intervals(region_dirs, id_dirs, shush=shush)


        #Get a list of the first interval from this orbit for each region.
        #(with placeholders for missing TIS results, see above).
        try:
            first_intervals = [at[orbit_ind][0] for at in all_all_time_intervals]
            copyintervals = copy.deepcopy(first_intervals)
        except IndexError:
            first_intervals=[]
            copyintervals = []
            if fixit:
                intervals = [at[orbit_ind] for at in all_all_time_intervals]
                for intr in intervals:
                    if intr:
                        first_intervals.append(intr[0])
                        copyintervals.append(intr[0])
                    else:
                        copyintervals.append('')
                        
        #Find the longest of the first intervals across all regions. Use that time interval for plotting.      
        durations = [(fi[1]-fi[0]).to(u.s).value for fi in first_intervals]
        maxint = np.argmax(durations)
        time_interval = first_intervals[maxint]
        region_dir = region_dirs[maxint]



    if method=='fit':
        all_time_intervals, all_time_intervals_list = find_all_intervals(working_dir, shush=shush)
        time_interval = all_time_intervals[orbit_ind][0]
        #irrelevant in this case:
        region_dir=''
    
    obsid=obsids[orbit_ind]
    #Plot NuSTAR over AIA, and save an aia regions file in the corresponding region+time interval directory
    dict, file = nu_aia_coalign(time_interval, working_dir, nushift, save_dict=True, input_aia=aiamaps[orbit_ind],
                            regionmethod=method, obsid=obsid, region_dir=region_dir)

    if (method in ['input', 'double']) and (len(region_dirs) > 1):
        #For all the other regions, copy the generated aia region file into their first time interval directories, as we
        #want the same shift for all regions (and the files contain all regions. 
        import subprocess
        for i in range(0, len(copyintervals)):
            time=copyintervals[i]
            if time:
                r = region_dirs[i]
                timestring = time[0].strftime('%H-%M-%S')
                stopstring = time[1].strftime('%H-%M-%S')
                timestring=timestring+'_'+stopstring
                regcopy = r+timestring+'/'+timestring+'_aia_region.pickle'
                status = subprocess.call('cp '+file+' '+regcopy, shell=True)






def coalign_based_on_prior(time_intervals, working_dir, reference_interval, dorotation=True,
                          input_aias=[]):

    """
    Take a list of lead time intervals (sequential suborbits) from observations of the 
    same AR within a period of time where solar rotation is not significant.

    Using one prior NuSTAR center of mass + AIA co-alignment shift for an earlier
    time interval (saved by a sucessful run of nu_aia_coalign), produce AIA
    co-alignment shifts for all later time intervals based on the assumption of a 
    stationary active region. This is to say, we assume that the new co-alignment shifts
    can be found by subtracting the change in NuSTAR COM from the OLD shift.

    Note this also neglects NuSTAR COM changes that are due to changes in source morphology.
    It assumes all changes in NuSTAR COM are due to changes in pointing. 

    As may be clear, this is to be used with caution. Check ALL saved images.

    Set dorotation=True to assume the target moves west at 10 arcseconds/hr (solar rotation).
    
    """
    import pickle

    #print(input_aias)

    #Make a new working directory for images if it doesn't yet exist
    save_path = pathlib.Path(working_dir+'/coalign_images/')
    if not save_path.exists():
        save_path.mkdir()

    time = reference_interval
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring

    file=working_dir+timestring+'/'+timestring+'_aia_region.pickle'

    with open(file, 'rb') as f:
        data = pickle.load(f)

    refCOM=data['nuCOM']
    nushift=data['nushift']

    i=0
    for t in time_intervals:
        print(i)

        if dorotation:
            tdiff = (t[0]-reference_interval[0]).to(u.hr)
            rotation = (10*u.arcsec/u.hr)*tdiff
            #print('time diff: ', tdiff, ' rotation: ', rotation)
            #print('')
        else:
            rotation=0*u.arcsec


        if input_aias:
            #print(input_aias[i])
            nunuCOM = nu_aia_coalign(t, working_dir, nushift, justCOM=True, input_aia=input_aias[i])
        else:
            nunuCOM = nu_aia_coalign(t, working_dir, nushift, justCOM=True)
        #print('OG NuSTAR COM: ', refCOM)
        #print('New NuSTAR COM: ', nunuCOM)
        xchange = nunuCOM[0]-refCOM[0]
        ychange = nunuCOM[1]-refCOM[1]

        #print('difference: ', xchange, ychange)
        
        
        nunushift = [(nushift[0]*u.arcsec-xchange+rotation).value, (nushift[1]*u.arcsec-ychange).value]
        #print('OG shift: ', nushift*u.arcsec)
        #print('New shift: ', nunushift*u.arcsec)


        time = t
        timestring = time[0].strftime('%H-%M-%S')
        stopstring = time[1].strftime('%H-%M-%S')
        timestring=timestring+'_'+stopstring
        file=working_dir+timestring+'/'+timestring+'_aia_region.pickle'
        try:
            with open(file, 'rb') as f:
                data = pickle.load(f)
            #print('datamap')
            #print(data['map'])
            dict = nu_aia_coalign(t, working_dir, nunushift, save_dict=True, input_aia = data['map'],
                                     savefigdir=working_dir+'/coalign_images/')
        except FileNotFoundError: 
            if input_aias:
                dict = nu_aia_coalign(t, working_dir, nunushift, save_dict=True, 
                                     savefigdir=working_dir+'/coalign_images/', input_aia=input_aias[i])
            else:
                dict = nu_aia_coalign(t, working_dir, nunushift, save_dict=True, 
                                     savefigdir=working_dir+'/coalign_images/')

        #refCOM=nunuCOM
        #nushift=nunushift
        i+=1



def make_all_aia_dicts(all_time_intervals, working_dir, key, additional_path=''):
    """
    Make AIA region files for ALL time intervals, using lead time interval regions as produced by 
    coalign_based_on_prior(). 

    Put them in nustar OBSID-specific directories. 
    
    """
    import pickle
    import pathlib

    aia_dict_dir=working_dir+'all_aia_dicts_'+key+'/'
    #Make a new working directory for prepped data/etc if it doesn't yet exist
    save_path = pathlib.Path(aia_dict_dir)
    if not save_path.exists():
        save_path.mkdir()

    if additional_path:
        other_aia_dict_dir=additional_path+'all_aia_dicts_'+key+'/'
        #Make a new working directory for prepped data/etc if it doesn't yet exist
        save_path = pathlib.Path(other_aia_dict_dir)
        if not save_path.exists():
            save_path.mkdir()

    suborbit_directories = []
    for at in range(0, len(all_time_intervals)):
        #Get lead interval, for which we've saved a region.
        lead_interval = all_time_intervals[at][0]
        time=lead_interval
        timestring = time[0].strftime('%H-%M-%S')
        stopstring = time[1].strftime('%H-%M-%S')
        timestring=timestring+'_'+stopstring
        file=working_dir+timestring+'/'+timestring+'_aia_region.pickle'
        try:
            with open(file, 'rb') as f:
                data = pickle.load(f)  
        except FileNotFoundError: 
            print('Something is wrong, no prepared region file found for ', timestring)
            print('Exiting.')
            return

        try:
            aiareg = data['aiaregdict']
        except KeyError:
            aiareg = data['aiaregdicts']
        #print(aiareg)
    
        #orbit-specific directories
        reffile = glob.glob(working_dir+timestring+'/'+'*.evt')[0]
        obsid = reffile.split('/')[-1][2:13]
        suborbit_dir=aia_dict_dir+'orbit_'+obsid
        suborbit_directories.append(suborbit_dir)
        #Make a new working directory for prepped data/etc if it doesn't yet exist
        save_path = pathlib.Path(suborbit_dir)
        if not save_path.exists():
            save_path.mkdir()

        if additional_path:
            other_suborbit_dir=other_aia_dict_dir+'orbit_'+obsid
            #suborbit_directories.append(suborbit_dir)
            #Make a new working directory for prepped data/etc if it doesn't yet exist
            save_path = pathlib.Path(other_suborbit_dir)
            if not save_path.exists():
                save_path.mkdir()
            
    
        # #alternately, can make suborbit-specific directories. 
        # suborbit_dir=working_dir+'suborbit_'+timestring
        # suborbit_directories.append(suborbit_dir)
        # #Make a new working directory for prepped data/etc if it doesn't yet exist
        # save_path = pathlib.Path(suborbit_dir)
        # if not save_path.exists():
        #     save_path.mkdir()
            
        make_interval_dicts(all_time_intervals[at], aiareg, where=suborbit_dir)
        if additional_path:
            make_interval_dicts(all_time_intervals[at], aiareg, where=other_suborbit_dir)
        

    return suborbit_directories

        


def make_interval_dicts(time_intervals, regiondict, where='./'):

    """
    Takes a list of time intervals + makes a bunch of pickle files
    containing the time, region (for input into NCCS).
    """
    import copy
    import pickle

    for t in time_intervals:
        if not isinstance(regiondict, list):
            regiondicts = [copy.deepcopy(regiondict)]
        else:
            regiondicts = copy.deepcopy(regiondict)

        dict_ = {}

        num=0
        for i in range(0, len(regiondicts)):
            r = regiondicts[i]
            r['time_interval'] = t
            dict_['region'+str(num)] = r
            num+=1

        time = t
        timestring = time[0].strftime('%H-%M-%S')
        stopstring = time[1].strftime('%H-%M-%S')
        timestring=timestring+'_'+stopstring

        filename = where+'/'+timestring+'_aia_prep.pickle'

        #print(dict_)
        
        with open(filename, 'wb') as f:
            pickle.dump(dict_, f, pickle.HIGHEST_PROTOCOL)
    

def read_interval_dicts(time_interval, where='./', bltr=False, common_string='_aia_prep',
                       xrt_region_input=True):

    """
    Takes a time interval for which you have made a 
    pickle file w/ time, region + maybe NCCS AIA inputs,
    and reads it in. 
    
    """
    import pickle

    #print(time_interval)

    time = time_interval
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    #print(timestring)
    timestring=timestring+'_'+stopstring

    filename = where+'/'+timestring+common_string+'.pickle'

    with open(filename, 'rb') as f:
        data = pickle.load(f)

    #print(data)

    if 'centerx' in data.keys():
        
        if bltr:
            offset = [data['centerx'], data['centery']]
            rad = data['radius']*u.arcsec
            
            xx = offset[0].value
            yy = offset[1].value
            
            #Set broad box for plotting (using region object)
            bl=[(xx-600)*u.arcsec, (yy-600)*u.arcsec]
            tr=[(xx+800)*u.arcsec,(yy+800)*u.arcsec]
    
            if xrt_region_input:
                region_input = {'center': offset,
                      'radius': rad}
    
                return data, bl, tr, region_input
            else:
                return data, bl, tr
    
        #print(data.keys())
        return data

    else:
        if bltr:
            offsets_x, offsets_y = [], []
            xrt_region_inputs = []
            for k in data.keys():
                regdata = data[k]
                if regdata is None:
                    return
                offsets_x.append(regdata['centerx'].value)
                offsets_y.append(regdata['centery'].value)
                
                if xrt_region_input:
                    region_input = {'center': [regdata['centerx'], regdata['centery']],
                          'radius': regdata['radius']*u.arcsec}
                    xrt_region_inputs.append(region_input)
    
            xx = np.mean(offsets_x)
            yy = np.mean(offsets_y)
    
            #Set broad box for plotting (using region object)
            bl=[(xx-600)*u.arcsec, (yy-600)*u.arcsec]
            tr=[(xx+800)*u.arcsec,(yy+800)*u.arcsec]

            if xrt_region_input:
                return data, bl, tr, xrt_region_inputs
            
            else:
                return data, bl, tr
        else:
            return data












            
            






