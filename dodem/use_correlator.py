"""
Code related to searching for pointing shifts, etc, using Reed's correlator methods. 

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
    import nustar_dem_prep as nu


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
    import time_interval_selection as tis

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








