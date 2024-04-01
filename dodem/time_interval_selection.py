import numpy as np
import matplotlib.pyplot as plt
import pickle
import lightcurves as lc
import nustar_dem_prep as nu

import matplotlib.dates as mdates
import astropy.time
import datetime
import glob



def real_count_lightcurves(datapath, timerange, working_dir, erange):
    
    """
    For inspection: making lightcurves for full-detector, in various grade combinations.
    
    
    Keywords
    ---------
    
    datapath - path to NuSTAR OBSID directory for orbit of interest.
    
    timerange - time interval of interest within NuSTAR orbit
            FORMAT LIKE, 
                time=(astropy.time.Time('2018-05-29T19:08:00', scale='utc'), 
                        astropy.time.Time('2018-05-29T19:14:00', scale='utc')
    
    working_dir - where to place/look for time-interval directories
    
    erange - tuple of energy bounds (range of events of interest)
    
    
    """

    evtA = glob.glob(datapath+'/event_cl/*A06_cl.evt')
    evtB = glob.glob(datapath+'/event_cl/*B06_cl.evt')
    hkA  = glob.glob(datapath+'/hk/*A_fpm.hk')
    hkB  = glob.glob(datapath+'/hk/*B_fpm.hk')

    res = prepare_nustar_grade_lightcurves(evtA, evtB, hkA, hkB, timebin=5, erange=erange, 
                                          save_dir=working_dir, return_lightcurves=True)

    times_convertedA, countratesA, lvtA, countsA_, times_convertedB, countratesB, lvtB, countsB_ = res

    #For the 6-10 keV emission, we have returned lightcurve values


    totals=[]
    grades=['0', '0-4', '21-24']
    for g in range(0,len(grades)): 

        fig, ax1 = plt.subplots(1, 1, figsize=(15, 5))

        total_counts = np.array(countsA_)[g,:] + np.array(countsB_)[g,:]

        ax1.stairs(np.array(countsA_)[g,0:-1], times_convertedA, label='FPMA')
        ax1.stairs(np.array(countsB_)[g,0:-1], times_convertedA, label='FPMB')
        ax1.stairs(total_counts[0:-1], times_convertedA, label='Sum')
        
        rangeinds = np.where(np.logical_and(times_convertedA > timerange[0], times_convertedA < timerange[1]))[0]
        ax1.set_ylim(0, np.max(total_counts[rangeinds])*1.5)
        ax1.set_xlim(timerange[0], timerange[1])
                             
        ax1.legend()
        ax1.set_title('Real Counts - Grade '+grades[g])
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        ax1.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))

        totals.append(total_counts)


    corr14 = np.array(totals)[0,:] - 0.25*np.array(totals)[2,:]
    corr14[np.where(corr14 < 0.)] = 0
    corr54 = np.array(totals)[1,:] - 1.25*np.array(totals)[2,:]
    corr54[np.where(corr54 < 0.)] = 0

    fig, ax1 = plt.subplots(1, 1, figsize=(15, 5))

    #ax1.plot(times_convertedA, corr14)
    ax1.stairs(corr14[0:-1], times_convertedA, label='Grade 0 - 0.25*Grades 21-24')
    ax1.stairs(corr54[0:-1], times_convertedA, label='Grade 0-4 - 1.25*Grades 21-24')

    rangeinds = np.where(np.logical_and(times_convertedA > timerange[0], times_convertedA < timerange[1]))[0]
    ax1.set_ylim(0, np.max(corr54[rangeinds])*1.5)
    ax1.set_xlim(timerange[0], timerange[1])
    ax1.set_title('Real Counts - Pile-Up Corr')
    ax1.legend()
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax1.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))




def find_time_intervals_plus(datapath, timerange, working_dir, countmin=10, erange=[6.,10], lctype='corr54',
                            nofit=True, fast_min_factor=1, minimum_seconds=[]):    
    """
    
    Component Methods:

    - find_interval_fast(): fast method: from .evt files, determine a time interval with sufficient counts 
    in the chosen GRADES and ENERGIES
        - this can be done with any grade condition/pile up correction scheme

    - check_interval_slow(): slow method: actually make NuSTAR spectral data products

    Workflow:

    - 1. Using FAST method, find first time interval with sufficient counts (target: actual desired 
    minimum X a factor to correct for some events occuring outside the chosen region).  
    
    - 2. Using SLOW method, check PROPOSED TIME INTERVAL
        - if it has sufficient counts, update start time to end of new interval, and start again at 1. 
        - if not, repeat FAST method with the same start point, but double the correction factor (so, if 
        originally we had factor=1.5 (1.5x target counts required for fast method), now we require 3x the 
        target counts.
            - return to 1.
            - if it still fails, quit (this is not expected, may indicate the assumptions made here are 
            not good – mutliple sources in FOV?)
            
    - 3. If the FAST method fails due to not enough counts remaining (AKA we've reached the end of the 
    orbit without getting 10 counts in the final interval), we concatenate whatever time is remaining with 
    the prior interval.


    
    Keywords
    ---------
    
    datapath - path to NuSTAR OBSID directory for orbit of interest.
    
    timerange - time interval of interest within NuSTAR orbit
            FORMAT LIKE, 
                time=(astropy.time.Time('2018-05-29T19:08:00', scale='utc'), 
                        astropy.time.Time('2018-05-29T19:14:00', scale='utc')
                        
    working_dir - where to place/look for time-interval directories
    
    erange - tuple of energy bounds (range of events of interest)
    
    countmin - desired minimum real counts in each interval
                    
    fast_min_factor - multiplied by the countmin to get the target counts as found via the faster method.
                        This is intended to be > 1, and account for emission falling outside the chosen 
                        region. 
    
    lctype - what grades/pile-up correction to use? Options:
    
            'grade0' – return grade 0 lightcurve (FPMA,B sum)
            'grade04' - return grade 0-4 lightcurve (FPMA,B sum)
            'corr14' - return grade 0 - (1/4)*grades 21-24 (FPMA, B sum)
            'corr54' - return grade 0-4 - (5/4)*grades 21-24 (FPMA, B sum)
    
    
    """
    
    
    #Things will break if your datapath does not end in the standard nustar obsid data directory
    obsid = datapath.split('/')[-2]
        
    #Find evt and hk files for the NuSTAR orbit
    #Note - if you have made more .evt files in the event_cl directory, this may break. Intended to work with
    #untouched post-pipeline directories. 
    evtA = glob.glob(datapath+'/event_cl/*A06_cl.evt')
    evtB = glob.glob(datapath+'/event_cl/*B06_cl.evt')
    hkA  = glob.glob(datapath+'/hk/*A_fpm.hk')
    hkB  = glob.glob(datapath+'/hk/*B_fpm.hk')
    
    #Make lightcurves for the orbit, in energy range, in different grades (0, 0-4, 21-24).
    res = prepare_nustar_grade_lightcurves(evtA, evtB, hkA, hkB, timebin=5, erange=erange, 
                              save_dir=working_dir, return_lightcurves=True)
    
    #Make array of counts for each grade, fpmA,B summed.
    times_convertedA, count_lc = make_count_lightcurve(res, plot_lcs=False, lctype=lctype)
    
    #Limit to time range of interest (trim time and count arrays)
    interval = np.where(np.logical_and(times_convertedA >= timerange[0], times_convertedA <= timerange[1]))
    intervaltimes = times_convertedA[interval]
    intervalcounts = count_lc[interval]
    if minimum_seconds:
        timestep=(intervaltimes[1]-intervaltimes[0]).total_seconds()
        minimum_steps = int(minimum_seconds/timestep)
        if minimum_steps < 3:
            print('Requiring at least 3 time bins...')
            minimum_steps = 3
        #print(minimum_steps, type(minimum_steps))
    
    
    #BEGINNING OF BIG LOOP
    
    
    new_intervals = []

    start_here=0
    stop_yet=False
    
    times_failed=0
    failed_intervals=[]
    
    #Defining 
    og_fast_min_factor=fast_min_factor

    while stop_yet == False:
        print('')
        #Make proposed interval
        res_ = find_interval_fast(intervalcounts, start_here, countmin*fast_min_factor)
        
        if not res_:
            #Indicates that we've reached the end of the interval of interest.
            print('Let us combine this last bit with the prior interval')
            proposed_interval = new_intervals[-1]
            proposed_interval[1] = astropy.time.Time(intervaltimes[-1])
            stop_yet=True
            continue
            
        else:
            int_counts, startdex, endex = res_
            print('Fast Method Counts: ', int_counts)
            #Make proposed interval longer by an amount equivalent to the assumption that (region_factor)% of
            #events will be in the region.
            #adjust_endex = startdex + round((endex-startdex)/region_factor)
            
            if minimum_seconds:
                dur_s = (intervaltimes[endex]-intervaltimes[startdex]).total_seconds()
                if dur_s < minimum_seconds:
                    print('Time interval shorter than chosen minimum of ', minimum_seconds, ' seconds.')
                    print('Extending to a ', minimum_seconds, ' second-long interval.')
                    endex = startdex+minimum_steps
                    if endex > (len(intervaltimes)-1):
                        print('Remainder of orbit < ', minimum_seconds,' seconds from current start.')
                        print('Combining this last bit with prior interval.')
                        proposed_interval = new_intervals[-1]
                        proposed_interval[1] = astropy.time.Time(intervaltimes[-1])
                        stop_yet=True
                        continue
                        
                             
            
            proposed_interval = astropy.time.Time([intervaltimes[startdex], intervaltimes[endex]], scale='utc')
            
        print(proposed_interval)

        check = check_interval_slow(proposed_interval, erange, datapath, obsid, working_dir, 
                                    nofit=nofit,lctype=lctype, countmin=countmin)

        if check[1]:
            print('Found Time Interval', proposed_interval[0].strftime('%H-%M-%S'), 
                  proposed_interval[1].strftime('%H-%M-%S'))
            print('Counts: ', check[0])
            #Append interval to our list, set a new start index for the next interval, and reset the fast method 
            #factor to the original factor (in case it's been changed).    
            new_intervals.append(proposed_interval)
            if endex == (len(intervaltimes)-1):
                stop_yet=True
                continue
                
            start_here = endex
            fast_min_factor=og_fast_min_factor

        else:
            print('Not Enough', proposed_interval[0].strftime('%H-%M-%S'), proposed_interval[1].strftime('%H-%M-%S'))
            print('Counts: ', check[0])
            #Make a note that we didn't get enough counts in the interval we tried.
            times_failed+=1
            failed_intervals.append(proposed_interval)
            #If it's the first time we've failed with this starting index, double the factor used for the fast
            #method count target (i.e. if we were looking for 1.5x our target counts, now looking for 3x our 
            #target counts with the initial fast method). 
            if fast_min_factor==og_fast_min_factor:
                print('Starting over with requirement for twice the counts in fast interval')
                fast_min_factor=og_fast_min_factor*2
            else:
                print('It STILL did not work - weird! Quitting.')
                return
            
    print('Finishing with ', len(new_intervals), ' new intervals, and ', times_failed, ' failed intervals.' )
    print('Failure %: ', times_failed/len(new_intervals))
    
    data = {'time_intervals': new_intervals}
    
    
    timestring = timerange[0].strftime('%H-%M-%S')
    stopstring = timerange[1].strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring
    
    filename = timestring+'_'+lctype+'_'+str(erange[0])+'-'+str(erange[1])+'keV_min'+str(countmin)+'time_intervals.pickle'
    
    with open(filename, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL) 
    
    return new_intervals, failed_intervals
    
def check_interval_slow(time_int, erange, datapath, obsid, nustar_path, nofit=True,
                       lctype='corr14', countmin=10):
    
    if lctype=='grade0':
        pile_up_corr=False
        adjacent_grades=False
    
    if lctype=='grade04':
        pile_up_corr=False
        adjacent_grades=True
    
    
    if lctype=='corr14':
        pile_up_corr=True
        adjacent_grades=False
    
    if lctype=='corr54':
        pile_up_corr=True
        adjacent_grades=True
    
    regfile='starter_region.reg' 
    gtifile=datapath+'event_cl/nu'+obsid+'A06_gti.fits'
    res = nu.combine_fpm(time_int, [erange], nustar_path, make_nustar=True,
                         pile_up_corr=pile_up_corr, adjacent_grades=adjacent_grades,
                                              gtifile=gtifile, datapath=datapath, regfile=regfile, 
                                                edit_regfile=True, actual_total_counts=True,
                                                nofit=nofit, clobber=False, countmin=countmin)
   
    return res

    
    
def find_interval_fast(counts, startindex, countmin):
    
    """
    Find indices (startindex, end) for interval with countmin counts. 
    """
    
    t=0+startindex
    int_counts=0
    while int_counts < countmin:
        try:
            int_counts+=counts[t]
        except IndexError:
            print('We have reached the end of the full time range, with still only ',int_counts,' counts in this interval')
            return []
        t+=1    
        
    return int_counts, startindex, t
    

def make_count_lightcurve(res, plot_lcs=False, lctype='corr14'):
    """
    Takes as input the output of prepare_nustar_grade_lightcurves (which makes lightcurves for both FPM, 
    all grades) and (based on keywords) returns a summed/pile-up-corrected lightcurve.
    
    Options for type of lightcurve (lctype):
    
    'grade0' – return grade 0 lightcurve (FPMA,B sum)
    'grade04' - return grade 0-4 lightcurve (FPMA,B sum)
    'corr14' - return grade 0 - (1/4)*grades 21-24 (FPMA, B sum)
    'corr54' - return grade 0-4 - (5/4)*grades 21-24 (FPMA, B sum)
    
    """
    times_convertedA, countratesA, lvtA, countsA_, times_convertedB, countratesB, lvtB, countsB_ = res
    
    totals=[]
    grades=['0', '0-4', '21-24']
    for g in range(0,len(grades)): 
        total_counts = np.array(countsA_)[g,:] + np.array(countsB_)[g,:]
        totals.append(total_counts)

        if plot_lcs:
            fig, ax1 = plt.subplots(1, 1, figsize=(15, 5))

            ax1.stairs(np.array(countsA_)[g,0:-1], times_convertedA, label='FPMA')
            ax1.stairs(np.array(countsB_)[g,0:-1], times_convertedA, label='FPMB')
            ax1.stairs(total_counts[0:-1], times_convertedA, label='Sum')
            ax1.set_xlim(timerange[0], timerange[1])
            ax1.legend()
            ax1.set_title('Normalized Lightcurves - Grade '+grades[g])
            ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            ax1.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))

    if lctype=='grade0':
        return times_convertedA, np.array(totals)[0,:]
    
    if lctype=='grade04':
        return times_convertedA, np.array(totals)[1,:]
    
    if lctype=='corr14':
        corr14 = np.array(totals)[0,:] - 0.25*np.array(totals)[2,:]
        #corr14[np.where(corr14 < 0.)] = 0
        return times_convertedA, corr14
    
    if lctype=='corr54':
        corr54 = np.array(totals)[1,:] - 1.25*np.array(totals)[2,:]
        #corr54[np.where(corr54 < 0.)] = 0
        return times_convertedA, corr54
        
        

def prepare_nustar_grade_lightcurves(evtA, evtB, hkA, hkB, timebin=10, erange=[2.,10.], livetime_corr=True, 
                               return_lightcurves=False, save_dir='./'):
    """
    Returns FPMA + B lightcurves. Wrapper for get_a_nustar_lightcurve() which does just one.
    
    Using some stuff from this example, but customizing for my use:
    https://github.com/ianan/nustar_sac/blob/master/python/example_nustar_lightcurve.ipynb
    
    Default behavior is to just save lightcurves to a file - set return_lightcurves=True to spit them all out
    directly as well. 
    """
    #Load in the evt file (has the list of photons)
    evtdataA, hdrA = lc.load_nufiles(evtA[0])
    # Load in the hk file (has the livetime info)
    lvdataA, lvhdrA = lc.load_nufiles(hkA[0])
    evtdataB, hdrB = lc.load_nufiles(evtB[0])
    lvdataB, lvhdrB = lc.load_nufiles(hkB[0])
    
    kevA = evtdataA['PI']*0.04+1.6
    erange_evtdataA = evtdataA[np.where(np.logical_and(kevA > erange[0],kevA < erange[1]))]
    kevB = evtdataB['PI']*0.04+1.6
    erange_evtdataB = evtdataB[np.where(np.logical_and(kevB > erange[0],kevB < erange[1]))]
    
    grades=['0', '0-4', '21-24']
    #grade=-1 is non-sensical, included to make same where logic work for grade 0 only case
    gradebounds = [[-1,0], [-1,4], [20,24]]
    
    countratesA=[]
    countratesB=[]  
    countsA_=[]
    countsB_=[] 
    
    for g in range(0,len(grades)): 
        gb=gradebounds[g]        
        grad = erange_evtdataA['GRADE']
        grad_erange_evtdataA = erange_evtdataA[np.where(np.logical_and(grad > gb[0], grad <= gb[1]))]
        grad = erange_evtdataB['GRADE']
        grad_erange_evtdataB = erange_evtdataB[np.where(np.logical_and(grad > gb[0], grad <= gb[1]))]    
    
        times_convertedA, countrateA, lvtA, countsA = lc.get_a_nustar_lightcurve(grad_erange_evtdataA, hdrA, lvdataA, lvhdrA, 
                                                                 timebin=timebin, livetime_corr=livetime_corr)
        times_convertedB, countrateB, lvtB, countsB = lc.get_a_nustar_lightcurve(grad_erange_evtdataB, hdrB, lvdataB, lvhdrB, 
                                                                 timebin=timebin, livetime_corr=livetime_corr)
        
    
        data = {'Livetime-Corrected?': livetime_corr,
                'Time Bin (s)': timebin,
                'Energy Range': erange,
                'file paths': [evtA, evtB, hkA, hkB],
                'FPMA_countrate': countrateA,
                'FPMB_countrate': countrateB,
                'FPMA_counts': countsA,
                'FPMB_counts': countsB,
                'FPMA_times': times_convertedA,
                'FPMB_times': times_convertedB,
                'FPMA_livetime': lvtA,
                'FPMB_livetime': lvtB,
                'Grade Expression': grades[g]
                }
        


        with open(save_dir+'NuSTAR_lightcurve_'+grades[g]+'_'+str(erange[0])+'_to_'+str(erange[1])+'_keV.pickle', 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
                
                
        countratesA.append(countrateA)
        countratesB.append(countrateB)
        countsA_.append(countsA)
        countsB_.append(countsB)


    return times_convertedA, countratesA, lvtA, countsA_, times_convertedB, countratesB, lvtB, countsB_


def get_saved_intervals(timerange, lctype='grade0', basedir='./', countmin=10, erange=[6.,10.], custom_file=[]):
    """
    Little wrapper, reads in a file of the type made by find_intervals() above, 
    containing DEM time intervals.
    
    Keywords
    ---------
    
    timerange - broad time interval of interest within NuSTAR orbit
            FORMAT LIKE, 
                time=(astropy.time.Time('2018-05-29T19:08:00', scale='utc'), 
                        astropy.time.Time('2018-05-29T20:07:00', scale='utc')
    
    lctype - what grades/pile-up correction were used to make intervals? Options:
    
            'grade0' – return grade 0 lightcurve (FPMA,B sum)
            'grade04' - return grade 0-4 lightcurve (FPMA,B sum)
            'corr14' - return grade 0 - (1/4)*grades 21-24 (FPMA, B sum)
            'corr54' - return grade 0-4 - (5/4)*grades 21-24 (FPMA, B sum)

    basedir - path to where the file is located. 
    
    countmin - minimum real counts/interval used to make intervals
    
    erange - energy range used to make intervals
    
    custom_file - Set to a specific file to use (still within the basedir).
    
    """
    

    if custom_file:
        filename=custom_file
    else:
        timestring = timerange[0].strftime('%H-%M-%S')
        stopstring = timerange[1].strftime('%H-%M-%S')
        timestring=timestring+'_'+stopstring
        filename = timestring+'_'+lctype+'_'+str(erange[0])+'-'+str(erange[1])+'keV_min'+str(countmin)+'time_intervals.pickle'
        
    with open(filename, 'rb') as f:
            data = pickle.load(f)
            
    return data['time_intervals']
    
    







