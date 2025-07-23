import drms
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.dates import *
from datetime import datetime as dt_obj
import pickle
import csv

import nustar_dem_prep as nu
import nustar_utilities as nuutil
import all_nu_analysis as ana
import visualize_dem_results as viz




def add_AR_ages(all_targets, agecsv):

    agedict = parse_age_csv(agecsv)
    
    for kk in all_targets.keys():
        sks = all_targets[kk]['sub_keys']
        agds=[]
        for sk in sks:
            if sk in agedict.keys():
                #print(sk, agedict[sk])
                agds.append(agedict[sk])
    
        all_targets[kk]['Region Ages'] = agds


    return all_targets

def parse_age_csv(agecsv, oldway=False):
    agedict = {}

    with open(agecsv, 'r', newline='') as file:
        csv_reader = csv.reader(file)
    
        first=0
        for row in csv_reader:
            #print(row)
            if first > 0:
                #Lower bound of possible range, upper bound of possible range, date of secondary new flux emergence (if observed).
                valstr = [row[2], row[3], row[4]]
                valstr_ = [float(v) if v else None for v in valstr]
                if oldway:
                    agedict[row[0]] = valstr_
                    
                elif np.any(valstr_):
    
                    if valstr_[2] is not None:
                        #For cases with secondary flux emergence, that age is used (+/- 1 day)
                        agedict[row[0]]=[valstr_[2], valstr_[2]-1, valstr_[2]+1]
        
                    elif valstr_[0] == valstr_[1]:
                        #print('samesies:', valstr_)
                        agedict[row[0]]=[valstr_[0], valstr_[0]-1, valstr_[0]+1]
        
                    else:
                       #print('diffsies:', valstr_)
                       agedict[row[0]]=[valstr_[0]+(valstr_[1]-valstr_[0])/2., valstr_[0], valstr_[1]] 
                    
                    
            else:
                first+=1

    return agedict


def get_drms_client():

    #"To be able to access the JSOC DRMS from Python, we first need to import the drms module and create an instance of the Client class"
    #https://docs.sunpy.org/projects/drms/en/stable/tutorial.html
    c = drms.Client() 
    #selecting SHARP series
    c.series(r'hmi\.sharp_')
    si = c.info('hmi.sharp_720s')
    #Checking for the parameters we are interested in
    for try_ in ['USFLUX', 'ERRVF', 'MEANALP', 'MEANJZH', 'TOTUSJZ']:
        if try_ in si.keywords.index.tolist():
            print(try_)

    return c



def get_harp_and_times(key, dict):
    
    id_dirs = dict[key]['datapaths']
    harps = dict[key]['HARP']
    
    orbits=[]
    for id in id_dirs:
        evt_data, hdr, obsid = nu.return_submap(datapath=id, fpm='A', return_evt_hdr=True, return_obsid=True)
        time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
        timerange = [time0.datetime, time1.datetime]
        #print(timerange[0].strftime('%H-%M-%S'), timerange[1].strftime('%H-%M-%S'))
        orbits.append(timerange)
        
    return harps, orbits


def parse_tai_string(tstr,datetime=True):
    year   = int(tstr[:4])
    month  = int(tstr[5:7])
    day    = int(tstr[8:10])
    hour   = int(tstr[11:13])
    minute = int(tstr[14:16])
    if datetime: return dt_obj(year,month,day,hour,minute)
    else: return year,month,day,hour,minute



def add_harp_params(all_targets, shush=True):

    keys_to_check = all_targets.keys()

    c = get_drms_client()
    
    for kk in keys_to_check:
        harps, orbits = get_harp_and_times(kk, all_targets)
        
        harplist = []
        for hh in harps:
            keys, segments = c.query('hmi.sharp_720s['+str(hh)+'][? (QUALITY<65536) ?]', 
                                     key='T_REC, USFLUX, ERRVF, MEANALP, MEANJZH, TOTUSJZ, LON_FWT, LAT_FWT', seg='Br')
            if keys.empty: 
                print('No result for HARP ', hh, ' key: ', kk)
                harplist.append([])
                continue
    
            t_rec = np.array([parse_tai_string(keys.T_REC[i],datetime=True) for i in range(keys.T_REC.size)])
            usflux = np.array(keys.USFLUX)
            lon_fwt = np.array(keys.LON_FWT)
    
    
            harpdict = {'harptimes': t_rec,
                          'usflux_los': usflux,
                          'usflux_rad': usflux/np.cos(np.deg2rad(lon_fwt)),
                          'meanalp': np.array(keys.MEANALP),
                          'meanjzh': np.array(keys.MEANJZH),
                          'totjz': np.array(keys.TOTUSJZ),
                          'lon_fwt': lon_fwt,
                          'lat_fwt': np.array(keys.LAT_FWT)
                           }
    
            olist=[]
            one=0
            for o in orbits:
                if interval_intersect([t_rec[0], t_rec[-1]], o):
                    if not shush:
                        print(kk+': HARP '+str(hh)+' in orbit: ', o[0].strftime('%H-%M-%S'), o[1].strftime('%H-%M-%S'))
                    olist.append(np.where(np.logical_and(t_rec > o[0], t_rec < o[1]))[0])
                    one=1
                                  
                else:
                    olist.append([])
                    
                    
            if not one==1:
                harplist.append([])
            else:
                harpdict['per_orbit_intersect_indices'] = olist
                harplist.append(harpdict)
    
        all_targets[kk]['HARP params'] = harplist


    return all_targets




def to_float(d, epoch):
    return (d - epoch)/ np.timedelta64(1, 's')


def make_interp_usflux(usflux, t_rec, seconds=20, plot=False):
    
    ref_time = t_rec[0]
    new_times = np.arange(ref_time, t_rec[-1], np.timedelta64(seconds, 's'), dtype='datetime64')

    new_times_float = [to_float(d, np.datetime64(ref_time)) for d in new_times]
    t_rec_float = [to_float(np.datetime64(d), np.datetime64(ref_time)) for d in t_rec]

    interp_usflux = np.interp(new_times_float, t_rec_float, usflux)

    if plot:
        import matplotlib.dates as mdates
        #from matplotlib.ticker import NullFormatter, ScalarFormatter
        fig, ax = plt.subplots(figsize=(8,6))
        plt.plot(new_times, interp_usflux)
        plt.scatter(t_rec, usflux)
        plt.xticks(rotation=45)
        #ax.xaxis.set_major_formatter(mdates.DateFormatter('%D %H:%M:%S'))
        # format the x-axis with universal time
        locator = AutoDateLocator()
        locator.intervald[HOURLY] = [24] # only show every 3 hours
        formatter = DateFormatter('%D')
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)

    return interp_usflux, new_times


def get_usflux(res_files, interval_type='all'):

    rfs=res_files

    if interval_type != 'all':
        flare_res = ana.get_saved_flares(add_stdv_flares=True)
        early_starts = flare_res[0]
        late_stops = flare_res[1]

    usflxs=[]
    em10s=[]
    em7s=[]
    em5s=[]
    powerlaws=[]
    peaks=[]

    for ff in rfs:
        data, timestring, time = viz.load_DEM(ff)

        if interval_type != 'all':
            flare = ana.check_for_flare(time, early_starts, late_stops)
            if flare:
                if interval_type != 'flare':
                    continue

            else:
                if interval_type == 'flare':
                    continue

        try:
            avgusflx = data['average_usflux']
            avgusflx_rad = data['average_usflux_rad']
            avgmalp = data['average_meanalp']
            avgmjzh = data['average_meanjzh']
            avgtotjz = data['average_totjz']
            if np.log10(avgtotjz) < 11:
                print('low helicity: ', ff)
        except KeyError:
            #print(ff)
            continue
                
        usflxs.append(np.array([avgusflx, avgusflx_rad, avgmalp, avgmjzh, avgtotjz]))
        em10s.append(data['above_10MK'][0])
        em7s.append(data['above_7MK'][0])
        em5s.append(data['above_5MK'][0])
        powerlaws.append(np.array([data['powerlaws'][0][0],data['powerlaws'][1][0]]))
        peaks.append(data['max_temp'])

    return_all=True
    if return_all:
        return files_with_usflx, usflxs, em10s, em7s, em5s, powerlaws, peaks
    else:
        return files_with_usflx


def save_usflux(res_files, harp, c, interp_seconds=20, interval_type='all', plot=True, label=''):
    """
    For a given HARP patch + array of DEM result files, gets the HARP unsigned flux values, 
    then opens each file, finds the average value during the relevant time interval (note:interpolated!)
    and saves it as a new entry in the result dictionary, 'average_usflux'.

    Interpolation: 
    ----------------
    HARP parameters are reported every 12 minutes. We do an interpolation of those values onto an array
    with times spaced every (interp_seconds) seconds (20s default: less than typical minimum 30s DEM time
    interval length, so ensures that there will be an enclosed value for any DEM time interval occuring during
    the time the HARP patch is defined. 
    
    """
    hh=harp
    rfs=res_files
    #print(res_files)

    files_with_usflx=[]
    usflxs=[]
    em10s=[]
    em7s=[]
    em5s=[]
    powerlaws=[]
    peaks=[]
    all_times=[]

    if interval_type != 'all':
        flare_res = ana.get_saved_flares(add_stdv_flares=True)
        early_starts = flare_res[0]
        late_stops = flare_res[1]

    keys, segments = c.query('hmi.sharp_cea_720s['+str(hh)+'][? (QUALITY<65536) ?]', 
                                 key='T_REC, USFLUX, ERRVF, MEANALP, MEANJZH, TOTUSJZ, LON_FWT, LAT_FWT', seg='Br')
    if keys.empty: 
        print('No result for HARP ', hh, ' key: ', kk)
        return
    t_rec = np.array([parse_tai_string(keys.T_REC[i],datetime=True) for i in range(keys.T_REC.size)])

    usflux = np.array(keys.USFLUX)
    meanalp = np.array(keys.MEANALP)
    meanjzh = np.array(keys.MEANJZH)
    totjz = np.array(keys.TOTUSJZ)
    lon_fwt = np.array(keys.LON_FWT)
    lat_fwt = np.array(keys.LAT_FWT)
    
    interp_usflux, new_times = make_interp_usflux(usflux, t_rec, seconds=interp_seconds, plot=False)
    interp_meanalp, new_times = make_interp_usflux(meanalp, t_rec, seconds=interp_seconds, plot=False)
    interp_meanjzh, new_times = make_interp_usflux(meanjzh, t_rec, seconds=interp_seconds, plot=False)
    interp_totjz, new_times = make_interp_usflux(totjz, t_rec, seconds=interp_seconds, plot=False)
    interp_lon_fwt, new_times = make_interp_usflux(lon_fwt, t_rec, seconds=interp_seconds, plot=False)
    interp_lat_fwt, new_times = make_interp_usflux(lat_fwt, t_rec, seconds=interp_seconds, plot=False)

    for ff in rfs:
        #print(ff)
        data, timestring, time = viz.load_DEM(ff)

        if interval_type != 'all':
            flare = ana.check_for_flare(time, early_starts, late_stops)
            if flare:
                if interval_type != 'flare':
                    continue

            else:
                if interval_type == 'flare':
                    continue
                
        times_ = np.where(np.logical_and(new_times > time[0].datetime, new_times < time[1].datetime))[0]
        if len(times_) > 0:
            avgusflx=np.mean(interp_usflux[times_])
            data['average_usflux'] = avgusflx

            avgusflx_rad=np.mean(interp_usflux[times_]/np.cos(np.deg2rad(interp_lon_fwt))[times_])
            data['average_usflux_rad'] = avgusflx_rad
            
            avgmalp=np.mean(interp_meanalp[times_])
            data['average_meanalp'] = avgmalp

            avgmjzh=np.mean(interp_meanjzh[times_])
            data['average_meanjzh'] = avgmjzh

            avgtotjz=np.mean(interp_totjz[times_])
            data['average_totjz'] = avgtotjz

            with open(ff, 'wb') as f:
                 # Pickle the 'data' dictionary using the highest protocol available.
                 pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)  

            files_with_usflx.append(ff)
            usflxs.append(np.array([avgusflx, avgusflx_rad, avgmalp, avgmjzh, avgtotjz]))
            em10s.append(data['above_10MK'][0])
            em7s.append(data['above_7MK'][0])
            em5s.append(data['above_5MK'][0])
            powerlaws.append(np.array([data['powerlaws'][0][0],data['powerlaws'][1][0]]))
            peaks.append(data['max_temp'])

            all_times.extend(new_times[times_])

    


    all_times.sort()
    
    if plot and len(all_times) > 0:
        marker_style = dict(linestyle='', markersize=8, fillstyle='full',color='darkgreen', markeredgecolor='darkgreen')
        
        fig, ax = plt.subplots(figsize=(12,6))
        the_brad = (keys.USFLUX)/np.cos(np.deg2rad(keys.LON_FWT))/1e22
        the_blos = (keys.USFLUX)/1e22
        ax.plot(t_rec, the_brad,'o',**marker_style, label='B_rad')
        ax.plot(t_rec, the_blos, 'o', color='blue', label='B_los')
        ax.set_ylabel('Maxwells x 1e22')
        ax.set_xlabel('Time')
        ax.axvspan(all_times[0], all_times[-1], alpha=0.7, color='Pink', label='NuSTAR Times')
        ax.legend()
        plt.savefig('./figures_etc/unsigned_flux_time_profile_with_nutimes_harp_'+str(harp)+'_'+label+'.png')
        plt.close()


            

    return_all=True
    if return_all:
        return files_with_usflx, usflxs, em10s, em7s, em5s, powerlaws, peaks
    else:
        return files_with_usflx




def interval_intersect(int1, int2):
    """for two intervals
    each a tuple of datetimes
    determine if they intersect.
    """

    intersection=False
    #If int_1 starts before int2 and stops after the begining of int2
    b4 = int1[0] < int2[0]
    ea = int1[1] > int2[0]
    es = np.where(np.logical_and(b4, ea))

    if es[0].size > 0:
        intersection=True

    #If int1 starts during int2
    b4 = int1[0] > int2[0]
    ea = int1[0] < int2[1]
    es = np.where(np.logical_and(b4, ea))

    if es[0].size > 0:
        intersection=True

    return intersection
    