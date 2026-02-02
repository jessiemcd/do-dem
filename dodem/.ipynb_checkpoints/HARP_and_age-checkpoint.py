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

from astropy import units as u



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
    for try_ in ['T_REC', 'USFLUX', 'ERRVF', 'MEANALP', 'MEANJZH', 'TOTUSJZ', 'TOTUSJH', 'LON_FWT', 'LAT_FWT', 'AREA']:
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
        print(kk)
        harps, orbits = get_harp_and_times(kk, all_targets)
        
        harplist = []
        for hh in harps:
            keys, segments = c.query('hmi.sharp_720s['+str(hh)+'][? (QUALITY<65536) ?]', 
                                     key='T_REC, USFLUX, ERRVF, MEANALP, MEANJZH, TOTUSJZ, TOTUSJH, LON_FWT, LAT_FWT, AREA', seg='Br')
            if keys.empty: 
                print('No result for HARP ', hh, ' key: ', kk)
                harplist.append([])
                continue
    
            t_rec = np.array([parse_tai_string(keys.T_REC[i],datetime=True) for i in range(keys.T_REC.size)])
            usflux = np.array(keys.USFLUX)
            lon_fwt = np.array(keys.LON_FWT)
    

            #HARP PARAMS documentation: http://jsoc.stanford.edu/doc/data/hmi/sharp/old/sharp.MB.htm
            #Also: https://solarb.mssl.ucl.ac.uk/JSPWiki/Wiki.jsp?page=JsocHarp
            harpdict = {'harptimes': t_rec,
                          'usflux_los': usflux, #Total unsigned flux in Maxwells
                          'usflux_rad': usflux/np.cos(np.deg2rad(lon_fwt)), 
                          'meanalp': np.array(keys.MEANALP), #Mean twist parameter, alpha, in 1/Mm
                          'meanjzh': np.array(keys.MEANJZH), #Mean current helicity in G^2/m
                          'totjz': np.array(keys.TOTUSJZ), #Total unsigned vertical current, in Amperes
                          'totjh': np.array(keys.TOTUSJH), #Total unsigned current helicity in G^2/m
                          'lon_fwt': lon_fwt, #Longitude (deg)
                          'lat_fwt': np.array(keys.LAT_FWT), #Lattitude (deg)
                           'harp_area': np.array(keys.AREA*3.044e16) #Area in cm^2; 
                                                                    #see: https://solarb.mssl.ucl.ac.uk/JSPWiki/Wiki.jsp?page=JsocHarp
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


# def save_usflux(res_files, harp, c, interp_seconds=20, interval_type='all', plot=True, label=''):
#     """

#     NOTE: WE HAVE NOW MOVED TO USING add_harp_params() INSTEAD (no interpolation). THIS WAS THE OLD
#     METHOD!!!!!!
    
#     For a given HARP patch + array of DEM result files, gets the HARP unsigned flux values, 
#     then opens each file, finds the average value during the relevant time interval (note:interpolated!)
#     and saves it as a new entry in the result dictionary, 'average_usflux'.

#     Interpolation: 
#     ----------------
#     HARP parameters are reported every 12 minutes. We do an interpolation of those values onto an array
#     with times spaced every (interp_seconds) seconds (20s default: less than typical minimum 30s DEM time
#     interval length, so ensures that there will be an enclosed value for any DEM time interval occuring during
#     the time the HARP patch is defined. 
    
#     """
#     hh=harp
#     rfs=res_files
#     #print(res_files)

#     files_with_usflx=[]
#     usflxs=[]
#     em10s=[]
#     em7s=[]
#     em5s=[]
#     powerlaws=[]
#     peaks=[]
#     all_times=[]

#     if interval_type != 'all':
#         flare_res = ana.get_saved_flares(add_stdv_flares=True)
#         early_starts = flare_res[0]
#         late_stops = flare_res[1]

#     keys, segments = c.query('hmi.sharp_cea_720s['+str(hh)+'][? (QUALITY<65536) ?]', 
#                                  key='T_REC, USFLUX, ERRVF, MEANALP, MEANJZH, TOTUSJZ, LON_FWT, LAT_FWT', seg='Br')
#     if keys.empty: 
#         print('No result for HARP ', hh, ' key: ', kk)
#         return
#     t_rec = np.array([parse_tai_string(keys.T_REC[i],datetime=True) for i in range(keys.T_REC.size)])

#     usflux = np.array(keys.USFLUX)
#     meanalp = np.array(keys.MEANALP)
#     meanjzh = np.array(keys.MEANJZH)
#     totjz = np.array(keys.TOTUSJZ)
#     lon_fwt = np.array(keys.LON_FWT)
#     lat_fwt = np.array(keys.LAT_FWT)
    
#     interp_usflux, new_times = make_interp_usflux(usflux, t_rec, seconds=interp_seconds, plot=False)
#     interp_meanalp, new_times = make_interp_usflux(meanalp, t_rec, seconds=interp_seconds, plot=False)
#     interp_meanjzh, new_times = make_interp_usflux(meanjzh, t_rec, seconds=interp_seconds, plot=False)
#     interp_totjz, new_times = make_interp_usflux(totjz, t_rec, seconds=interp_seconds, plot=False)
#     interp_lon_fwt, new_times = make_interp_usflux(lon_fwt, t_rec, seconds=interp_seconds, plot=False)
#     interp_lat_fwt, new_times = make_interp_usflux(lat_fwt, t_rec, seconds=interp_seconds, plot=False)

#     for ff in rfs:
#         #print(ff)
#         data, timestring, time = viz.load_DEM(ff)

#         if interval_type != 'all':
#             flare = ana.check_for_flare(time, early_starts, late_stops)
#             if flare:
#                 if interval_type != 'flare':
#                     continue

#             else:
#                 if interval_type == 'flare':
#                     continue
                
#         times_ = np.where(np.logical_and(new_times > time[0].datetime, new_times < time[1].datetime))[0]
#         if len(times_) > 0:
#             avgusflx=np.mean(interp_usflux[times_])
#             data['average_usflux'] = avgusflx

#             avgusflx_rad=np.mean(interp_usflux[times_]/np.cos(np.deg2rad(interp_lon_fwt))[times_])
#             data['average_usflux_rad'] = avgusflx_rad
            
#             avgmalp=np.mean(interp_meanalp[times_])
#             data['average_meanalp'] = avgmalp

#             avgmjzh=np.mean(interp_meanjzh[times_])
#             data['average_meanjzh'] = avgmjzh

#             avgtotjz=np.mean(interp_totjz[times_])
#             data['average_totjz'] = avgtotjz

#             with open(ff, 'wb') as f:
#                  # Pickle the 'data' dictionary using the highest protocol available.
#                  pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)  

#             files_with_usflx.append(ff)
#             usflxs.append(np.array([avgusflx, avgusflx_rad, avgmalp, avgmjzh, avgtotjz]))
#             em10s.append(data['above_10MK'][0])
#             em7s.append(data['above_7MK'][0])
#             em5s.append(data['above_5MK'][0])
#             powerlaws.append(np.array([data['powerlaws'][0][0],data['powerlaws'][1][0]]))
#             peaks.append(data['max_temp'])

#             all_times.extend(new_times[times_])

    


#     all_times.sort()
    
#     if plot and len(all_times) > 0:
#         marker_style = dict(linestyle='', markersize=8, fillstyle='full',color='darkgreen', markeredgecolor='darkgreen')
        
#         fig, ax = plt.subplots(figsize=(12,6))
#         the_brad = (keys.USFLUX)/np.cos(np.deg2rad(keys.LON_FWT))/1e22
#         the_blos = (keys.USFLUX)/1e22
#         ax.plot(t_rec, the_brad,'o',**marker_style, label='B_rad')
#         ax.plot(t_rec, the_blos, 'o', color='blue', label='B_los')
#         ax.set_ylabel('Maxwells x 1e22')
#         ax.set_xlabel('Time')
#         ax.axvspan(all_times[0], all_times[-1], alpha=0.7, color='Pink', label='NuSTAR Times')
#         ax.legend()
#         plt.savefig('./figures_etc/unsigned_flux_time_profile_with_nutimes_harp_'+str(harp)+'_'+label+'.png')
#         plt.close()


            

#     return_all=True
#     if return_all:
#         return files_with_usflx, usflxs, em10s, em7s, em5s, powerlaws, peaks
#     else:
#         return files_with_usflx




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


def get_all_orbit_avg(indslist, parray, lonarray, lonthresh=60):

    """
    Keywords
    ---------

    indlist: list of lists of indices (specifying times of interest)
    parayy: parameter array (to be indexed by the indices above, function returns mean and stdv of 
                input parameter at those indices).
    lonarray: array of corresponding longitude values
    lonthresh: mean, stdv will only include values from times where the longitude is within the chosen
                    symmetric bounds (e.g. for lonthreshold=60, values must be within +/- 60 from center).

                    
    """

    allvals = []
    alllons = []
    #alltimes = []
    for o in indslist:
        #print('indices: ', o)
        #print(parray[o])
        allvals.extend(parray[o])
        alllons.extend(lonarray[o])
        #alltimes.extend(timearray[o])

    alllons = np.array(alllons)

    #print('longitudes: ', alllons)
    #print('times: ', alltimes)

    goodloninds = np.where(np.logical_and(alllons < lonthresh, alllons > -1*lonthresh))[0]

    if len(goodloninds) > 0:
        allvals = np.array(allvals)
        vv = np.mean(allvals[goodloninds])
        sv = np.std(allvals[goodloninds])
        if len(allvals[goodloninds]) < 3:
            sv == 0.
    
        return np.array([vv, sv])
        
    else:
        return


lonthresh=30

def fetch_age_and_usflux(all_targets, lonthresh=30, plot=True,
                        return_only_goodlon=True, check_quiet=False):

    params_of_interest = []
    
    uflxr = [] #total radial unsigned flux in Maxwells
    uflxl = [] #total LOS unsigned flux in Maxwells
    tuvc = [] #total unsigned vertical current in Amperes
    malp = [] #mean twist parameter, alpha, in 1/Mm
    mheli = [] #mean current helicity in G^2/m
    tuheli = [] #total unsigned current helicity in G^2/m
    areas = [] #HARP area (cm^2)

    
    ages = [] #age in days
    regions = []
    
    for kk in all_targets.keys():
        #print(kk)
        #print(all_targets[kk])
        if kk != '07-may-21':
    
            #for all HARP patches observed by NuSTAR in the pointing/day associated with this key...
            for i in range(0, len(all_targets[kk]['HARP'])):
                #print(all_targets[kk]['HARP'][i])
                #If there are HARP parameters saved, as well as region ages.
                if all_targets[kk]['HARP params'][i] and all_targets[kk]['Region Ages']:
                    #If set, check if there are any quiet files:
                    if check_quiet:
                        files = all_targets[kk]['res_file_dict(s)'][i]['quiet files all-inst']
                        if len(files) == 0:
                            continue
    
                    #Returns: mean unsigned radial flux for this HARP over all times where values fall into NuSTAR orbits, 
                    #         mean +/- stdv of the same.          
                    uflxr_ = get_all_orbit_avg(all_targets[kk]['HARP params'][i]['per_orbit_intersect_indices'], 
                                      all_targets[kk]['HARP params'][i]['usflux_rad'], all_targets[kk]['HARP params'][i]['lon_fwt'],
                                               lonthresh=lonthresh)
                    if uflxr_ is not None:
                    
                        #Add the age + HARP to the lists. 
                        ages.append(np.array(all_targets[kk]['Region Ages'][i]))
                        regions.append(all_targets[kk]['sub_keys'][i])
                        uflxr.append(uflxr_)
                        
                        #Same thing but for the LOS unsigned flux.
                        uflxl.append(get_all_orbit_avg(all_targets[kk]['HARP params'][i]['per_orbit_intersect_indices'], 
                                      all_targets[kk]['HARP params'][i]['usflux_los'], all_targets[kk]['HARP params'][i]['lon_fwt'],
                                                  lonthresh=lonthresh))
                        
                        #Same thing but for the total unsigned vertical current
                        tuvc.append(get_all_orbit_avg(all_targets[kk]['HARP params'][i]['per_orbit_intersect_indices'], 
                                      all_targets[kk]['HARP params'][i]['totjz'], all_targets[kk]['HARP params'][i]['lon_fwt'],
                                                  lonthresh=lonthresh))

                        #Same thing but for the unsigned helicity.
                        tuheli.append(get_all_orbit_avg(all_targets[kk]['HARP params'][i]['per_orbit_intersect_indices'], 
                                      all_targets[kk]['HARP params'][i]['totjh'], all_targets[kk]['HARP params'][i]['lon_fwt'],
                                                  lonthresh=lonthresh))

                        #Same thing but for the mean current helicity.
                        mheli.append(get_all_orbit_avg(all_targets[kk]['HARP params'][i]['per_orbit_intersect_indices'], 
                                      all_targets[kk]['HARP params'][i]['meanjzh'], all_targets[kk]['HARP params'][i]['lon_fwt'],
                                                  lonthresh=lonthresh))  

                        #Same thing but for the mean twist parameter.
                        malp.append(get_all_orbit_avg(all_targets[kk]['HARP params'][i]['per_orbit_intersect_indices'], 
                                      all_targets[kk]['HARP params'][i]['meanalp'], all_targets[kk]['HARP params'][i]['lon_fwt'],
                                                  lonthresh=lonthresh))   

                        #Same thing but for the area.
                        areas.append(get_all_orbit_avg(all_targets[kk]['HARP params'][i]['per_orbit_intersect_indices'], 
                                      all_targets[kk]['HARP params'][i]['harp_area'], all_targets[kk]['HARP params'][i]['lon_fwt'],
                                                  lonthresh=lonthresh)) 
    
                else:
                    if not return_only_goodlon:
                        ages.append(np.array(all_targets[kk]['Region Ages'][i]))
                        uflxr.append(0.)
                        uflxl.append(0.)
                        tuvc.append(0.)
                        tuheli.append(0.)
                        mheli.append(0.)
                        malp.append(0.)
                        areas.append(0.)
                        regions.append(all_targets[kk]['sub_keys'][i])
                    else:
                        pass

    
    ages = np.array(ages)

    #For cases with stdv too ill-defined (< 3 total data points), make the uncertainty the mean of all other uncertainties.
    hparr = replace_bad_stdv(np.array(uflxr), shush=True)
    hparl = replace_bad_stdv(np.array(uflxl))
    # tuvc = replace_bad_stdv(np.array(tuvc))
    # tuheli = replace_bad_stdv(np.array(tuheli))
    # mheli = replace_bad_stdv(np.array(mheli))
    # malp = replace_bad_stdv(np.array(malp))

    data = {'usflx': {'name': 'Total Radial Unsigned Flux',
                     'vals': hparr,
                     'unit': u.Mx,
                     'norm': 1e22,
                     'regions': regions},
            'usflx_los': {'name': 'Total LOS Unsigned Flux',
                     'vals': hparl,
                     'unit': u.Mx,
                     'norm': 1e22,
                     'regions': regions},
            'tuvc': {'name': 'Total Unsigned Vertical Current',
                     'vals': replace_bad_stdv(np.array(tuvc)),
                     'unit': u.A,
                     'norm': 1e13,
                     'regions': regions},
            'tuheli': {'name': 'Total Unsigned Current Helicity',
                     'vals': replace_bad_stdv(np.array(tuheli)),
                     'unit': (u.Gauss)**2/u.m,
                     'norm': 1,
                     'regions': regions},
            'mheli': {'name': 'Mean Current Helicity',
                     'vals': replace_bad_stdv(np.array(mheli)),
                     'unit': (u.Gauss)**2/u.m,
                     'norm': 1,
                     'regions': regions},
            'malp': {'name': 'Mean Twist Parameter (alpha)',
                     'vals': replace_bad_stdv(np.array(malp)),
                     'unit': 1/u.Mm,
                     'norm': 1,
                     'regions': regions},
            'areas': {'name': 'HARP Area',
                     'vals': replace_bad_stdv(np.array(areas)),
                     'unit': u.cm**2,
                     'norm': 1e20,
                     'regions': regions},
            'ages': {'name': 'AR Age',
                     'vals': ages,
                     'unit': u.day,
                     'norm': 1,
                     'regions': regions},
            'regions': {'name': 'Region IDs',
                     'vals': regions,
                     'units': False,
                     'norm': 1,
                     'regions': regions}
           }
                     


    
    if plot:


        from scipy import odr
        #ODR documentation: https://docs.scipy.org/doc/scipy-1.16.0/reference/odr.html
            
        linear = odr.Model(flin)
        
        ageerr = ages[:,0]-ages[:,1]
        ww = [1/he for he in ageerr]
        herr = hparr[:,1]/1e22
        herl = hparl[:,1]/1e22

        # print(ageerr)
        # print(ages[:,0])
        # print(hparr[:,0]/1e22)
        # print(1/herr)
        
        fig = plt.figure(figsize=(8, 5), tight_layout = {'pad': 1})
        
        
        plt.errorbar(hparr[:,0]/1e22, ages[:,0], yerr=ageerr, xerr=herr, 
                     linestyle='none', color='red', label='Usigned Flux (Radial)')#, marker='o')
        
        m, b = np.polyfit(hparr[:,0]/1e22, ages[:,0], 1, w=ww)
        fity = np.array(hparr[:,0]/1e22)*m + b
        plt.plot(hparr[:,0]/1e22, fity, linestyle='dotted', label='linear fit, slope: '+str(round(m,3)), color='red')

        # mydata = odr.Data(hparr[:,0]/1e22, ages[:,0], wd=1/herr, we=ww)
        # myodr = odr.ODR(mydata, linear, beta0=[0.,1.]) # how to put beta0 values here
        # myoutput = myodr.run()
        # #myoutput.pprint()
        # m, b = myoutput.beta
        # fity = np.array(hparr[:,0]/1e22)*m + b
        # plt.plot(hparr[:,0]/1e22, fity, linestyle='dotted', label='TLS linear fit, slope: '+str(round(m,3)), color='red')
        
        
        plt.errorbar(hparl[:,0]/1e22, ages[:,0], yerr=ageerr, xerr=herl, 
                     linestyle='none', color='blue', label='Usigned Flux (LOS)')#, marker='o')
        
        m, b = np.polyfit(hparl[:,0]/1e22, ages[:,0], 1, w=ww)
        fity = np.array(hparl[:,0]/1e22)*m + b
        plt.plot(hparl[:,0]/1e22, fity, linestyle='dotted', label='linear fit, slope: '+str(round(m,3)), color='blue')

        # mydata = odr.Data(hparl[:,0]/1e22, ages[:,0], wd=1/herr, we=ww)
        # myodr = odr.ODR(mydata, linear, beta0=[0.,1.]) # how to put beta0 values here
        # myoutput = myodr.run()
        # #myoutput.pprint()
        # m, b = myoutput.beta
        # fity = np.array(hparl[:,0]/1e22)*m + b
        # plt.plot(hparl[:,0]/1e22, fity, linestyle='dotted', label='TLS linear fit, slope: '+str(round(m,3)), color='blue')
        

        plt.xlabel('Unsigned Flux (x1e22)')
        plt.ylabel('AR Age (Days)')
        plt.title('Age vs. Unsigned Flux, longitude threshold = '+str(lonthresh))
        plt.legend()


    return data

def replace_bad_stdv(arr, shush=True):

    jrstds = np.mean(arr[:,1][np.where(arr[:,1] != 0.)[0]])
    if not shush:
        print(jrstds)
        print(1/jrstds)
        print(arr[:,1]/1e22)
        print(1/arr[:,1]*1e22)
    arr[:,1][np.where(arr[:,1] == 0.)[0]] = jrstds

    return arr

def flin(B, x):
    #Linear function
    return B[0]*x + B[1]


def get_ages_vs_params(all_targets, param='above10', return_stdvs=False,
                      filetype='quiet files all-inst', skiphuh=False):
    
    ages_ = []
    params_ = []
    stdvs_ = []
    sks = []

    huh_list = ['12-sep-17', '03-may-21_1', '03-may-21_2', '13-sep-17']

    for kk in all_targets.keys():
        #print('')
        #print(kk)
        #print(len(all_targets[kk]['sub_keys']))
        #print(len(all_targets[kk]['Region Ages']))

        if kk != '07-may-21' and all_targets[kk]['Region Ages']:
            for sk in range(0, len(all_targets[kk]['sub_keys'])):
                ages = np.array(all_targets[kk]['Region Ages'][sk])
                files = all_targets[kk]['res_file_dict(s)'][sk][filetype]
                #print(kk, sk, len(files))
                if param=='dn_in':
                    params = np.zeros((len(files), 12))
                else:
                    params = np.zeros((len(files), 3))
    
                ind=0
                for f_ in files:
                    with open(f_, 'rb') as f:
                        data = pickle.load(f)
    
                    if param=='spexkT':
                        params[ind] = data['SPEX_dict']['kT_m_it'][0].value
                    if param=='above10':
                        #It needs to be this way to make the magnetic parameters plot nice
                        params[ind] = np.log10(data['above_10MK'])
                        #At some point it was changed to the below, presumably for a reason!
                        #params[ind] = data['above_10MK']
                    if param=='above7':
                        params[ind] = np.log10(data['above_7MK'])
                    if param=='above5':
                        params[ind] = np.log10(data['above_5MK'])
    
                    if param=='peak':
                        params[ind] = 10**data['max_temp']/1e6
    
                    if param=='upperpower':
                        params[ind] = list(data['powerlaws'][1])
    
                    if param=='lowerpower':
                        params[ind] = list(data['powerlaws'][0])
                        #if kk in huh_list:
                            #print(data['time_interval'][0])
                            #print(kk, sk, ' lp1s:', data['powerlaws'][0])

                    if param=='dn_in':
                        params[ind][0:len(data['dn_in'])] = data['dn_in']

                        
                    if param=='upperpower2':
                        params[ind] = list(data['powerlaws2'][1])
    
                    if param=='lowerpower2':
                        params[ind] = list(data['powerlaws2'][0])
                        #if kk in huh_list:
                            #print(data['time_interval'][0])
                            #print(kk, sk, ' lp2s:', data['powerlaws2'][0])

                    
                    ind+=1
                    
                if kk in huh_list and skiphuh:
                    continue

                if param=='dn_in':
                    mns = np.mean(params, axis=0)
                    if np.all(np.isfinite(mns)):
                        params_.append(mns)
                        ages_.append(ages)
                        sks.append(all_targets[kk]['sub_keys'][sk])

                else:      
                    stdv = np.std(params[:,0])
                    mns = np.mean(params[:,0])
    
                    if len(params[:,1]) < 3:
                        stdv = 0.
    
                    #propogating power law (or other built-in parameter) uncertainty
                    sqsum_err = (np.sum([pe**2 for pe in params[:,1]]))**(1/2)/len(params[:,1])
            
                    if np.isfinite(mns):
                        params_.append([mns, sqsum_err])
                        ages_.append(ages)
                        stdvs_.append(stdv)
                        #print(kk, sk, stdv)
                        sks.append(all_targets[kk]['sub_keys'][sk])



    
    params_ = np.array(params_)
    ages_ = np.array(ages_)
    stdvs_ = np.array(stdvs_)

    #print(stdvs_)
    jrstds = np.mean(stdvs_[np.where(stdvs_ != 0.)[0]])
    stdvs_[np.where(stdvs_ == 0.)[0]] = jrstds
    #print(stdvs_)

    #print('')

    if return_stdvs:
        return params_, ages_, stdvs_, sks
    else:
        return params_, ages_






    