from aiapy.calibrate import degradation
from aiapy.calibrate.util import get_correction_table, get_pointing_table
from aiapy.calibrate import register, update_pointing, degradation, estimate_error

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors

#https://docs.python.org/3/library/glob.html
import glob
import astropy.time
from astropy.coordinates import SkyCoord
from astropy import units as u
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
import regions


import hissw
import pathlib
import pickle
import copy
import scipy.io as io

#Jessie Code
import lightcurves as lc
import region_fitting as rf



"""

Preparing AIA Data for DEM: IDL/other helpers 
---------------------------------------------

        AIA
        ----
        
        NOTE: Self-contained (data is downloaded and prepped from scratch by functions defined in this file).
        
        aia_response_hissw_wrapper.pro - for getting AIA response (uses aia_get_response.pro from SSWIDL aia
                                         library). This is set up to run automatically via hissw if the 
                                         correctly-named file containing its output (the aia temperature responses)
                                         is not found in the current directory. 
                                         
                                         Note 1: we apply a time-dependent correction to the aia response using the
                                         degradation correction, so the response file only needs to be made once on
                                         a given system and can be used for any time interval.
                                         
        Note 2: getting aia uncertainties via aiapy.calibrate.estimate_error was failing on an issue downloading the
        error tables at the time this was written. There is a line here that hard-codes to a specific 
        existing error table file on Jessie's UMN machine; edit if that is not where you are using this code. 
        
  
"""

#AIA Error table - set path to location in your system.
errortab='/Users/jessieduncan/ssw/sdo/aia/response/aia_V3_error_table.txt'



def load_aia(time, bl, tr, plot=True, aia_exclude=[], aia_path='./', method='Middle',
            input_region=[], input_aia_region_dict=[], real_aia_err=False, aia_clobber=False,
             path_to_dodem='./',
            sunpy_dir='/Users/jessieduncan/sunpy/', 
             errortab='/Users/jessieduncan/ssw/sdo/aia/response/aia_V3_error_table.txt'):
    """
    -Looks for prepped (level 1.5) AIA files, calls aia_for_DEM to make them if not there. 
    -Extracts mean data values in region of choice.
    -Retrieves AIA degradation factors and corrects for degradation.
    -Returns DN/s/px for each channel
    
    Keywords
    ---------
	
    time - start time of interval for DEM (astropy Time)
         FORMAT LIKE, time=astropy.time.Time('2018-05-29T19:08:00', scale='utc')
         
    bl, tr - define rectangular region for DEM (bottom left, top right in arcsec)
        FORMAT LIKE, bl=[-200*u.arcsec,150*u.arcsec]
                      tr=[-50*u.arcsec,300*u.arcsec]
	
	plot - save AIA images for later reference (slow)	
	
	aia_exclude - wavelength channels to throw out (if not set, uses 94, 131, 171, 193, 211, 335)
	
	aia_path - location where we will place (or find) a directory specific to input time interval
	
	method
    	Middle - takes single AIA file per channel from near the center of the time interval
    	Average - averages results from all AIA files in a given channel during time interval    
    	
	input_region - type of region object used to select data. Currently supports: 
			'rectangle' (RectangleSkyObject)
			'circle' (CircleSkyObject)
			[] - if not set, uses rectangular region described by bl,tr
			
			if you want to add more, go edit map_to_dn_s_px()
			
	input_region_dict - dictionary of inputs needed to define chosen region object. See map_to_dn_s_px()
						for expected contents for each type of region.
						
						
	real_aia_err - set True to use aiapy.estimate_error + 10% in quadrature (else, return no error)
	
	aia_clobber - set True to start over prepping data from scratch (not use saved prepped map files)				
			
    	
    """
    
    clobber=aia_clobber
    
    if bool(bl) and bool(input_region):
        print("You provided both a region box and a specific region.")
        print("Bounding box will be used only for initial data prep (will save submap instead of full disk);")
        print("specific region will be used for DEM")
        print("")

        
    if bool(bl)==False and bool(input_region)==False:
        print("You need to supply either a region box (bl, tr) or a specific region.")
        return
    
    if bool(input_aia_region_dict)==False and bool(input_region)==True:
        print("To use this type of region (", input_region, "), you need to supply a dictionary of inputs.")
        return
    
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring
    #print(timestring)
    
    #Check for time-interval-specific directory, make one if we don't have one already
    save_path = pathlib.Path(aia_path) / timestring
    if not save_path.exists():
        save_path.mkdir()
         
    
    #====================================================================================== 
    if method=='Middle':
        
        #=========================
        ##Check if there are prepped files named like this method expects.
        #=========================
        
        ffp=sorted(glob.glob(aia_path+timestring+'/'+'maps_prep_*.fits'))
        #print(ffp)
        if len(ffp)!=6 or clobber==True:
            print('No prepped AIA data (or clobber==True), fetching some new AIA files at DEM time and converting to lev 1.5')
            aia_for_DEM(time, bl, tr, plot=plot, aia_path=aia_path, method='Middle', clobber=clobber, sunpy_dir=sunpy_dir)
            ffp=sorted(glob.glob(aia_path+timestring+'/'+'maps_prep_*.fits'))
            if len(ffp) > 6:
                print('More than six files found! Please resolve.')
                return

        aprep=sunpy.map.Map(ffp)
        
        #=========================

        #=========================
        #Get AIA degredation factors
        #=========================
        
        wavs = [m.meta['wavelnth'] for m in aprep]
        channels = wavs * u.angstrom
        deg_dict = get_degs(aia_path, timestring, channels, time)
                
        #=========================
        
        #=========================
        #Retrieve AIA data from each map
        #=========================

        aia_dn_s_px=[]
        aia_err_dn_s_px=[]
        for i in range(0,len(aprep)):
            m = aprep[i]
            wav = wavs[i]
            wstring=str(wav)
            if wav == 94:
                wstring=str(0)+str(wav)
            deg = deg_dict[wstring]
            if real_aia_err:
                aia_dn_s_px_, err = map_to_dn_s_px(m, deg, bl=bl, tr=tr, input_region=input_region, 
                                                  input_aia_region_dict=input_aia_region_dict, plot=plot,
                                                  timestring=timestring, aia_path=aia_path, 
                                                  real_aia_err=real_aia_err, errortab=errortab)
                aia_dn_s_px.append(aia_dn_s_px_)
                aia_err_dn_s_px.append(err)
            else:
                aia_dn_s_px.append(map_to_dn_s_px(m, deg, bl=bl, tr=tr, input_region=input_region, 
                                                  input_aia_region_dict=input_aia_region_dict, plot=plot,
                                                  timestring=timestring, aia_path=aia_path))
            
        #=========================
            
    #======================================================================================         

            
    #======================================================================================       
    if method=='Average':

        
        #=========================
        #Get AIA degredation factors
        #=========================
        
        waves=[94, 131, 171, 193, 211, 335]
        channels = waves * u.angstrom
        deg_dict = get_degs(aia_path, timestring, channels, time)
        
        #=========================
        
        #print('Before exclude:', waves)
        if bool(aia_exclude):
            main_list = list(set(waves) - set(aia_exclude))
            waves = sorted([waves[waves.index(x)] for x in main_list])
            #print('After:', waves)
            
            
        #=========================
        #Check if there are prepped files named like this method expects...
        #  ...if not, fetch them. Then make maps, get values, and average over the interval.
        #=========================
        
        aia_dn_s_px=[]
        aia_err_dn_s_px=[]
        #For each channel...
        for w in range(0,len(waves)):
            #Look for all prepped files...
            wstring=str(waves[w])
            if waves[w] == 94:
                wstring=str(0)+str(waves[w])
            ffp=sorted(glob.glob(aia_path+timestring+'/'+'map_t*_'+wstring+'.fits'))
            #If under two files, get more
            if len(ffp)<2 or clobber==True:
                print('Less than two files ready to average (or clobber set); we will go prep more.')
                aia_for_DEM(time, bl, tr, wav=[waves[w]], plot=plot, aia_path=aia_path, method='Average', 
                            clobber=clobber, sunpy_dir=sunpy_dir)
                ffp=sorted(glob.glob(aia_path+timestring+'/'+'map_t*_'+wstring+'.fits'))
                #If still less than two files, quit.
                if len(ffp)<2:
                    print('Still less than two files - quitting (something wrong).')
                    print('If your time interval is very short, perhaps use the Middle method option')
                    return

            #Maps of all files
            aprep=sunpy.map.Map(ffp)
            
            #degradation for this wavelength
            deg = deg_dict[wstring]

            
            #Get data from each map
            wav_dn_s_px=[]
            wav_err_dn_s_px=[]
            for i in range(0,len(aprep)):
                m = aprep[i]
                if m.exposure_time == 0*u.s:
                    print('File with no exposure time found - excluding.')
                    print('File was:', ffp[i])
                    continue
                if real_aia_err:
                    wav_dn_s_px_, err = map_to_dn_s_px(m, deg, bl=bl, tr=tr, input_region=input_region, 
                                                      input_aia_region_dict=input_aia_region_dict, plot=plot,
                                                      timestring=timestring, aia_path=aia_path, 
                                                      real_aia_err=real_aia_err, errortab=errortab)
                    wav_dn_s_px.append(wav_dn_s_px_)
                    wav_err_dn_s_px.append(err)
                else:
                    wav_dn_s_px.append(map_to_dn_s_px(m, deg, bl=bl, tr=tr, input_region=input_region,
                                                  input_aia_region_dict=input_aia_region_dict, plot=plot,
                                                  timestring=timestring, aia_path=aia_path))

                
             
            #Take the mean of all the files
            aia_dn_s_px.append(np.mean(wav_dn_s_px))
            if real_aia_err:
                aia_err_dn_s_px.append(np.mean(wav_err_dn_s_px))
            
         #=========================
    #====================================================================================== 
            
    aia_dn_s_px = np.array(aia_dn_s_px)
    
    # #  Load in the AIA responses from sswidl make_aiaresp_forpy.pro
    # Note, this involves a call to aia_get_response with keyword /evenorm set, which involves a 
    # normalization to agree with the SDO/EVE instrument.
    
    #Note this is NOT for a specific time interval, so it just goes in the working directory.
    try:
        aia_tresp=io.readsav('aia_tresp_en.dat')
    except FileNotFoundError: 
        print('No AIA response file found, so using HISSW to make one using SSWIDL aia_get_response.')
        ssw = hissw.Environment(ssw_packages=['sdo/aia', 'hessi'], ssw_paths=['aia', 'hessi'])
        agr_path = path_to_dodem+'/hissw_idl/aia_response_hissw_wrapper.pro'
        try:
            ssw_resp = ssw.run(agr_path)
            
            aia_tresp=io.readsav('aia_tresp_en.dat')
        except Exception:
            import traceback
            print(traceback.print_exc())
            print('Something is wrong with the SSWIDL run - make sure the following IDL script exists:')
            print(agr_path)
            print('')
            return

    #For each AIA channel, make a nicer formatted channel name:
    for i in np.arange(len(aia_tresp['channels'])):
        aia_tresp['channels'][i]=aia_tresp['channels'][i].decode("utf-8")
    chans=np.array(aia_tresp['channels']) 

    # Get the temperature response functions in the correct form for dem input
    #Define the list of temperatures for which we have AIA response values
    tresp_logt=np.array(aia_tresp['logt'])
    
    aia_tr = aia_tresp['tr']
    
    #If we are excluding any channels from the default list...
    if bool(aia_exclude):
        print('Excluding AIA: ', aia_exclude)
        main_list = list(set(waves) - set(aia_exclude))
        c = sorted([waves.index(x) for x in main_list])
        if method=='Middle':
            #print(c)
            aia_dn_s_px = aia_dn_s_px[c]
        chans = chans[c]
        aia_tr = aia_tr[c, :]
        
    if real_aia_err:
        print('Adding 10% error in quadrature with aiapy.estimate_error output.')
        tens = 0.1*np.copy(aia_dn_s_px)
        #print(aia_err_dn_s_px)
        newerr = [(aia_err_dn_s_px[i]**2+tens[i]**2)**(1/2) for i in range(0, len(aia_err_dn_s_px))]
        #print(newerr)
        
        return aia_dn_s_px, newerr, chans, aia_tr, tresp_logt
    else:
        return aia_dn_s_px, chans, aia_tr, tresp_logt
    
    

def get_degs(aia_path, timestring, channels, time):
    """
    Wrapper to look for saved aia degradation factors, or if not extract them using
    aiapy.calibrate. Returns list of degradation factors (by channel).
    """
 
    #print(channels)
    nc=len(channels)
   
    try:
        with open(aia_path+timestring+'/'+timestring+'degradation.pickle', 'rb') as f:
            deg_dict = pickle.load(f)
    except FileNotFoundError:    
 
        deg_dict = {'Time': time[0]}
        for i in np.arange(nc):
            wave = int(channels[i].value)
            wstring=str(wave)
            if wave == 94:
                wstring=str(0)+str(wave)
            #degradation is from aiapy.calibrate
            deg_dict[wstring]=degradation(channels[i],time[0])#,calibration_version=10)
 
        with open(aia_path+timestring+'/'+timestring+'degradation.pickle', 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(deg_dict, f, pickle.HIGHEST_PROTOCOL)
 
    return deg_dict


def aia_for_DEM(time, bl, tr, wav=[], plot=True, aia_path='./', method='Middle', clobber=False,
                sunpy_dir='/Users/jessieduncan/sunpy/'):
    """
    Finds and downloads an AIA image in each of six channels shortly
    after the chosen time. Doesn't return anything– just downloads files.
    
    Keywords
    ---------
    time - start time of interval for DEM (astropy Time)
         FORMAT LIKE, time=astropy.time.Time('2018-05-29T19:08:00', scale='utc')
         
    bl, tr - define rectangular region for DEM (bottom left, top right in arcsec)
        FORMAT LIKE, bl=[-200*u.arcsec,150*u.arcsec]
                      tr=[-50*u.arcsec,300*u.arcsec]
	
	wav - if you only one a specific wavelength channel (if not set, default is to fetch
			94, 131, 171, 193, 211, 335)
			
	method - depending on whether we are taking one file (Middle) or averaging (Average), 
			we will sample differently in time.
			
	aia_path - location of time interval directory (or where one will be placed)
	
	plot - set True to plot image maps for later reference
	
	clobber - set True to overwrite previously prepped files
	
	Returns
	--------
	Nothing (saves map files)
    
    """
    
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring
    #print(timestring)
    
    midway = time[0] + (time[1]-time[0]).to(u.s).value/2*u.s
    
    save_path = pathlib.Path(aia_path) / timestring
    if not save_path.exists():
        save_path.mkdir()

    waves=[94, 131, 171, 193, 211, 335]
    print(wav)
    if bool(wav):
        waves=wav
    in_dir = sunpy_dir
    fulldisk=True
    
    if method=='Middle':
        time_range = (midway-30*u.s, midway+30*u.s)
        one_of_each=[]
        sample_every=1*u.s
        
    if method=='Average':
        time_range = time
        sample_every=45*u.s
        files=1

    checker=0
    print('Waves:', waves)
    for w in waves:
        if method=='Middle':
            #Check for ANY files in the short time range selected around the midpoint of the full input time range
            files = lc.gather_aia_files(
                in_dir,
                time_range, 
                w,
                fulldisk,
                )
            #print(files)
            if files != []:
                one_of_each.append(in_dir+files[0])
        if files == [] or method=='Average':
            #If there are no files (Middle method) or if we're using the 'Average' method, look for all files in
            #the short time range (in Middle case) or full time range (in Average case).
            query = Fido.search(
                a.Instrument.aia,
                a.Physobs.intensity,
                a.Wavelength(w*u.angstrom),
                a.Time(time_range[0],time_range[1]),
                a.Sample(sample_every)
            )
            print(query)
            if method=='Average' or query.file_num == 1:
                print('LOOKING FOR + DOWNLOADING FILES')
                #If we're averaging (or if the query returns only one file), download them all.
                files = Fido.fetch(query, max_conn=1)
                count=1
                while files.errors != [] and count < 4:
                    files = Fido.fetch(query[0], max_conn=1)
                    count+=1
                if count == 4:
                    print('Failed to download this AIA file 4 times:')
                    print(query[0])
                    return
                if method=='Middle':
                    one_of_each.append(files[0])
            
            if method=='Middle' and query.file_num > 1:
                #If we're not averaging, and there's more than one file, just download the first one.
                files = Fido.fetch(query[0][0], max_conn=1)
                count=1
                while files.errors != [] and count < 4:
                    files = Fido.fetch(query[0][0], max_conn=1)
                    count+=1
                if count == 4:
                    print('Failed to download this AIA file 4 times:')
                    print(query[0][0])
                    return

                one_of_each.append(files[0])
                
            if method=='Average':
                #aiaprep all the maps for this wavelength, and save them for later!
                amaps=sunpy.map.Map(files)
                
                #Only load the pointing table once
                if checker == 0:
                    ptab = get_pointing_table(amaps[0].date - 12 * u.h, amaps[0].date + 12 * u.h)
                    checker+=1
                    
                # aiaprep the images, may take a while to run
                aprep=[]
                for m in amaps:
                    #Make saved versions of the maps for input into DEM code
                    wvn="{0:d}".format(1000+m.meta['wavelnth'])
                    wvn=wvn[1:]
                    filename=aia_path+timestring+'/'+'map_t'+m.date.strftime('%y-%m-%d_%H-%M-%S')+'_prep_'+wvn+'.fits'
                    checkfile = glob.glob(filename)
                    if bool(checkfile)==False or clobber==True:
                        print('PREPPING MAP')
                        mm = prep_this_map(m, bl, tr, ptab, filename, save=True)

    if method=='Middle':
        
        #sort into order based on file name (e.g. 94A will be last)
        ffa=sorted(one_of_each)
        amaps=sunpy.map.Map(ffa)

        # Get the wavelengths of the maps, get index of sort for this list of maps and reorder
        #(New order is increasing wavelength, ie 94A is first)
        wvn0 = [m.meta['wavelnth'] for m in amaps]
        srt_id = sorted(range(len(wvn0)), key=wvn0.__getitem__)
        amaps = [amaps[i] for i in srt_id]

        ptab = get_pointing_table(amaps[0].date - 12 * u.h, amaps[0].date + 12 * u.h)

        # aiaprep the images, may take a while to run
        for m in amaps:
            #Make saved versions of the maps for input into DEM code
            wvn="{0:d}".format(1000+m.meta['wavelnth'])
            wvn=wvn[1:]
            filename=aia_path+timestring+'/'+'maps_prep_'+wvn+'.fits'
            mm = prep_this_map(m, bl, tr, ptab, filename, save=True, plot=plot)



def prep_this_map(m, bl, tr, ptab, filename, save=True, plot=False):
    """
    Takes in a map and a pointing table. AIAPREPs the map, and saves a fits version 
    if save==True. Plots the prepped map if plot==True. Returns the prepped map.
    """
    
    #Update the pointing information for each map
    m_temp = update_pointing(m, pointing_table=ptab)
    #converts lev1 map to lev1.5 map 
    #(see: https://aiapy.readthedocs.io/en/stable/api/aiapy.calibrate.register.html?highlight=register)
    m = register(m_temp)
    #Make submap for each wavelength
    bottom_left = SkyCoord(bl[0]-100*u.arcsec, bl[1]-100*u.arcsec, frame=m.coordinate_frame)
    top_right = SkyCoord(tr[0]+100*u.arcsec,tr[1]+100*u.arcsec, frame=m.coordinate_frame)
    mm = m.submap(bottom_left=bottom_left, top_right=top_right)
    
    if plot:
        #plot them to see
        fig = plt.figure(figsize=(9, 7))
        mm.plot(cmap=cm.get_cmap('Spectral_r'))
        
        coords = SkyCoord(
            Tx=(bl[0], tr[0])*u.arcsec,
            Ty=(bl[1], tr[1])*u.arcsec,
            frame=mm.coordinate_frame,
        )


        mm.draw_quadrangle(
            coords,
            edgecolor="blue",
            linestyle="-",
            linewidth=2
        )
        
        plt.savefig(str(m.meta['wavelnth'])+'_aia_image.png')
    
    if save:
        #Make saved versions of the maps for input into DEM code
        wvn="{0:d}".format(1000+mm.meta['wavelnth'])
        wvn=wvn[1:]
        mm.save(filename, overwrite='True')
    
    return mm
    

def map_to_dn_s_px(m, deg, bl=[], tr=[], input_region=[], input_aia_region_dict=[], plot=False, 
                  timestring='', aia_path='', real_aia_err=False, 
                   errortab='/Users/jessieduncan/ssw/sdo/aia/response/aia_V3_error_table.txt'):
    """
    
    Inputs
    -------
    m - AIA sunpy map
    deg - degradation factor for m
    
    Keywords
    -------
    
    bl, tr OR input_region and input_aia_region_dict - define the region we care about
    timestring - also the name of the time interval directory
    aia_path - where time interval directory is located
	plot - set True to plot image map for later reference
	real_aia_err - set True to use aiapy.estimate_error + 10% in quadrature (else, return no error)
	errortab - point to AIA error table file (used if real_aia_err==True)  
	
	Outputs
	---------
	DN/px/s value
	
	OR (if real_aia_err==True):
	
	DN/px/s value, error
	
    """
    
    
    
    dur = m.meta['exptime']
    wav = m.meta['wavelnth']
    
    if bool(bl)==False and bool(input_region)==False:
        print('Need either bl, tr for submap, or input region!')
        return
    
    if bool(bl) and bool(input_region)==False:
        bottom_left = SkyCoord(bl[0],bl[1], frame=m.coordinate_frame)
        top_right = SkyCoord(tr[0],tr[1], frame=m.coordinate_frame)
        sub_temp = m.submap(bottom_left=bottom_left, top_right=top_right)
        data_mean = np.mean(sub_temp.data)
        num_pix = sub_temp.data.size
        
    if bool(input_region):
        subm = copy.deepcopy(m)
                
        if input_region=='rectangle':
            region_data=input_aia_region_dict
            region = regions.RectangleSkyRegion(
                SkyCoord(*region_data['center'], frame=subm.coordinate_frame ),
                width=region_data['width'],
                height=region_data['height'],
                angle=region_data['angle']
            )
            
        if input_region=='circle':
            region_data=input_aia_region_dict
            region = regions.CircleSkyRegion(
                SkyCoord(*region_data['center'], frame=subm.coordinate_frame ),
                region_data['radius']
            )            

        data = rf.get_region_data(subm, region, b_full_size=True)
        data_mean = np.mean(data[np.where(data > 0)])
        
        if plot:
            fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(projection=subm)
            subm.plot(axes=ax)
            (region.to_pixel(subm.wcs)).plot(ax=ax, color='red')
            norm = colors.PowerNorm(0.5, 0, 1e3) # Define a different normalization to make it easier to see
            plt.colorbar(norm=norm)
            wvn="{0:d}".format(1000+m.meta['wavelnth'])
            plt.savefig(aia_path+timestring+'/'+str(m.meta['wavelnth'])+'_'+wvn+'_input_region_aia_image.png')
    
    channel = wav * u.angstrom
    
    #  Correct the AIA data for the degradation
    cor_data=data_mean/deg
    
    if real_aia_err:
        
        
        #Alternate method (note: slow as hell) - use ssw methods via hissw
        
#         ssw = hissw.Environment(ssw_packages=['sdo/aia'], ssw_paths=['aia'])
#         agr_path = '/Users/jessieduncan/dems/aia_error_hissw_wrapper.pro'
#         inputs = {'channel': [wav], 'data': [cor_data]}
#         try:
#             ssw_resp = ssw.run(agr_path, args=inputs)
#         except Exception:
#             import traceback
#             print(traceback.print_exc())     
#         err = ssw_resp['err'][0]
        
        
        #cheat for quick if debugging non-aia
        #err=0.1*cor_data
        
        err = estimate_error(cor_data*(u.ct / u.pix), channel, error_table=errortab).value[0]
        
        #print('Degredation-corrected AIA DN/px: ', cor_data)
        #print('Per second: ', cor_data/dur)
        #print('Error: ', err)
        #print('Per second: ', err/dur)
        
        return cor_data/dur, err/dur 
    
    # Get into DN/s/px for the DEM stuff
    aia_dn_s_px=cor_data/dur
    
    #print(aia_dn_s_px)

    return aia_dn_s_px

    