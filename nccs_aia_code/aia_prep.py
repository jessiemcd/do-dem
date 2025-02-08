#Some needed code
#Extracted from do-dem: https://github.com/jessiemcd/do-dem

errortab = '/home/jmdunca2/aia_V3_error_table.txt'
#sunpy_dir='/Users/jessieduncan/sunpy/'

from aiapy.calibrate.util import get_correction_table, get_pointing_table
from aiapy.calibrate import register, update_pointing, degradation, estimate_error
from aiapy.calibrate import util

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

import pathlib
import pickle
import copy
import scipy.io as io

#Do-DEM Code
#import lightcurves as lc
#import region_fitting as rf

def aia_for_DEM_NCCS(time, bl, tr, wav=[], plot=True, aia_path='./', method='Middle', clobber=False,
                data_dir='./'):
    """
    Finds an AIA image on the NCCS in each of the input channels shortly
    after the chosen time. Doesn't return anything.

    Calls 
    ------
    prep_this_map
    
    Keywords
    ---------
    time - start time of interval for DEM (astropy Time)
         FORMAT LIKE, time=astropy.time.Time('2018-05-29T19:08:00', scale='utc')
         
    bl, tr - define rectangular region for DEM (bottom left, top right in arcsec)
        FORMAT LIKE, bl=[-200*u.arcsec,150*u.arcsec]
                      tr=[-50*u.arcsec,300*u.arcsec]
	
	wav - if you only want a specific wavelength channel (if not set, default is to fetch
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
    import os
    from astropy.utils.exceptions import AstropyUserWarning
    import glob
    
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring
    #print(timestring)
    
    midway = time[0] + (time[1]-time[0]).to(u.s).value/2*u.s
    
    save_path = pathlib.Path(aia_path) / timestring
    if not save_path.exists():
        save_path.mkdir()

    waves=[94, 131, 171, 193, 211, 335]
    #print(wav)
    if bool(wav):
        waves=wav
    in_dir = data_dir
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
    #print('Waves:', waves)
    one_of_each=[]              
    for w in waves:
        #print(w)
        wave_dir = data_dir+str(w)#+'/'
        #Note this code is modified from:
        #https://github.com/masek014/aia_lightcurves/blob/main/file_io.py
        
        
        dir_files = [f for f in os.listdir(
            wave_dir)]# if os.path.isfile(os.path.join(wave_dir, f))]
        
        fits_dir = os.path.abspath(wave_dir)
        fits_paths = [os.path.join(fits_dir, f) for f in dir_files]

        times = [] # Used for sorting the file names
        paths = []
        for p in fits_paths:
            try:
                t = checkfile(
                    path=p,
                    time_range=time_range,
                    wavelength=w*u.angstrom
                )
                if t is not None:
                    times.append(t)
                    paths.append(p)
            except OSError as e:  # Catch empty or corrupted fits files and non-fits files
                print(f'OSError with file {p}: {e}')
            except AstropyUserWarning as e:
                print(f'AstropyUserWarning with file {p}: {e}')

        if method=='Middle':
            #print(paths)
            if paths != []:
                one_of_each.append(paths[0])
            else:
                print('No '+str(w)+'files found in the time range:')
                print(time_range)
                print('in dir', wave_dir)
                print('exiting.')
                return           
                
        if method=='Average':
            if paths == []:
                print('No '+str(w)+'files found in the time range:')
                print(time_range)
                print('in dir', wave_dir)
                print('exiting.')
                return
            else:
                amaps=sunpy.map.Map(paths)
                if checker == 0:
                    if len(paths) == 1:
                        ref_map = amaps
                    else:
                        ref_map=amaps[0]
                    ptab = get_pointing_table(ref_map.date - 12 * u.h, ref_map.date + 12 * u.h)
                    checker+=1
                aprep=[]
                if len(paths) == 1:
                    m=amaps
                    prep_map_wrap(m, aia_path, timestring, clobber, bl, tr, ptab)
                else:    
                    for m in amaps:
                        prep_map_wrap(m, aia_path, timestring, clobber, bl, tr, ptab)
                    # #Make saved versions of the maps for input into DEM code
                    # wvn="{0:d}".format(1000+m.meta['wavelnth'])
                    # wvn=wvn[1:]
                    # #print(wvn)
                    # filename=aia_path+timestring+ \
                    #     '/'+'map_t'+m.date.strftime('%y-%m-%d_%H-%M-%S')+'_prep_'+wvn+'.fits'
                    # check_prepfile = glob.glob(filename)
                    # if bool(check_prepfile)==False or clobber==True:
                    #     print('PREPPING MAP')
                    #     mm = prep_this_map(m, bl, tr, ptab, filename, save=True)
    
    if method=='Middle':
        #print(one_of_each)
        ffa=sorted(one_of_each)
        #print(ffa)
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


def prep_map_wrap(m, aia_path, timestring, clobber, bl, tr, ptab):
    
    #Make saved versions of the maps for input into DEM code
    wvn="{0:d}".format(1000+m.meta['wavelnth'])
    wvn=wvn[1:]
    #print(wvn)
    filename=aia_path+timestring+ \
        '/'+'map_t'+m.date.strftime('%y-%m-%d_%H-%M-%S')+'_prep_'+wvn+'.fits'
    check_prepfile = glob.glob(filename)
    if bool(check_prepfile)==False or clobber==True:
        print('PREPPING MAP')
        mm = prep_this_map(m, bl, tr, ptab, filename, save=True)



def checkfile(path, time_range, wavelength):

    from astropy.io import fits
    
    with fits.open(path, output_verify='warn') as hdu:
        hdr = hdu[1].header
        obs_time = astropy.time.Time(
            hdr['DATE-OBS'], scale='utc', format='isot')
        same_time = (obs_time >= time_range[0]) and (
            obs_time <= time_range[1])
        #sanity check
        same_wavelength = (
            wavelength == hdr['WAVELNTH'] * u.Unit(hdr['WAVEUNIT']))
            
        if same_time and same_wavelength:
            return obs_time


def prep_this_map(m, bl, tr, ptab, filename, save=True, plot=False):
    """
    Takes in a map and a pointing table. AIAPREPs the map, and saves a fits version 
    if save==True. Plots the prepped map if plot==True. Returns the prepped map.
    """
    
    #Update the pointing information for each map
    #m_temp = update_pointing(m, pointing_table=ptab)
    
    try:
        m_temp = update_pointing(m, pointing_table=ptab)
    except TypeError:
        m.meta.pop('crpix1')
        m.meta.pop('crpix2')
        print('CRPIX issue on ', m.date, m.wavelength)
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
        mm.plot() #cmap=cm.get_cmap('Spectral_r'))
        
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
                   errortab=errortab):
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

        data = get_region_data(subm, region, b_full_size=True)
        data_mean = np.mean(data[np.where(data > 0)])
        
        #plot=True
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


def get_degs(aia_path, timestring, channels, time):
    """
    Wrapper to look for saved aia degradation factors, or if not extract them using
    aiapy.calibrate. Returns list of degradation factors (by channel). 
    """
    from aiapy.calibrate import degradation
    import pickle
    import numpy as np

    print(channels)
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

            #new edited version for while JSOC is down
            #aia_res_table = '/home/jmdunca2/aia_V10_20201119_190000_response_table.txt'
            #print('Using degradation factors from ', aia_res_table, ', assuming JSOC is down.')
            #correction_table = util.get_correction_table(correction_table=aia_res_table)
            #deg_dict[wstring]=degradation(channels[i],time[0], correction_table=correction_table)

        with open(aia_path+timestring+'/'+timestring+'degradation.pickle', 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(deg_dict, f, pickle.HIGHEST_PROTOCOL)

    return deg_dict

def load_aia(time, bl, tr, plot=True, NCCS=False, aia_exclude=[], aia_path='./', method='Middle',
             one_avg=False,
            input_region=[], input_aia_region_dict=[], real_aia_err=False, aia_clobber=False,
             path_to_dodem='./', NCCS_aia_resp_path='./',
             NCCS_save_path='./saved_AIA_dem_inputs/',
            data_dir='./', 
             errortab='/Users/jessieduncan/ssw/sdo/aia/response/aia_V3_error_table.txt'):
    """
    -Looks for prepped (level 1.5) AIA files, calls aia_for_DEM to make them if not there. 
    -Extracts mean data values in region of choice.
    -Retrieves AIA degradation factors and corrects for degradation.
    -Returns DN/s/px for each channel

    NOTE: sunpy_dir (which by pointed to e.g.'/Users/jessieduncan/sunpy/', where Fido 
    search/download places data) has been replaced by data_dir, which more accurately 
    reflects that this can be any path. 
    
    Keywords
    ---------
	
    time - start time of interval for DEM (astropy Time)
         FORMAT LIKE, time=astropy.time.Time('2018-05-29T19:08:00', scale='utc')
         
    bl, tr - define rectangular region for DEM (bottom left, top right in arcsec)
        FORMAT LIKE, bl=[-200*u.arcsec,150*u.arcsec]
                      tr=[-50*u.arcsec,300*u.arcsec]
	
	plot - save AIA images for later reference (slow)	

    NCCS - set True to use an alternate version of aia_for_DEM designed to find files on the NCCS
            instead of Fido search/download (or find them in your sunpy data directory). To be used
            when running on NCCS/ADAPT/PRISM. 
	
	aia_exclude - wavelength channels to throw out (if not set, uses 94, 131, 171, 193, 211, 335)
	
	aia_path - location where we will place (or find) a directory specific to input time interval
	
	method
    	Middle - takes single AIA file per channel from near the center of the time interval
    	Average - averages results from all AIA files in a given channel during time interval    

    one_avg 
        Set True to proceed with preparing data products while using the "average" method, even if
        only one file was found in the time interval (i.e. just that file will be used, technically not
        an average for that case). 
        
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

    import glob
    import numpy as np
    import pickle
    import scipy.io as io
                 
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
            print('No prepped AIA data (or clobber==True), ')
            print('fetching some new AIA files at DEM time and converting to lev 1.5')
            if NCCS:
                aia_for_DEM_NCCS(time, bl, tr, plot=plot, aia_path=aia_path, 
                                    method='Middle', clobber=clobber, data_dir=data_dir)
            else:
                aia_for_DEM(time, bl, tr, plot=plot, aia_path=aia_path, method='Middle', 
                            clobber=clobber, data_dir=data_dir)
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
        #Get AIA degredation factors - 
        #just do for all 6 channels, even if excluding some later.
        #=========================
        
        waves=[94, 131, 171, 193, 211, 335]
        channels = waves * u.angstrom
        deg_dict = get_degs(aia_path, timestring, channels, time)
        
        #=========================
        
        if bool(aia_exclude):
            print('Before exclude:', waves)
            main_list = list(set(waves) - set(aia_exclude))
            waves = sorted([waves[waves.index(x)] for x in main_list])
            print('After exclude:', waves)
            
            
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
                if NCCS:
                    aia_for_DEM_NCCS(time, bl, tr, wav=[waves[w]], plot=plot, aia_path=aia_path, 
                                   method='Average', clobber=clobber, data_dir=data_dir)
                else:
                    aia_for_DEM(time, bl, tr, wav=[waves[w]], plot=plot, aia_path=aia_path, 
                                method='Average', clobber=clobber, sunpy_dir=sunpy_dir)
                ffp=sorted(glob.glob(aia_path+timestring+'/'+'map_t*_'+wstring+'.fits'))
                #If still less than two files, quit.
                if len(ffp)<2:
                    print('Still less than two files.')
                    print('If your time intervals are very short, perhaps use the Middle method option')
                    if one_avg and len(ffp) == 1:
                        print('Continuing with one file.')
                    else:
                        return

            #Maps of all files
            aprep=sunpy.map.Map(ffp)
            
            #degradation for this wavelength
            deg = deg_dict[wstring]

            
            #Get data from each map
            wav_dn_s_px=[]
            wav_err_dn_s_px=[]
            if len(ffp) == 1:
                m = aprep
                if m.exposure_time == 0*u.s:
                    print('File with no exposure time found - excluding.')
                    print('File was:', ffp[i])
                    continue
                res = map_to_dn_s_px(m, deg, bl=bl, tr=tr, input_region=input_region, 
                                                      input_aia_region_dict=input_aia_region_dict, plot=plot,
                                                      timestring=timestring, aia_path=aia_path, 
                                                      real_aia_err=real_aia_err, errortab=errortab)

                if real_aia_err:
                    wav_dn_s_px_, err = res
                    wav_dn_s_px.append(wav_dn_s_px_)
                    wav_err_dn_s_px.append(err)
                else:
                    wav_dn_s_px.append(res)

            else:
                for i in range(0,len(aprep)):
                    m = aprep[i]
                    if m.exposure_time == 0*u.s:
                        print('File with no exposure time found - excluding.')
                        print('File was:', ffp[i])
                        continue
                    res = map_to_dn_s_px(m, deg, bl=bl, tr=tr, input_region=input_region, 
                                                          input_aia_region_dict=input_aia_region_dict, plot=plot,
                                                          timestring=timestring, aia_path=aia_path, 
                                                          real_aia_err=real_aia_err, errortab=errortab)
                    if real_aia_err:
                        wav_dn_s_px_, err = res
                        wav_dn_s_px.append(wav_dn_s_px_)
                        wav_err_dn_s_px.append(err)
                    else:
                        wav_dn_s_px.append(res)

                
             
            #Take the mean of all the files
            aia_dn_s_px.append(np.mean(wav_dn_s_px))
            if real_aia_err:
                aia_err_dn_s_px.append(np.mean(wav_err_dn_s_px))
            
         #=========================
    #====================================================================================== 
    #Regardless of method, moving on to get responses
    #====================================================================================== 
            
    aia_dn_s_px = np.array(aia_dn_s_px)

    if NCCS:
        aia_resp_path = NCCS_aia_resp_path
    else:
        aia_resp_path = 'aia_tresp_en.dat'
    
    # #  Load in the AIA responses from sswidl make_aiaresp_forpy.pro
    # Note, this involves a call to aia_get_response with keyword /evenorm set, which involves a 
    # normalization to agree with the SDO/EVE instrument.
    
    #Note this is NOT for a specific time interval, so it just goes in the working directory.
    try:
        aia_tresp=io.readsav(aia_resp_path)
    except FileNotFoundError: 
        if NCCS:
            print('Check your aia response path – no file found!')
            return
            
        print('No AIA response file found, so using HISSW to make one using SSWIDL aia_get_response.')
        ssw = hissw.Environment(ssw_packages=['sdo/aia', 'hessi'], ssw_paths=['aia', 'hessi'])
        agr_path = path_to_dodem+'/hissw_idl/aia_response_hissw_wrapper.pro'
        try:
            ssw_resp = ssw.run(agr_path)
            
            aia_tresp=io.readsav(aia_resp_path)
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

    deminput_dict = {'aia_dn_s_px': aia_dn_s_px,
                     'chans': chans,
                     'aia_tr': aia_tr,
                     'tresp_logt': tresp_logt}
                            
    if real_aia_err:
        print('Adding 10% error in quadrature with aiapy.estimate_error output.')
        tens = 0.1*np.copy(aia_dn_s_px)
        #print(aia_err_dn_s_px)
        newerr = [(aia_err_dn_s_px[i]**2+tens[i]**2)**(1/2) for i in range(0, len(aia_err_dn_s_px))]
        #print(newerr)

        res = aia_dn_s_px, newerr, chans, aia_tr, tresp_logt
        if NCCS:
            deminput_dict['newerr']=newerr

            
    else:
        res = aia_dn_s_px, chans, aia_tr, tresp_logt

    if NCCS:
        save_path = pathlib.Path(NCCS_save_path)
        if not save_path.exists():
            save_path.mkdir()
        
        with open(NCCS_save_path+'/'+timestring+'_saved_AIA_DEM_Inputs.pickle', 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(deminput_dict, f, pickle.HIGHEST_PROTOCOL)

        return deminput_dict
        
    return res




def get_region_data(map_obj: sunpy.map.Map,
    region: regions.SkyRegion,
    fill_value: float = 0,
    b_full_size: bool = False
) -> np.ndarray:
    """
    Shared by Reed Masek for use in doing occulted flare DEMs with novel regions.
    
    Get the map data contained within the provided region.

    Parameters
    ----------
    map_obj : sunpy.map.Map
        The map containing the region of interest.
    region : regions.SkyRegion
        The bounding region. Can be any SkyRegion defined in the regions package.
    fill_value : float
        The default null value in indices outside the region.
    b_full_size : bool
        Specifies whether the returned array, region_data,
        is the same shape as the input array, data.
        The default is False since it is wasteful in memory.

    Returns
    -------
    region_data : np.ndarray
        An array containing only the pixel information within
        the provided reg.
    """

    map_data = map_obj.data
    #print(map_obj.wcs)
    reg_mask = (region.to_pixel(map_obj.wcs))
    #print(dir(reg_mask))
    reg_mask=reg_mask.to_mask()
    #print(dir(reg_mask))
    #print(reg_mask.bbox)
    xmin, xmax = reg_mask.bbox.ixmin, reg_mask.bbox.ixmax
    ymin, ymax = reg_mask.bbox.iymin, reg_mask.bbox.iymax
    #print('bound values:', xmin, xmax, ymin, ymax)
    #print('Shape of reg_mask data:', reg_mask.data.shape)
    #print('Shape of map data, indexed with bound values:', map_data[ymin:ymax, xmin:xmax].shape)
    #print('Shape of map data, no change:', map_data.shape)
    #print('Y bound max-min, X bound max-min:', ymax-ymin, xmax-xmin)
    region_data = np.where(reg_mask.data==1, map_data[ymin:ymax, xmin:xmax], fill_value)

    if b_full_size:
        a = np.full(
            shape=map_data.shape,
            fill_value=fill_value,
            dtype=region_data.dtype
        )
        a[ymin:ymax, xmin:xmax] = region_data
        region_data = a

    return region_data  


def read_deminputs(file):
    """
    Wrapper to read in AIA DEM inputs. 
    
    """
    import pickle
    
    with open(file, 'rb') as f:
            deminput_dict = pickle.load(f)

    if 'newerr' in deminput_dict:
        res = deminput_dict['aia_dn_s_px'], deminput_dict['newerr'], deminput_dict['chans'], \
            deminput_dict['aia_tr'], deminput_dict['tresp_logt']
    else:
        res = deminput_dict['aia_dn_s_px'], deminput_dict['chans'], \
            deminput_dict['aia_tr'], deminput_dict['tresp_logt']

    return res 

def circle_region_wrapper(offset, rad):
    """
    Makes a box (1200 arcseconds square) and a circle region, for use as DEM inputs).
    """
    xx = offset[0].value
    yy = offset[1].value

    #Set broad box for plotting (using region object)
    bl=[(xx-700)*u.arcsec, (yy-700)*u.arcsec]
    tr=[(xx+700)*u.arcsec,(yy+700)*u.arcsec]
    #print(tr[0]-bl[0], tr[1]-bl[1])
    
    input_region='circle'
    input_aia_region_dict={'center': (xx,  yy)*u.arcsec,
                      'radius': rad}

    return bl, tr, input_region, input_aia_region_dict


def prep_interval_dir(interval_dir, data_dir, NCCS_save_path, map_save_path,
                      NCCS_aia_resp_path, errortab, resprint=False,
                     clobber=False, aia_clobber=False):
    
    
    files = glob.glob(interval_dir+'/*.pickle')
    files.sort()

    for f in files:
        print(f)
        file_prep(f, data_dir, NCCS_save_path, map_save_path,
                NCCS_aia_resp_path, errortab,
                clobber=clobber, aia_clobber=clobber)
    
    if resprint:
        for f in files:
            with open(f, 'rb') as f_:
                    data = pickle.load(f_)
        
            #print(data.keys())
            print(data['time_interval'])

def file_prep(f, data_dir, NCCS_save_path, map_save_path,
                NCCS_aia_resp_path, errortab,
                clobber=False, aia_clobber=False):
    
    with open(f, 'rb') as f_:
        data = pickle.load(f_)

        if not clobber:
            if 'aia_dn_s_px' in data:
                print('already prepped, and clobber=False, skipping.')
                return
    
        rad = data['radius']*u.arcsec
        offset = [data['centerx'], data['centery']]
        time = data['time_interval']
        
        bl, tr, input_region, input_aia_region_dict = circle_region_wrapper(offset, rad)
    
        deminputs = load_aia(time, bl, tr, plot=False, NCCS=True, 
                             aia_exclude=[], aia_path=map_save_path, 
                 method='Average', one_avg=True,
                 input_region=input_region, input_aia_region_dict=input_aia_region_dict, 
                 real_aia_err=True, errortab=errortab,
                 aia_clobber=aia_clobber, 
                 NCCS_aia_resp_path=NCCS_aia_resp_path,
                 NCCS_save_path=NCCS_save_path,
                 path_to_dodem='./',
                 data_dir=data_dir)
        
        data.update(deminputs)
        print(data.keys())
    
        with open(f, 'wb') as f_:
            pickle.dump(data, f_, pickle.HIGHEST_PROTOCOL)
 

def aia_prep_orbit(data_dir, regions_dir, map_save_path):
    
    where='./'
    
    obsid = regions_dir.split('/')[4][-11:]
    
    files = glob.glob(regions_dir+'/*.pickle')
    files.sort()
    
    pystrings = []
    for ff in files:
        timestring=ff.split('/')[5][0:17]
        pyfile = 'aia_prep_'+timestring+'.py'
        pystring = 'python '+pyfile+' > '+' prep_out_'+timestring+'.txt &'
        pystrings.append(pystring)
        
        with open(where+'aia_prep_template.py', 'r') as f:
            lines = f.read()
            llist = lines.split('\n')
            #print(llist)
            llist[6] = 'map_save_path = "'+map_save_path+'"'
            llist[14] = 'data_dir = "'+data_dir+'"'
            llist[15] = 'f = "'+ff+'"'
            #print(llist)
            newlist = '\n'.join(llist)
        
            with open(pyfile, 'w') as file_out:
                file_out.seek(0)
                file_out.write(newlist)
                file_out.truncate()
    
    with open(where+'aia_prep_all_template.sh', 'r') as f:
        lines = f.read()
        llist = lines.split('\n')
        llist.extend(['','',])
        llist.extend(pystrings)
        llist.extend(['wait', '', 'echo "all orbit scripts finished"'])
        pylist = '\n'.join(llist)
    
        with open(where+'aia_prep_all_'+obsid+'.sh', 'w') as file_out:
            file_out.seek(0)
            file_out.write(pylist)
            file_out.truncate()
    
    print(pylist)


                       
