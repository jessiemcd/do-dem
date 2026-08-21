import nustar_dem_prep as nu
import all_nu_analysis as ana
import nustar_utilities as nuutil
import region_fitting as rf

import pickle
import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u

import sunpy.map
import glob
from astropy.io import fits
import nustar_pysolar as nustar

from astropy.coordinates import SkyCoord
from astropy import coordinates as coord
from regions import CircleSkyRegion




def make_ghost_region(key, file, show=False, gorbit=0, gregion={}, 
                                 eng_tr = [[2.5, 3.5], [3.5, 6.], [6., 10.]],
                                 go_ahead=False, write_new_file=True, 
                             path_to_dodem = '/Users/jmdunca2/do-dem/'):
    """
    For a given observation, makes a ghost ray plot for every orbit. When you've picked which orbit and which
    region will be the ghost ray reference region, use the gorbit and gregion commands to select it. The region 
    will be plotted and saved for later use. 
    
    """
    contourlevs = [1,2,3,5,10,20,25,50,75,99]
    colors = ['orange', 'blue', 'teal', 'indianred', 'springgreen', 'pink', 'purple', 'gold', 'yellow', 'green']

    
    with open(file, 'rb') as f:
        data = pickle.load(f)
    
    ARDict = data[key]
    
    id_dirs = ARDict['datapaths']
    obsids = ARDict['obsids']
    working_dir = ARDict['working_dir']
    #prepped_aia_dir = ARDict['prepped_aia']
    method=ARDict['method']

    
    
    orbit=0
    for id in id_dirs:
        sunpos = glob.glob(id+'/event_cl/*sunpos*.evt')
        count=1

        fig = plt.figure(figsize=(20,8), tight_layout = {'pad': 1})

        for s in sunpos:
            with fits.open(s, ignore_missing_simple=True) as hdu:
                hdr = hdu[1].header
            time0 = nuutil.convert_nustar_time(hdr['TSTART'])
            time1 = nuutil.convert_nustar_time(hdr['TSTOP'])

            with fits.open(s) as hdu:
                evt_data = hdu[1].data
                hdr = hdu[1].header

            cleanevt = nustar.filter.event_filter(evt_data, energy_low=2.5, energy_high=15.,
                                             no_bad_pix_filter=True, no_grade_filter=True)
        
            nustar_map = nustar.map.make_sunpy(cleanevt, hdr)
            bl = SkyCoord( *(-1650, -1250)*u.arcsec, frame=nustar_map.coordinate_frame)
            tr = SkyCoord( *(1650, 1250)*u.arcsec, frame=nustar_map.coordinate_frame)
            submap = nustar_map.submap(bottom_left=bl, top_right=tr)
            ax = fig.add_subplot(1,2,count, projection=submap, aspect='auto')#, sharex=True, sharey=True)
            submap.plot(axes=ax, cmap='Greys')#
            submap.draw_contours(contourlevs, colors=colors[0:len(contourlevs)], axes=ax) #axs[count])
            submap.draw_limb(color='red', axes=ax) #axs[count])
            count+=1

            fpm=s.split('.')[0][-13]
            obsid=s.split('.')[0][-24:-13]
            region_dir=working_dir
            if method == 'fit':
                regfile = region_dir+'gauss_cen_'+obsid+'_'+fpm+'_'+'single.reg'
            if method == 'input':
                print('ONLY DID IT FOR FIRST REGION')
                regfile = region_dir+'gauss_cen_'+obsid+'_'+fpm+'_user_input_0.reg'
            if method == 'double':
                print('ONLY DID IT FOR FIRST REGION')
                regfile = region_dir+'gauss_cen_'+obsid+'_'+fpm+'_0.reg'

            print(regfile)
            #read in regfile (DEM region) and make region object. For multi-region pointings, just the 
            #first region will be plotted.
            offset, rad = rf.read_regfile(regfile, time0, time1, 'hourangle')
            center = SkyCoord( *(offset[0].value, offset[1].value)*u.arcsec, frame=submap.coordinate_frame)
            region = CircleSkyRegion(
                    center = center,
                    radius = rad
                )
        
            #print(offset, rad)
            
            og_region = region.to_pixel(submap.wcs)
            og_region.plot(axes=ax, color='white', ls='--', lw=3, zorder=5)

            #If this is the orbit we are using for extracting the ghost ray background, and an input ghost ray region 
            #has been entered:
            if orbit==gorbit:
                #If we want to write a new file containing the ghost ray region we entered as an input:
                if write_new_file and gregion:
                    center = SkyCoord( *(gregion['centerx'], gregion['centery'])*u.arcsec, frame=submap.coordinate_frame)
                    region = CircleSkyRegion(
                            center = center,
                            radius = gregion['radius']*u.arcsec
                        )
                
                    midway = time0 + (time1-time0).to(u.s).value/2*u.s
                    newfile=working_dir+'ghost_ref_'+obsid+'_'+fpm
                    rf.write_regfile('starter_region.reg', midway, region, newfile=newfile) 

                else:
                    newfile=working_dir+'ghost_ref_'+obsid+'_'+fpm

                    offset, rad = rf.read_regfile(newfile+'.reg', time0, time1, 'hourangle')
                    
                    center = SkyCoord( *(offset[0].value, offset[1].value)*u.arcsec, frame=submap.coordinate_frame)
                    region = CircleSkyRegion(
                            center = center,
                            radius = rad
                        )
                    print('here')
                    print('offset')
                    
                og_region = region.to_pixel(submap.wcs)
                og_region.plot(axes=ax, color='red', ls='--', lw=3, zorder=5)

                #If we want to move forward with calculating the ghost ray background rate (time consuming, 
                #only do when sure of inputs). 
                if go_ahead:
                    time = [time0, time1]
                    nustar_path = working_dir
                    gtifile=id+'event_cl/nu'+obsid+'A06_gti.fits'
                    
                    res = nu.load_nustar(time, eng_tr, nustar_path, fpm, 
                                   make_nustar=True, gtifile=gtifile, datapath=id, regfile=newfile+'.reg', 
                                    edit_regfile=False, compare_fpm=False, combine_fpm=False, actual_total_counts=True, centroid_region=False,
                                       use_fit_regfile=False, clobber=False, default_err=0.2, pile_up_corr=True, special_pha='',
                                       adjacent_grades=True, nuradius=150, path_to_dodem=path_to_dodem, shush=False,
                                       twogauss=False, onegauss=False, direction='', guess=[], guess2=[], return_for_pile_up_figure=False,
                                       energy_percents=False)
                    
                    rate, erate, nutrs, tresp, logt, fpm, atc = res
                    
                    ARDict['Ghost_Ray_Rate_'+fpm] = rate
                

        if orbit==gorbit:  
            plt.savefig(working_dir+'ghost_plot_efilter2.5-15_'+obsids[orbit]+'ghost_ref.png')
        else:
            plt.savefig(working_dir+'ghost_plot_efilter2.5-15_'+obsids[orbit]+'.png')
        orbit+=1


    data[key].update(ARDict)
    with open(file, 'wb') as f:
         # Pickle the 'data' dictionary using the highest protocol available.
         pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)  
