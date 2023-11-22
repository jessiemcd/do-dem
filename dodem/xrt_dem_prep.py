import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import hissw

#https://docs.python.org/3/library/glob.html
import glob
import astropy.time
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import sunpy.map
import scipy.io as io
import os
import copy
import regions


"""

Preparing XRT Data for DEM: IDL/other helpers + other needed prep 
----------------------------------------------------------------------
        
        NOTE: Actual download + prep of XRT data is not done automatically (yet, feel free to add). 
                IDL FUNCTION lets_get_that_xrt can be used for this (input time interval + pointing, 
                and nustar obsid or other label to define a folder where data will be placed. 
                
                Generally, XRT data for a wide interval can be downloaded at once. The python code will
                go in (to the xrt_path directory) and select data relevant to each DEM interval. 
                
                Sample run: 
                
                IDL> nutimes='2018-may-29 '+['1909','1955'] 
                IDL> nupointing=[-110., 250]
                IDL> lgtx = lets_get_that_xrt(nutimes=nutimes, nupointing=nupointing, $
                IDL>    download=1, prep=1, obsid=80410203001)
                
                This will download level 0 files, and place them in './80410203001_coobs/'
                It will then prep files, and place level 1 files + grade maps in './80410203001_coobs/XRT_for_DEM/'
                
                

        
        xrt_tresp_hissw_wrapper_aiaspec.pro - for getting XRT response, called automatically using HISSW if there
                                        is not a proper response file in the xrt_path. 
                                        
                                        The response is date-specific. So this should be re-done for data 
                                        on different days.
                                        
                                        IDL-only version (to run in IDL session, if prefered): 
                                        make_xrtresp_forpy.pro

"""

exposure_dict={'Be_thin': [],
                'Be_thick': [],
              'Al_poly': []}

def load_xrt(xrt_path, time, bl, tr, xrt_exclude=[], plot=True, method='First', exposure_dict=exposure_dict,
            input_xrt_region_dict=[], input_xrt_region=[], real_xrt_err=False, 
             path_to_dodem='./', working_dir='./'):
    """
    Prep steps to do a DEM with XRT data.
    
    Keywords
    ---------
    
    xrt_path - Where XRT data and response files should be.
    
    time     - time interval for DEM (tuple of astropy Time objects)
                FORMAT LIKE, 
                time=(astropy.time.Time('2018-05-29T19:08:00', scale='utc'), 
                        astropy.time.Time('2018-05-29T19:14:00', scale='utc')
                
    bl, tr   - To make XRT sub-map (only using data from a certain region
                to do the DEM). 
                
    xrt_exclude - List of XRT filter combination IDs to NOT use, even if
                    there is data availible. 
    
    plot - set True to plot an XRT image of the sub-map used for DEM.
    
    method – set method='First' to use the first file in each XRT filter 
                combination during the time interval. Set method='Average'
                to average all of the files with that filter combination.
    
    exposure_dict - Dictionary of maximum exposure times for each XRT filter combination. Used for excluding files
                    with longer exposure times in conditions where we have noticed those images are saturated.
    
    """
    
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring
    
    #======================================================
    #GET XRT DATA
    
    #print('XRT EXCLUDE:', xrt_exclude)
    
    if plot:
        if not os.path.exists(working_dir+'/xrt_images/'):
            os.makedirs(working_dir+'/xrt_images/')
    
    data_method=1
    
    if data_method==0:
        #Look for XRT data - assumes desired files are the only thing in the xrt_path.
        #Not recommended
        xdata=sorted(glob.glob(xrt_path+'XRT_*.fits'))
        xgm=sorted(glob.glob(xrt_path+'gm_XRT_*.fits'))
        if len(xdata)==0:
            print('Did not find any XRT data at '+xrt_path)
            print('Not using XRT.')
            return
        
        #print(xdata, xgm)
        
    if data_method==1:
        
        #Look for XRT data
        #dur=120
        filters = ['Be_thin', 'Be_med', 'Al_poly']
        time_range = time
        
        xdata=[]
        xgm=[]
        #Only making the actual data/filter list here if using the averaging method - if not, it's done below.
        xdnspxs=[]
        if real_xrt_err:
            xrt_errs=[]
        filters_=[]
        for w in filters:
            #no exposure limit (will be updated if limit IS set)
            exposure_lim=[]
            files, gmfiles = gather_xrt_files(
                        xrt_path,
                        astropy.time.Time(time_range),
                        w, True
            )   
            
            #print(xrt_path)
            
            if bool(files) == False:
                print('No ', w, ' files.')
                continue

            files = np.sort(files)
            gmfiles = np.sort(gmfiles)   
            
            if w in exposure_dict.keys():
                exposure_lim = exposure_dict[w]
            
            if len(files) >= 1 and method=='First':
                #If using first file in each interval, append that file to the list (as well as grade map)
                xdata.append(files[0])
                xgm.append(gmfiles[0])
                       
                
            if len(files) >= 1 and method=='Average':
                xdnspx_all = []
                if real_xrt_err:
                    errs_all=[]
                files = [xrt_path+f for f in files]
                gmfiles = [xrt_path+f for f in gmfiles]
                for f in range(0, len(files)):
                    filename = files[f].split('/')[-1]
                    filename_stem = working_dir+'/xrt_images/'+filename.split('.')[0]
                    thefilter=filename.split('_')[3]+'-'+filename.split('_')[4]
                    #print('XRT EXCLUDE:', xrt_exclude, 'the filter', thefilter)
                    if thefilter in xrt_exclude:
                        #print('Availible filter excluded: ', thefilter)
                        continue
                    else:                        
                        if real_xrt_err:
                            #Find DN/s/px from file in region
                            res = load_xrt_filter(files[f], gmfiles[f], bl, tr, saveimage=filename_stem, plot=plot,
                                                          real_xrt_err=real_xrt_err, exposure_lim=exposure_lim,
                                                          input_xrt_region_dict=input_xrt_region_dict,
                                                          input_region=input_xrt_region)
                            if res is not None:
                                xdnspx, err = res
                                xdnspx_all.append(xdnspx)
                                errs_all.append(err)
                        else:
                        
                            res = load_xrt_filter(files[f], gmfiles[f], bl, tr, saveimage=filename_stem, plot=plot,
                                              exposure_lim=exposure_lim, input_xrt_region_dict=input_xrt_region_dict,
                                              input_region=input_xrt_region)
                            if res is not None:
                                xdnspx_all.append(res)
                
                #If the current filter has at least one file fitting our criteria, and we haven't decided to exclude it...
                if thefilter not in xrt_exclude and bool(xdnspx_all):
                    xdnspx = np.mean(xdnspx_all)
                    xdnspxs.append(xdnspx)
                    if real_xrt_err:
                        errmean = np.mean(errs_all)
                        xrt_errs.append(errmean)
                    filters_.append(thefilter)


    #======================================================
    
    #======================================================
    #PREPARE DATA FROM EACH FILTER - non-averaging case
    
    if method=='First':
        #this is a list of files from ALL filters (one each)
        xdata = sorted([xrt_path+f for f in xdata])
        xgm = sorted([xrt_path+f for f in xgm])
        #If we're using xrt (and found some data), prep for use in DEM
        #For each file/filter
        for i in range(0, len(xdata)):
            #Extract filter name and add to list
            filename = xdata[i].split('/')[-1]
            thefilter=filename.split('_')[3]+'-'+filename.split('_')[4]
            if thefilter in xrt_exclude:
                #print('Availible filter excluded: ', thefilter)
                continue
            #if thefilter in filters:
            #    continue
            filters_.append(thefilter)
            
            if real_xrt_err:
                #Find DN/s/px from file in region
                xdnspx, err = load_xrt_filter(xdata[i], xgm[i], bl, tr, saveimage=thefilter, plot=plot,
                                              real_xrt_err=real_xrt_err, input_xrt_region_dict=input_xrt_region_dict,
                                              input_region=input_xrt_region, exposure_lim=exposure_lim)
                xrt_errs.append(err)
            else:
                #Find DN/s/px from file in region
                xdnspx = load_xrt_filter(xdata[i], xgm[i], bl, tr, saveimage=thefilter, plot=plot,
                                     input_xrt_region_dict=input_xrt_region_dict, input_region=input_xrt_region,
                                        exposure_lim=exposure_lim)
            print(xdnspx, ' DN/s/px in '+thefilter)
            xdnspxs.append(xdnspx)
        
    
    if bool(xdnspxs) == False:
        print('No XRT files found, exiting.')
        return
    
    filters=filters_
    xdnspxs=np.array(xdnspxs)
    if real_xrt_err:
        xrt_errs=np.array(xrt_errs)
    
    #======================================================
    
    #======================================================
    #GET RESPONSE
    
    #Look for xrt response file 
    #tr_xrt_file=glob.glob('xrt_tresp_*.dat')
    datestring=time[0].value[0:10]
    respfile = 'xrt_tresp_'+datestring.split('-')[0]+datestring.split('-')[1]+datestring.split('-')[2]+'_aiaspec.dat'
    print(respfile)
    
    try:
        tr_xrt=io.readsav(xrt_path+respfile)
    except FileNotFoundError: 
        print('No response file found, now trying to make the response file using HISSW to run IDL code.')
        
        ssw = hissw.Environment(ssw_packages=['hinode/xrt'], ssw_paths=['xrt', 'aia'])
        #Default emission model
        #agr_path = '/Users/jessieduncan/dems/xrt_tresp_hissw_wrapper.pro'
        #AIA emission model (more consistent with assumptions made in calculating AIA response)
        agr_path = path_to_dodem+'/hissw_idl/xrt_tresp_hissw_wrapper_aiaspec.pro'
        inputs = {'filters': filters, 'time': [datestring], 'xrt_path': [xrt_path]}
        try:
            ssw_resp = ssw.run(agr_path, args=inputs)
            tr_xrt=io.readsav(xrt_path+respfile)
        except Exception:
            import traceback
            print(traceback.print_exc())
            print('Something is wrong with the SSWIDL run - make sure the following IDL script exists:')
            print(agr_path)
            print('And make sure search_network.pro is in the working directory, or in your IDL path!')
            print('')
            print('It will ALSO fail if there are no XRT files in xrt_path in the time_range.')
            return

    #re-format filter names
    filters_res=np.array(tr_xrt['filters'])
    for i in range(0,len(filters_res)):
        filters_res[i]=filters_res[i].decode('utf-8')
        
    print(filters)
    print(filters_res)

    #Check that we have a response for every filter we have data from
    check=all(item in filters_res for item in list(filters))
    
    print('check:',check)
    if check==False:
        print('Not all XRT filters from data files are in response file.')
        print('Data filters: ', filters)
        print('Response filters: ', list(filters_res))
        print("Going to make a nice new xrt response file using HISSW to run IDL code.")
        
        ssw = hissw.Environment(ssw_packages=['hinode/xrt'], ssw_paths=['xrt'])
        agr_path = path_to_dodem+'/hissw_idl/xrt_tresp_hissw_wrapper_aiaspec.pro'
        inputs = {'filters': filters, 'time': [datestring], 'xrt_path': [xrt_path]}
        try:
            ssw_resp = ssw.run(agr_path, args=inputs)
            tr_xrt=io.readsav(xrt_path+respfile)
        except Exception:
            import traceback
            print(traceback.print_exc())
            print('Something is wrong with the SSWIDL run - make sure the following IDL script exists:')
            print(agr_path)
            return
        
        #re-format filter names
        filters_res=np.array(tr_xrt['filters'])
        for i in range(0,len(filters_res)):
            filters_res[i]=filters_res[i].decode('utf-8')
            
        #Check that we have a response for every filter we have data from
        check=all(item in filters_res for item in list(filters))
        if check==False:
            print('Huh, filter list is still inconsistent. Better go check on that. Not using XRT')
            print('Data filters: ', filters)
            print('Response filters: ', list(filters_res))
            return

    #Need to make sure that we include the correct response for each filter in the data, 
    #in the correct order.
    response=np.zeros((len(filters), len(tr_xrt['logt'])))
    for i in range(len(filters)):
        index = list(filters_res).index(filters[i])
        response[i,:]=tr_xrt['tr'][index,:]
        
    #======================================================
    
    
    
    if real_xrt_err:
        print('Adding 10% error in quadrature with real error.')
        tens = 0.1*np.copy(xdnspxs)
        #print(aia_err_dn_s_px)
        newerr = [(xrt_errs[i]**2+tens[i]**2)**(1/2) for i in range(0, len(xrt_errs))]
        #print(newerr)
        
        return xdnspxs, newerr, filters, response, tr_xrt['logt'] 
    
    
    return xdnspxs, filters, response, tr_xrt['logt']





def load_xrt_filter(data, gm, bl, tr, plot, saveimage='test', exposure_lim=[],
                   input_xrt_region_dict=[], input_region=[], real_xrt_err=False,
                   grade_mask=True):
    """
    For a single XRT filter, load in the data, plot an image (save it),
    and return DN/s/px (non-binned pixels) for the region.
    
    data and grade map files are outputs of make_xrt_for_python.pro
    
    Keywords
    ---------
    data - data file 
    gm - grade map file
    bl, tr - define regtangular region for DEM (bottom left, top right in arcsec)
        FORMAT LIKE, bl=[-200*u.arcsec,150*u.arcsec]
                      tr=[-50*u.arcsec,300*u.arcsec]
    plot - set to plot + save.              
    saveimage - string for name of saved image
    
    exposure_lim - upper and lower limits for exposure time of images in this image's filter combination
                    If set, will be used to check for compliance. If the exposure time is out of the
                    chosen range, return None
    
    input_region - If set, will expect for input_xrt_region_dict to be a dictionary
                    of parameters to supply in order to make a region object.
                    Currently set up for rectangular or circular regions only. 
                    This region will be used instead of the sub-map defined by bl, tr
                    to extract data for DEM.
                    
    input_xrt_region_dict - See above                
    
    real_xrt_err - If set, return an uncertainty estimate (using expression from Lee et al. 2017)
                    in addition to DN/s/px. 
                    
    grade_mask - If set, use grade map (xrt_prep output) to mask out bad pixels, dust, bleeding (in all images)
                    and also contamination spots in Al-poly images. 
                    
    """
    #======================================================
    #MAKE MAP (& PLOT?)
    
    xmap=sunpy.map.Map(data)
    #(it's level 1)
    #print(xmap.processing_level)
    xgmmap=sunpy.map.Map(gm)
    
    #Make submap using bl, tr
    bottom_left = SkyCoord(bl[0],bl[1], frame = xmap.coordinate_frame)
    top_right = SkyCoord(tr[0],tr[1], frame=xmap.coordinate_frame)
    regxmap = xmap.submap(bottom_left=bottom_left, top_right=top_right)
    reggmmap = xgmmap.submap(bottom_left=bottom_left, top_right=top_right)
    
    #This prints the filters and exposure time
    #print(xmap.fits_header['EC_FW1_'], xmap.fits_header['EC_FW2_'], xmap.exposure_time)
    
    #If we are using an exposure time condition, check the image satisfies it or else exit
    if bool(exposure_lim):
        if regxmap.exposure_time > exposure_lim[1] or regxmap.exposure_time < exposure_lim[0]:
            if plot==True:
                fig = plt.figure(figsize=(9, 7))
                regxmap.plot()
                plt.savefig(saveimage+'_exposure_excluded_xrt_image.png')
                fig = plt.figure(figsize=(9, 7))
                reggmmap.plot()
                plt.savefig(saveimage+'_exposure_excluded_xrt_grade_map.png')                
                
            return None
    
    #Grade map data array
    rd = reggmmap.data

    #zeros = len(np.where(rd == 0)[0])
    ones = len(np.where(rd == 1)[0])
    #twos = len(np.where(rd == 2)[0])
    #fours = len(np.where(rd == 4)[0])
    #eights = len(np.where(rd == 8)[0])
    #sixteens = len(np.where(rd == 16)[0])
    #thirtytwos = len(np.where(rd == 32)[0])
    
#     print('OK 0s:', zeros)
#     print('Saturated 1s:', ones)
#     print('Bloom/Bleed 2s:', twos)
#     print('Contamination Spot 4s:', fours)
#     print('Dust Speck 8s:', eights)
#     print('Hot Pixel 16s:', sixteens)
#     print('Dust Growth 32s:', thirtytwos)
#     print('')    
    
    #Check for saturation, exit if there is any (exposure time limits should have taken care of this)
    if ones > 0:
        print('Saturation in ', regxmap.meta['ec_fw1_'], ' at ', regxmap.meta['date_obs'], ', skipping.')
        print(ones, ' total saturated pixels.')
        if plot==True:
            fig = plt.figure(figsize=(9, 7))
            regxmap.plot()
            plt.savefig(saveimage+'_saturation_excluded_xrt_image.png')
            fig = plt.figure(figsize=(9, 7))
            reggmmap.plot()
            plt.savefig(saveimage+'_saturation_excluded_xrt_grade_map.png')   
            
        return None
    
    
    #If we are removing bad event grades...
    if grade_mask:
        if regxmap.meta['ec_fw1_'] == 'Al_poly':
            #For AL-poly files, mask out contamination spots along with hot pixels, bleeding, dust
            mask = np.where(rd > 1, 0, 1)
            #number_of_pixels = np.sum(mask)
            clean_image = mask*regxmap.data
        else:
            #For Be-thin, Be-med, mask out only hot pixels, bleeding, + dust
            mask = np.where(np.logical_or(rd == 2, rd > 4), 0, 1)
            #number_of_pixels = np.sum(mask)
            clean_image = mask*regxmap.data
    
        
    #======================================================
    
    #======================================================
    #GET EMISSION FROM MAP
    
    #  What is the DN/s/px from the region???
    
    #What is the exposure time?
    dur=regxmap.exposure_time.value
    # What pixel binning per dimension
    chipsum=regxmap.meta['chip_sum']
    
    if bool(input_region)==False:
        #  Get a DN/s/px (non-binned pixels) for the sub-map (not using region object)
        
        if grade_mask:
            positive_pix = clean_image[np.where(clean_image > 0)]      
        else:
            positive_pix = regxmap.data[np.where(regxmap.data > 0)]

        xdnspx=np.mean(positive_pix)/dur/chipsum**2

        if real_xrt_err:
            #Using expression from Lee et al (2017) - Kathy recommended
            uncertainty_list = (1 + np.sqrt(positive_pix + 0.75))/dur/chipsum**2
            err_xrt = np.mean(uncertainty_list)
        
        if xdnspx > 0 and plot:
            fig = plt.figure(figsize=(9, 7))
            regxmap.plot()
            plt.savefig(saveimage+'_xrt_image.png')
            fig = plt.figure(figsize=(9, 7))
            reggmmap.plot()
            plt.savefig(saveimage+'_xrt_grade_map.png')                
        
    else:
        #i.e.: if we ARE using an input region object. 
        region_data=input_xrt_region_dict
        subm = copy.deepcopy(xmap)

        #Make the region object if it's a rectangle
        if input_region=='rectangle':
            region = regions.RectangleSkyRegion(
                SkyCoord(*region_data['center'], frame=subm.coordinate_frame ),
                width=region_data['width'],
                height=region_data['height'],
                angle=region_data['angle']
            )        
        
        #Make the region object if it's a circle
        if input_region=='circle':
            region = regions.CircleSkyRegion(
                SkyCoord(*region_data['center'], frame=subm.coordinate_frame ),
                region_data['radius']
            )     
        
        data = get_region_data(subm, region, b_full_size=True, fill_value=-20)            
        positive_pix = data[np.where(data > 0)]
        xdnspx=np.mean(positive_pix)/dur/chipsum**2
        uncertainty_list = (1 + np.sqrt(positive_pix + 0.75))/dur
        
        #print('Number of pixels with positive values (pre-mask):', len(positive_pix))
        #print('Number of pixels with negative values (pre-mask):', len(data[np.where(data < 0)]))
        
        #print('Pre-mask value in', xmap.meta['ec_fw1_'] , xdnspx)
            
        if grade_mask:
            subgm = copy.deepcopy(xgmmap)
            gm_data = get_region_data(subgm, region, b_full_size=True, fill_value=-20)
            
            rd=gm_data
            
#             if xmap.meta['ec_fw1_'] == 'Al_poly':
#                 fig = plt.figure(figsize=(9, 7))
#                 plt.imshow(gm_data)
#                 plt.title('Grade Mask Data - in region')
            
            ones = len(np.where(rd == 1)[0])
            
            if ones > 0:
                print('Saturation in ', regxmap.meta['ec_fw1_'], ' at ', regxmap.meta['date_obs'], ', skipping.')
                print(ones, ' total saturated pixels.')
                return None
            
            if xmap.meta['ec_fw1_'] == 'Al_poly':
                #For AL-poly files, mask out contamination spots along with hot pixels, bleeding, dust
                #First term: array that is 1 at locations of good pixels and excluded pixels
                #Second term: array that is 1 at locations outside the region (subtracted out)
                mask = np.where(gm_data <= 1, 1, 0) - np.where(gm_data == -20, 1, 0)
#                 print(np.min(mask), np.max(mask))
#                 fig = plt.figure(figsize=(9, 7))
#                 plt.imshow(mask)
#                 plt.title('Mask')
                #number_of_pixels = np.sum(mask)
                clean_image = mask*data
                
#                 fig = plt.figure(figsize=(9, 7))
#                 plt.imshow(mask)
#                 plt.title('Clean Image')
                
            else:
                #For Be-thin, Be-med, mask out only hot pixels, bleeding, + dust
                #First term: array that is 0 at locations of pixels to remove, and 1 otherwise
                #Second term: array that is 1 at locations outside the region (subtracted out)
                mask = np.where(np.logical_or(gm_data == 2, gm_data > 4), 0, 1) - np.where(gm_data == -20, 1, 0)
                #print(np.min(mask), np.max(mask))
                #print('Number of pixels in mask with value 1:', np.sum(mask)) 
#                 fig = plt.figure(figsize=(9, 7))
#                 plt.imshow(mask)
#                 plt.title('Mask')
#                 plt.colorbar()
                clean_image = mask*data 
              
            
            positive_pix = clean_image[np.where(clean_image > 0)]
            xdnspx=np.mean(positive_pix)/dur/chipsum**2
            uncertainty_list = (1 + np.sqrt(positive_pix + 0.75))/dur
            
            
        if real_xrt_err:
            #Using expression from Lee et al (2017) - Kathy recommended
            uncertainty_list = (1 + np.sqrt(positive_pix + 0.75))/dur/chipsum**2
            err_xrt = np.mean(uncertainty_list)
        
        if xdnspx > 0 and plot:
            fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(projection=subm)
            subm.plot(axes=ax)
            (region.to_pixel(subm.wcs)).plot(ax=ax, color='blue')
            norm = colors.PowerNorm(0.5, 0, 1e3) # Define a different normalization to make it easier to see           
            plt.colorbar(norm=norm)
            plt.savefig(saveimage+'_input_region_xrt_image.png')
    
    #======================================================
    
    if xdnspx < 0:
        print('Negative XRT emission - excluding file.')
        if plot==True:
                fig = plt.figure(figsize=(9, 7))
                regxmap.plot()
                plt.savefig(saveimage+'_exposure_excluded_xrt_image.png')
        return        
      
    if real_xrt_err:
        #print('load filter returns:', xdnspx, err_xrt)
        return xdnspx, err_xrt
    else:
        return xdnspx
    
    
    

def gather_xrt_files(
    in_dir: str,
    time_range: tuple[astropy.time.Time], 
    filter_: str,
    gm: bool,
) -> list[str]:
    """
    Checks in_dir for XRT files in a folder (in_dir) that fall within the specified time_range.
    Returns a list of files names sorted by time.
    From a specific filter.
    
    set gm to return list of grade map files as well 
    """

    times = []
    files = []
    gmfiles = []
    
    dir_files = [f for f in os.listdir(in_dir) if os.path.isfile(os.path.join(in_dir, f))]
    data_files = [f for f in dir_files if f[0] == 'X']
    if gm:
        gm_files = [f for f in dir_files if f[0] == 'g']
        
        
    for f in data_files:
        try:
            with fits.open(f'{in_dir}/{f}') as hdu:
                hdr = hdu[0].header
                obs_time = astropy.time.Time(hdr['DATE_OBS'], format='isot')
                if obs_time >= time_range[0] and obs_time <= time_range[1]:
                    if hdr['EC_FW1_'] == filter_:
                        times.append(obs_time)
                        files.append(f)
        except OSError as e: # Catch empty or corrupted fits files
            print(f'OSError with file {f}: {e}')
       
    files = [f for _, f in sorted(zip(times, files))]
    
    if gm:
        for f in gm_files:
            try:
                with fits.open(f'{in_dir}/{f}') as hdu:
                    hdr = hdu[0].header
                    obs_time = astropy.time.Time(hdr['DATE_OBS'], format='isot')
                    #print(hdr)
                    if obs_time >= time_range[0] and obs_time <= time_range[1]:
                        if hdr['EC_FW1_'] == filter_:
                            times.append(obs_time)
                            gmfiles.append(f)
            except OSError as e: # Catch empty or corrupted fits files
                print(f'OSError with file {f}: {e}')

        gmfiles = [f for _, f in sorted(zip(times, gmfiles))]
        
        return files, gmfiles
    
        
    return files    
    
    
    
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
    
    print(reg_mask.shape, map_data.shape)
    
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
    