#code for making lightcurves (and associated tasks)

import numpy as np
import astropy.time
from astropy.io import fits
import astropy.units as u
import pickle
import os
import sunpy.map
from sunpy.net import attrs as a
from sunpy.net import Fido
from sunpy import timeseries as ts
from datetime import timezone
import datetime

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import regions
from astropy.coordinates import SkyCoord

import region_fitting as rf




default_kwargs = {
        'linestyle': 'dashed',
        'linewidth': 1,
        'marker': 'o',
        'markersize': 2,
    }

wavelengths = [94, 131, 171, 193, 211, 335]
filters = ['Be_thin', 'Be_med', 'Al_poly']


def gather_aia_files(
    in_dir: str,
    time_range: tuple[astropy.time.Time], 
    wave: int,
    fulldisk: bool,
) -> list[str]:
    """
    Checks in_dir for AIA fits files that fall within the specified time_range.
    Returns a list of files names sorted by time.
    From a specific channel.
    Set fulldisk=True to return files that are full-disc only (False if you want cutouts)
    """

    times = []
    files = []
    dir_files = [f for f in os.listdir(in_dir) if os.path.isfile(os.path.join(in_dir, f))]
    for f in dir_files:
        #print(f)
        try:
            with fits.open(f'{in_dir}/{f}') as hdu:
                if len(hdu) == 1:
                    hdr = hdu[0].header
                else:
                    hdr = hdu[1].header
                obs_time = astropy.time.Time(hdr['T_OBS'], format='isot')
                if obs_time >= time_range[0] and obs_time <= time_range[1]:
                    if hdr['WAVELNTH'] == wave:
                        full = hdr['NAXIS1'] == 4096 and hdr['NAXIS2'] == 4096
                        #print(hdr['NAXIS1'])
                        if full == fulldisk:
                            times.append(obs_time)
                            files.append(f)
        except OSError as e: # Catch empty or corrupted fits files
            #print(f'OSError with file {f}: {e}')
            print(f'OSError with file {f}')

    files = [f for _, f in sorted(zip(times, files))]

    return files

def make_channel_lightcurve(in_dir, time_range, wave, fulldisk=False, aiaregion=[]):
    """
    time range should be: tuple[astropy.time.Time]
    """
    
    files = gather_aia_files(
        in_dir,
        astropy.time.Time(time_range),
        wave,
        fulldisk
        )
    #print(files)
    paths = [in_dir+ff for ff in files]
    amaps_=sunpy.map.Map(paths)

    amaps = [m/m.exposure_time for m in amaps_]
    
    
    data_mean = []
    data_total = []
    data_time = []
    exp_time = []
    for i in range(0,len(amaps)):
        m=amaps[i]

        if aiaregion:
            region_data=aiaregion
            regionobj = regions.CircleSkyRegion(
                        SkyCoord(*region_data['center'], frame=m.coordinate_frame ),
                        region_data['radius']
                    )    
            data = rf.get_region_data(m, regionobj, b_full_size=True)
            data_mean.append(np.mean(data[np.where(data > 0)]))
            data_total.append(np.sum(data[np.where(data > 0)]))
        else:
            data_mean.append(np.mean(m.data))
            data_total.append(np.sum(m.data))
        data_time.append(m.date)
        exp_time.append(m.exposure_time)
        #exp_time.append(m.meta['exptime'])
        
    times_converted = [t.datetime for t in data_time]
        
    return data_total, times_converted, exp_time

exposure_dict={'Be_thin': [],
                'Be_thick': [],
              'Al_poly': []}

def prepare_lightcurves(in_dir, channels, time_range, instrument, fulldisk=False,
                        plot=True, exposure_dict=exposure_dict, save_dir='./', 
                        aiaregion=[]):
    """
    
    Makes lightcurves from AIA or XRT data files in in_dir from the time range (and with
    the fulldisk condition selected). Pickles the lightcurves for later use.
    
    Doesn't return anything. 
    
    Keywords:
    ---------
    
    in_dir - location of data files
    
    channels - names of each channel/filter
    
    time_range - start, end times, between which to search for data
                like: time_range = ('2018-05-29 22:24:00', '2018-05-29 22:50:00')
    
    instrument - 'AIA' or 'XRT'
    
    fulldisk - if using AIA data, chose whether we're searching for fulldisk or cutout data
    
    plot - set True to plot lightcurves. 
    
    
    Sample call - AIA: 
    -------------------
    
    channels = [94, 131, 171, 193, 211, 335]
    in_dir = '/Users/jessieduncan/sunpy/data/'
    time_range = ('2018-05-29 22:24:00', '2018-05-29 22:50:00')
    fulldisk=False
    instrument = 'AIA'
    plot=True
    
    prepare_lightcurves(in_dir, channels, time_range, instrument, fulldisk=fulldisk, plot=plot)
    
    
    
    Sample call - XRT: 
    -------------------
    
    filters = ['Be_thin', 'Be_med', 'Al_poly']
    in_dir='./other_xrt/'
    time_range = ('2018-05-29 22:24:00', '2018-05-29 22:50:00')
    instrument='XRT'
    plot=True
    
    prepare_lightcurves(in_dir, channels, time_range, instrument, plot=plot)

    """
    


    for w in channels:

        if plot:
            fig, ax = plt.subplots(figsize=(15, 5))
        
        if instrument == 'AIA':
            data_total, times_converted, exp_time = make_channel_lightcurve(
                in_dir,
                astropy.time.Time(time_range),
                w,
                fulldisk=fulldisk,
                aiaregion=aiaregion
            )

            chanlabel='AIA'+str(w)
            
        if instrument == 'XRT':
            data_total, times_converted, exp_time = make_xrt_filter_lightcurve(
            in_dir,
            astropy.time.Time(time_range),
            w,
            exposure_dict)

            chanlabel='XRT_'+w

        if plot:
            exp_time_vals = np.array([e.value for e in exp_time])     
            #exp_time_vals = np.array(exp_time)   
            #print(exp_time)
            #print('')
            #plt.plot(times_converted, np.array(data_total)/exp_time_vals, label=chanlabel, **default_kwargs)
            plt.plot(times_converted, np.array(data_total), label=chanlabel, **default_kwargs)

        data = {
            'times_'+chanlabel: times_converted,
            'data_total_'+chanlabel: data_total,
            'exp_time_'+chanlabel: exp_time,
        }


        with open(save_dir+chanlabel+'_lightcurve.pickle', 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)



        if plot:
            #ax.set_ylim(880000, 900000)   
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))
            plt.legend()
        
    return 

def load_lightcurves(instrument, wavelengths=wavelengths, erange=[2.,10.], fexviii=False, lc_dir='./'):
    """
    Assuming you've already run prepare_lightcurves for your chosen instrument, 
    loads in the pickled results, which should be in files named after each 
    wavelength in the working directory.
    
    Makes one big dictionary with lightcurves, times, exposure times for each 
    channel included (and fexviii, if using with AIA and you have all the needed
    channels.)
    
    To plot, see plot_multi_lightcurves.
    
    """
    
    if instrument == 'GOES':
        try:
            with open(lc_dir+'GOES_lightcurve.pickle', 'rb') as f:
                data = pickle.load(f)
            return data
        except FileNotFoundError:
            print('Not finding prepared GOES data in this directory: ', lc_dir)
            print('Please use get_goes() to do this.')
            return            
    
    if instrument == 'NuSTAR':
        try:
            with open(lc_dir+'NuSTAR_lightcurve_'+str(erange[0])+'_to_'+str(erange[1])+'_keV.pickle', 'rb') as f:
                data = pickle.load(f)
            return data
        except FileNotFoundError:
            print('Not finding prepared NuSTAR data in the range: ', erange, ' in this directory: ', lc_dir)
            print('Please use prepare_nustar_lightcurves() to do this.')
            return 
    
    #print('Using: ', instrument, wavelengths)
    
    all_dict = {
        'chan_list': wavelengths}
    
    if fexviii==True:
        if 94 in wavelengths and 171 in wavelengths and 211 in wavelengths:
            print('Adding Fe-XVIII')
        else:
            print("You need to be including AIA 94, 171, and 211 \AA to make FeXVIII lightcurves.")
            print("Ignoring FeXVIII=True")
            fexviii=False

    clrs=make_colors(len(wavelengths)+1)
    ind=0
    for w in wavelengths:
        
        if instrument=='AIA':
            chanlabel='AIA'+str(w)
        if instrument=='XRT':
            chanlabel='XRT_'+w
            
        try:
            with open(lc_dir+chanlabel+'_lightcurve.pickle', 'rb') as f:
                data = pickle.load(f)
        except FileNotFoundError:
            print('Not finding the following prepared lightcurve file: ', lc_dir+chanlabel+'_lightcurve.pickle')
            print('Use prepare_lightcurves() to prepare data for this instrument.')
            return
            
        all_dict.update(data)
     
    if fexviii==True:
        #============================================================================
        #Checking all arrays are the same length, and correcting if they are one unit off
        #============================================================================
        l94 = len(all_dict['data_total_AIA94'])
        l171 = len(all_dict['data_total_AIA171'])
        
        all_dict['times_AIA94'][-1] < all_dict['times_AIA171'][-2]
        
        if abs(l94-l171) == 1:
            print('Adjusting 171 length as it is one off from 94 \AA ', l94, l171)
            if all_dict['times_AIA94'][-1] > all_dict['times_AIA171'][-2]:
                new171 = all_dict['data_total_AIA171'][0:-1]
                
            elif all_dict['times_AIA94'][0] < all_dict['times_AIA171'][1]:
                new171 = all_dict['data_total_AIA171'][1:]
        else:
            new171 = all_dict['data_total_AIA171']
            
        l211 = len(all_dict['data_total_AIA211'])
        if abs(l94-l211) == 1:
            print('Adjusting 211 length as it is one off from 94 \AA ', l94, l211)
            if all_dict['times_AIA94'][-1] > all_dict['times_AIA211'][-2]:
                new211 = all_dict['data_total_AIA211'][0:-1]
                
            elif all_dict['times_AIA94'][0] < all_dict['times_AIA211'][1]:
                new211 = all_dict['data_total_AIA211'][1:]
        else:
            new211 = all_dict['data_total_AIA211']
            
        #============================================================================
        
        #====================
        #Making FEXVIII array
        #====================
        
        fexviii = np.array(all_dict['data_total_AIA94']) - np.array(new171)/450. - np.array(new211)/120.
        all_dict.update({'data_total_fexviii': fexviii})
        
        
    #print('')    


    return all_dict


def plot_nustar_lightcurves(save_dir='./', timerange=[],
                            eranges = [[2.,4.],[4.,6.],[6.,10.]]):

    """
    To be run after prepare_nustar_lightcurves has already been run for each energy range for the obsid in question.
    
    Single orbit/single obsid functionality.
    
    
    """
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(15, 15))
    instrument='NuSTAR'
    clrs=make_colors(26)
    ind=8
    for er in eranges:
        erange=er

        data = load_lightcurves(instrument, erange=erange, lc_dir=save_dir)
        
        if data is None:
            print('Missing prepared data - exiting.')
            return

        times_convertedA = data['FPMA_times']
        countrateA = data['FPMA_countrate']
        lvtA = data['FPMA_livetime']
        times_convertedB = data['FPMB_times']
        countrateB = data['FPMB_countrate']
        lvtB = data['FPMB_livetime']

        maxA = max(countrateA[np.isfinite(countrateA)])
        maxB = max(countrateB[np.isfinite(countrateB)])

        ax1.plot(times_convertedA, countrateA/maxA, 
                 label='NuSTAR FPMA Counts '+str(erange[0])+' to '+str(erange[1])+' keV (norm)',
                 **default_kwargs, color=clrs[ind])
        ax1.plot(times_convertedB, countrateB/maxB, 
                 label='NuSTAR FPMB Counts '+str(erange[0])+' to '+str(erange[1])+' keV (norm)', 
                 **default_kwargs, color=clrs[ind+1])
        
        ax2.plot(times_convertedA, countrateA, 
                 label='NuSTAR FPMA Counts '+str(erange[0])+' to '+str(erange[1])+' keV',
                 **default_kwargs, color=clrs[ind])
        ax2.plot(times_convertedB, countrateB, 
                 label='NuSTAR FPMB Counts '+str(erange[0])+' to '+str(erange[1])+' keV', 
                 **default_kwargs, color=clrs[ind+1])        
        
        #Normalized boxcar
        n_bx = 5
        arr_lc = np.array(countrateA/maxA)
        avg_lc = boxcar_average(arr_lc, n_bx)
        avg_lc[0:3]=arr_lc[0:3]
        avg_lc[-3:]=arr_lc[-3:]
        ax1.plot(times_convertedA, avg_lc, color=clrs[ind])
        
        #Non-normalized boxcar
        arr_lc = np.array(countrateA)
        avg_lc = boxcar_average(arr_lc, n_bx)
        avg_lc[0:3]=arr_lc[0:3]
        avg_lc[-3:]=arr_lc[-3:]
        ax2.plot(times_convertedA, avg_lc, color=clrs[ind])

        arr_lc = np.array(countrateB/maxB)
        avg_lc = boxcar_average(arr_lc, n_bx)
        avg_lc[0:3]=arr_lc[0:3]
        avg_lc[-3:]=arr_lc[-3:]

        ax1.plot(times_convertedB, avg_lc, color=clrs[ind+1])
                
        #Non-normalized boxcar
        arr_lc = np.array(countrateB)
        avg_lc = boxcar_average(arr_lc, n_bx)
        avg_lc[0:3]=arr_lc[0:3]
        avg_lc[-3:]=arr_lc[-3:]
        ax2.plot(times_convertedB, avg_lc, color=clrs[ind+1])        
        
        ind+=2  
        
        
    ax3.plot(times_convertedA, lvtA, 
                 label='NuSTAR FPMA Livetime',
                 **default_kwargs, color='Black')

    #print('times A, B: ', len(times_convertedA), len(times_convertedB)) 
    #print('lvt A, B: ', len(lvtA), len(lvtB))
    ax3.plot(times_convertedB, lvtB, 
                 label='NuSTAR FPMA Livetime',
                 **default_kwargs, color='Red')

    
    if bool(timerange)==False:
        timerange= [ times_convertedA[0], times_convertedA[-1]]
    print('Using time limits:')
    print(timerange)

    range3 = [0, 2*np.max(lvtA[np.where(np.logical_and(times_convertedA > timerange[0], times_convertedA < timerange[1]))])]


    #ax1.set_ylim(range1[0], range1[1])
    ax1.set_xlim(timerange[0], timerange[1])
    ax1.legend(ncol=2)
    ax1.set_title('Normalized Lightcurves')
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax1.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))

    ax2.legend(ncol=2)
    #ax2.set_ylim(range2[0], range2[1])
    ax2.set_xlim(timerange[0], timerange[1])
    ax2.set_title('Lightcurves')
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax2.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))
    ax2.set_yscale('log')

    ax3.set_ylim(range3[0],range3[1])
    ax3.set_xlim(timerange[0], timerange[1])
    ax3.set_title('Livetimes')
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax3.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))
    ax3.legend()  
    
    plt.savefig(save_dir+'NuSTAR_lightcurves.png')







def plot_multi_lightcurves(plotaia=True, markxrt=True, plotnustar=True, plotGOESA=True, plotGOESB=True, 
                           plotnustar_stats=False, plotnustar_nonnorm=True, 
                           wavelengths=wavelengths, filters=filters, fexviii=True, aianorm=True,
                           nccs_aia=False, nccs_aia_dir='',
                           nustar_file_list=[],
                          rangeaia=[0.,1.], rangenustar=[0.,1.], rangenustar_s=[0.01,1.],
                           rangegA=[5e-9,1e-7], rangegB=[5e-9,1e-7], #GOES RANGES
                          eranges = [[2.,4.],[4.,6.],[6.,8.],[8.,10.]],
                          timerange=[datetime.datetime(2018, 5, 29, 18, 45), datetime.datetime(2018, 5, 29, 19, 45)],
                          plotname='', neupert=False):
    """
    Plot wrapper.
    
    Note that you need to be plotting NuSTAR data for NuSTAR stats (accepted event % vs. rejected) to also be plotted.
    (Set both plotnustar=True, plotnustar_stats=True). rangenustar_s is the y-range for the event statistic ratio. 
    
    """
    if plotaia==False and fexviii==True:
        print('Not set up to do fexviii without other AIA, not doing it.')
        fexviii=False

    instruments = (plotaia, plotnustar, plotnustar_nonnorm, plotnustar_stats, fexviii, plotGOESA, plotGOESB)
    panels = len(np.nonzero(instruments)[0])

    #print(panels)
    #print(instruments)

    fig, axes = plt.subplots(panels, 1, figsize=(9, panels*2))
    count=0
    ranges=[]

    
    if plotGOESA:
        goaxA=axes[count]
        count+=1
        ranges.append(rangegA)
        clr='cornflowerblue'
        instrument='GOES'
        data = load_lightcurves(instrument)
        ylabel = data['GOES flux label']
        goestimes = data['GOES Times']
        xrsacounts = data['XRSA counts']
        xrsalabel = data['XRSA label']

        gts = [t.datetime for t in goestimes]
        maxA = max(xrsacounts[np.isfinite(xrsacounts)])
        goaxA.semilogy(gts, xrsacounts, 
                     label=xrsalabel,
                     **default_kwargs, color=clr)

        n_bx = 10
        arr_lc = np.array(xrsacounts)
        avg_lc = boxcar_average(arr_lc, n_bx)
        avg_lc[0:3]=arr_lc[0:3]
        avg_lc[-3:]=arr_lc[-3:]

        goaxA.semilogy(gts, avg_lc, color=clr)
        
    if plotGOESB:
        goaxB=axes[count]
        count+=1
        ranges.append(rangegB)
        clr='red'
        instrument='GOES'
        data = load_lightcurves(instrument)
        
        ylabel = data['GOES flux label']
        goestimes = data['GOES Times']
        xrsbcounts = data['XRSB counts']
        xrsblabel = data['XRSB label']
        
        gts = [t.datetime for t in goestimes]
        

        goaxB.semilogy(gts, xrsbcounts, 
                     label=xrsblabel, 
                     **default_kwargs, color=clr)

        n_bx = 10
        arr_lc = np.array(xrsbcounts)
        avg_lc = boxcar_average(arr_lc, n_bx)
        avg_lc[0:3]=arr_lc[0:3]
        avg_lc[-3:]=arr_lc[-3:]

        goaxB.semilogy(gts, avg_lc, color=clr)
        goaxB.set_yticks([1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9], labels=['X','M','C','B','A','A0.1'])
        
    
    if plotnustar:
        nuax=axes[count]
        #nuax.set_title('Normalized Lightcurves')
        count+=1
        ranges.append(rangenustar)
        if plotnustar_stats:
            nuax_s=axes[count]
            #nuax.set_title('Normalized Lightcurves')
            count+=1
            ranges.append(rangenustar_s)

        if plotnustar_nonnorm:
            nuax_n=axes[count]
            count+=1

            minset=[]
            maxset=[]
            
            
        instrument='NuSTAR'
        clrs=make_colors(26)
        ind=8
        enum=len(eranges)

        if nustar_file_list:
            num = len(nustar_file_list)
            dn_ins = np.zeros((enum, num))
            times = []
            for i in range(0, num):
                f_ = nustar_file_list[i]
                with open(f_, 'rb') as f:
                    data = pickle.load(f)

                dn_ins[:,i] = np.array(data['dn_in'][-1*enum:])
                times.append(data['time_interval'][0])
                
            times_converted=[t.datetime for t in times]

        for er in range(0, len(eranges)):
            erange=eranges[er]

            if nustar_file_list:
                countrate = dn_ins[er,:]
                max_ = max(countrate[np.isfinite(countrate)])
                nuax.plot(times_converted, countrate/max_, 
                         label='NuSTAR Counts '+str(erange[0])+' to '+str(erange[1])+' keV (norm, FPM sum)',
                         **default_kwargs, color=clrs[ind])

                n_bx = 3
                arr_lc = np.array(countrate/max_)
                avg_lc = boxcar_average(arr_lc, n_bx)
                avg_lc[0:3]=arr_lc[0:3]
                avg_lc[-3:]=arr_lc[-3:]
    
                nuax.plot(times_converted, avg_lc, color=clrs[ind])

            else:
                
                data = load_lightcurves(instrument, erange=erange)
    
                times_convertedA = data['FPMA_times']
                countrateA = data['FPMA_countrate']
                lvtA = data['FPMA_livetime']
                times_convertedB = data['FPMB_times']
                countrateB = data['FPMB_countrate']
                lvtB = data['FPMB_livetime']

                maxA = max(countrateA[np.isfinite(countrateA)])
                maxB = max(countrateB[np.isfinite(countrateB)])
    
                nuax.plot(times_convertedA, countrateA/maxA, 
                         label='NuSTAR FPMA Counts '+str(erange[0])+' to '+str(erange[1])+' keV (norm)',
                         **default_kwargs, color=clrs[ind])
                nuax.plot(times_convertedB, countrateB/maxB, 
                         label='NuSTAR FPMB Counts '+str(erange[0])+' to '+str(erange[1])+' keV (norm)', 
                         **default_kwargs, color=clrs[ind+1])
                
                n_bx = 5
                arr_lc = np.array(countrateA/maxA)
                avg_lc = boxcar_average(arr_lc, n_bx)
                avg_lc[0:3]=arr_lc[0:3]
                avg_lc[-3:]=arr_lc[-3:]
    
                nuax.plot(times_convertedA, avg_lc, color=clrs[ind])
                
                arr_lc = np.array(countrateB/maxB)
                avg_lc = boxcar_average(arr_lc, n_bx)
                avg_lc[0:3]=arr_lc[0:3]
                avg_lc[-3:]=arr_lc[-3:]
    
                nuax.plot(times_convertedB, avg_lc, color=clrs[ind+1])
                #ind+=2

            if plotnustar_nonnorm:

                if nustar_file_list:
                    nuax_n.plot(times_converted, countrate, 
                     label='NuSTAR Counts '+str(erange[0])+' to '+str(erange[1])+' keV (FPM sum)',
                     **default_kwargs, color=clrs[ind])

                    minset.append(np.nanmin(countrate)*0.8)
                    maxset.append(np.nanmax(countrate)*1.2)
                else:
                    
                    nuax_n.plot(times_convertedA, countrateA, 
                         label='NuSTAR FPMA Counts '+str(erange[0])+' to '+str(erange[1])+' keV',
                         **default_kwargs, color=clrs[ind])
                    nuax_n.plot(times_convertedB, countrateB, 
                         label='NuSTAR FPMB Counts '+str(erange[0])+' to '+str(erange[1])+' keV', 
                         **default_kwargs, color=clrs[ind+1])

                    minset.append(np.nanmin(countrateA)*0.8)
                    maxset.append(np.nanmax(countrateA)*1.2)

            ind+=2

        if plotnustar_nonnorm:
            ranges.append([np.min(minset), np.max(maxset)])
                


        if plotnustar_stats:
            if nustar_file_list:
                data = load_lightcurves(instrument, erange=erange)
                times_convertedA = data['FPMA_times']
                times_convertedB = data['FPMB_times']
            
            acc_A=data['FPMA_accepted']
            acc_B=data['FPMB_accepted']
            rej_A=data['FPMA_rejected']
            rej_B=data['FPMB_rejected']  

            nuax_s.plot(times_convertedA, acc_A/(acc_A+rej_A), 
                     label='NuSTAR FPMA % Events Accepted',
                     **default_kwargs, color=clrs[ind+1])
            nuax_s.plot(times_convertedB, acc_B/(acc_B+rej_B), 
                     label='NuSTAR FPMB % Events Accepted',
                     **default_kwargs, color=clrs[ind+2])
            nuax_s.set_yscale('log')

#             ax2.plot(times_convertedA, countrateA, 
#                      label='NuSTAR FPMA Counts '+str(erange[0])+' to '+str(erange[1])+' keV', **default_kwargs)
#             ax2.plot(times_convertedB, countrateB, 
#                      label='NuSTAR FPMB Counts '+str(erange[0])+' to '+str(erange[1])+' keV', **default_kwargs)


    
    
    if markxrt:
        instrument='XRT'
        all_xrt = load_lightcurves(instrument, wavelengths=filters)

        xrtshades=['black', 'grey', 'brown']
        ind=0
        for f in filters:
            chanlabel='XRT_'+f
            v=0
            for t in all_xrt['times_'+chanlabel]:
                if v==0:
                    for ax in axes:
                        ax.axvline(t, color=xrtshades[ind], label=chanlabel+' times')
                    #ax2.axvline(t, color=xrtshades[ind], label=chanlabel+' times')
                    #ax3.axvline(t, color=xrtshades[ind], label=chanlabel+' times')
                    #if plotGOES:
                    #    ax4.axvline(t, color=xrtshades[ind], label=chanlabel+' times')
                else:
                    for ax in axes:
                        ax.axvline(t, color=xrtshades[ind])
                    #ax2.axvline(t, color=xrtshades[ind])
                    #ax3.axvline(t, color=xrtshades[ind])
                    #if plotGOES:
                    #    ax4.axvline(t, color=xrtshades[ind])
                v+=1
            ind+=1

    
    
    
    if plotaia:
        aiaax=axes[count]
        count+=1
        #aiaax.set_title('Normalized Lightcurves')
        ranges.append(rangeaia)
        if fexviii:
            dzax=axes[count]
            #dzax.set_title('Normalized Lightcurves')
            count+=1
        clrs=make_colors(len(wavelengths)+1)
        instrument='AIA'
        if nccs_aia:
            import glob
            files = glob.glob(nccs_aia_dir+'/*')
            files.sort()
            aia_dn_s_pxs=[]
            times=[]
            for f_ in files:
                with open(f_, 'rb') as f:
                    data = pickle.load(f)
            
                aia_dn_s_pxs.append(np.array(data['aia_dn_s_px']))
                times.append(data['time_interval'][0])
            wavelengths = [int(ww[1:]) for ww in data['chans']]
            ts=[t.datetime for t in times]
            aia_dn_s_pxs = np.array(aia_dn_s_pxs)
            
        else:
            all_aia = load_lightcurves(instrument, wavelengths=wavelengths, fexviii=fexviii)
        
        ind=0
        for w in wavelengths:
            chanlabel='AIA'+str(w)
            if nccs_aia:
                times_converted=ts
                corr_totals=aia_dn_s_pxs[:,ind]
                
            else:
                times_converted = all_aia['times_'+chanlabel]
                data_total = all_aia['data_total_'+chanlabel]
                exp_time = all_aia['exp_time_'+chanlabel]

                exp_times = np.array([e.value for e in exp_time])

                corr_totals= data_total/(exp_times)

            n_bx = 5
            if aianorm:
                arr_lc = np.array(corr_totals/max(corr_totals))
                corrplot=corr_totals/max(corr_totals)
            else:
                arr_lc = np.array(corr_totals)
                corrplot=corr_totals
            avg_lc = boxcar_average(arr_lc, n_bx)
            avg_lc[0:3]=arr_lc[0:3]
            avg_lc[-3:]=arr_lc[-3:]

            #print(max(corr_totals))

            aiaax.plot(times_converted, avg_lc, color=clrs[ind])
            aiaax.plot(times_converted, corrplot, label=chanlabel, **default_kwargs, color=clrs[ind])
            
            #if w == 94 or w==171 or w==211:
                #ax3.plot(times_converted, data_total/max(data_total), label=chanlabel, **default_kwargs, color=clrs[ind])
                

            #ax2.plot(times_converted, corr_totals, label=chanlabel, **default_kwargs, color=clrs[ind])
            ind+=1

#         ax1.plot(all_aia['times_AIA94'], all_aia['data_total_fexviii']/max(all_aia['data_total_fexviii']), 
#                  **default_kwargs, color=clrs[ind], label='Del Zanna Fe18')
#         ax2.plot(all_aia['times_AIA94'], all_aia['data_total_fexviii']/max(all_aia['data_total_fexviii']),
#                  **default_kwargs, color=clrs[ind], 
#                  label='Del Zanna Fe18')

        if fexviii:
            if nccs_aia:
                fexviii = np.array(aia_dn_s_pxs[:,0] - aia_dn_s_pxs[:,2]/450. - aia_dn_s_pxs[:,4]/120.)
                max18 = np.max(fexviii)
                dzax.plot(times_converted, fexviii, #/max18, 
                     **default_kwargs, color=clrs[ind], 
                     label='Del Zanna Fe18')
                #print(fexviii)
                #print(max18)

                ranges.append([np.nanmin(fexviii), np.nanmax(fexviii)])
                
            else:
                max18 = max(all_aia['data_total_fexviii'][np.isfinite(all_aia['data_total_fexviii'])])
                dzax.plot(all_aia['times_AIA94'], all_aia['data_total_fexviii']/max18, 
                     **default_kwargs, color=clrs[ind], 
                     label='Del Zanna Fe18')

            if neupert and nccs_aia:
                neuax = dzax.twinx()
                neuax.plot(times_converted, np.gradient(fexviii),
                           **default_kwargs, #color=clrs[ind+1], 
                             label='Del Zanna Fe18 - derivative')
                neuax.legend(fontsize=10, loc='lower right')
                


    #Formatting:
    
    #print(type(timerange[0]), type(timerange[1]))


    for i in range(0,len(axes)):
        ax=axes[i]
        range_=ranges[i]
        ax.set_ylim(range_[0], range_[1])
        ax.set_xlim(timerange[0], timerange[1])
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))
        ax.legend(fontsize=10, loc='upper right')
        #2024 AGU figure neupert effect peak
        #ax.axvline(datetime.datetime(2024, 7, 17, 5, 7), linestyle='dotted')
        
        

    # ax1.set_ylim(range1[0], range1[1])
    # ax1.set_xlim(timerange[0], timerange[1])
    # ax1.legend(ncol=3)
    # ax1.set_title('Normalized Lightcurves')
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    # ax1.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))

    # ax2.legend(ncol=3)
    # ax2.set_ylim(range2[0], range2[1])
    # ax2.set_xlim(timerange[0], timerange[1])
    # ax2.set_title('Lightcurves')
    # ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    # ax2.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))

    # ax3.set_ylim(range3[0],range3[1])
    # ax3.set_xlim(timerange[0], timerange[1])
    # ax3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    # ax3.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))
    # ax3.legend()
    
    # if plotGOES:
    #     ax4.set_ylim(range4[0],range4[1])
    #     ax4.set_xlim(timerange[0], timerange[1])
    #     ax4.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    #     ax4.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))
    #     ax4.legend(ncol=2)
    
    plt.savefig(plotname+'multi_lightcurves.png', transparent=True)

    return


def get_goes(time_range, satellite=14, peek=False):
    """
    Download GOES file from the correct day, select interval of interest, and save data.
    
    """

    if satellite >= 16:
        print('Using 1s GOES flux ("flx1s")')
        result_goes = Fido.search(a.Time(time_range[0], time_range[1]), a.Instrument("XRS"), a.goes.SatelliteNumber(satellite),
                             a.Resolution("flx1s"))
    else:
        result_goes = Fido.search(a.Time(time_range[0], time_range[1]), a.Instrument("XRS"), a.goes.SatelliteNumber(satellite))
        
    # print(result_goes)
    # print(dir(result_goes))
    # print(result_goes.all_colnames)
    # print(type(result_goes))
    # print(result_goes['Resolution'])
    # print('')
    file_goes = Fido.fetch(result_goes)
    goes_ = ts.TimeSeries(file_goes)
    # print(type(goes_))
    # print(goes_.header)
    # print(dir(goes_[0]))
    if peek:
        goes_.peek()

    if not isinstance(goes_, list):

        goes_interval = goes_.truncate(time_range[0],  time_range[1])   
        
        data = {'GOES flux label': "Flux (Wm$^{-2}$$s^{-1}$)",
                'XRSA label': " 0.5-4 $\AA$",
                'XRSB label': " 1-8 $\AA$",
                'GOES Times': goes_interval.time,
                'XRSA counts': goes_interval.quantity("xrsa"),
                'XRSB counts': goes_interval.quantity("xrsb")
                }
        
        with open('GOES_lightcurve.pickle', 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
    else:
        #IF the time range extends over multiple days, append together data from 
        #the different GOES files.
        #print(goes_)
        times=[]
        xrsa=[]
        xrsb=[]
        for g in goes_:
            if g == goes_[0]:
                tr2 = g.time_range.end.datetime
                goes_interval = g.truncate(time_range[0],  tr2)
            elif g == goes_[-1]:
                tr1 = g.time_range.start.datetime
                goes_interval = g.truncate(tr1,  time_range[1])
            else:
                tr1 = g.time_range.start.datetime
                tr2 = g.time_range.end.datetime
                goes_interval = g.truncate(tr1, tr2)

            times.extend(goes_interval.time)
            xrsa.extend(goes_interval.quantity("xrsa").value)
            xrsb.extend(goes_interval.quantity("xrsb").value)

        times = astropy.time.Time(times)
        xrsa = xrsa*goes_interval.quantity("xrsa")[0].unit
        xrsb = xrsb*goes_interval.quantity("xrsb")[0].unit

        
        data = {'GOES flux label': "Flux (Wm$^{-2}$$s^{-1}$)",
                'XRSA label': " 0.5-4 $\AA$",
                'XRSB label': " 1-8 $\AA$",
                'GOES Times': times,
                'XRSA counts': xrsa,
                'XRSB counts': xrsb
                }
        print(times[0], times[-1])
        with open('GOES_lightcurve.pickle', 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
            
    return



def load_nufiles(f):
    """
    """
    hdulist = fits.open(f)
    dat = hdulist[1].data
    hdr = hdulist[1].header
    hdulist.close()
        
    return dat, hdr

def get_a_nustar_lightcurve(evtdata, hdr, lvdata, lvhdr, timebin=10, livetime_corr=True, 
                            event_stats=False):
    """
    
    """
    # Convert times from sec relative to ref time into normal times
    mjd_ref_time=astropy.time.Time(hdr['mjdrefi'],format='mjd')
    alltimes=astropy.time.Time(mjd_ref_time+evtdata['time']*u.s,format='mjd')
    
    # This is the time binning of the livetimes not the actual livetimes - which is lvdata['livetime']
    livetime_tbins=astropy.time.Time(mjd_ref_time+lvdata['time']*u.s,format='mjd')
    
    # Use the 1sec time binning of the livetime for the binning of the counts
    # And just work everything relative to the first livetime time bin

    #time delta between every event time and the first livetime bin
    td=(alltimes-livetime_tbins[0]).sec

    #time delta between every livetime bin and the first livetime bin
    tdedgs=(livetime_tbins-livetime_tbins[0]).sec
    #print(tdedgs)
    tdedgs=tdedgs[0::timebin]
    #print(tdedgs)
    
    # Use histogram to bin events to get number per 1s time bins
    counts, bed=np.histogram(td, bins=tdedgs)

    # #helpful for debugging issues with the livetime and time arrays not ending up the same size
    # print('firstlast: ', lvdata['time'][0], lvdata['time'][-1])
    # print('length: ', len(lvdata['livetime'])) #2
    # print(timebin)
    # print('length lvt: ', len(lvdata['livetime'][0::timebin]))
    # print('length livetime tbins: ', len(livetime_tbins)) #2
    
    lvt=lvdata['livetime'][0::timebin]
    if livetime_corr:
        countrate=counts/lvt[:-1]
    else:
        countrate=counts

    #Previously we got rid of the last entry here, but that makes this fail for obsid 21012007001
    #(magixs 2 observation, day 2). So switched off for now.
    #times_converted = [t.datetime for t in livetime_tbins[:-1]]
    times_converted = [t.datetime for t in livetime_tbins]
    times_converted = [t.replace(tzinfo=timezone.utc) for t in times_converted]
    #print('length times_converted: ', len(times_converted)) #1
    times_converted = np.array(times_converted)[0::timebin]

    #print('length times_converted (resample): ', len(times_converted))
    #print('')
    if abs(len(times_converted)-len(countrate)) == 1:
        times_converted=times_converted[:-1]
        lvt=lvt[:-1]

    if event_stats:
        acc_sample = lvdata['NACCEPT'][0::timebin]
        rej_sample = lvdata['NREJECT'][0::timebin]
        all_sample = lvdata['NEV'][0::timebin]

        return times_converted, countrate, lvt, counts, acc_sample[:-1], rej_sample[:-1], all_sample[:-1] 
        
    else:
        return times_converted, countrate, lvt, counts


def prepare_nustar_lightcurves(evtA, evtB, hkA, hkB, timebin=10, erange=[2.,10.], livetime_corr=True, 
                               return_lightcurves=False, save_dir='./', event_stats=False):
    """
    Returns FPMA + B lightcurves. Wrapper for get_a_nustar_lightcurve() which does just one.
    
    Using some stuff from this example, but customizing for my use:
    https://github.com/ianan/nustar_sac/blob/master/python/example_nustar_lightcurve.ipynb
    
    Default behavior is to just save lightcurves to a file - set return_lightcurves=True to spit them all out
    directly as well. 
    """
    #Load in the evt file (has the list of photons)
    evtdataA, hdrA = load_nufiles(evtA[0])
    # Load in the hk file (has the livetime info)
    lvdataA, lvhdrA = load_nufiles(hkA[0])
    evtdataB, hdrB = load_nufiles(evtB[0])
    lvdataB, lvhdrB = load_nufiles(hkB[0])

    # lvtAtimes = np.array(lvdataA['time'])
    # lvtBtimes = np.array(lvdataB['time'])
    # mjd_ref_time=astropy.time.Time(hdrA['mjdrefi'],format='mjd')
    # print(mjd_ref_time)
    # mjd_ref_time=astropy.time.Time(hdrB['mjdrefi'],format='mjd')
    # print(mjd_ref_time)
    
    # adjust=0
    # for i in range(0, len(lvtAtimes)): 
    #     test = lvtAtimes[i]-lvtBtimes[i]-adjust
    #     if np.abs(test) > 0:
    #         timeA = astropy.time.Time(mjd_ref_time+lvtAtimes[i]*u.s,format='mjd')
    #         timeB = astropy.time.Time(mjd_ref_time+lvtBtimes[i]*u.s,format='mjd')
    #         print(i, test, lvtAtimes[i], lvtBtimes[i])
    #         print(timeA.datetime)
    #         adjust+=test

    
    kevA = evtdataA['PI']*0.04+1.6
    erange_evtdataA = evtdataA[np.where(np.logical_and(kevA > erange[0],kevA < erange[1]))]
    kevB = evtdataB['PI']*0.04+1.6
    erange_evtdataB = evtdataB[np.where(np.logical_and(kevB > erange[0],kevB < erange[1]))]
    
    #print('A')
    resA = get_a_nustar_lightcurve(erange_evtdataA, hdrA, lvdataA, lvhdrA, 
                                    timebin=timebin, livetime_corr=livetime_corr,
                                  event_stats=event_stats)
    #print('B')
    resB = get_a_nustar_lightcurve(erange_evtdataB, hdrB, lvdataB, lvhdrB, 
                                    timebin=timebin, livetime_corr=livetime_corr,
                                   event_stats=event_stats)

    data = {'Livetime-Corrected?': livetime_corr,
            'Time Bin (s)': timebin,
            'Energy Range': erange,
            'file paths': [evtA, evtB, hkA, hkB],
            'FPMA_countrate': resA[1],
            'FPMB_countrate': resB[1],
            'FPMA_counts': resA[3],
            'FPMB_counts': resB[3],
            'FPMA_times': resA[0],
            'FPMB_times': resB[0],
            'FPMA_livetime': resA[2],
            'FPMB_livetime': resB[2]
            }

    if event_stats:
        data['FPMA_accepted']=resA[4]
        data['FPMB_accepted']=resB[4]
        data['FPMA_rejected']=resA[5]
        data['FPMB_rejected']=resB[5]  
        data['FPMA_all_ev']=resA[6]
        data['FPMB_all_ev']=resB[6]

    with open(save_dir+'NuSTAR_lightcurve_'+str(erange[0])+'_to_'+str(erange[1])+'_keV.pickle', 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

    if return_lightcurves==True:
        return times_convertedA, countrateA, lvtA, countsA, times_convertedB, countrateB, lvtB, countsB
    else:
        return


def gather_xrt_files(
    in_dir: str,
    time_range: tuple[astropy.time.Time], 
    filter_: str,
    gm: bool,
) -> list[str]:
    """
    Checks in_dir for AIA fits files that fall within the specified time_range.
    Returns a list of files names sorted by time.
    From a specific channel.
    Set fulldisk=True to return files that are full-disc only (False if you want cutouts)
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

exposure_dict={'Be_thin': [],
                'Be_thick': [],
              'Al_poly': []}

def make_xrt_filter_lightcurve(
    in_dir: str,
    time_range: tuple[astropy.time.Time], 
    filter_: str,
    exposure_dict: dict,
) -> list[str]:
    
    if filter_ in exposure_dict.keys():
        exposure_lim = exposure_dict[filter_]
    
        files = gather_xrt_files(
            in_dir,
            astropy.time.Time(time_range),
            filter_,
            False
            )
        paths = [in_dir+ff for ff in files]
        #xmaps=sunpy.map.Map(paths)
        xmaps=[]
        for p in paths:
            #print('P:', p)
            xmap=sunpy.map.Map(p)
            #print(xmap.exposure_time, filter_)
            #if bool(exposure_lim):
            #    print(exposure_lim[0], exposure_lim[1])
            if bool(exposure_lim) and exposure_lim[0] < xmap.exposure_time < exposure_lim[1]:
                xmaps.append(xmap)
            elif bool(exposure_lim) == False:
                xmaps.append(xmap)

    
            #print('')
        
        data_mean = []
        data_total = []
        data_time = []
        exp_time = []
        for i in range(0,len(xmaps)):
            m=xmaps[i]
            if np.sum(m.data) < 0:
                continue
            data_mean.append(np.mean(m.data))
            data_total.append(np.sum(m.data))
            data_time.append(m.date)
            exp_time.append(m.exposure_time)

        times_converted = [t.datetime for t in data_time]

        return data_total, times_converted, exp_time
    else:
        print('Filter ', filter_, ' not in exposure dictionary')
        return

def plot_with_stdv(aia_inputs=[94], fexviii=True, nustar_inputs=[[2.,4.],[4.,6.],[6.,10.]], 
                   goes_inputs=['xrsa','xrsb'], smooth=24, 
                   plot_each=True, plot_logic=True, remove_fexviii_max=False, 
                   analyze_transients=True, transient_number=3,
                   timerange=[datetime.datetime(2018, 5, 29, 22, 22), datetime.datetime(2018, 5, 29, 23, 19)],
                  excluded_range=[], save_dir='./', savestring='test', show=False):
    """
    For an input set of types of lightcurve (prepared using prepare_lightcurves - see above), 
    over an input time interval:
        - plots each lightcurve in a separate plot (optional)
        - adds hotizontal line for mean in interval, with error bars defining 1 STDV
        - makes two arrays: 
                one which is 1 where lightcurve exceeds 1STDV region (0 elsewhere)
                one which is -1 where lightcurve is less than the 1STDV region (0 elsewhere)
                
    Additionally, plots total number of lightcurves/instrument exceeding + less than 1STDV as function 
    of time, all on the same axes so instruments can be compared. This summary plot is saved.
    
    To exclude an instrument entirely, set its inputs to an empty list (e.g. aia_inputs=[]).
    
    Keywords
    ---------
    
    Instrument Inputs
    ------------------
    
    aia_inputs - list of integer AIA channels, or empty list 
   
    fexviii -            set True to include FeXVIII curve (made from other AIA channels)
    remove_fexviii_max - set True to remove whatever is the maximum value in the fexviii array (helpful for
                          weird outliers surrounding AIA outages). 
    
    nustar_inputs - list of tuples (NuSTAR energy ranges), or empty list
    
    goes_inputs - list of GOES strings ('xrsa', 'xrsb'), or empty list
    
    Method Inputs
    --------------
    
    plot_each - Set True to plot individual lightcurves (with 1-stdv range).
    
    plot_logic - Set True to plot -1, 0, 1 values re whether data is out of the 1-stdv range at each point.
    
    timerange - Must be a subset of the time for which we have prepped data in the files made by prepare_lightcurves.
    
    excluded_range - If there is a pointing shift or flare we want to exclude (so far, only one interval allowed).
    
    analyze_transients - Set True to use transient_analysis function to automatically categorize where multiple
                        channels are above the stdv level. 
    smooth – interval to smooth each lightcurve over (in seconds). Timestep naturally found for each instrument will 
                be used together with this to define boxcar width (in # bins) for averaging. 
    
    """
    
    clrs=make_colors(26)
    ind=0


    number=len(aia_inputs)+2*len(nustar_inputs)+len(goes_inputs)
    if plot_logic:
        number+=1
    if fexviii:
        number+=1
    fig, axes = plt.subplots(number, 1, figsize=(15, number*2), sharex=True)
    plt.subplots_adjust(hspace=0)

    nn=0

    if nustar_inputs:
        instrument='NuSTAR'
        
        n_all_outupsA=[]
        n_all_outdownsA=[]
        n_all_outupsB=[]
        n_all_outdownsB=[]
        first=0
        for er in nustar_inputs:
            erange=er
            data = load_lightcurves(instrument, erange=erange, lc_dir=save_dir)

            times_convertedA = data['FPMA_times']
            countrateA = data['FPMA_countrate']
            lvtA = data['FPMA_livetime']
            times_convertedB = data['FPMB_times']
            countrateB = data['FPMB_countrate']
            lvtB = data['FPMB_livetime']

            tstep=(times_convertedA[1]-times_convertedA[0]).seconds

            accA = data['FPMA_accepted']
            accB = data['FPMB_accepted']
            rejA = data['FPMA_rejected']
            rejB = data['FPMB_rejected'] 
            allA = data['FPMA_all_ev']
            allB = data['FPMB_all_ev']

            if first==0:
                #BECAUSE these things (accepted event percentage, time range, non-zero count rates)
                #should not be energy-range dependent, just doing for the first range and then keeping 
                #selected indices constant.

                accA = data['FPMA_accepted']
                accB = data['FPMB_accepted']
                rejA = data['FPMA_rejected']
                rejB = data['FPMB_rejected'] 
                allA = data['FPMA_all_ev']
                allB = data['FPMB_all_ev'] 
                
                evsumA=accA/(accA+rejA)
                evsumB=accB/(accB+rejB)
                evsum = (np.array(evsumA)+np.array(evsumB))/2.
                #Condition that accepted event percent is above 95
                accthresh_inds = np.where(evsum > 0.9)
                #print('acc', len(accthresh_inds[0]))
                #Condition that both countrates are not 0 (this is usually SAA)
                nonzero_inds = np.where(np.logical_and(countrateA > 0, countrateB > 0))
                #print('nz', len(nonzero_inds[0]))
                #Intersection of the prior two conditions
                allinds = np.intersect1d(accthresh_inds, nonzero_inds)
                
                if bool(excluded_range): 
                    times1 = np.where(np.logical_and(times_convertedA > timerange[0], times_convertedA < excluded_range[0]))
                    times2 = np.where(np.logical_and(times_convertedA > excluded_range[1], times_convertedA < timerange[1]))
                    times = np.concatenate((times1[0], times2[0]))
                else:
                    #Trim with timerange:
                    times = np.where(np.logical_and(times_convertedA > timerange[0], times_convertedA < timerange[1]))
                    #print('nt', len(times[0]))
                #Intersection of the prior two conditions, plus the timerange
                allinds = np.intersect1d(allinds, times)

                first+=1

            times_convertedA = times_convertedA[allinds]
            times_convertedB = times_convertedB[allinds]
            countrateA = countrateA[allinds]            
            countrateB = countrateB[allinds]  
            lvtA = lvtA[allinds]
            lvtB = lvtB[allinds]


            # print(countrateB)
            
            # if bool(excluded_range): 
            #     times1 = np.where(np.logical_and(times_convertedA > timerange[0], times_convertedA < excluded_range[0]))
            #     times2 = np.where(np.logical_and(times_convertedA > excluded_range[1], times_convertedA < timerange[1]))
            #     times = np.concatenate((times1[0], times2[0]))
            # else:
            #     #Trim with timerange:
            #     times = np.where(np.logical_and(times_convertedA > timerange[0], times_convertedA < timerange[1]))


            # times_convertedA = times_convertedA[times]
            # countrateA = countrateA[times]
            # lvtA = lvtA[times]
            # times_convertedB = times_convertedB[times]
            # countrateB = countrateB[times]
            # lvtB = lvtB[times]

            std = np.std(countrateA[np.isfinite(countrateA)])
            mean_val = np.mean(countrateA[np.isfinite(countrateA)])
            means = np.full(len(countrateA), mean_val)

            #tstep=(times_convertedA[1]-times_convertedA[0]).seconds
            #print('NuSTAR Timestep (s): ', tstep)

            n_bx = round(smooth/tstep)
            
            arr_lc = np.array(countrateA)
            avg_lc = boxcar_average(arr_lc, n_bx)
            avg_lc[0:3]=arr_lc[0:3]
            avg_lc[-3:]=arr_lc[-3:]
            
            #outs = np.where(np.logical_or(avg_lc < mean_val-std, avg_lc > mean_val+std), 1, 0)
            outs_up = np.where(avg_lc > mean_val+std, 1, 0)
            outs_down = np.where(avg_lc < mean_val-std, -1, 0)
            outs = outs_up+outs_down
            
            
            n_all_outupsA.append(outs_up)
            n_all_outdownsA.append(outs_down)
            
            
            if plot_each:
                ax=axes[nn]
                ax.plot(times_convertedA, countrateA, 
                         label='NuSTAR FPMA Counts '+str(erange[0])+' to '+str(erange[1])+' keV',
                         **default_kwargs, color=clrs[ind])

                ax.plot(times_convertedA, avg_lc, color=clrs[ind])
                ax.errorbar(times_convertedA, means, yerr=std, color='Black')
                ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
                ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=10))
                ax.legend()

                nn+=1
                
            if plot_logic:
            
                ax = axes[-1]
                ax.plot(times_convertedA, outs)
                ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
                ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=10))


            std = np.std(countrateB[np.isfinite(countrateB)])
            mean_val = np.mean(countrateB[np.isfinite(countrateB)])
            means = np.full(len(countrateB), mean_val)
            
            #tstep=(times_convertedB[1]-times_convertedB[0]).seconds
            #print('NuSTAR Timestep (s): ', tstep)
            #n_bx = round(smooth/tstep)

            arr_lc = np.array(countrateB)
            avg_lc = boxcar_average(arr_lc, n_bx)
            avg_lc[0:3]=arr_lc[0:3]
            avg_lc[-3:]=arr_lc[-3:]
            
            #outs = np.where(np.logical_or(avg_lc < mean_val-std, avg_lc > mean_val+std), 1, 0)
            outs_up = np.where(avg_lc > mean_val+std, 1, 0)
            outs_down = np.where(avg_lc < mean_val-std, -1, 0)
            outs = outs_up+outs_down
            
            n_all_outupsB.append(outs_up)
            n_all_outdownsB.append(outs_down)

            if plot_each:
                ax=axes[nn]

                ax.plot(times_convertedB, countrateB, 
                         label='NuSTAR FPMB Counts '+str(erange[0])+' to '+str(erange[1])+' keV', 
                         **default_kwargs, color=clrs[ind+1])
                ax.plot(times_convertedB, avg_lc, color=clrs[ind+1])
                ax.errorbar(times_convertedB, means, yerr=std, color='Black')
                ind+=2

                ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
                ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=10))

                ax.legend()

                nn+=1
                
            if plot_logic:

                ax = axes[-1]
                ax.plot(times_convertedA, outs)
                ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
                ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=10))

        
        
        
    
    if aia_inputs or fexviii:
        instrument='AIA'
        #ADDS ALL FEXVIII COMPONENTS IF FEXVIII IS TRUE
        if fexviii:
            if 94 in aia_inputs and 171 in aia_inputs and 211 in aia_inputs:
                wavelengths=aia_inputs
            else:
                wavelengths=aia_inputs.copy()
                wavelengths.append(94)
                wavelengths.append(171)
                wavelengths.append(211)
                wavelengths = list(set(wavelengths))
                
        all_aia=load_lightcurves(instrument, wavelengths=wavelengths, fexviii=fexviii)
        
        if fexviii:
            w=94
            chanlabel='AIA'+str(w)

            times_converted = np.array(all_aia['times_'+chanlabel])
            data_total = np.array(all_aia['data_total_fexviii'])
            exp_time = all_aia['exp_time_'+chanlabel]
            exp_times = np.array([e.value for e in exp_time])
            
            #Getting rid of "0 exposure" times
            fintot = np.where(exp_times > 0.)
            
            times_converted = times_converted[fintot]
            data_total = data_total[fintot]
            exp_times = exp_times[fintot]
            
            if remove_fexviii_max:
                #Getting rid of "0 exposure" times
                dmax = np.where(data_total != max(data_total))

                times_converted = times_converted[dmax]
                data_total = data_total[dmax]
                exp_times = exp_times[dmax]
     
            
            #Trim with timerange:
            times = np.where(np.logical_and(times_converted > timerange[0], times_converted < timerange[1]))

            times_converted = times_converted[times]
            data_total = data_total[times]
            exp_times = exp_times[times]
            
            
            #Save THIS times array for use when we plot fexviii deviations in summary (will have different
            #size than other AIA if we have removed the max value).
            times_converted_ = times_converted
            

            corr_totals= data_total/(exp_times)
            

            std = np.std(corr_totals)
            mean_val = np.mean(corr_totals)
            means = np.full(len(corr_totals), mean_val)
            
            tstep=(times_converted[1]-times_converted[0]).seconds
            #print('AIA Timestep (s): ', tstep)

            n_bx = round(smooth/tstep)

            arr_lc = np.array(corr_totals)
            avg_lc = boxcar_average(arr_lc, n_bx)
            avg_lc[0:3]=arr_lc[0:3]
            avg_lc[-3:]=arr_lc[-3:]
            
            #outs = np.where(np.logical_or(avg_lc < mean_val-std, avg_lc > mean_val+std), 1, 0)
            fexviii_outs_up = np.where(avg_lc > mean_val+std, 1, 0)
            fexviii_outs_down = np.where(avg_lc < mean_val-std, -1, 0)
            outs = fexviii_outs_up+fexviii_outs_down
            
            if plot_each:
                ax=axes[nn]
                ax.plot(times_converted, avg_lc, label='Fe-XVIII boxcar', color=clrs[ind])
                ax.plot(times_converted, corr_totals, label='Fe-XVIII', **default_kwargs, color=clrs[ind])
                ax.errorbar(times_converted, means, yerr=std, color='Black', label='1 stdv from mean')
                
                ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
                ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=10))
                ax.legend()

                nn+=1
                
            if plot_logic:

                ax=axes[-1]
                ax.plot(times_converted, outs)
                ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
                ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=10))
            
        
        if aia_inputs:
            all_outups =[]
            all_outdowns=[]
            all_times=[]
            #print(all_aia.keys())
            for k in range(0, len(aia_inputs)):
                w = aia_inputs[k]
                
                
                chanlabel='AIA'+str(w)

                times_converted = np.array(all_aia['times_'+chanlabel])
                data_total = np.array(all_aia['data_total_'+chanlabel])
                exp_time = all_aia['exp_time_'+chanlabel]
                exp_times = np.array([e.value for e in exp_time])
                
                #print('AIA '+str(w)+' First Times:', times_converted[0:2])
                
                #Getting rid of "0 exposure" times
                fintot = np.where(exp_times > 0.)
            
                times_converted = times_converted[fintot]
                data_total = data_total[fintot]
                exp_times = exp_times[fintot]

                #Trim with timerange:
                times = np.where(np.logical_and(times_converted > timerange[0], times_converted < timerange[1]))

                times_converted = times_converted[times]
                data_total = data_total[times]
                exp_times = exp_times[times]
                
                
                #This whole section is a response to the highly annoying problem that AIA files are from slightly 
                #different times and thus sometimes the lightcurves have different lengths within the same time 
                #interval. Later, we want to sum the above/below arrays to find times when multiple channels are 
                #anomalous at once, so we need to fix this.
                
                #Setting up a "reference length" from the first wavelength used
                if k==0:
                    reflen = len(times_converted)-1
                    #print(reflen)
                    times_converted = times_converted[0:-1]
                    data_total = data_total[0:-1]
                    exp_times = exp_times[0:-1]
                    #print(len(exp_times))
                    refstart = times_converted[0]
                    refend = times_converted[-1]
                    #print(refstart, refend)
                    
                #For all other wavelengths, trimming to match reference length (in a smart-ish way to retain coverage
                #of more of the same time interval). 
                if k > 0:
                    vlen = len(times_converted)
                    while vlen > reflen:
                        start = times_converted[0]
                        end = times_converted[-1]
                        
                        #if start/end values outside reference window, trim until
                        #array is the same length as the reference array
                        if start < refstart and vlen > reflen:
                            times_converted = times_converted[1:]
                            data_total = data_total[1:]
                            exp_times = exp_times[1:]
                            start = times_converted[0]
                            vlen = len(times_converted)
                            
                        if end > refend and vlen > reflen:
                            times_converted = times_converted[0:-1]
                            data_total = data_total[0:-1]
                            exp_times = exp_times[0:-1]
                            end = times_converted[-1]
                            vlen = len(times_converted)
                          
                    
                    

                corr_totals= data_total/(exp_times)

                std = np.std(corr_totals)
                mean_val = np.mean(corr_totals)
                
                tstep=(times_converted[1]-times_converted[0]).seconds
                #print('AIA Timestep (s): ', tstep)

                n_bx = round(smooth/tstep)


                arr_lc = np.array(corr_totals)
                avg_lc = boxcar_average(arr_lc, n_bx)
                avg_lc[0:3]=arr_lc[0:3]
                avg_lc[-3:]=arr_lc[-3:]

                means = np.full(len(corr_totals), mean_val)

                #outs = np.where(np.logical_or(avg_lc < mean_val-std, avg_lc > mean_val+std), 1, 0)
                outs_up = np.where(avg_lc > mean_val+std, 1, 0)
                outs_down = np.where(avg_lc < mean_val-std, -1, 0)
                outs = outs_up+outs_down

                all_outups.append(outs_up)
                all_outdowns.append(outs_down)
                all_times.append(times_converted)

                if plot_each:

                    ax=axes[nn]

                    ax.plot(times_converted, avg_lc, label=chanlabel+' boxcar', color=clrs[ind])
                    ax.plot(times_converted, corr_totals, label=chanlabel, **default_kwargs, color=clrs[ind])
                    ax.errorbar(times_converted, means, yerr=std, color='Black', label='1 stdv from mean')

                    ax.legend()

                    nn+=1
                    
                if plot_logic:

                    ax=axes[-1]
                    ax.plot(times_converted, outs)


                ind+=1
            
            
    if goes_inputs:
        instrument='GOES'
        data = load_lightcurves(instrument)
        
        ylabel = data['GOES flux label']
        goestimes = data['GOES Times']
        gts = np.array([t.datetime for t in goestimes])
        
        #Trim with timerange:
        times = np.where(np.logical_and(gts > timerange[0], gts < timerange[1]))
        
        gts = gts[times]
        
        goesoutups=[]
        goesoutdowns=[]
        
        if 'xrsa' in goes_inputs:
            xrsalabel = data['XRSA label']
            xrsacounts = data['XRSA counts']
            xrsacounts = xrsacounts[times]
            xrsaclr=0
            
            std = np.std(xrsacounts).value
            mean_val = np.mean(xrsacounts).value
            means = np.full(len(xrsacounts), mean_val)
            
            tstep=(gts[1]-gts[0]).seconds
            #print('GOES Timestep (s): ', tstep)

            n_bx = round(smooth/tstep)

            arr_lc = np.array(xrsacounts)
            avg_lc = boxcar_average(arr_lc, n_bx)
            avg_lc[0:3]=arr_lc[0:3]
            avg_lc[-3:]=arr_lc[-3:]
            
            #outs = np.where(np.logical_or(avg_lc < mean_val-std, avg_lc > mean_val+std), 1, 0)
            outs_up = np.where(avg_lc > mean_val+std, 1, 0)
            outs_down = np.where(avg_lc < mean_val-std, -1, 0)
            outs = outs_up+outs_down
            
            goesoutups.append(outs_up)
            goesoutdowns.append(outs_down)
            
            if plot_each:
                #FIGURE A
                ax=axes[nn]

                ax.semilogy(gts, xrsacounts, 
                             label=xrsalabel,
                             **default_kwargs, color=clrs[xrsaclr])

                ax.plot(gts, avg_lc, color=clrs[xrsaclr], label='XRSA Boxcar')
                #plt.errorbar(gts, means, yerr=std, color='Black', label='1 stdv from mean')
                ax.axhline(mean_val, color='Black', label='1 stdv from mean')
                ax.axhline(mean_val+std, color='Black')
                ax.axhline(mean_val-std, color='Black')
                ax.legend()

                nn+=1
                
            if plot_logic:
                
                ax=axes[-1]
                ax.plot(gts, outs)

            
        if 'xrsb' in goes_inputs:
            xrsblabel = data['XRSB label']
            xrsbcounts = data['XRSB counts']
            xrsbcounts = xrsbcounts[times]
            xrsbclr=0
            
            std = np.std(xrsbcounts).value
            mean_val = np.mean(xrsbcounts).value
            means = np.full(len(xrsbcounts), mean_val)
            
            tstep=(gts[1]-gts[0]).seconds
            #print('GOES Timestep (s): ', tstep)

            n_bx = round(smooth/tstep)
            
            arr_lc = np.array(xrsbcounts)
            avg_lc = boxcar_average(arr_lc, n_bx)
            avg_lc[0:3]=arr_lc[0:3]
            avg_lc[-3:]=arr_lc[-3:]
            
            #outs = np.where(np.logical_or(avg_lc < mean_val-std, avg_lc > mean_val+std), 1, 0)
            outs_up = np.where(avg_lc > mean_val+std, 1, 0)
            outs_down = np.where(avg_lc < mean_val-std, -1, 0)
            outs = outs_up+outs_down
            
            goesoutups.append(outs_up)
            goesoutdowns.append(outs_down)
            
            if plot_each:
                #FIGURE B
                ax=axes[nn]

                ax.semilogy(gts, xrsbcounts, 
                             label=xrsblabel, 
                             **default_kwargs, color=clrs[xrsbclr])

                ax.plot(gts, avg_lc, color=clrs[xrsbclr])
                ax.axhline(mean_val, color='Black', label='1 stdv from mean')
                ax.axhline(mean_val+std, color='Black')
                ax.axhline(mean_val-std, color='Black')
                ax.legend()

                nn+=1
                
            if plot_logic:
                
                ax=axes[-1]
                ax.plot(gts, outs)


    plt.savefig(save_dir+'/'+savestring+'_stdv.png')
    if not show:
        plt.close()
            
    fig = plt.figure(figsize=(15, 5))
    
    lsty='dashed'
    
    aia_res=[]
    fexviii_res=[]
    nustar_res=[]
    goes_res=[]

    interp=False
    
    if goes_inputs:
        totalgoesoutups = np.sum(goesoutups, axis=0)
        totalgoesoutdowns = np.sum(goesoutdowns, axis=0)
        plt.plot(gts, totalgoesoutups, label='GOES Above 1-Sigma', color='green', linestyle=lsty)
        plt.plot(gts, totalgoesoutdowns, label='GOES Below 1-Sigma', color='forestgreen', linestyle=lsty)
        
        goes_res = {'out_ups': goesoutups, 'out_downs': goesoutdowns, 'times': gts, 'channels': goes_inputs}
        
        #If using GOES, interpolate the other result arrays to use the GOES timebins
        interp=True
        GTS = [g.timestamp() for g in list(gts)]
        interp_sum = totalgoesoutups
        interp_string='Sum of: GOES, '
        
        #print('GOES Time Step:', gts[1]-gts[0])

    
    if aia_inputs:
        totaloutups = np.sum(all_outups, axis=0)
        totaloutdowns = np.sum(all_outdowns, axis=0)
        
        if interp==False:
            plt.plot(times_converted, totaloutups, label='AIA Above 1-Sigma', color='orange', linestyle=lsty)
            plt.plot(times_converted, totaloutdowns, label='AIA Below 1-Sigma', color='brown', linestyle=lsty)
            
        if interp==True:
            TC = [g.timestamp() for g in list(times_converted)]
            new_aia_up_arr = np.interp(GTS,TC,totaloutups)
            new_aia_down_arr = np.interp(GTS,TC,totaloutdowns)
            plt.plot(gts, new_aia_up_arr, label='AIA Above 1-Sigma', color='orange', linestyle=lsty)
            plt.plot(gts, new_aia_down_arr, label='AIA Below 1-Sigma', color='brown', linestyle=lsty) 
            
            interp_sum = interp_sum + new_aia_up_arr
            interp_string=interp_string+'AIA, '
        
        aia_res={'out_ups': all_outups, 'out_downs': all_outdowns, 'times': all_times, 'waves': aia_inputs}
        
        #print('AIA '+str(aia_inputs[-1])+' Time Step:', times_converted[1]-times_converted[0])
        
    if fexviii:
        if interp==False:
            plt.plot(times_converted_, fexviii_outs_up, label='FeXVIII Above 1-Sigma', color='red', linestyle=lsty)
            plt.plot(times_converted_, fexviii_outs_down, label='FeXVIII Below 1-Sigma', color='indianred', linestyle=lsty)
        
        if interp==True:
            TC = [g.timestamp() for g in list(times_converted_)]
            new_fe_up_arr = np.interp(GTS,TC,fexviii_outs_up)
            new_fe_down_arr = np.interp(GTS,TC,fexviii_outs_down)
            plt.plot(gts, new_fe_up_arr, label='FeXVIII Above 1-Sigma', color='red', linestyle=lsty)
            plt.plot(gts, new_fe_down_arr, label='FeXVIII Below 1-Sigma', color='indianred', linestyle=lsty) 
            
            interp_sum = interp_sum + new_fe_up_arr
            interp_string=interp_string+'Fe-XVIII, '
        
        fexviii_res={'out_ups': fexviii_outs_up, 'out_downs': fexviii_outs_down, 'times': times_converted_}
        
        #print('AIA 94 Time Step:', times_converted_[1]-times_converted_[0])
    
    if nustar_inputs:
        n_totaloutupsA = np.sum(n_all_outupsA, axis=0)
        n_totaloutdownsA = np.sum(n_all_outdownsA, axis=0)
        n_totaloutupsB = np.sum(n_all_outupsB, axis=0)
        n_totaloutdownsB = np.sum(n_all_outdownsB, axis=0)
        
        if interp==False:
            plt.plot(times_convertedA, n_totaloutupsA, label='NuSTAR A Above 1-Sigma', color='dodgerblue', linestyle=lsty)
            plt.plot(times_convertedA, n_totaloutdownsA, label='NuSTAR A Below 1-Sigma', color='cornflowerblue',
                     linestyle=lsty)
            plt.plot(times_convertedA, n_totaloutupsB, label='NuSTAR B Above 1-Sigma', color='blue', linestyle=lsty)
            plt.plot(times_convertedA, n_totaloutdownsB, label='NuSTAR B Below 1-Sigma', color='darkblue', linestyle=lsty)
            
        if interp==True:   
            TC = [g.timestamp() for g in list(times_convertedA)]
            new_A_up_arr = np.interp(GTS,TC,n_totaloutupsA)
            new_A_down_arr = np.interp(GTS,TC,n_totaloutdownsA)
            plt.plot(gts, new_A_up_arr, label='NuSTAR A Above 1-Sigma', color='dodgerblue', linestyle=lsty)
            plt.plot(gts, new_A_down_arr, label='NuSTAR A Below 1-Sigma', color='cornflowerblue',
                     linestyle=lsty)
            TC = [g.timestamp() for g in list(times_convertedB)]
            new_B_up_arr = np.interp(GTS,TC,n_totaloutupsB)
            new_B_down_arr = np.interp(GTS,TC,n_totaloutdownsB)    
            plt.plot(gts, new_B_up_arr, label='NuSTAR B Above 1-Sigma', color='blue', linestyle=lsty)
            plt.plot(gts, new_B_down_arr, label='NuSTAR B Below 1-Sigma', color='darkblue', linestyle=lsty)
            
            interp_sum = interp_sum + new_A_up_arr
            interp_sum = interp_sum + new_B_up_arr
            interp_string=interp_string+'NuSTAR (A+B)'
            
        
        nustar_res={'out_upsA': n_all_outupsA, 'out_downsA': n_all_outdownsA, 'timesA': times_convertedA, 
                        'out_upsB': n_all_outupsB, 'out_downsB': n_all_outdownsB, 'timesB': times_convertedB,
                       'eranges': nustar_inputs}
        
        #print('NuSTAR Time Step:', times_convertedA[1]-times_convertedA[0])
    
    
#         if interp:
#             plt.plot(gts, interp_sum, label='Interpolated Sum', color='Black')
  
    
    plt.legend(ncol=3)
    plt.xlim(timerange)
    plt.savefig(save_dir+'/'+savestring+'_quiescence_summary.png')
    if not show:
        plt.close()
    

    if analyze_transients:
        res = transient_analysis(aia_res, goes_res, nustar_res, fexviii_res, timerange, transient_number=transient_number,
                                 show=show, save_dir=save_dir, savestring=savestring)
    
    
    return res

def transient_analysis(aia_res, goes_res, nustar_res, fexviii_res, timerange, interp=True,
                      transient_number=2, show=True, save_dir='./', savestring='test'):
    """
    Takes in arrays corresponding to times where each instrument is above/below the mean+-stdv window during
    the observation, and quantifies intervals where there are multiple instruments above/below. 
    
    Transient number: minimum number of channels which must be above mean+stdv for the time interval to
    count as a transient.
    
    """
    
    #================================================================================================
    #For each instrument with an existing input, make a dictionary containing IDs for each sub-instrument
    #channel, as well as arrays of time windows where lightcurve is above and below the mean+-stdv levels.
    #================================================================================================
    
    aia_windows=[]
    goes_windows=[]
    nu_windows=[]
    fexviii_windows=[]
    
    dicts=[]
    
    if aia_res:
        outups=np.array(aia_res['out_ups'])
        outdowns=np.array(aia_res['out_downs'])
        times=np.array(aia_res['times'])
        
        aid=[]
        for aa in aia_res['waves']:
            aid.append(str(aa))
        
        aia_windows={'ID': aid}
        
        for i in range(0, len(aia_res['waves'])):
            w = aia_res['waves'][i]            
            uwindows, dwindows = windows(times[i,:], outups[i,:], outdowns[i,:])
            aia_windows[str(w)+'_uwins'] = uwindows
            aia_windows[str(w)+'_dwins'] = dwindows
            
            
        dicts.append(aia_windows)
            
        
    if goes_res:
        outups=np.array(goes_res['out_ups'])
        outdowns=np.array(goes_res['out_downs'])
        times=np.array(goes_res['times'])
        
        goes_windows={'ID': goes_res['channels']}
        
        for i in range(0, len(goes_res['channels'])):
            ch = goes_res['channels'][i]
            uwindows, dwindows = windows(times, outups[i,:], outdowns[i,:])
            goes_windows[ch+'_uwins'] = uwindows
            goes_windows[ch+'_dwins'] = dwindows

        dicts.append(goes_windows)
        

        
        
    if nustar_res:
        
        nid=[]
        for er in nustar_res['eranges']:
            nid.append(str(er[0])+'_'+str(er[1])+'A')
            nid.append(str(er[0])+'_'+str(er[1])+'B')
            
        nu_windows={'ID': nid}
        
        #FPMA
        outups=np.array(nustar_res['out_upsA'])
        outdowns=np.array(nustar_res['out_downsA'])
        times=np.array(nustar_res['timesA'])
        
        for i in range(0, len(nustar_res['eranges'])):
            er = nustar_res['eranges'][i]
            uwindows, dwindows = windows(times, outups[i,:], outdowns[i,:])
            nu_windows[str(er[0])+'_'+str(er[1])+'A'+'_uwins'] = uwindows
            nu_windows[str(er[0])+'_'+str(er[1])+'A'+'_dwins'] = dwindows
            
        #FPMB
        outups=np.array(nustar_res['out_upsB'])
        outdowns=np.array(nustar_res['out_downsB'])
        times=np.array(nustar_res['timesB'])
        
        for i in range(0, len(nustar_res['eranges'])):
            er = nustar_res['eranges'][i]
            uwindows, dwindows = windows(times, outups[i,:], outdowns[i,:])
            nu_windows[str(er[0])+'_'+str(er[1])+'B'+'_uwins'] = uwindows
            nu_windows[str(er[0])+'_'+str(er[1])+'B'+'_dwins'] = dwindows
        
        dicts.append(nu_windows)
        
    if fexviii_res:
        
        outups=np.array(fexviii_res['out_ups'])
        outdowns=np.array(fexviii_res['out_downs'])
        times=np.array(fexviii_res['times'])
        
        fexviii_windows={'ID': ['Fe-XVIII']}
        uwindows, dwindows = windows(times, outups, outdowns)
        fexviii_windows['Fe-XVIII_uwins'] = uwindows
        fexviii_windows['Fe-XVIII_dwins'] = dwindows
        
        dicts.append(fexviii_windows)
        
        
    #================================================================================================
    #================================================================================================

    #Make big lists of all up + down windows, and a big list of all the sub-instrument IDs
    
    uwins = []
    dwins = []
    IDs = []
    for d in dicts:
        ids=d['ID']
        for i_d in ids:
            uwins.append(d[i_d+'_uwins'])
            dwins.append(d[i_d+'_dwins'])
            IDs.append(i_d)
    count=0
    each_window_has=[]
    each_window=[]
    #For every sub-instrument list of windows...
    for i in range(0, len(uwins)):
        wins=uwins[i]
        #For every window in that list...
        for j in range(0, len(wins)):
            win=wins[j]
            this_window_has=[IDs[i]+'-'+str(j)]
            thiswin=win
            #Want to check if this window intersects with the windows from any other sub-instruments.
            
            #With replacement (gives repeats):
            others = uwins[:i] + uwins[i+1:]
            oIDs = IDs[:i] + IDs[i+1:]
            
            #Without replacement
            #others = uwins[i+1:]
            #oIDs = IDs[i+1:]
            
            #For every OTHER sub-instrument list of windows...
            for o in range(0, len(others)):
                oo = others[o]
                #For every window in that list...
                for k in range(0, len(oo)):
                    kk = oo[k]
                    #If the other window starts or ends in the middle of the window of interest.
                    if win[0] < kk[0] < win[1] or win[0] < kk[1] < win[1]:
                        allwin=thiswin+kk
                        #Update the size of the window to now include the full extent defined by
                        #both win and kk 
                        thiswin=[min(allwin), max(allwin)]
                        #print(thiswin)

                        this_window_has.append(oIDs[o]+'-'+str(k))
                        count+=1
            
            #print('Every sucessive overlap with '+IDs[i]+'-'+str(j)+': ', this_window_has)
            #print('Range: ', thiswin[0].strftime('%H:%M:%S'), thiswin[1].strftime('%H:%M:%S'))
            each_window_has.append(this_window_has)
            each_window.append(thiswin)
            #print('')
    
    #print(each_window_has)
    #print('Overlaps:', count)
    print('')
    
    #We know that some of the windows may have have overlap with eachother. We want to simplify to combine these. 
    newwindows=[]
    newwindow_labels=[]
    for i in range(0, len(each_window)):
        ew = each_window[i]
        
        flip=0
        for nw in newwindows:
            #print(nw)
            #print(ew)
            if nw[0] < ew[0] < nw[1] or nw[0] < ew[1] < nw[1]:
                flip=1
                
        if flip==1:
            #print('stopping 1')
            continue
        
        #Define lists of all OTHER windows
        others =  each_window[:i] + each_window[i+1:]
        oh =  each_window_has[:i] + each_window_has[i+1:]
        #I
        allin = ew
        allstr = each_window_has[i]
        for o in range(0, len(others)):
            oo=others[o]
            if oo[0] < ew[0] < oo[1] or oo[0] < ew[1] < oo[1]:
                #print(each_window_has[i])
                #print(oh[o])
                #print('')
                allin.extend(oo)
                allstr.extend(oh[o])

        #print(list(set(allstr)))
        #print(allstr)
        newwin = [min(allin), max(allin)]
        #print('New win: ', newwin)
        #print('New win: ', newwin[0].strftime('%H:%M:%S'), newwin[1].strftime('%H:%M:%S'))
        flip=0
        for nw in newwindows:
            #print(nw[0].strftime('%H:%M:%S'), nw[1].strftime('%H:%M:%S'))
            #print(ew)
            if nw[0] <= newwin[0] < nw[1] or nw[0] < newwin[1] <= nw[1]:
                #print('not adding')
                flip=1
                
        if flip==1:
            #print('stopping 2')
            continue
                        

        #print(allstr)
        newwindows.append(newwin)
        newwindow_labels.append(list(set(allstr)))
        #print('')
        
    print('Number of Windows (>= 1 channel): ', len(newwindows))
    
    #Removing windows with less than our required number of transients
    if transient_number > 1:
        counts = [len(nl) for nl in newwindow_labels]
        above_num = np.where(np.array(counts) >= transient_number)
        newwindows = [newwindows[i] for i in above_num[0]]
        newwindow_labels = [newwindow_labels[i] for i in above_num[0]]
        
    print('Number of Windows (>= '+str(transient_number)+' channels): ', len(newwindows))
    
    newwindows_sort = sorted(newwindows)
    #print(newwindows_sort)
    
    firsts = [w[0] for w in newwindows]
    inds = np.argsort(np.array(firsts))
    newwindows = [newwindows[i] for i in inds]
    newwindow_labels = [newwindow_labels[i] for i in inds]
    
    # for i in range(0, len(newwindows_sort)):
    #     print('')
    #     print('Window ', i)
    #     print('Channels: ', sorted( newwindow_labels[i]))
    #     print('Time Range: ', newwindows[i][0].strftime('%H:%M:%S'), newwindows[i][1].strftime('%H:%M:%S'))

    fig, ax = plt.subplots(1, figsize=(15, 3))
    
    ind=0
    clrs=make_colors(len(newwindows))
    for ew in newwindows:
        #print(newwindow_labels[ind])
        plt.axvline(ew[0], label=str(len(newwindow_labels[ind]))+' total', color=clrs[ind])
        plt.axvline(ew[1], color=clrs[ind])
        ind+=1
        
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    #ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))
        
    
    plt.legend(ncol=3)
    plt.xlim(timerange)
    plt.savefig(save_dir+'/'+savestring+'_transients.png')
    if not show:
        plt.close()
        
        
    return newwindows

def windows(time, outup, outdown):
    """
    Take arrays of times, above-stdv intervals + below-stdv intervals + returns two lists
    of start-stop times for the windows where lightcurve is above and below.
   
    """
    ai = np.where(outup == 1)[0]
    ca = consecutive(ai)
    awindows=[]
    for c in ca:
        if len(c) > 0:
            c=np.array(c)
            awindows.append([time[c[0]], time[c[-1]]])
    
    bi = np.where(outdown == -1)[0]
    cb = consecutive(bi)
    bwindows=[]
    for c in cb:
        if len(c) > 0:
            c=np.array(c)
            bwindows.append([time[c[0]], time[c[-1]]])
    
    return awindows, bwindows
    
    
def consecutive(data, stepsize=1):
    #https://stackoverflow.com/questions/7352684/how-to-find-the-groups-of-consecutive-elements-in-a-numpy-array
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)  
    


def boxcar_average(
    arr: np.ndarray,
    N: int,
    insert_val: int | float = np.nan
) -> np.ndarray:
    """
    Perform a boxcar average (AKA moving average) on the provided data.
    
    Parameters
    ----------
    arr : np.ndarray
        The input array on which the boxcar average will be performed.
    N : int
        The boxcar width.
    insert_val : some quantity
        This quantity is padded on the left and right sides of the
        averaged data so that the size of the output array
        matches the size of the input array.
    
    Returns
    -------
    bc : np.ndarray
        The array containing the averaged data.
    """

    # print(f'Applying boxcar average. Provided N={N}')
    # N = adjust_n(len(arr), N)
    if N > len(arr):
        raise ValueError(f'Provided N={N} greater than the '\
            f'number of available data points, {len(arr)}')
    
    bc = np.convolve(arr, np.ones(N)/N, mode='same')
    #for _ in range((N-1)//2):
    #    bc = np.insert(bc, 0, insert_val)
    #    bc = np.insert(bc, len(bc), insert_val)

    return bc

def make_colors(number):
    """
    Makes a list of colors depending on which instruments are being used (so color-channel correspondence remains 
    consistent across DEMs using different instruemnt combinations.
    
    Color table assumes set numbers of inputs per included instrument:
    
    six AIA channels
    two XRT filters
    three NuSTAR energy ranges
    ten EIS lines.
    
    """
    #normal
    aiaclrs=['darkgreen','darkcyan','gold','sienna','indianred','darkorange']
    xrtclrs=['darkslateblue', 'dodgerblue', 'cornflowerblue']
    nuclrs=['purple', 'mediumpurple', 'plum']
    eisclrs=['seagreen', 'mediumseagreen', 'springgreen', 'green', 'mediumspringgreen', 'mediumaquamarine', 
          'aquamarine', 'turquoise', 'lightseagreen', 'mediumturquoise', 'lawngreen', 'cadetblue', 
             'slategray', 'darkslateblue']
    
    allcolors = aiaclrs
    allcolors.extend(xrtclrs)
    allcolors.extend(nuclrs)
    allcolors.extend(eisclrs)
    
    if number > 26:
        print('You are asking for too many colors, will now repeat my max list of 26')
        while len(allcolors) < number:
            allcolors.extend(allcolors)
        return allcolors
    else:
        return allcolors[0:number]
      
    return allcolors
