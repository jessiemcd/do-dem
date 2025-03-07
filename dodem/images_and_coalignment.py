"""
Making images, doing coalignment, etc. 
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

#dodem code
import region_fitting as rf
import nustar_dem_prep as nu
import nustar_utilities as nuutil
import nustar_pysolar as nustar


import shutil
import glob
from astropy.io import fits
import matplotlib.colors as colors
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from regions import CircleSkyRegion
from sunpy.coordinates import Helioprojective

import copy
import pickle
import os



def make_region(regiondict, map):
    """
    wrapper - make circular region out of dictionary

    """

    regcenter=SkyCoord(regiondict['centerx'],regiondict['centery'], frame=map.coordinate_frame)

    region = CircleSkyRegion(
            center = regcenter,
            radius = regiondict['radius'] * u.arcsec
        )

    return region


def regfile_to_regdict(regfile, time):

    """
    Read in region file (circular region) and return a dictionary of its contents.
    """

    offset, rad = rf.read_regfile(regfile, time[0], time[1], regRAunit='hourangle')

    regiondict = {'radius': rad.value,
              'centerx': offset[0],
              'centery': offset[1]}

    return regiondict




def check_region_emission(all_time_intervals, working_dir, grade='0_4', plot=True, efilter=[],
                         keep_problem_intervals=True):

    """
    Look at the evolution of the % of the whole-FOV emission that is in the chosen region.
    Input: all_time_intervals, i.e. output of time_interval_selection.find_all_intervals()
    or similar. 
    """


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





def plot_rectangle(coordx, coordy, width, height, angle, map_):

    """
    With the output this gives you, you can draw a rectangle like:
    
    map_.draw_quadrangle(rectangle,
            axes=ax, edgecolor="red", linestyle="--", linewidth=2)

    """
    
    rotation_angle = angle * u.deg
    center_coord = SkyCoord(coordx * u.arcsec, coordy * u.arcsec, frame=map_.coordinate_frame)
    width_ = width * u.arcsec
    height_ = height * u.arcsec
    offset_frame = SkyOffsetFrame(origin=center_coord, rotation=rotation_angle)
    rectangle = SkyCoord(lon=[-1/2, 1/2] * width_, lat=[-1/2, 1/2] * height_, frame=offset_frame)

    return rectangle




def get_orbit_aiamaps(aia_dir, id_dirs, wave=94):

    import sunpy.map
    from aiapy.calibrate.util import get_correction_table, get_pointing_table
    from aiapy.calibrate import register, update_pointing, degradation, estimate_error

    aiamaps = []
    for id in id_dirs:
        print(id)
        evt_data, hdr = nu.return_submap(datapath=id, fpm='A', return_evt_hdr=True)
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



def aia_search_fido(midtime):
    
    
    from aiapy.calibrate.util import get_correction_table, get_pointing_table
    from aiapy.calibrate import register, update_pointing, degradation, estimate_error
    from sunpy.net import Fido
    from sunpy.net import attrs as a

    query = Fido.search(
            a.Instrument.aia,
            a.Physobs.intensity,
            a.Wavelength(94*u.angstrom),
            a.Time(midtime-12*u.s, midtime))
            #a.Time(time_range[0],time_range[1]))#,
            #a.Sample(sample_every)
            #)
    print(query)
    files = Fido.fetch(query, max_conn=1)
    print(files)
    try:
        amap=sunpy.map.Map(files)
        #If things keep failing on the pointing table line (and the JSOC is not down), 
        #it may be an issue with the actual AIA map – try another time.
        ptab = get_pointing_table(amap.date - 12 * u.h, amap.date + 12 * u.h)
    except AttributeError:
        amap=sunpy.map.Map(files[0])
        #If things keep failing on the pointing table line (and the JSOC is not down), 
        #it may be an issue with the actual AIA map – try another time.
        ptab = get_pointing_table(amap.date - 12 * u.h, amap.date + 12 * u.h)
        
    try:
        m_temp = update_pointing(amap, pointing_table=ptab)
    except TypeError:
        amap.meta.pop('crpix1')
        amap.meta.pop('crpix2')
        print('CRPIX issue on ', files)
        m_temp = update_pointing(amap, pointing_table=ptab)
        
    #converts lev1 map to lev1.5 map 
    #(see: https://aiapy.readthedocs.io/en/stable/api/aiapy.calibrate.register.html?highlight=register)
    m = register(m_temp)

    return m


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
    import time_interval_selection as tis

    if method=='input' or method =='double':
        #Find the top level directory for each region
        if method=='input':
            region_dirs = find_region_dirs(working_dir)
        if method=='double':
            region_dirs = find_direction_dirs(working_dir, sep_axis)
            #print(region_dirs)
        all_all_time_intervals, fixit = tis.region_time_intervals(region_dirs, id_dirs, shush=shush)

        #return all_all_time_intervals

        if len(all_all_time_intervals) > 1:
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
        else:
            time_interval = all_all_time_intervals[0][0][0]
            region_dir = region_dirs[0]

    #print(time_interval[0], time_interval[1])

    if method=='fit':
        all_time_intervals, all_time_intervals_list = tis.find_all_intervals(working_dir, shush=shush)
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
        print(copyintervals)
        for i in range(0, len(copyintervals)):
            time=copyintervals[i]
            if time:
                r = region_dirs[i]
                timestring = time[0].strftime('%H-%M-%S')
                stopstring = time[1].strftime('%H-%M-%S')
                timestring=timestring+'_'+stopstring
                regcopy = r+timestring+'/'+timestring+'_aia_region.pickle'
                status = subprocess.call('cp '+file+' '+regcopy, shell=True)



def nu_aia_coalign(time_interval, working_dir, nushift, regionmethod='fit',
                   obsid='', region_dir='',
                   input_aia=[], 
                   grade='0_4', justCOM=False, save_dict=False,
                  savefigdir=[]):
    """
    nushift in x, y arcseconds.
    
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
        m, nu_smap, aiareg, COMxy = nuevtplot(evtA=evtA, evtB=evtB,
              savefigdir=savefigdir, AIA94=True, input_aia=input_aia,
              regiondictA=regiondictA, regiondictB=regiondictB,
             regionsave=False, #regionsavename=regionsavename, 
                             overlimb=True, nushift=nushift) 
    else:
        m, nu_smap, aiareg, COMxy = nuevtplot(evtA=evtA, evtB=evtB,
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



path_to_dodem = '/Users/jmdunca2/do-dem/'
def nuevtplot(evtA=[], evtB=[], datapath='./',
              AIA94=False, nushift=[], input_aia=[],
              savefigdir='./',
              regiondictA=[], regiondictB=[], 
              regionsave=False, regionsavename='region',
              starter_region=path_to_dodem+'starter_region.reg',
             overlimb=False):
    """
    Previously called "orbitplot".
    
    For a given path to obsid (set "datapath"), make all-orbit FPMA, B plots.

    OR, do the same for any specific evt files by selecting (evtA, evtB).

    OR, by setting AIA94 = True, do either of those things but as contours 
        overplotted on (you guessed it) an AIA 94 \AA image from the midpoint 
        of the observed time. 

    Set overlimb=True to handle the case where the NuSTAR data extends off the solar limb.
    We have to plot differently in such cases, as CompositeMap does not handle this. 
        
    """

    fig = plt.figure(figsize=(14,22))


    evt_data, hdr = nu.return_submap(datapath=datapath, fpm='A', specific_evt=evtA, return_evt_hdr=True)
    time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
    midtime = time0 + (time1-time0).to(u.s).value/2*u.s

    reffile = evtA
    obsid = reffile.split('/')[-1][2:13]
    #print(obsid)

    regionsavename = regionsavename+'_'+obsid

    from scipy import ndimage
    nustar_map_for_com = nustar.map.make_sunpy(evt_data, hdr)
    #Take center of mass and get it into world coordinates
    com = ndimage.measurements.center_of_mass(nustar_map_for_com.data)
    #print(np.sum(nustar_map_for_com.data))
    com_world = nustar_map_for_com.pixel_to_world(com[1]*u.pix, com[0]*u.pix)
    COMxy=[com_world.Tx, com_world.Ty]
    #print(com_world.Tx)
    #print(type(com_world))

    aia_regiondict=[]
    

    if regionsave:
        regionsavetime=midtime

    timestring = time0.strftime('%H-%M-%S')
    stopstring = time1.strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring

    if AIA94:
        import sunpy.map

        if input_aia:
            m = input_aia
        else:
            m = aia_search_fido(midtime)

        fpms = ['A','B']
        regiondicts = [regiondictA, regiondictB]
        #print(regiondicts)
        evts=[evtA, evtB]
        for i in [0,1]:
            fpm=fpms[i]
            regiondict=regiondicts[i]
            evt=evts[i]

            # if regiondict:
            #     input_region='circle'
            #     input_aia_region_dicts=[]
            #     for r in regiondict:
            #         rx, ry, rr = r['centerx'], r['centery'], r['radius']
            #         input_aia_region_dicts.append({'center': (rx.value,  ry.value)*u.arcsec,
            #                                           'radius': rr})
                

            #From COM, above.
            xx = COMxy[0].value
            yy = COMxy[1].value
            
            #Set broad box for plotting (using region object)
            bl=[(xx-600)*u.arcsec, (yy-600)*u.arcsec]
            tr=[(xx+800)*u.arcsec,(yy+800)*u.arcsec]
            #print(tr[0]-bl[0], tr[1]-bl[1])

  
            bottom_left = SkyCoord(bl[0]-100*u.arcsec, bl[1]-100*u.arcsec, frame=m.coordinate_frame)
            top_right = SkyCoord(tr[0]+100*u.arcsec,tr[1]+100*u.arcsec, frame=m.coordinate_frame)
            mm = m.submap(bottom_left=bottom_left, top_right=top_right)
    
            #Load in NuSTAR map, specifically submaps covering [-1250,1250] arcseconds in each dimension.
            nu_smap = nu.return_submap(datapath=datapath, fpm=fpm, specific_evt=evt, bl=bl, tr=tr)

            
            if overlimb:
                #Need to reproject a specific way if the NuSTAR FOV extends over the limb.
                
                #Make a new observer object, using the observer from the AIA map
                new_observer = mm.observer_coordinate
                
                #Note: for some reason, using the same data shape with a reprojection can cause most of the data to 
                #be cropped out and set to NaN values.
                #To be safe, we will use a VERY large output shape.
                #out_shape = nu_smap.data.shape
                out_shape = (1500,1500)
                
                out_ref_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime=new_observer.obstime,
                                         frame='helioprojective', observer=new_observer,
                                         rsun=nu_smap.coordinate_frame.rsun)
                
                out_header = sunpy.map.make_fitswcs_header(
                    out_shape,
                    out_ref_coord,
                    scale=u.Quantity(nu_smap.scale)
                    )

                with Helioprojective.assume_spherical_screen(nu_smap.observer_coordinate):
                    nu_reproject = nu_smap.reproject_to(out_header)

            else:
                #Combining into composite map. Change 'levels' keyword to change levels of NuSTAR contours
                comp_map = sunpy.map.Map(mm, nu_smap, composite=True)
                comp_map.set_levels(index=1, levels=[5,10,30,50,70,90], percent=True)
                #comp_map.peek()

            #return m, nu_smap

            #Making plot limits
            world_coords = SkyCoord(Tx=[bl[0], tr[0]], Ty=[bl[1],tr[1]], frame=mm.coordinate_frame)
            apixel_coords_x, apixel_coords_y = mm.wcs.world_to_pixel(world_coords)

            if nushift:      
                from scipy.ndimage.interpolation import shift
                
                #Note: axes have same scale, used axis1 for both coordinates
                xshift=nushift[0]/nu_smap.scale.axis1.value
                yshift=nushift[1]/nu_smap.scale.axis1.value
                
                #Making shifted NuSTAR submap
                shifted_nu_smap = shift(nu_smap.data, [yshift, xshift], mode='constant')
                shifted_nu_smap=sunpy.map.Map(shifted_nu_smap, nu_smap.meta)

                ax = fig.add_subplot(3,2,(i+3), projection=mm)

                if overlimb:
                    mm.plot(axes=ax, zorder=0)
                    
                    with Helioprojective.assume_spherical_screen(shifted_nu_smap.observer_coordinate):
                        nu_reproject_shift = shifted_nu_smap.reproject_to(out_header)                    
                    levels = np.array([1, 5, 10, 30, 50, 70, 90, 95])*u.percent 
                    nu_reproject_shift.draw_contours(levels, axes=ax, alpha=1, zorder=1)
    
                else:
                    shift_comp_map = sunpy.map.Map(mm, shifted_nu_smap, composite=True)
                    shift_comp_map.set_levels(index=1, levels=[5, 10,30,50,70,90], percent=True)
                
                    shift_comp_map.plot()
                    shift_comp_map.draw_limb()
                    
                ax.set_xlim(apixel_coords_x)
                ax.set_ylim(apixel_coords_y)    
                ax.text(0.1,0.9, 'FPM'+fpm+'  shift:'+str(nushift)+' arcsec', color='white',transform=ax.transAxes, fontsize=20)
                bloc = evt.find(fpm, -20)
                ax.text(0.1, 0.85, evt[bloc:-4], color='pink',transform=ax.transAxes, fontsize=20)
                ax.text(0.1, 0.75, 'NuSTAR shift to match AIA.', color='blue',transform=ax.transAxes, fontsize=20)

                #magixs2 orbit 6
                #rectangle = plot_rectangle(525, -275, 800, 800, -60, mm)
                #magixs2 orbit 1
                #rectangle = plot_rectangle(1175, -150, 800, 800, -60, mm)
                #mm.draw_quadrangle(rectangle, axes=ax, edgecolor="red", linestyle="--", linewidth=2)

                ax = fig.add_subplot(3,2,(i+5), projection=mm)

                if overlimb:                
                    mm.plot(axes=ax, zorder=0)                    
                    levels = np.array([1, 5, 10, 30, 50, 70, 90, 95])*u.percent 
                    nu_reproject.draw_contours(levels, axes=ax, alpha=1, zorder=1)

                else:                
                    comp_map.plot()
                    comp_map.draw_limb()
                    
                ax.set_xlim(apixel_coords_x)
                ax.set_ylim(apixel_coords_y)    
                ax.text(0.1,0.9, 'FPM'+fpm, color='Red',transform=ax.transAxes, fontsize=20)
                bloc = evt.find(fpm, -20)
                ax.text(0.1, 0.85, evt[bloc:-4], color='pink',transform=ax.transAxes, fontsize=20)
                #ax.text(0.1, 0.75, 'Both regions, for analysis.', color='blue',transform=ax.transAxes, fontsize=20)


                if regiondict:
                    ax.text(0.1, 0.75, 'AIA+NuSTAR regions, for analysis.', color='blue',transform=ax.transAxes, fontsize=20)
                    num=0
                    aia_regiondicts=[]
                    for r in regiondict:
                        region=make_region(r, mm)
                        og_region = region.to_pixel(mm.wcs)                    
                        og_region.plot(axes=ax, color='blue', ls='--', lw=3, label='Chosen NuSTAR Region')

                        aia_regiondict = r.copy()
                        aia_regiondict['centerx'] = (r['centerx'].value + nushift[0])*u.arcsec
                        aia_regiondict['centery'] = (r['centery'].value + nushift[1])*u.arcsec
                        aia_regiondicts.append(aia_regiondict)

                        region=make_region(aia_regiondict, mm)
                        og_region = region.to_pixel(mm.wcs)                    
                        og_region.plot(axes=ax, color='pink', ls='--', lw=3, label='Chosen AIA Region')

                        if regionsave:
                            if len(regiondict) > 1:
                                rf.write_regfile(starter_region, regionsavetime, region, newfile=regionsavename+fpm+'_AIA_'+num)
                            else:
                                rf.write_regfile(starter_region, regionsavetime, region, newfile=regionsavename+fpm+'_AIA')

                
    
            ax = fig.add_subplot(3,2,(i+1), projection=mm)

            if overlimb:
                mm.plot(axes=ax, zorder=0)                
                levels = np.array([0, 1, 5, 10, 30, 50, 70, 90, 95])*u.percent 
                nu_reproject.draw_contours(levels, axes=ax, alpha=1, zorder=1)

            else:
                comp_map.plot()
                #comp_map.draw_limb()
            
            ax.set_xlim(apixel_coords_x)
            ax.set_ylim(apixel_coords_y)
            ax.text(0.1,0.9, 'FPM'+fpm, color='Red',transform=ax.transAxes, fontsize=20)
            bloc = evt.find(fpm, -20)
            ax.text(0.1, 0.85, evt[bloc:-4], color='pink',transform=ax.transAxes, fontsize=20)

            if regiondict:
                ax.text(0.1, 0.75, 'Initial data + input region.', color='blue',transform=ax.transAxes, fontsize=20)
                
                for r in regiondict:
                    region=make_region(r, mm)
                    og_region = region.to_pixel(mm.wcs)                    
                    og_region.plot(axes=ax, color='blue', ls='--', lw=3, label='Chosen Region')
                
                if regionsave:
                    if len(regiondict) > 1:
                        rf.write_regfile(starter_region, regionsavetime, region, newfile=regionsavename+fpm+'_'+num)
                    else:
                        rf.write_regfile(starter_region, regionsavetime, region, newfile=regionsavename+fpm)
                        
            else:
               ax.text(0.1, 0.75, 'Initial data.', color='blue',transform=ax.transAxes, fontsize=20) 


    else:
        cmap = plt.cm.get_cmap('plasma')

        fpms = ['A','B']
        regiondicts = [regiondictA, regiondictB]
        evts=[evtA, evtB]
        for i in [0,1]:
            fpm=fpms[i]
            regiondict=regiondicts[i]
            evt=evts[i]
            submap = nu.return_submap(datapath=datapath, fpm=fpm, specific_evt=evt)
            norm = colors.Normalize(0, np.max(submap.data))
        
            ax = fig.add_subplot(121, projection=submap)
            submap.plot(axes=ax, norm=norm, cmap=cmap)
            submap.draw_limb()
            ax.text(0.1,0.9, 'FPM'+fpm, color='Red',transform=ax.transAxes, fontsize=30)
            if evt:
                aloc = evt.find(fpm, -20)
                ax.text(0.1, 0.85, evt[aloc:-4], color='pink', transform=ax.transAxes, fontsize=20)
    
            if regiondict:
                for r in regiondict:
                    region=make_region(r, submap)
                    og_region = region.to_pixel(submap.wcs)
                    og_region.plot(axes=ax, color='blue', ls='--', lw=3, label='Chosen Region')
                    if regionsave:
                        if len(regiondict) > 1:
                            rf.write_regfile(starter_region, regionsavetime, region, newfile=regionsavename+fpm+'_'+num)
                        else:
                            rf.write_regfile(starter_region, regionsavetime, region, newfile=regionsavename+fpm)
    


    plt.savefig(savefigdir+'/Both_FPM_images_'+timestring+evt[bloc+3:-4]+'.png')

    if AIA94:
        if nushift:
            return m, nu_smap, aia_regiondicts, COMxy
        else:            
            return m, nu_smap, regiondict, COMxy
        #return m










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

    #Make a new working directory for prepped data/etc if it doesn't yet exist
    aia_dict_dir=working_dir+'all_aia_dicts_'+key+'/'
    save_path = pathlib.Path(aia_dict_dir)
    if not save_path.exists():
        save_path.mkdir()

    #This is used to make an additional directory containing aia region files for all regions
    #in a given observation.
    if additional_path:
        other_aia_dict_dir=additional_path+'all_aia_dicts_'+key+'/'
        #Make a new working directory for prepped data/etc if it doesn't yet exist
        save_path = pathlib.Path(other_aia_dict_dir)
        if not save_path.exists():
            save_path.mkdir()

    
    suborbit_directories = []
    #For each orbit...
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

    #print(time_interval)

    time = time_interval
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    #print(timestring)
    timestring=timestring+'_'+stopstring

    filename = where+'/'+timestring+common_string+'.pickle'
    #print('filename: ', filename)
    
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







# path_to_dodem = '/Users/jmdunca2/do-dem/'
# def nuevtplot(evtA=[], evtB=[], datapath='./',
#               AIA94=False, nushift=[], input_aia=[],
#               savefigdir='./',
#               regiondictA=[], regiondictB=[], 
#               regionsave=False, regionsavename='region',
#               starter_region=path_to_dodem+'starter_region.reg',
#              overlimb=False):
#     """
#     Previously called "orbitplot".
    
#     For a given path to obsid (set "datapath"), make all-orbit FPMA, B plots.

#     OR, do the same for any specific evt files by selecting (evtA, evtB).

#     OR, by setting AIA94 = True, do either of those things but as contours 
#         overplotted on (you guessed it) an AIA 94 \AA image from the midpoint 
#         of the observed time. 

#     Set overlimb=True to handle the case where the NuSTAR data extends off the solar limb.
#     We have to plot differently in such cases, as CompositeMap does not handle this. 
        
#     """

#     fig = plt.figure(figsize=(14,22))


#     evt_data, hdr = nu.return_submap(datapath=datapath, fpm='A', specific_evt=evtA, return_evt_hdr=True)
#     time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
#     midtime = time0 + (time1-time0).to(u.s).value/2*u.s

#     reffile = evtA
#     obsid = reffile.split('/')[-1][2:13]
#     #print(obsid)

#     regionsavename = regionsavename+'_'+obsid

#     if not regiondictA:
#         from scipy import ndimage
#         print('No FPMA region, so we will base the submap box on the FPMA COM')
#         nustar_map_for_com = nustar.map.make_sunpy(evt_data, hdr)
#         #Take center of mass and get it into world coordinates
#         com = ndimage.measurements.center_of_mass(nustar_map_for_com.data)
#         print(np.sum(nustar_map_for_com.data))
#         com_world = nustar_map_for_com.pixel_to_world(com[1]*u.pix, com[0]*u.pix)
#         offset=[com_world.Tx, com_world.Ty]
#         #print(com_world.Tx)
#         #print(type(com_world))

#         aia_regiondict=[]
    

#     if regionsave:
#         regionsavetime=midtime

#     timestring = time0.strftime('%H-%M-%S')
#     stopstring = time1.strftime('%H-%M-%S')
#     timestring=timestring+'_'+stopstring

#     if AIA94:
#         import sunpy.map

#         if input_aia:
#             m = input_aia
#         else:
        
#             from aiapy.calibrate.util import get_correction_table, get_pointing_table
#             from aiapy.calibrate import register, update_pointing, degradation, estimate_error
#             from sunpy.net import Fido
#             from sunpy.net import attrs as a

#             query = Fido.search(
#                     a.Instrument.aia,
#                     a.Physobs.intensity,
#                     a.Wavelength(94*u.angstrom),
#                     a.Time(midtime-12*u.s, midtime))
#                     #a.Time(time_range[0],time_range[1]))#,
#                     #a.Sample(sample_every)
#                     #)
#             print(query)
#             files = Fido.fetch(query, max_conn=1)
#             print(files)
#             try:
#                 amap=sunpy.map.Map(files)
#                 #If things keep failing on the pointing table line, it may be an issue with the actual AIA map – try another time.
#                 ptab = get_pointing_table(amap.date - 12 * u.h, amap.date + 12 * u.h)
#             except AttributeError:
#                 amap=sunpy.map.Map(files[0])
#                 #If things keep failing on the pointing table line, it may be an issue with the actual AIA map – try another time.
#                 ptab = get_pointing_table(amap.date - 12 * u.h, amap.date + 12 * u.h)
                

#             try:
#                 m_temp = update_pointing(amap, pointing_table=ptab)
#             except TypeError:
#                 amap.meta.pop('crpix1')
#                 amap.meta.pop('crpix2')
#                 print('CRPIX issue on ', files)
#                 m_temp = update_pointing(amap, pointing_table=ptab)
                
#             #converts lev1 map to lev1.5 map 
#             #(see: https://aiapy.readthedocs.io/en/stable/api/aiapy.calibrate.register.html?highlight=register)
#             m = register(m_temp)

#         fpms = ['A','B']
#         regiondicts = [regiondictA, regiondictB]
#         evts=[evtA, evtB]
#         for i in [0,1]:
#             fpm=fpms[i]
#             regiondict=regiondicts[i]
#             evt=evts[i]

#             if regiondict:
#                 offset = [regiondict['centerx'], regiondict['centery']]
#                 rad=regiondict['radius']*u.arcsec

        
#             xx = offset[0].value
#             yy = offset[1].value
            
#             #Set broad box for plotting (using region object)
#             bl=[(xx-600)*u.arcsec, (yy-600)*u.arcsec]
#             tr=[(xx+800)*u.arcsec,(yy+800)*u.arcsec]
#             #print(tr[0]-bl[0], tr[1]-bl[1])

#             if regiondict:
#                 input_region='circle'
#                 input_aia_region_dict={'center': (xx,  yy)*u.arcsec,
#                                   'radius': rad}
    
#             bottom_left = SkyCoord(bl[0]-100*u.arcsec, bl[1]-100*u.arcsec, frame=m.coordinate_frame)
#             top_right = SkyCoord(tr[0]+100*u.arcsec,tr[1]+100*u.arcsec, frame=m.coordinate_frame)
#             mm = m.submap(bottom_left=bottom_left, top_right=top_right)
    
#             #Load in NuSTAR map, specifically submaps covering [-1250,1250] arcseconds in each dimension.
#             nu_smap = nu.return_submap(datapath=datapath, fpm=fpm, specific_evt=evt, bl=bl, tr=tr)

            
            
#             if overlimb:
#                 #Need to reproject a specific way if the NuSTAR FOV extends over the limb.
                
#                 #Make a new observer object, using the observer from the AIA map
#                 new_observer = mm.observer_coordinate
                
#                 #Note: for some reason, using the same data shape with a reprojection can cause most of the data to 
#                 #be cropped out and set to NaN values.
#                 #To be safe, we will use a VERY large output shape.
#                 #out_shape = nu_smap.data.shape
#                 out_shape = (1500,1500)
                
#                 out_ref_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime=new_observer.obstime,
#                                          frame='helioprojective', observer=new_observer,
#                                          rsun=nu_smap.coordinate_frame.rsun)
                
#                 out_header = sunpy.map.make_fitswcs_header(
#                     out_shape,
#                     out_ref_coord,
#                     scale=u.Quantity(nu_smap.scale)
#                     )

#                 with Helioprojective.assume_spherical_screen(nu_smap.observer_coordinate):
#                     nu_reproject = nu_smap.reproject_to(out_header)

#             else:
#                 #Combining into composite map. Change 'levels' keyword to change levels of NuSTAR contours
#                 comp_map = sunpy.map.Map(mm, nu_smap, composite=True)
#                 comp_map.set_levels(index=1, levels=[5,10,30,50,70,90], percent=True)
#                 #comp_map.peek()

#             #return m, nu_smap

#             #Making plot limits
#             world_coords = SkyCoord(Tx=[bl[0], tr[0]], Ty=[bl[1],tr[1]], frame=mm.coordinate_frame)
#             apixel_coords_x, apixel_coords_y = mm.wcs.world_to_pixel(world_coords)

#             if nushift:      
#                 from scipy.ndimage.interpolation import shift
                
#                 #Note: axes have same scale, used axis1 for both coordinates
#                 xshift=nushift[0]/nu_smap.scale.axis1.value
#                 yshift=nushift[1]/nu_smap.scale.axis1.value
                
#                 #Making shifted NuSTAR submap
#                 shifted_nu_smap = shift(nu_smap.data, [yshift, xshift], mode='constant')
#                 shifted_nu_smap=sunpy.map.Map(shifted_nu_smap, nu_smap.meta)

#                 ax = fig.add_subplot(3,2,(i+3), projection=mm)

#                 if overlimb:
#                     mm.plot(axes=ax, zorder=0)
                    
#                     with Helioprojective.assume_spherical_screen(shifted_nu_smap.observer_coordinate):
#                         nu_reproject_shift = shifted_nu_smap.reproject_to(out_header)                    
#                     levels = np.array([1, 5, 10, 30, 50, 70, 90, 95])*u.percent 
#                     nu_reproject_shift.draw_contours(levels, axes=ax, alpha=1, zorder=1)
    
#                 else:
#                     shift_comp_map = sunpy.map.Map(mm, shifted_nu_smap, composite=True)
#                     shift_comp_map.set_levels(index=1, levels=[10,30,50,70,90], percent=True)
                
#                     shift_comp_map.plot()
#                     shift_comp_map.draw_limb()
                    
#                 ax.set_xlim(apixel_coords_x)
#                 ax.set_ylim(apixel_coords_y)    
#                 ax.text(0.1,0.9, 'FPM'+fpm+'  shift:'+str(nushift)+' arcsec', color='white',transform=ax.transAxes, fontsize=20)
#                 bloc = evt.find(fpm, -20)
#                 ax.text(0.1, 0.85, evt[bloc:-4], color='pink',transform=ax.transAxes, fontsize=20)
#                 ax.text(0.1, 0.75, 'NuSTAR shift to match AIA.', color='blue',transform=ax.transAxes, fontsize=20)

#                 #magixs2 orbit 6
#                 #rectangle = plot_rectangle(525, -275, 800, 800, -60, mm)
#                 #magixs2 orbit 1
#                 #rectangle = plot_rectangle(1175, -150, 800, 800, -60, mm)
#                 #mm.draw_quadrangle(rectangle, axes=ax, edgecolor="red", linestyle="--", linewidth=2)

#                 ax = fig.add_subplot(3,2,(i+5), projection=mm)

#                 if overlimb:                
#                     mm.plot(axes=ax, zorder=0)                    
#                     levels = np.array([5, 10, 30, 50, 70, 90, 95])*u.percent 
#                     nu_reproject.draw_contours(levels, axes=ax, alpha=1, zorder=1)

#                 else:                
#                     comp_map.plot()
#                     comp_map.draw_limb()
                    
#                 ax.set_xlim(apixel_coords_x)
#                 ax.set_ylim(apixel_coords_y)    
#                 ax.text(0.1,0.9, 'FPM'+fpm, color='Red',transform=ax.transAxes, fontsize=20)
#                 bloc = evt.find(fpm, -20)
#                 ax.text(0.1, 0.85, evt[bloc:-4], color='pink',transform=ax.transAxes, fontsize=20)
#                 #ax.text(0.1, 0.75, 'Both regions, for analysis.', color='blue',transform=ax.transAxes, fontsize=20)


#                 if regiondict:
#                     ax.text(0.1, 0.75, 'Both regions, for analysis.', color='blue',transform=ax.transAxes, fontsize=20)
#                     region=make_region(regiondict, mm)
#                     og_region = region.to_pixel(mm.wcs)                    
#                     og_region.plot(axes=ax, color='blue', ls='--', lw=3, label='Chosen NuSTAR Region')

#                     aia_regiondict = regiondict.copy()
#                     aia_regiondict['centerx'] = (regiondict['centerx'].value + nushift[0])*u.arcsec
#                     aia_regiondict['centery'] = (regiondict['centery'].value + nushift[1])*u.arcsec

#                     region=make_region(aia_regiondict, mm)
#                     og_region = region.to_pixel(mm.wcs)                    
#                     og_region.plot(axes=ax, color='pink', ls='--', lw=3, label='Chosen AIA Region')

#                     if regionsave:
#                         rf.write_regfile(starter_region, regionsavetime, region, newfile=regionsavename+fpm+'_AIA')

                
    
#             ax = fig.add_subplot(3,2,(i+1), projection=mm)

#             if overlimb:
#                 mm.plot(axes=ax, zorder=0)                
#                 levels = np.array([0, 1, 5, 10, 30, 50, 70, 90, 95])*u.percent 
#                 nu_reproject.draw_contours(levels, axes=ax, alpha=1, zorder=1)

#             else:
#                 comp_map.plot()
#                 #comp_map.draw_limb()
            
#             ax.set_xlim(apixel_coords_x)
#             ax.set_ylim(apixel_coords_y)
#             ax.text(0.1,0.9, 'FPM'+fpm, color='Red',transform=ax.transAxes, fontsize=20)
#             bloc = evt.find(fpm, -20)
#             ax.text(0.1, 0.85, evt[bloc:-4], color='pink',transform=ax.transAxes, fontsize=20)

#             if regiondict:
#                 ax.text(0.1, 0.75, 'Initial data + input region.', color='blue',transform=ax.transAxes, fontsize=20)
#                 region=make_region(regiondict, mm)
#                 og_region = region.to_pixel(mm.wcs)
                
#                 og_region.plot(axes=ax, color='blue', ls='--', lw=3, label='Chosen Region')
#                 if regionsave:
#                     rf.write_regfile(starter_region, regionsavetime, region, newfile=regionsavename+fpm)
#             else:
#                ax.text(0.1, 0.75, 'Initial data.', color='blue',transform=ax.transAxes, fontsize=20) 


#     else:

#         #Load in NuSTAR maps, specifically submaps covering [-1250,1250] arcseconds in each dimension.
#         fpm='A'
#         submapA = nu.return_submap(datapath=datapath, fpm=fpm, specific_evt=evtA)
#         fpm='B'
#         submapB = nu.return_submap(datapath=datapath, fpm=fpm, specific_evt=evtB)
        

#         #Plot FPMA
#         cmap = plt.cm.get_cmap('plasma')
#         norm = colors.Normalize(0, np.max(submapA.data))
        
#         ax = fig.add_subplot(121, projection=submapA)
#         submapA.plot(axes=ax, norm=norm, cmap=cmap)
#         submapA.draw_limb()
#         ax.text(0.1,0.9, 'FPMA', color='Red',transform=ax.transAxes, fontsize=30)
#         if evtA:
#             aloc = evtA.find('A', -20)
#             ax.text(0.1, 0.85, evtA[aloc:-4], color='pink',transform=ax.transAxes, fontsize=20)
    
#         if regiondictA:
#             region=make_region(regiondictA, submapA)
#             og_region = region.to_pixel(submapA.wcs)
#             og_region.plot(axes=ax, color='blue', ls='--', lw=3, label='Chosen Region')
#             if regionsave:
#                 rf.write_regfile(starter_region, regionsavetime, region, newfile=regionsavename+'A')
    
        
#         #Plot FPMB
#         cmap = plt.cm.get_cmap('plasma')
#         norm = colors.Normalize(0, np.max(submapB.data))
        
#         ax = fig.add_subplot(122, projection=submapB)
#         submapB.plot(axes=ax, norm=norm, cmap=cmap)
#         submapB.draw_limb()
#         ax.text(0.1,0.9, 'FPMB', color='Red',transform=ax.transAxes, fontsize=30)
#         if evtB:
#             bloc = evtB.find('B', -20)
#             ax.text(0.1, 0.85, evtB[bloc:-4], color='pink',transform=ax.transAxes, fontsize=20)
    
#         if regiondictB:
#             region=make_region(regiondictB, submapB)
#             og_region = region.to_pixel(submapB.wcs)
#             og_region.plot(axes=ax, color='blue', ls='--', lw=3, label='Chosen Region')
#             if regionsave:
#                 rf.write_regfile(starter_region, regionsavetime, region, newfile=regionsavename+'B')


#     plt.savefig(savefigdir+'/Both_FPM_images_'+timestring+evtB[bloc+3:-4]+'.png')

#     if AIA94:
#         if nushift:
#             return m, nu_smap, aia_regiondict
#         else:            
#             return m, nu_smap, regiondict
#         #return m

    