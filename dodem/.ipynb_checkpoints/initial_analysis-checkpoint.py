"""
Should this be in do-dem ? Maybe not!

Code useful for doing initial analysis of NuSTAR data – making images + spectra. Depends on nustar_dem_prep + others. 
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



def return_submap(datapath='./', fpm='A', specific_evt=[], 
                  bl=[], tr=[], nusmooth=True,
                  return_evt_hdr=False, plot=False):
    """
    wrapper - convert to solar coordinates and make submap for nice plot
    """
    
    #Get evt file
    if specific_evt:
        evt_file=specific_evt
    else:
        evt_file = glob.glob(datapath+'/event_cl/*'+fpm+'06_cl.evt')[0]

    #print(evt_file)
    #Convert evt file to solar coordinates         
    nu.convert_wrapper(evt_file)
    
    #Get solar coordinates file
    if specific_evt:
        sun_file = specific_evt[:-4]+'_sunpos.evt'
    else:
        sun_file = glob.glob(datapath+'/event_cl/*'+fpm+'06_cl_sunpos.evt')[0]
    #print('Using solar coodinates file:', sun_file)
    
    with fits.open(sun_file) as hdu:
        evt_data = hdu[1].data
        hdr = hdu[1].header

    

    if return_evt_hdr:
        return evt_data, hdr

    nustar_map = nustar.map.make_sunpy(evt_data, hdr, norm_map=True)
    
    if nusmooth:
       from scipy import ndimage
       import sunpy.map
       #Smoothing the data; change sigma to smooth more or less
       dd=ndimage.gaussian_filter(nustar_map.data, sigma=2, mode='nearest')
       nustar_map=sunpy.map.Map(dd, nustar_map.meta)
        
    if bl:
        submap = nustar_map.submap(bottom_left=SkyCoord(*bl, frame=nustar_map.coordinate_frame),
                          top_right=SkyCoord(*tr, frame=nustar_map.coordinate_frame))
    else:    
        bl = SkyCoord( *(-1250, -1250)*u.arcsec, frame=nustar_map.coordinate_frame)
        tr = SkyCoord( *(1250, 1250)*u.arcsec, frame=nustar_map.coordinate_frame)
        submap = nustar_map.submap(bottom_left=bl, top_right=tr)

    if plot:
        print(submap.unit)
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(1,1,1, projection=submap)
        submap.plot(axes=ax)
        submap.draw_contours(5, axes=ax)#, index=1, percent=True, fill=True)

        #print(submap.contour(0.05))
        #ax.plot_coord(submap.contour(0.05))
        #print(submap.contour(0.05))

    #print(np.max(submap.data))
    #print(np.min(submap.data))
    #print(np.mean(submap.data))

    return submap

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


    evt_data, hdr = return_submap(datapath=datapath, fpm='A', specific_evt=evtA, return_evt_hdr=True)
    time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
    midtime = time0 + (time1-time0).to(u.s).value/2*u.s

    reffile = evtA
    obsid = reffile.split('/')[-1][2:13]
    #print(obsid)

    regionsavename = regionsavename+'_'+obsid

    if not regiondictA:
        from scipy import ndimage
        print('No FPMA region, so we will base the submap box on the FPMA COM')
        nustar_map_for_com = nustar.map.make_sunpy(evt_data, hdr)
        #Take center of mass and get it into world coordinates
        com = ndimage.measurements.center_of_mass(nustar_map_for_com.data)
        print(np.sum(nustar_map_for_com.data))
        com_world = nustar_map_for_com.pixel_to_world(com[1]*u.pix, com[0]*u.pix)
        offset=[com_world.Tx, com_world.Ty]
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
                #If things keep failing on the pointing table line, it may be an issue with the actual AIA map – try another time.
                ptab = get_pointing_table(amap.date - 12 * u.h, amap.date + 12 * u.h)
            except AttributeError:
                amap=sunpy.map.Map(files[0])
                #If things keep failing on the pointing table line, it may be an issue with the actual AIA map – try another time.
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

        fpms = ['A','B']
        regiondicts = [regiondictA, regiondictB]
        evts=[evtA, evtB]
        for i in [0,1]:
            fpm=fpms[i]
            regiondict=regiondicts[i]
            evt=evts[i]

            if regiondict:
                offset = [regiondict['centerx'], regiondict['centery']]
                rad=regiondict['radius']*u.arcsec

        
            xx = offset[0].value
            yy = offset[1].value
            
            #Set broad box for plotting (using region object)
            bl=[(xx-600)*u.arcsec, (yy-600)*u.arcsec]
            tr=[(xx+800)*u.arcsec,(yy+800)*u.arcsec]
            #print(tr[0]-bl[0], tr[1]-bl[1])

            if regiondict:
                input_region='circle'
                input_aia_region_dict={'center': (xx,  yy)*u.arcsec,
                                  'radius': rad}
    
            bottom_left = SkyCoord(bl[0]-100*u.arcsec, bl[1]-100*u.arcsec, frame=m.coordinate_frame)
            top_right = SkyCoord(tr[0]+100*u.arcsec,tr[1]+100*u.arcsec, frame=m.coordinate_frame)
            mm = m.submap(bottom_left=bottom_left, top_right=top_right)
    
            #Load in NuSTAR map, specifically submaps covering [-1250,1250] arcseconds in each dimension.
            nu_smap = return_submap(datapath=datapath, fpm=fpm, specific_evt=evt, bl=bl, tr=tr)

            
            
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
                    shift_comp_map.set_levels(index=1, levels=[10,30,50,70,90], percent=True)
                
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
                    levels = np.array([5, 10, 30, 50, 70, 90, 95])*u.percent 
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
                    ax.text(0.1, 0.75, 'Both regions, for analysis.', color='blue',transform=ax.transAxes, fontsize=20)
                    region=make_region(regiondict, mm)
                    og_region = region.to_pixel(mm.wcs)                    
                    og_region.plot(axes=ax, color='blue', ls='--', lw=3, label='Chosen NuSTAR Region')

                    aia_regiondict = regiondict.copy()
                    aia_regiondict['centerx'] = (regiondict['centerx'].value + nushift[0])*u.arcsec
                    aia_regiondict['centery'] = (regiondict['centery'].value + nushift[1])*u.arcsec

                    region=make_region(aia_regiondict, mm)
                    og_region = region.to_pixel(mm.wcs)                    
                    og_region.plot(axes=ax, color='pink', ls='--', lw=3, label='Chosen AIA Region')

                    if regionsave:
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
                region=make_region(regiondict, mm)
                og_region = region.to_pixel(mm.wcs)
                
                og_region.plot(axes=ax, color='blue', ls='--', lw=3, label='Chosen Region')
                if regionsave:
                    rf.write_regfile(starter_region, regionsavetime, region, newfile=regionsavename+fpm)
            else:
               ax.text(0.1, 0.75, 'Initial data.', color='blue',transform=ax.transAxes, fontsize=20) 


    else:

        #Load in NuSTAR maps, specifically submaps covering [-1250,1250] arcseconds in each dimension.
        fpm='A'
        submapA = return_submap(datapath=datapath, fpm=fpm, specific_evt=evtA)
        fpm='B'
        submapB = return_submap(datapath=datapath, fpm=fpm, specific_evt=evtB)
        

        #Plot FPMA
        cmap = plt.cm.get_cmap('plasma')
        norm = colors.Normalize(0, np.max(submapA.data))
        
        ax = fig.add_subplot(121, projection=submapA)
        submapA.plot(axes=ax, norm=norm, cmap=cmap)
        submapA.draw_limb()
        ax.text(0.1,0.9, 'FPMA', color='Red',transform=ax.transAxes, fontsize=30)
        if evtA:
            aloc = evtA.find('A', -20)
            ax.text(0.1, 0.85, evtA[aloc:-4], color='pink',transform=ax.transAxes, fontsize=20)
    
        if regiondictA:
            region=make_region(regiondictA, submapA)
            og_region = region.to_pixel(submapA.wcs)
            og_region.plot(axes=ax, color='blue', ls='--', lw=3, label='Chosen Region')
            if regionsave:
                rf.write_regfile(starter_region, regionsavetime, region, newfile=regionsavename+'A')
    
        
        #Plot FPMB
        cmap = plt.cm.get_cmap('plasma')
        norm = colors.Normalize(0, np.max(submapB.data))
        
        ax = fig.add_subplot(122, projection=submapB)
        submapB.plot(axes=ax, norm=norm, cmap=cmap)
        submapB.draw_limb()
        ax.text(0.1,0.9, 'FPMB', color='Red',transform=ax.transAxes, fontsize=30)
        if evtB:
            bloc = evtB.find('B', -20)
            ax.text(0.1, 0.85, evtB[bloc:-4], color='pink',transform=ax.transAxes, fontsize=20)
    
        if regiondictB:
            region=make_region(regiondictB, submapB)
            og_region = region.to_pixel(submapB.wcs)
            og_region.plot(axes=ax, color='blue', ls='--', lw=3, label='Chosen Region')
            if regionsave:
                rf.write_regfile(starter_region, regionsavetime, region, newfile=regionsavename+'B')


    plt.savefig(savefigdir+'/Both_FPM_images_'+timestring+evtB[bloc+3:-4]+'.png')

    if AIA94:
        if nushift:
            return m, nu_smap, aia_regiondict
        else:            
            return m, nu_smap, regiondict
        #return m

    

def plot_grade_spectra(working_dir, timestring, fpm):

    """
    Wrapper to plot spectra.
    """
    
    #unphysical grades
    arf_files_unphys, rmf_files_unphys, pha_files_unphys = nu.find_nuproducts(working_dir, timestring, fpm,
                                                                               special_pha=False, grade='21_24')
    engs,cnts_u,lvtm,ontim=nuutil.read_pha(pha_files_unphys[0]) 

    #Grade 0 files: 

    arf_files, rmf_files, pha_files = nu.find_nuproducts(working_dir, timestring, fpm, grade='0')
    engs,cnts,lvtm,ontim=nuutil.read_pha(pha_files[0])
    
    fig=plt.figure(figsize=(10,5))
    plt.plot(engs, (cnts-0.25*cnts_u), label='corr')
    plt.plot(engs, cnts, label='og grade 0')
    plt.plot(engs, cnts_u, label='unphys')
    plt.yscale('log')
    plt.xlim([0,13])
    plt.legend()
    plt.savefig(working_dir+timestring+'/'+timestring+fpm+'pile_up.png')
    plt.close(fig)


    #Grade 0-4 files: 
    
    arf_files, rmf_files, pha_files = nu.find_nuproducts(working_dir, timestring, fpm, grade='0')
    engs,cnts,lvtm,ontim=nuutil.read_pha(pha_files[0])
    
    fig=plt.figure(figsize=(10,5))
    plt.plot(engs, (cnts-(5./4)*cnts_u), label='corr')
    plt.plot(engs, cnts, label='og grades 0-4')
    plt.plot(engs, cnts_u, label='unphys')
    plt.yscale('log')
    plt.xlim([0,13])
    plt.legend()
    plt.savefig(working_dir+timestring+'/'+timestring+fpm+'_adjacent_pile_up.png')   
    plt.close(fig)

    