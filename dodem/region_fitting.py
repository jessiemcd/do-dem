
import nustar_utilities as nuutil

import numpy as np
import astropy.units as u
import nustar_pysolar as nustar
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import coordinates as coord
from regions import CircleSkyRegion
from scipy import ndimage
from pathlib import Path
import sunpy

import os
import copy

import regions
import sunpy.map
import glob


"""
Various code related to NuSTAR data regions:
-Options for fitting the NuSTAR data regions to avoid the chip gap, where instrument response is less well-known. 
-Some code from Reed Masek (March 2023); Jessie-modified.
-Used in nustar_dem_prep.py

"""

def get_file_region(evt_file, time0, time1, regfile, plotfile=False, regRAunit='hourangle',
                   nofit=False, radius=150, working_dir='./'):
    """
    Takes in a file and auto-generates a region for making spectral data products. Returns the new region file name.
    
    Requires both a NuSTAR event file and a region file – this will be used as a template to
    make a new, ideal region file for the event file. 
    
    Set plotfile=True to plot the original region from the template file (default False).
    
    If plotting from region file, set regRAunit to specify RA units (sometimes given in degrees,
    sometimes in hourangle).
    
    """
    #Find midpoint of time interval
    midway = time0 + (time1-time0).to(u.s).value/2*u.s
    
    #Open NuSTAR .evt file
    with fits.open(evt_file) as hdu:
        evt_data = hdu[1].data
        hdr = hdu[1].header
        
    #Make NuSTAR map and submap
    nustar_map = nustar.map.make_sunpy(evt_data, hdr)
    bl = SkyCoord( *(-1250, -1250)*u.arcsec, frame=nustar_map.coordinate_frame)
    tr = SkyCoord( *(1250, 1250)*u.arcsec, frame=nustar_map.coordinate_frame)
    submap = nustar_map.submap(bottom_left=bl, top_right=tr)
    cmap = plt.cm.get_cmap('plasma')
    norm = colors.Normalize(0, np.max(submap.data))
    
    
    #Take center of mass and get it into world coordinates
    com = ndimage.measurements.center_of_mass(submap.data)
    com_world = submap.pixel_to_world(com[1]*u.pix, com[0]*u.pix)
    
    if nofit:
        
        print("No region fitting: using "+str(radius)+" arcsec circle around COM")
        #Default region: 100", centered at COM
        region = CircleSkyRegion(
                center = com_world,
                radius = radius * u.arcsec
            )
        fig = plt.figure(figsize=(16,10))
        ax = fig.add_subplot(121, projection=submap)
        submap.plot(axes=ax, norm=norm, cmap=cmap)
        
        #Diameter of plot window (pixels) - to make things easier to read
        d=300
        ax.set_xlim(com[1]-d/2, com[1]+d/2)
        ax.set_ylim(com[0]-d/2, com[0]+d/2)
        
        og_region = region.to_pixel(submap.wcs)
        og_region.plot(axes=ax, color='green', ls='--', lw=3, label='Chosen Region')
        
        regdata = get_region_data(nustar_map, region, 0)
        
        percent = np.sum(regdata)/np.sum(nustar_map.data)
        
        newregfile = write_regfile(regfile, midway, region, 
                             Path(evt_file).parent.as_posix()+'/'+Path(evt_file).parts[-1][0:-4]+'_COM_region')
        
    else:
        #Default region: 100", centered at COM
        region = CircleSkyRegion(
                center = com_world,
                radius = 100 * u.arcsec
            )
        fit_region = copy.deepcopy(region) 
        #Set scootch=True to find an optimal region by moving (rather than shrinking) the original region.
        fitted_region, percent = fit_region_within_chipgap(evt_data, hdr, fit_region, scooch=True)

        fig = plt.figure(figsize=(16,10))
        ax = fig.add_subplot(121, projection=submap)

        cmap = plt.cm.get_cmap('plasma')
        norm = colors.Normalize(0, np.max(submap.data))
        submap.plot(axes=ax, norm=norm, cmap=cmap)

        #Diameter of plot window (pixels) - to make things easier to read
        d=300
        ax.set_xlim(com[1]-d/2, com[1]+d/2)
        ax.set_ylim(com[0]-d/2, com[0]+d/2)


        og_region = region.to_pixel(submap.wcs)
        og_region.plot(axes=ax, color='pink', ls='--', lw=3, label='OG Region')

        fitted_pix_region = fitted_region.to_pixel(submap.wcs)
        fitted_pix_region.plot(axes=ax, color='lightgreen', ls='--', lw=3, label='Fitted Region')
        
        newregfile = write_regfile(regfile, midway, fitted_region, 
                             './'+Path(evt_file).parent.as_posix()+'/'+Path(evt_file).parts[-1][0:-4]+'_fit_region')
    
    if plotfile:
        res = check_region(evt_file, time0, time1, regfile=True, file=regfile, regRAunit=regRAunit)
        offset, positions, com = res[1:]
    
        ax.plot_coord(coord.SkyCoord(offset[0], offset[1], frame=submap.coordinate_frame), "+", color='Black',
                         label='Region File')
        ax.plot_coord(coord.SkyCoord(positions[:,0]*u.arcsec, positions[:,1]*u.arcsec, 
                                         frame=submap.coordinate_frame), "o", color='Black', markersize=0.5)
    
    plt.legend()
    
    #Add second plot panel showing detID of each event
    det_map = make_det_map(evt_data, hdr)
    det_submap = det_map.submap(bottom_left=bl, top_right=tr)

    ax = fig.add_subplot(122, projection=det_submap)

    det_submap.plot(axes=ax)
    ax.set_xlim(com[1]-d/2, com[1]+d/2)
    ax.set_ylim(com[0]-d/2, com[0]+d/2)    
    
    if nofit:
        plt.savefig(Path(evt_file).parent.as_posix()+'/'+Path(evt_file).parts[-1][0:-4]+'_COM_region.png')
    else:
        plt.savefig(Path(evt_file).parent.as_posix()+'/'+Path(evt_file).parts[-1][0:-4]+'_fit_region.png')
    
    debug=0
    if debug==1:
        return newregfile, percent, fitted_region
    
    return newregfile, percent 


def write_regfile(regfile, time, region, newfile='sample'):
    
    """
    
    Read in a region file + change the region specified.
    
    Times 0,1 (astropy.time.Time objects) for coordinate conversion needed. (Defining window, we use midtime).
    
    Expects region file made in ds9 GUI, circular region, in fk5 coordinates, like:
        (RA, DEC, RAD) in (hourangle, degrees, arcsec).
        
    Returns name of new region file.
    
    Keywords
    --------
    
    regfile - existing circular region file (to be used as a template for our new one).
    region - expects circular region object
    time - data time interval
    newfile - name of new region file to save

    """

    #Open the old file, put contents into string
    f = open(regfile, "r")
    regstring = f.read()
    cs = regstring.split('\n')[-2]
    cs = cs.split('(')[-1]
    cs = cs.split(')')[0]
    cs = cs.split(',')
    
    newcs = copy.deepcopy(cs)
    newcs[2]=str(region.radius.value)+'"'
    
    #print([region.center.Tx.value, region.center.Ty.value]*u.arcsec)

    #Get RA, DEC from region in heliocentric coordinates
    RA, DEC = nuutil.get_sky_position(time, [region.center.Tx.value, region.center.Ty.value]*u.arcsec)
    #print(RA,DEC)
    RA = coord.Angle(RA, u.deg)
    DEC = coord.Angle(DEC, u.deg)
    newcs[0] = RA.to_string(unit=u.hour, sep=':')[0:-4]
    newcs[1] = '+'+DEC.to_string(unit=u.deg, sep=':')[0:-5]

    #Edit copy of region file contents string
    newcs_string = 'circle('+newcs[0]+','+newcs[1]+','+newcs[2]+')'
    cs = regstring.split('\n')[-2]
    split_text = regstring.split('\n')
    new_split_text = copy.deepcopy(split_text)
    new_split_text[-2] = newcs_string
    new_regstring = '\n'.join(new_split_text)
    #print(new_regstring)
    
    #Open the new region file + write to it
    text_file = open(newfile+".reg", "w")
    n = text_file.write(new_regstring)
    text_file.close()
    
    f = open(newfile+".reg", "r")
    regstring = f.read()

    
    return newfile+".reg"





def read_regfile(regfile, time0, time1, regRAunit):
    """
    
    Read in a circular region file, return offset and radius (in arcsec from solar center).
    
    Times 0,1 (astropy.time.Time objects) for coordinate conversion needed. (Defining window)
    
    Expects region file made in ds9 GUI, circular region, in fk5 coordinates, like:
        (RA, DEC, RAD) in (hourangle, degrees, arcsec).
        
    Returns coordinates of circle center + radius in arcseconds from center of solar disk.

    """
    
    f = open(regfile, "r")
    regstring = f.read()
    #print(regstring)
    cs = regstring.split('\n')[-2]
    cs = cs.split('(')[-1]
    cs = cs.split(')')[0]
    cs = cs.split(',')
    #print(cs)
    
    if regRAunit=='degrees':
        ra = coord.Angle(cs[0], unit=u.deg)
    if regRAunit=='hourangle':
        ra = coord.Angle(cs[0], unit=u.hr)
        
    dec = coord.Angle(cs[1], unit=u.deg)
    c = coord.SkyCoord(ra=ra, dec=dec, frame='fk5')


    #print('RA + DEC:', ra.to(u.deg), dec.to(u.deg))
    rad = float(cs[2][0:-1])*u.arcsec
    #print('RADIUS:', rad)
    
    #time0=nuutil.convert_nustar_time(time0, leap=5, astropy_time=True, from_datetime=True)
    #time1=nuutil.convert_nustar_time(time1, leap=5, astropy_time=True, from_datetime=True)
    #offset = nuutil.RADEC_to_solar(time0, time1, c.ra.deg, c.dec.deg)
    midway = time0 + (time1-time0).to(u.s).value/2*u.s
    offset = nuutil.sky_to_offset(midway, ra.to(u.deg), dec.to(u.deg))
    
    return offset, rad


    
def just_plot_region(nufile, time0, time1, regRAunit='hourangle', file=''):  
    """
    Makes a plot of NuSTAR file data and superimposes a region from a region file.
    
    Returns the NuSTAR map. 
    
    """
    
    with fits.open(nufile) as hdu:
        evt_data = hdu[1].data
        hdr = hdu[1].header

    nustar_map = nustar.map.make_sunpy(evt_data, hdr)

    bl = SkyCoord( *(-1100, -1100)*u.arcsec, frame=nustar_map.coordinate_frame)
    tr = SkyCoord( *(1100, 1100)*u.arcsec, frame=nustar_map.coordinate_frame)
    submap = nustar_map.submap(bottom_left=bl, top_right=tr)
    
    #Take center of mass and get it into world coordinates
    com = ndimage.measurements.center_of_mass(submap.data)

    fig = plt.figure(figsize=(16,10))
    ax = fig.add_subplot(121, projection=submap)

    cmap = plt.cm.get_cmap('plasma')
    norm = colors.Normalize(0, np.max(submap.data))
    submap.plot(axes=ax, norm=norm, cmap=cmap)
    submap.draw_limb()

    #Diameter of plot window (pixels) - to make things easier to read
    d=300
    ax.set_xlim(com[1]-d/2, com[1]+d/2)
    ax.set_ylim(com[0]-d/2, com[0]+d/2)
    
    res = check_region(nufile, time0, time1, regfile=True, file=file, regRAunit=regRAunit, shush=True)
    offset, positions, com = res[1:]
    print(offset)
    print(com)

    ax.plot_coord(coord.SkyCoord(offset[0], offset[1], frame=submap.coordinate_frame), "+", color='Black',
                     label='Region File')
    ax.plot_coord(coord.SkyCoord(positions[:,0]*u.arcsec, positions[:,1]*u.arcsec, 
                                     frame=submap.coordinate_frame), "o", color='Black', markersize=0.5)
    
    return submap
    

def check_region(nufile, time0, time1, regfile=True, file='', regobj=False, region='', 
                 regRAunit='hourangle', shush=True):
    """
    For a given NuSTAR file (in solar coordinates), time interval, and region file, 
    make a nice plot that compares the data to the region in the region file. 
    
    Returns the map, the circular region offset, a set of positions around the region circle, and the NuSTAR COM.
    
    Expects region file made in ds9 GUI, circular region, in fk5 coordinates, like:
        (RA, DEC, RAD) in (hourangle, degrees, arcsec).
        
    """
    if regfile==False and reobj==False:
        print('You need to input EITHER a region file or a region object.')
    
   
    if regfile:
        if shush==False:
            print('Reading region file...')
        offset, rad = read_regfile(file, time0, time1, regRAunit)
        positions = np.array(circle_coords(offset, rad))
        
    if regobj:
        if shush==False:
            print('Using region object...')
        offset = [region.center.Tx, region.center.Ty]
        rad = region.radius
    
    # Load in the evt
    hdulist = fits.open(nufile)
    evtdata=hdulist[1].data
    hdr = hdulist[1].header
    hdulist.close()

    nustar_map = nustar.map.make_sunpy(evtdata, hdr, norm_map=True)
    #Smoothing the data; change sigma to smooth more or less
    dd=ndimage.gaussian_filter(nustar_map.data, sigma=2, mode='nearest')
    nustar_map=sunpy.map.Map(dd, nustar_map.meta)

    com = ndimage.measurements.center_of_mass(nustar_map.data)


    return nustar_map, offset, positions, com


def circle_coords(offset, radius):
    """
    Takes in offset ([x,y]*u.arcsec) and radius (also in arcsec).
    Generates list of points on circle for plotting. 
    
    """

    #The lower this value the higher quality the circle is with more points generated
    stepSize = 0.01

    #Generated vertices
    positions = []
    t = 0
    while t < 2 * np.pi:
        positions.append((radius.value * np.cos(t) + offset[0].value, radius.value * np.sin(t) + offset[1].value))
        t += stepSize

    return positions

  
    

def make_det_map_array(evt_data):
    """
    REED WROTE
    Makes a 2999x2999 array (the size of a Sunpy NuSTAR map)
    tracking the most recent DET ID triggered in each pixel.

    Parameters
    ----------
    evt_data : FITS record
        The photon event list containing the data.
        Must have at least the 'X', 'Y', and 'DET_ID' columns.

    Returns
    -------
    arr : np array
        The array containing the det information.
        It is structured as follows: arr[y,x] = det_id
        for each pixel (x,y) that appears in the grid.
        The null value (no photon event) is -1.
    """

    arr = np.full((2999,2999), -1)

    x, y = evt_data['X'], evt_data['Y']
    dets = evt_data['DET_ID']
    arr[y,x] = dets

    return arr


def make_det_map(evt_data, hdr):
    """
    REED WROTE
    Makes a map consisting of the det information for each pixel.
    Note that if photons from different detectors hit the same pixel,
    the most recent value is what will be contained within the pixel.

    Parameters
    ----------
    evt_data : FITS record
        The photon list containing the det information.
    hdr : FITS header
        The header corresponding to evt_data.

    Returns
    -------
    det_map : Sunpy map
        The map containing the det information for each pixel.
    """

    # Replace the map data with the DET data.
    det_map = nustar.map.make_sunpy(evt_data, hdr)
    det_map.data[:] = make_det_map_array(evt_data)

    return det_map


def get_region_data(map_obj, region, fill_value):
    """
    REED WROTE
    Get the map data contained within the provided region.
    Parameters
    ----------
    map_obj : sunpy.map.Map
        The map containing the region of interest.
    region : regions.SkyRegion
        The bounding region.
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
    reg_mask = (region.to_pixel(map_obj.wcs)).to_mask()
    xmin, xmax = reg_mask.bbox.ixmin, reg_mask.bbox.ixmax
    ymin, ymax = reg_mask.bbox.iymin, reg_mask.bbox.iymax
    
    #print('rmask:', reg_mask.data.shape)
    #print('md:', map_data[ymin:ymax, xmin:xmax].shape)
    
    region_data = np.where(reg_mask.data==1, map_data[ymin:ymax, xmin:xmax], fill_value)

    return region_data


def find_dets_in_region(det_map, region):
    """
    REED WROTE
    Determines which detectors appear within the given region.

    Parameters
    ----------
    det_map : Sunpy map
        The map containing the detector data.
    region : CircleSkyRegion
        The region of interest.

    Returns
    -------
    dets : list
        A list containing the dets that are within the region.
    """

    det_data = get_region_data(det_map, region, fill_value=-1)
    dets = list(np.unique(det_data))
    dets.sort()

    return dets


def fit_region_within_chipgap(evt_data, hdr, region, scooch=False):
    """
    REED WROTE
    Fit the inital region within the chip gap by reducing the radius
    if the region contains events from more than one detector.
    
    JESSIE added version to prioritize moving the region rather than shrinking
    (increase total active region coverage) – set scooch=True to use this option.
    
    """

    det_map = make_det_map(evt_data, hdr)
    dets = find_dets_in_region(det_map, region)

    if scooch==False:
        while len(dets) > 1 and region.radius.value > 25:
            # Decrement by 2.5 arcseconds per iteration since each detector pixel is about 2.5 arcseconds
            region.radius = (region.radius.value - 2.5) * u.arcsec
            dets = find_dets_in_region(det_map, region)
            if -1 in dets:
                dets.remove(-1)
        return region
    
    if scooch==True:
        
        nustar_map = nustar.map.make_sunpy(evt_data, hdr)
        
#         -For a given shift, (x,y), determine whether there is only one detector represented.
#         -List the shift coodinates (x,y) where there is only one detector represented.
#         -Find the shift where the region includes the maximum amount of NuSTAR emission.

        holdvalx=region.center.Tx.value
        holdvaly=region.center.Ty.value
        winners=[]
        for x in range(-70,40):
            xshift=2.5*x
            for y in range(-40,60):
                yshift=2.5*y
                
                #print(holdvalx + xshift, holdvaly + yshift)
                
                region.center = SkyCoord( *((holdvalx + xshift), (holdvaly + yshift))*u.arcsec,
                                     frame=nustar_map.coordinate_frame)
                dets = find_dets_in_region(det_map, region)
                if -1 in dets:
                    dets.remove(-1)
                if len(dets) == 1:
                    #print('hi')
                    real_data = get_region_data(nustar_map, region, fill_value=0)
                    winners.append([x,y,np.sum(real_data)])
                    
        #mags = [(w[0]**2+w[1]**2)**(1/2) for w in winners]
        #print('shift magnitudes: ', mags)
        contents = [w[2] for w in winners]
        #print(contents)
        #print('max contents: ', np.max(contents))
        
        windex = np.where(np.array(contents) == np.max(contents))
        #print(windex)
        ans = np.array(winners)[windex]
        #print(ans)
        #print(winners)
        
        bl = SkyCoord( *(-700, -700)*u.arcsec, frame=nustar_map.coordinate_frame)
        tr = SkyCoord( *(700, 700)*u.arcsec, frame=nustar_map.coordinate_frame)
        submap = nustar_map.submap(bottom_left=bl, top_right=tr)
        
        region.center = SkyCoord( *((holdvalx + 2.5*ans[0]), (holdvaly + 2.5*ans[1]))*u.arcsec,
                                     frame=submap.coordinate_frame)
        
        print('Counts in region:', ans[2])
        print('Counts in map:', np.sum(nustar_map.data))
        print('Ratio:', ans[2]/np.sum(nustar_map.data))
        #print(ans)
        
        return region, ans[2]/np.sum(nustar_map.data)
    
    
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
    
    
def get_region_percentage(time_interval, fpm, nustar_path='./', nofit=False):
    """
    Find percentage of total NuSTAR emission which is contained in the saved optimized region
    for a given time interval + fpm. 
    """
    
    time=time_interval
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring
    
    sun_file = glob.glob(nustar_path+timestring+'/*'+fpm+'06_cl_sunpos.evt')[0]
    if nofit:
        file = glob.glob(nustar_path+timestring+'/*'+fpm+'06_cl_sunpos_COM_region.reg')[0]
    else:
        file = glob.glob(nustar_path+timestring+'/*'+fpm+'06_cl_sunpos_fit_region.reg')[0]
    
    with fits.open(sun_file) as hdu:
        evt_data = hdu[1].data
        hdr = hdu[1].header
        
    map_obj = nustar.map.make_sunpy(evt_data, hdr)
    
    
    res = check_region(sun_file, time[0], time[1], regfile=True, file=file)
    off_set = SkyCoord( *res[1], frame=map_obj.coordinate_frame)
    
    region = CircleSkyRegion(
            center = off_set,
            radius = 100 * u.arcsec
        )
    real_data = get_region_data(map_obj, region, 0)
    
    return np.sum(real_data)/np.sum(map_obj.data)
    