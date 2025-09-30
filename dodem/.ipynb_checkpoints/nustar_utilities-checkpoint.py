import sys
from os.path import *
import os
import numpy as np

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.io import fits
from sunpy.time import parse_time
from sunpy.coordinates import sun

from skyfield.api import load, EarthSatellite
from astropy import units as u
from astropy import coordinates as coord


# Useful functions for making the solar temperature response from NuSTAR data...
# loading in rmf, arf etc
# Mostly adapted from https://github.com/KriSun95/krispy/blob/master/krispy/nu_spec.py

#Also, coordinate conversion routines,
#Adapted from: https://github.com/ianan/nustar_pysolar/blob/master/nustar_pysolar/utils.py


def read_pha(file, return_dat_hdr=False):
    ''' Takes a .pha file and extracts useful information from it.
    
    Parameters
    ----------
    file : Str
            String for the .pha file of the spectrum under investigation.
            
    Returns
    -------
    The energy [kev], counts, livetime [s] and ontime [s] for the observation. 
    '''

    hdul = fits.open(file)
    data = hdul[1].data
    hdr = hdul[0].header
    hdul.close()
    engs=1.6+0.04*data['channel']
    
    if return_dat_hdr:
    	return hdul
    else:
    	return engs, data['counts'], hdr['LIVETIME'],hdr['ONTIME']

def read_arf(file):
    ''' Takes a .arf file and extracts useful information from it.
    
    Parameters
    ----------
    file : Str
            String for the .arf file of the spectrum under investigation.
            
    Returns
    -------
    The low and high boundary of energy bins, and the ancillary response [cm^2] (data['specresp']).  
    '''

    hdul = fits.open(file)
    data = hdul[1].data
    hdul.close()
    
    return data['energ_lo'], data['energ_hi'], data['specresp']


def read_rmf(file):
    ''' Takes a .rmf file and returns the actual rmf matrix (modified from Kris' version).
    
    Parameters
    ----------
    file : Str
            String for the .rmf file of the spectrum under investigation.
            
    Returns
    -------
    The low and high boundary of energy bins (data['energ_lo'], data['energ_hi']), 2D redistribution matrix [counts per photon]. 
    '''

    hdul = fits.open(file)
    data = hdul[2].data
    hdul.close()
    
# Taken from 
# https://github.com/KriSun95/nustarFittingExample/blob/master/nustarFittingExample/NuSTAR%20Spectrum.ipynb
    fchan_array = col2arr_py(data['f_chan'])
    nchan_array = col2arr_py(data['n_chan'])
    
    mat_array = vrmf2arr_py(data=data['matrix'],  
                                n_grp_list=data['n_grp'],
                                f_chan_array=fchan_array, 
                                n_chan_array=nchan_array)
    
    return data['energ_lo'], data['energ_hi'], mat_array


def col2arr_py(data):
    ''' Takes a list of parameters for each energy channel from a .rmf file and returns it in the correct format.

    From: https://lost-contact.mit.edu/afs/physics.wisc.edu/home/craigm/lib/idl/util/vcol2arr.pro
    
    Parameters
    ----------
    data : array/list-like object
            One parameter's array/list from the .rmf file.
            
    Returns
    -------
    A 2D numpy array of the correctly ordered input data.
    
    Example
    -------
    data = FITS_rec([(  1.6 ,   1.64,   1, [0]   , [18]  , [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]), 
                     (  1.64,   1.68,   1, [0]   , [20]  , [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),
                     (  1.68,   1.72,   2, [0,22], [20,1], [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]), 
                     dtype=(numpy.record, [('ENERG_LO', '>f4'), ('ENERG_HI', '>f4'), ('N_GRP', '>i2'), 
                                           ('F_CHAN', '>i4', (2,)), ('N_CHAN', '>i4', (2,)), ('MATRIX', '>i4', (2,))]))
                          
    >>> col2arr_py(data['F_CHAN'])
    array([[  0.,   0.],
           [  0.,   0.],
           [  0.,  22.]])
    ## max row length of 2 so 2 columns, each row is an energy channel. 
    '''
    
    chan = np.array(data)

    nc = np.array([len(n) for n in data]) # number of entries in each row
    accum_nc_almost = [nc[i]+sum(nc[0:i]) for i in range(len(nc))] # running total in each row
    
    # need 0 as start with 0 arrays
    accum_nc = np.array([0] + accum_nc_almost) # this acts as the index as if the array has been unraveled

    ## number of columns is the length of the row with the max number of entries (nc)
    ncol = np.max(nc)
    ## number of rows is just the number of rows chan just has
    nrow = len(chan)

    chan_array = np.zeros(shape=(nrow, ncol))

    for c in range(ncol):
        # indices where the number of entries in the row are greater than the column
        where = (nc > c).nonzero()[0] 

        # cycle through the rows to be filled in:
        ## if this row is one that has more values in it than the current column number then use the appropriate chan 
        ## number else make it zero
        chan_array[:,c] = [chan[n][c] if (n in where) else 0 for n in range(nrow)] 
        
    return chan_array

def vrmf2arr_py(data=None, n_grp_list=None, f_chan_array=None, n_chan_array=None):
    ''' Takes redistribution parameters for each energy channel from a .rmf file and returns it in the correct format.

    From: https://lost-contact.mit.edu/afs/physics.wisc.edu/home/craigm/lib/idl/spectral/vrmf2arr.pro
    
    Parameters
    ----------
    data : array/list-like object
            Redistribution matrix parameter array/list from the .rmf file. Units are counts per photon.
            Default : None
            
    no_of_channels : int
            Number of channels/ energy bins.
            Default : None
            
    f_chan_array : numpy.array
            The index of each sub-set channel from each energy bin from the .rmf file run through col2arr_py().
            Default : None
            
    n_chan_array : numpy.array
            The number of sub-set channels in each index for each energy bin from the .rmf file run through col2arr_py().
            Default : None
            
    Returns
    -------
    A 2D numpy array of the correctly ordered input data with dimensions of energy in the rows and channels in 
    the columns.
    
    Example
    -------
    data = FITS_rec([(  1.6 ,   1.64,   1, [0]   , [18]  , [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]), 
                     (  1.64,   1.68,   1, [0]   , [20]  , [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),
                     (  1.68,   1.72,   2, [0,22], [20,1], [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]), ...,
                     dtype=(numpy.record, [('ENERG_LO', '>f4'), ('ENERG_HI', '>f4'), ('N_GRP', '>i2'), 
                                           ('F_CHAN', '>i4', (2,)), ('N_CHAN', '>i4', (2,)), ('MATRIX', '>i4', (2,))]))
                          
    >>> vrmf2arr_py(data['MATRIX'])
    array([[0.00033627, 0.0007369 , 0.00113175, ..., 0.        , 0.        , 0.        ],
           [0.00039195, 0.00079259, 0.00138341, ..., 0.        , 0.        , 0.        ],
           [0.00042811, 0.00083381, 0.00157794, ..., 0.        , 0.        , 0.        ],
                                                ...,
           [0.        , 0.        , 0.        , ..., 0.00408081, 0.00409889, 0.00403308],
           [0.        , 0.        , 0.        , ..., 0.00405333, 0.00413722, 0.00413216],
           [0.        , 0.        , 0.        , ..., 0.        , 0.        , 0.        ]])
    ## rows = energy, columns = channels 
    '''
    
    # unravel matrix array, can't use numpy.ravel as this has variable length rows
    ## now can index the start of each row with the running total
    unravel_dmat = []
    for n in data:
        for nn in n:
            unravel_dmat.append(nn)

    no_of_channels = len(n_grp_list)

    nrows = len(data)
    ncols = no_of_channels
    nc = np.array([len(n) for n in data])
    accum_nc_almost = [nc[i]+sum(nc[0:i]) for i in range(len(nc))]
    accum_nc = np.array([0] + accum_nc_almost) 
    # sorted wobble of diagonal lines, the indices were off by one left and right
    ## i.e. this is the running index so should start at zero

    mat_array = np.zeros(shape=(nrows, ncols))

    for r in range(nrows):
        if nc[r] > 0:
            # in IDL code the second index is -1 but that's because IDL's index boundaries 
            ## are both inclusive sod rop the -1, i.e. was accum_nc[r+1]-1
            row = unravel_dmat[accum_nc[r]:accum_nc[r+1]] 

            c=0

            # for number of sub-set channels in each energy channel groups
            for ng in range(n_grp_list[r]):
                # want redist. prob. for number of sub-set channels 
                ## if c+m is larger than len(row)-1 then only want what we can get
                wanted_r = [row[int(c+m)] for m in np.arange(n_chan_array[r,ng]) if c+m <= len(row)-1 ]

                # now fill in the entries in mat_array from the starting number of the sub-set channel, 
                ## the fchan_array[r, ng]
                for z,wr in enumerate(wanted_r):
                    mat_array[r, int(f_chan_array[r, ng])+z] = wr

                # move the place that the that the index for row starts from along 
                c = c + n_chan_array[r,ng]

            # if dgrp[r] == 0 then above won't do anything, need this as not to miss out the 0th energy channel
            if n_grp_list[r] == 0:
                wanted_r = [row[int(c+m)] for m in np.arange(n_chan_array[r,0]) if c+m <= len(row)-1 ]
                for z,wr in enumerate(wanted_r):
                    mat_array[r, int(f_chan_array[r, 0])+z] = wr
                    
    return mat_array



def convert_nustar_time(t, leap=5, astropy_time=True, from_datetime=False):
    '''Converts MET seconds to a datetime object, or does the opposite.
    
    (Expanded from nustar_pysolar.utils.convert_nustar_time by Jessie Duncan, 3/8/2023)
    Now it goes both ways.
    Original Found: https://github.com/ianan/nustar_pysolar/blob/master/nustar_pysolar/utils.py
    
    Uses the astropy.time.Time method to do the conversion since you have
    to go to MJD first.
    
    Default is to subtract off 5 leap seconds.
    '''
    mjdref = 55197*u.d
    from astropy.time import Time
    
    if from_datetime == False:
        #Converting from NuSTAR time (MET seconds) to astropy Time object or datetime object
        met = (t - leap)* u.s + mjdref
        #print(met.to(u.d))
        met_datetime = Time(met.to(u.d), format = 'mjd').datetime

        if astropy_time is True:
            met_astro = Time(met_datetime)
            return met_astro
        else:
            return met_datetime
        
    else:
        #Converting back to NuSTAR time (MET seconds)
        if astropy_time is False:
            t = Time(t)

        tls = (t.mjd*u.d - mjdref + leap*u.s).to(u.s)
        tv = round(tls.value, 1)
        return tv
    


def RADEC_to_solar(time0, time1, RA, DEC, print_motion=True):
    '''
    For a certain NuSTAR orbit with a certain pointing (in celestial 
    coordinates), returns equivalent coordinates on the solar disk as 
    x,y offsets from solar center. 
    
    Keywords:
    
    time0 – observation start (TSTART in NuSTAR .evt header)
    time1 – obsertation stop (TSTOP in NuSTAR .evt header)
    RA/DEC – celestial coordinates (RA_OBJ, DEC_OBJ in NuSTAR .evt header)
    
    The aimtime associated with the coordinates is taken to be the center
    of the interval between time0, time1. 
    
    (optional)
    Set print_motion=True to find the exact rate of apparent solar motion in
    each direction (RA/DEC) at the aim time.
    
    IMPORTANT NOTE: 
    
    The apparent motion of the Sun is significant (>100" in RA 
    and ~20" in DEC per hour). Thus, the accuracy of the resulting pointing 
    (solar coordinates) will be sensitive to differences between the exact 
    aimtime used to make the RA/DEC coordinates, and the aimtime used here. 
    
    
    Returns: 
    
    Offset – Coordinates on the solar disk as x,y offsets from solar center 
    (in arcsec).
    
    '''
    midtime = (time1-time0)*0.5+time0
    #print(midtime)
    utc = convert_nustar_time(midtime)
    #print(utc)
    
    #Generate "skyfield objects", indicating the locations of the Earth and Sun in a 
    #format that can be computed at different times.
    #load_path = './'
    #load=Loader(load_path)

    ts = load.timescale()
    planets = load('de436.bsp')
    earth = planets['Earth']
    Sun = planets['Sun']
    
    #Turning a "Time Object" into a time usable by at(), and then finding the position 
    #of the Sun as viewed from the position of the Earth.
    tcheck = ts.from_astropy(utc)
    #print('Tcheck:', tcheck)
    geocentric = earth.at(tcheck).observe(Sun).apparent()
    #RA/DEC of Sun at tcheck
    this_ra_geo, this_dec_geo, dist = geocentric.radec()
    #print('solar RA, DEC:', this_ra_geo.to(u.deg), this_dec_geo.to(u.deg))

    #Apparent motion of the Sun (to know how far off you'll be spatially given a certain 
    #difference in aim time)
#     from skyfield import framelib
#     eeod = framelib.true_equator_and_equinox_of_date
#     dec, ra, distance, dec_rate, ra_rate, range_rate = (
#         geocentric.frame_latlon_and_rates(eeod)
#     )
#     if print_motion==True:
#         print(f'RA:  {ra_rate.arcseconds.per_hour:+.2f} asec/hr')
#         print(f'Dec: {dec_rate.arcseconds.per_hour:+.2f} asec/hr')


    
    #here we can find the offset in RA/DEC between the sun position and the NuSTAR pointing. 
    #Note thatthese offsets are not along the solar NS and EW directions, but rather along 
    #lines of RA/DEC. A coordinate conversion is needed.
    ra_offset = (RA-this_ra_geo._degrees)*u.deg
    dec_offset = (DEC-this_dec_geo._degrees)*u.deg
    
    print('Offsets, input - solar:', ra_offset, dec_offset)
    
    #Angle between solar north and geocentric north, measured eastward from geocentric north
    sun_np=sun.P(utc).cgs
    #print('solar vs. geometric north angle:', sun_np)
    #Rotation matrix to rotate offset into coordinate system aligned w/ solar north
    rotMatrix = np.array([[np.cos(sun_np), np.sin(sun_np)],
                          [-np.sin(sun_np), np.cos(sun_np)]])
    
    offset = np.array([ra_offset.to(u.arcsec).value, dec_offset.to(u.arcsec).value]) * u.arcsec
    
    
    
    print('Offset (not projected):', offset)
    
    # Project the offset onto the Sun
    delta_offset = np.dot(offset, rotMatrix)
    
    # Account for East->West conversion for +X direction in heliophysics coords
    delta_offset = delta_offset*[-1., 1.]
    
    #print(utc, delta_offset)
    
    return delta_offset
    
# def get_skyfield_position(time, offset, load_path=None, parallax_correction=False):
    
#     """Code for converting solar coordinates to astrometric (J200) RA/Dec coordinates.
#     Parameters
#     ----------
#     time: Date that is parsable by sunpy.time.parse_time()
#     i.e.,
#     time='2016-07-26T19:53:15.00'
#     offset: Offset from the center of the Sun. Must have units from astropy:
#     i.e.: offset = np.array([1000, 150]) * u.arcsec
#     load_path (optional): Relative path from currently location to store bsp files
#     parallax_correction: Use the NuSTAR TLE to correct for orbital parallax
#     Returns
#     ----------
#     sky_position: Two-element array giving the [RA, Dec] coordinates of the
#     target location. Note this is given in astrometric (J2000) RA/Dec, which is what
#     we need for the NuSTAR planning system.
#     Notes
#     ----------
#     Syntax:
#     skyfield_position = get_skyfield_position(time, offset)
#     """

#     from astropy.time import Time
# #     Replaced with newer sunpy v1 function
# #     from sunpy import sun
#     from sunpy.coordinates import sun
#     from nustar_pysolar.utils import skyfield_ephem
#     start_date = parse_time(time)
#     utc = Time(start_date)

#     observer, sunephem, ts = skyfield_ephem(load_path=load_path,
#                                         parallax_correction=parallax_correction,
#                                         utc=utc)

#     tcheck = ts.from_astropy(utc)
#     geocentric = observer.at(tcheck).observe(sunephem)
#     this_ra_geo, this_dec_geo, dist = geocentric.radec()


#     # Get the solar north pole angle. cgs --> radians
# #     sun_np = sunpy.sun.solar_north(t=time).cgs
#     #     Update for sunpy v1.0+
#     sun_np=sun.P(time).cgs

#     # Get the center of the Sun, and assign it degrees.
#     # Doing it this was is necessary to do the vector math below.
#     sun_pos = np.array([this_ra_geo.to(u.deg).value, this_dec_geo.to(u.deg).value])*u.deg


#     # Rotation matrix for a counter-clockwise rotation since we're going
#     # back to celestial north from solar north
#     rotMatrix = np.array([[np.cos(sun_np), np.sin(sun_np)],
#                           [-np.sin(sun_np), np.cos(sun_np)]])

#     # Project the offset onto the Sun
#     delta_offset = np.dot(offset, rotMatrix)

#     # Scale to RA based on the declination.
#     delta_offset = delta_offset * np.array([1. / np.cos(sun_pos[1]), 1.])

#     # Account for the fact that +Ra == East and we have defined +X = West
#     delta_offset = delta_offset * [-1.0, 1.0]

#     # Apply the offset and return the sky position.
#     sky_position = sun_pos + delta_offset

#     return sky_position



def get_sky_position(time, offset):
    """Code for converting solar offsets to pointing position.
    Parameters
    ----------
    time: Date that is parsable by sunpy.time.parse_time()
    i.e.,
    time='2016-07-26T19:53:15.00'
    offset: Offset from the center of the Sun. Must have units from astropy:
    i.e.: offset = np.array([1000, 150]) * u.arcsec
    Returns
    ----------
    sky_position: Two-element array giving the [RA, Dec] coordinates of the
    Notes
    ----------
    Syntax:
    sky_position = get_sky_position(time, offset)
    """

    from astropy.coordinates import get_sun
    from astropy.time import Time
#     Replaced with newer sunpy v1 function
#     from sunpy import sun
    from sunpy.coordinates import sun

    # Convert the date into something that's usable by astropy.


    start_date = parse_time(time)
    astro_time = Time(start_date)

    # Use astropy get_sun for Sun sky position.
    # sunpy has a similar function, but it may be giving a different
    # epoch for the RA and dec. We need them in J2000 RA and dec.

    astro_sun_pos = get_sun(astro_time)
    #print('astropy sun position:', astro_sun_pos)

    # Get the solar north pole angle. cgs --> radians
#     Update for sunpy v1.0+
#     sun_np=sun.solar_north(t=time).cgs
    sun_np=sun.P(time).cgs
    
    #print('solar angle:', sun_np)

    # Get the center of the Sun, and assign it degrees.
    # Doing it this was is necessary to do the vector math below.
    sun_pos = np.array([astro_sun_pos.ra.deg, astro_sun_pos.dec.deg])* u.deg
    #print(sun_pos)
    


    # Rotation matrix for a counter-clockwise rotation since we're going
    # back to celestial north from solar north
    rotMatrix = np.array([[np.cos(sun_np), np.sin(sun_np)],
                         [-np.sin(sun_np),  np.cos(sun_np)]])

    # Project the offset onto the Sun
    delta_offset = np.dot(offset, rotMatrix)
    #print('Input Offset:', offset)
    #print('Projected Offset (rotate offset by solar north pole angle:', delta_offset)

    # Scale to RA based on the declination.
    delta_offset2 = delta_offset * np.array([1. / np.cos(sun_pos[1]), 1.])
    #print('Scaled Offset:', delta_offset2)
    #print('Back to projected offset:', delta_offset2 / np.array([1. / np.cos(sun_pos[1]), 1.]))

    # Account for the fact that +Ra == East and we have defined +X = West
    delta_offset = delta_offset2 * [-1.0, 1.0]
    

    # Apply the offset and return the sky position.
    sky_position = sun_pos + delta_offset
    #print('Sun Position, RA/DEC', sun_pos)
    #print('Offset, RA/DEC', delta_offset) 
    #print('sky_position:', sky_position)

    return sky_position



def sky_to_offset(time, RA, DEC):
    
    """
    GOAL: perfect opposite of get_sky_position()
    
    Parameters
    ----------
    time: Date that is parsable by sunpy.time.parse_time()
    i.e.,
    time='2016-07-26T19:53:15.00'
    [RA, Dec] coordinates to be converted into offsets from the center of the Sun. 
                Must be Astropy Quanitities (units u.deg). 

    ----------
    offset: Offset from the center of the Sun. Must have units from astropy:
    i.e.: offset = np.array([1000, 150]) * u.arcsec
    ----------
    Syntax:
    offset = sky_to_offset(time, RA, DEC)
    
    """
    from astropy.coordinates import get_sun
    from astropy.time import Time
#     Replaced with newer sunpy v1 function
#     from sunpy import sun
    from sunpy.coordinates import sun

    
    sky_position = [RA.value, DEC.value]*u.deg
    #print(sky_position)
    
    start_date = parse_time(time)
    astro_time = Time(start_date)

    # Use astropy get_sun for Sun sky position.
    # sunpy has a similar function, but it may be giving a different
    # epoch for the RA and dec. We need them in J2000 RA and dec.

    astro_sun_pos = get_sun(astro_time)
    #print('astropy sun position:', astro_sun_pos)

    # Get the solar north pole angle. cgs --> radians
    #Update for sunpy v1.0+
    #sun_np=sun.solar_north(t=time).cgs
    #Adding negative factor to rotate in the opposite direction
    sun_np=-1*sun.P(time).cgs
    
    #print('solar angle:', sun_np)

    # Get the center of the Sun, and assign it degrees.
    # Doing it this was is necessary to do the vector math below.
    sun_pos = np.array([astro_sun_pos.ra.deg, astro_sun_pos.dec.deg])* u.deg
    #print(sun_pos)
    
    delta_offset = sky_position - sun_pos
    delta_offset = delta_offset.to(u.arcsec)
    
    # Account for the fact that +Ra == East and we have defined +X = West
    delta_offset = delta_offset * [-1.0, 1.0]
    
    #Remove the scaling to RA based on the declination.
    delta_offset = delta_offset / np.array([1. / np.cos(sun_pos[1]), 1.])
    #print(delta_offset)
    
    # Rotation matrix for a counter-clockwise rotation since we're going
    # back to celestial north from solar north
    rotMatrix = np.array([[np.cos(sun_np), np.sin(sun_np)],
                         [-np.sin(sun_np),  np.cos(sun_np)]])
    
    # Rotate onto solar disk
    offset = np.dot(delta_offset, rotMatrix)
    
    return offset
    
    





