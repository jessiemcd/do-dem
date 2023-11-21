import numpy as np
import pandas as pd
from astropy.time import Time
from sunpy.time import parse_time
from sunpy.coordinates import sun

from skyfield.api import load, EarthSatellite
from astropy import units as u
from astropy import coordinates as coord


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
    
    




