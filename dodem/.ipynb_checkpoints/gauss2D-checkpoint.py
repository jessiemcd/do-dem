import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy
import nustar_pysolar as nustar
import glob
from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion
from astropy import coordinates as coord
import copy

import nustar_utilities as nuutil
import nustar_dem_prep as nu




def abs_dif_cord(cord):
    return ((cord[0].Tx.value-cord[1].Tx.value)**2 + (cord[0].Ty.value-cord[1].Ty.value)**2)**(1/2)


def nu_fit_gauss(file, twogaussians=True, boxsize=200, plot=False, plotmoments=False, plotfile='',
                guess=[], guess2=[], plotregion=[], write_input_regions=False,
                 plotgaussregions=False, 
                 write_regions=False, region_dir='./'):
    """Takes in a nustar .evt file and fits one (or two) 2D gaussians 
    to the distribution of data once made into a sunpy map.
    """

    obsid = file.split('/')[-1][2:13]
    fpm = file.split('/')[-1][13]
    #print(fpm)
    
    #Read in .evt file
    with fits.open(file) as hdu:
        evt_data = hdu[1].data
        hdr = hdu[1].header

    #Make sunpy map
    nustar_map = nustar.map.make_sunpy(evt_data, hdr)
    nudata=nustar_map.data

    #for when there are ~issues~
    #print(nudata.shape)
    #plt.imshow(nudata)
    #print(np.min(nudata), np.max(nudata))
    #print(nudata)

    #NuSTAR start time, stoptime: 
    time0 = nuutil.convert_nustar_time(hdr['TSTART'])
    time1 = nuutil.convert_nustar_time(hdr['TSTOP'])

    #Take moments of the data to find centroid, window to 2*boxsize square around it
    m = moments(nudata)
    #print(m)
    cencoords=[round(m[1]),round(m[2])]  
    coords=cencoords
    #print('centroid data coordinates: ', coords)
    xbox = [coords[1]-boxsize, coords[1]+boxsize]
    ybox = [coords[0]-boxsize, coords[0]+boxsize]    
    #boxdata = nudata[ybox[0]:ybox[1], xbox[0]:xbox[1]]
    boxdata = nudata[ybox[0]:ybox[1], xbox[0]:xbox[1]]
    #print(boxdata.shape)

    params = fitgaussian(boxdata, twogaussians=twogaussians, guess=guess, guess2=guess2)
    paramsl = list(params)
    if twogaussians:
        fit = two_gaussians(*params)
    else:
        fit = gaussian(*params)

    if plot:
        fig = plt.figure(figsize=(12,6))
        ax = fig.add_subplot(121)
        plt.imshow(boxdata)
        ss=boxdata.shape
        xs=np.arange(0, ss[0], 1)
        ys=np.arange(0, ss[1], 1)

        if plotmoments:
            #Take moments     
            m = moments(boxdata)
            print(m)
            mmt=gaussian(*m)
            plt.contour(xs, ys, mmt(*np.indices(boxdata.shape)), cmap=plt.cm.Reds)

        plt.contour(xs, ys, fit(*np.indices(boxdata.shape)), cmap=plt.cm.Greens)
        if twogaussians:
            plt.scatter([paramsl[7], paramsl[2]], [paramsl[6], paramsl[1]], color='Red')


    if twogaussians:
        from astropy import coordinates as coord
        
        cen1 = paramsl[1:3]
        cen2 = paramsl[6:8]
        #print(cen1, cen2)

        #x and y flipped vs. the way the fit output has it.
        cen1_ = [xbox[0]+cen1[1], ybox[0]+cen1[0]]
        cen2_ = [xbox[0]+cen2[1], ybox[0]+cen2[0]]

        #cen1_ = cen1
        #cen2_ = cen2

        if plot:
            #fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(122, projection=nustar_map)
            nustar_map.plot(axes=ax)
            nustar_map.draw_limb()
            levels = np.array([5, 10, 30, 50, 70, 90, 95])*u.percent 
            nustar_map.draw_contours(levels, axes=ax, alpha=1, zorder=1)
            if plotregion:
                num=0
                if write_input_regions:
                    inputcens=[]
                for r in plotregion:
                    center = SkyCoord( *(r['centerx'], r['centery'])*u.arcsec, frame=nustar_map.coordinate_frame)
                    region = CircleSkyRegion(
                            center = center,
                            radius = r['radius'] * u.arcsec
                        )
                    og_region = region.to_pixel(nustar_map.wcs)
                    og_region.plot(axes=ax, color='green', ls='--', lw=3)

                    if write_input_regions:
                        midway = time0 + (time1-time0).to(u.s).value/2*u.s
                        write_regfile('starter_region.reg', midway, region, newfile=region_dir+'gauss_cen_'+obsid+'_'+fpm+'_user_input_'+str(num))
                        num+=1
                        inputcens.append(center)

        num=0
        worldcens=[]
        for cen in [cen1_, cen2_]:

            cen1_world = nustar_map.pixel_to_world(cen[0]*u.pix, cen[1]*u.pix)
            print(cen1_world.Tx, cen1_world.Ty)
            print('')
            if plot:
                ax.plot_coord(coord.SkyCoord(cen1_world.Tx, cen1_world.Ty, frame=nustar_map.coordinate_frame), "o", color='Red',
                         label='Center')

                if plotgaussregions:
                    region = CircleSkyRegion(
                            center = cen1_world,
                            radius = 150 * u.arcsec
                        )
                    og_region = region.to_pixel(nustar_map.wcs)
                    og_region.plot(axes=ax, color='red', ls='--', lw=3)

                if write_regions:
                    midway = time0 + (time1-time0).to(u.s).value/2*u.s
                    region = CircleSkyRegion(
                            center = cen1_world,
                            radius = 150 * u.arcsec
                        )
                    
                    write_regfile('starter_region.reg', midway, region, newfile=region_dir+'gauss_cen_'+obsid+'_'+fpm+'_'+str(num))
                    num+=1
                    
            #print(cen1_world.Tx, cen1_world.Ty)
            
            worldcens.append(cen1_world)
            
        if plot:    
            ax.set_xlim((cen1_[0]-boxsize), (cen1_[0]+boxsize))
            ax.set_ylim((cen1_[1]-boxsize), (cen1_[1]+boxsize))

            if plotfile:
                plt.savefig(plotfile)
                plt.close()

        if write_input_regions:
            return params, inputcens, nustar_map, time0, time1
        else:
            return params, worldcens, nustar_map, time0, time1
        
    else:
        from astropy import coordinates as coord
        cen1 = paramsl[1:3]
        #x and y flipped vs. the way the fit output has it.
        cen1_ = [xbox[0]+cen1[1], ybox[0]+cen1[0]]

        if plot:
            #fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(122, projection=nustar_map)
            nustar_map.plot(axes=ax)
            if plotregion:
                for r in plotregion:
                    center = SkyCoord( *(r['centerx'], r['centery'])*u.arcsec, frame=nustar_map.coordinate_frame)
                    region = CircleSkyRegion(
                            center = center,
                            radius = r['radius'] * u.arcsec
                        )
                    og_region = region.to_pixel(nustar_map.wcs)
                    og_region.plot(axes=ax, color='green', ls='--', lw=3)

        cen1_world = nustar_map.pixel_to_world(cen1_[0]*u.pix, cen1_[1]*u.pix)
        if plot:
            ax.plot_coord(coord.SkyCoord(cen1_world.Tx, cen1_world.Ty, frame=nustar_map.coordinate_frame), "o", color='Red',
                     label='Center')
            print(cen1_world.Tx, cen1_world.Ty)
            print('')
            if plotgaussregions:
                region = CircleSkyRegion(
                        center = cen1_world,
                        radius = 150 * u.arcsec
                    )
                og_region = region.to_pixel(nustar_map.wcs)
                og_region.plot(axes=ax, color='red', ls='--', lw=3)

            
            ax.set_xlim((cen1_[0]-boxsize), (cen1_[0]+boxsize))
            ax.set_ylim((cen1_[1]-boxsize), (cen1_[1]+boxsize))
            if plotfile:
                plt.savefig(plotfile)
                plt.close()
    
    return params, cen1_world, nustar_map, time0, time1


def write_regfile(regfile, time, region, newfile='sample'):
    
    """
    
    Read in a region file + change the region specified.
    
    Time (astropy.time.Time object) for coordinate conversion needed.
    
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
    #print(regstring)
    cs = regstring.split('\n')[-2]
    cs = cs.split('(')[-1]
    cs = cs.split(')')[0]
    cs = cs.split(',')
    #print(cs)
    
    newcs = copy.deepcopy(cs)
    newcs[2]=str(region.radius.value)+'"'
    
    #print([region.center.Tx.value, region.center.Ty.value]*u.arcsec)

    #Get RA, DEC from region in heliocentric coordinates
    RA, DEC = nuutil.get_sky_position(time, [region.center.Tx.value, region.center.Ty.value]*u.arcsec)
    #print(RA,DEC)
    RA = coord.Angle(RA, u.deg)
    DEC = coord.Angle(DEC, u.deg)
    #print(DEC)
    newcs[0] = RA.to_string(unit=u.hour, sep=':')[0:-4]
    decstring=DEC.to_string(unit=u.deg, sep=':')[0:-5]
    if decstring[0] == '-':
        newcs[1] = DEC.to_string(unit=u.deg, sep=':')[0:-5]
    else:
        newcs[1] = '+'+DEC.to_string(unit=u.deg, sep=':')[0:-5]
        
    #print(newcs[1])
    

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







#Functions to do 2D fitting of dot locations. 
#First 3 functions adapted from https://scipy-cookbook.readthedocs.io/items/FittingData.html

def fitgaussian(data, twogaussians=False, guess=[], guess2=[]):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit

    twogaussians: if True, fit two 2D gaussians to the data array. If False, just one. 
    guess: set to a set of x,y coordinates that's a guess for the second gaussian center 
            (only does anything with twogaussians==True). 
    
    """
    
    if twogaussians:
        params = moments(data)*2
        paramsl=list(params)

        #bounds on the two gaussians:
        #–heights must be greater than zero (no negative gaussian components)
        #-widths must be greater than zero
        #-x,y must be within the following width from the centroid as found with the other moments
        #(for both gaussians)
        boundxwidth=150
        boundywidth=150
        #height, x, y, width_x, width_y, height, x, y, width_x, width_y
        bounds=([0, paramsl[6]-boundxwidth, paramsl[7]-boundywidth, 0, 0,
               0, paramsl[6]-boundxwidth, paramsl[7]-boundywidth, 0, 0], 
                [np.inf, paramsl[6]+boundxwidth, paramsl[7]+boundywidth, np.inf, np.inf,
               np.inf, paramsl[6]+boundxwidth, paramsl[7]+boundywidth, np.inf, np.inf])

        if guess:
            paramsl[6] = guess[0]
            paramsl[7] = guess[1]

        if guess2:
            paramsl[1] = guess2[0]
            paramsl[2] = guess2[1]
        
        params = tuple(paramsl)
            
            
        #errorfunction: generates a gaussian using a set of parameters (height, x, y, width_x, width_y, see gaussian)
        #on the same grid as the input data, then subtracts the data. 
        errorfunction = lambda p: np.ravel(two_gaussians(*p)(*np.indices(data.shape)) -
                                 data)
        pp = scipy.optimize.least_squares(errorfunction, params, bounds=bounds)
        if not pp.success:
            print(pp)
        p = pp.x
        
    else:
        params = moments(data)
        paramsl=list(params)
        #bounds on the gaussian:
        #–height must be greater than zero (no negative gaussian components)
        #-width must be greater than zero
        #-x,y must be within the following width from the centroid
        boundxwidth=150
        boundywidth=150
        #height, x, y, width_x, width_y, height, x, y, width_x, width_y
        bounds=([0, paramsl[1]-boundxwidth, paramsl[2]-boundywidth, 0, 0], 
                [np.inf, paramsl[1]+boundxwidth, paramsl[2]+boundywidth, np.inf, np.inf])
        if guess:
            paramsl[1] = guess[0]
            paramsl[2] = guess[1]

        params = tuple(paramsl)
        
        errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
        pp = scipy.optimize.least_squares(errorfunction, params, bounds=bounds)
        if not pp.success:
            print(pp)
        p = pp.x
        
        #pp = scipy.optimize.leastsq(errorfunction, params, full_output=True)
        #p = pp[0]
    #The first entry is the optimized list of parameters, the second is the covarience matrix (the diagonal of
    #which is the varience values, or the 1-sigma uncertainties of each parameter.)
    #print(pp[0])
    #print(pp)
    #perr = np.sqrt(np.diag(pp[1]))
    #print(perr)

    return p

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-x)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-y)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)


def two_gaussians(height, center_x, center_y, width_x, width_y,
                  height2, center_x2, center_y2, width_x2, width_y2):
    """Returns a sum of two gaussian functions, with the respective given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    width_x2 = float(width_x2)
    width_y2 = float(width_y2)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2) + height2*np.exp(
                -(((center_x2-x)/width_x2)**2+((center_y2-y)/width_y2)**2)/2)


def per_orbit_manual_params(in_dir, guess=[], guess2=[], plot=True, 
                            plotregion=[], write_input_regions=True, 
                              plotgaussregions=True, 
                                region_dir='./'):

    """

    guess parameters are for the two-gaussian fit we'll be doing, as a helpful reference. We are NOT going
    to write the gaussian-fit regions to files (for that functionality, see per_orbit_twogauss_params below).
    
    """

    files = glob.glob(in_dir+'/event_cl/nu*06_cl.evt')
    for f in files:
        nu.convert_wrapper(f)
    sunfiles=glob.glob(in_dir+'/event_cl/nu*06_cl_sunpos.evt')

    percents_=[]
    for s in sunfiles:
        #If you don't like your two gaussians, try again with the "guess" keyword – enter a coordinate around where
        #your missing gaussian should be – in pixel coordindates as shown on the left plot.
        res = nu_fit_gauss(s, twogaussians=True, boxsize=200, plot=plot, plotmoments=False, guess=guess, guess2=guess2,
                          plotregion=plotregion, write_input_regions=write_input_regions, 
                           plotgaussregions=plotgaussregions, 
                          region_dir=region_dir)
        
        
        params, inputcens, nustar_map, time0, time1 = res

        percents=[]
        num=0
        for w in inputcens:
            from regions import CircleSkyRegion
            region = CircleSkyRegion(
                center = w,
                radius = 150*u.arcsec
                )
        
            regdata = get_region_data(nustar_map, region, 0)
            percent = np.sum(regdata)/np.sum(nustar_map.data)
            percents.append(percent)
            print('Percent of data in region '+str(num)+': ', percent)
            num+=1

        percents_.append(percents)

    percents_ = np.array(percents_)
    
    fast_min_factors=[]
    for i in range(0, len(inputcens)):
        pe = np.mean(percents_[:,i])
        fast_min_factors.append(round(1/pe*2))

    #print(fast_min_factors)

    return fast_min_factors




def per_orbit_twogauss_params(in_dir, sep_axis='SN', guess=[], guess2=[], plot=True, plotregion=[],
                             write_input_regions=False, 
                              plotgaussregions=False, write_regions=False, region_dir='./'):

    """
    Did you need a guess for one of the centers to make it work? Set it here (data coordinates - look
    at the left output image). 
    
    Run once, then determine if the regions are better separated along the north-south (NS) 
    or east-west (EW) axes:

    sep_axis = 'EW'
    or
    sep_axis = 'SN'

    If halfway between, either will work.
    
    """

    files = glob.glob(in_dir+'/event_cl/nu*06_cl.evt')
    for f in files:
        nu.convert_wrapper(f)
    sunfiles=glob.glob(in_dir+'/event_cl/nu*06_cl_sunpos.evt')

    seps=[]
    percents=[]
    for s in sunfiles:
        #If you don't like your two gaussians, try again with the "guess" keyword – enter a coordinate around where
        #your missing gaussian should be – in pixel coordindates as shown on the left plot.
        res = nu_fit_gauss(s, twogaussians=True, boxsize=200, plot=plot, plotmoments=False, guess=guess, guess2=guess2,
                          plotregion=plotregion, write_input_regions=write_input_regions, 
                           plotgaussregions=plotgaussregions, write_regions=write_regions,
                          region_dir=region_dir)
        params, worldcens, nustar_map, time0, time1 = res
        separation = abs_dif_cord(worldcens)
        print('Separation between double centers: ', separation)
        seps.append(separation)

        relevant_centers=[]
        percents_=[]
        for w in worldcens:
            if sep_axis=='EW':
                relevant_centers.append(w.Tx.value)
            elif sep_axis=='SN':
                relevant_centers.append(w.Ty.value)
            from regions import CircleSkyRegion
            region = CircleSkyRegion(
                    center = w,
                    radius = 150*u.arcsec
                )
            regdata = get_region_data(nustar_map, region, 0)
            percent = np.sum(regdata)/np.sum(nustar_map.data)
            percents_.append(percent)
            print('Percent of data in region: ', percent)

        #print(relevant_centers)
        #print(np.argsort(relevant_centers))
        percents.append(np.array(percents_)[np.argsort(relevant_centers)])
    
    region_radius = int(np.min(seps)/2)
    print('Region radius (generally non-overlapping): ', region_radius)
    
    percents = np.array(percents)
    percents1 = [percents[0,0],percents[1,0]]
    percents2 = [percents[0,1], percents[1,1]]
    percent_estimate = [np.mean(percents1), np.mean(percents2)]
    fast_min_factors = [round(1/pe*2) for pe in percent_estimate]
    
    return sep_axis, guess, guess2, fast_min_factors


def per_orbit_onegauss_params(in_dir, guess=[], plot=True, plotregion=[], plotgaussregions=False):

    """
    Test the fitting on the whole orbit, and set a guess if needed. 
    
    """

    files = glob.glob(in_dir+'/event_cl/nu*06_cl.evt')
    for f in files:
        nu.convert_wrapper(f)
    sunfiles=glob.glob(in_dir+'/event_cl/nu*06_cl_sunpos.evt')

    percents=[]
    for s in sunfiles:
        #If you don't like your gaussian, try again with the "guess" keyword – enter a coordinate around where
        #your missing gaussian should be – in pixel coorindates as shown on the left plot.
        res = nu_fit_gauss(s, twogaussians=False, boxsize=200, plot=plot, plotmoments=False, guess=guess,
                           plotregion=plotregion, plotgaussregions=plotgaussregions)
        params, worldcen, nustar_map, time0, time1 = res
        
        from regions import CircleSkyRegion
        region = CircleSkyRegion(
                center = worldcen,
                radius = 150*u.arcsec
            )
        
        regdata = get_region_data(nustar_map, region, 0)
        percent = np.sum(regdata)/np.sum(nustar_map.data)
        percents.append(percent)
        print('Percent of data in region: ', percent)

    percent_estimate = np.mean(percents)
    fast_min_factor = round(1/percent_estimate*2)
    
    return guess, fast_min_factor







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








