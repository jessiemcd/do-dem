import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy
import nustar_pysolar as nustar

import nustar_utilities as nuutil




def abs_dif_cord(cord):
    return ((cord[0].Tx.value-cord[1].Tx.value)**2 + (cord[0].Ty.value-cord[1].Ty.value)**2)**(1/2)


def nu_fit_gauss(file, twogaussians=True, boxsize=200, plot=False, plotmoments=False,
                guess=[]):
    """Takes in a nustar .evt file and fits one (or two) 2D gaussians 
    to the distribution of data once made into a sunpy map.
    """
    
    #Read in .evt file
    with fits.open(file) as hdu:
        evt_data = hdu[1].data
        hdr = hdu[1].header

    #Make sunpy map
    nustar_map = nustar.map.make_sunpy(evt_data, hdr)
    nudata=nustar_map.data

    #NuSTAR start time, stoptime: 
    time0 = nuutil.convert_nustar_time(hdr['TSTART'])
    time1 = nuutil.convert_nustar_time(hdr['TSTOP'])

    #Take moments of the data to find centroid, window to 2*boxsize square around it
    m = moments(nudata)
    cencoords=[round(m[1]),round(m[2])]  
    coords=cencoords
    #print('centroid data coordinates: ', coords)
    xbox = [coords[1]-boxsize, coords[1]+boxsize]
    ybox = [coords[0]-boxsize, coords[0]+boxsize]    
    #boxdata = nudata[ybox[0]:ybox[1], xbox[0]:xbox[1]]
    boxdata = nudata[ybox[0]:ybox[1], xbox[0]:xbox[1]]

    params = fitgaussian(boxdata, twogaussians=twogaussians, guess=guess)
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

        worldcens=[]
        for cen in [cen1_, cen2_]:

            cen1_world = nustar_map.pixel_to_world(cen[0]*u.pix, cen[1]*u.pix)
            #print(cen1_world)
            if plot:
                ax.plot_coord(coord.SkyCoord(cen1_world.Tx, cen1_world.Ty, frame=nustar_map.coordinate_frame), "o", color='Red',
                         label='Center')
            #print(cen1_world.Tx, cen1_world.Ty)
            
            worldcens.append(cen1_world)
        if plot:    
            ax.set_xlim((cen1_[0]-boxsize), (cen1_[0]+boxsize))
            ax.set_ylim((cen1_[1]-boxsize), (cen1_[1]+boxsize))


        return params, worldcens, nustar_map, time0, time1
    
    return params

#Functions to do 2D fitting of dot locations. 
#First 3 functions adapted from https://scipy-cookbook.readthedocs.io/items/FittingData.html

def fitgaussian(data, twogaussians=False, guess=[]):
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
        boundxwidth=500
        boundywidth=500
        #height, x, y, width_x, width_y, height, x, y, width_x, width_y
        bounds=([0, paramsl[6]-boundxwidth, paramsl[7]-boundywidth, 0, 0,
               0, paramsl[6]-boundxwidth, paramsl[7]-boundywidth, 0, 0], 
                [np.inf, paramsl[6]+boundxwidth, paramsl[7]+boundywidth, np.inf, np.inf,
               np.inf, paramsl[6]+boundxwidth, paramsl[7]+boundywidth, np.inf, np.inf])

        if guess:
            paramsl[6] = guess[0]
            paramsl[7] = guess[1]
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
        errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
        pp = scipy.optimize.leastsq(errorfunction, params, full_output=True)
        p = pp[0]
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
