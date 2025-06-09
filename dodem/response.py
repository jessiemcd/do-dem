import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import shutil

from astropy.io import fits
from dataclasses import dataclass


def update_header_response_files(
    infile: str,
    respfile: str,
    ancrfile: str = 'none',
    outfile: str = None,
    overwrite: bool = False
):
    """
    Updates RESPFILE and ANCRFILE keywords in the header of 'infile' to the
    provided respfile and optional ancrfile keywords. If outfile is provided,
    then infile will be copied to outfile, and the header of outfile will be modified.
    """
    
    respfile = respfile.split('/')[-1]
    ancrfile = ancrfile.split('/')[-1]
    if outfile is None:
        outfile = infile
    else:
        shutil.copyfile(infile, outfile)
    with fits.open(outfile) as hdu:
        hdu[1].header.set('RESPFILE', respfile)
        hdu[1].header.set('ANCRFILE', ancrfile)
        hdu.writeto(outfile, overwrite=overwrite)


def make_srm_file(
    out_file: str,
    rmf_file: str,
    arf_file: str = None,
    data_file: str = None,
    new_data_file: str = None
):
    """
    Generates a file, with a name provided by out_file, for the
    SRM using the RMF and optional ARF files. The created SRM file
    is the same structure as the RMF file.

    The data file corresponds to the provided RMF and ARF files
    (e.g. a PHA file) and will be updated to point to the newly
    generated SRM file if provided. If new_data_file is None,
    then data_file will be overwritten.
    """

    with fits.open(rmf_file) as hdu:
        low_threshold = hdu[2].header['LO_THRES']

    handler = ResponseHandler(rmf_file, arf_file)
    rmf = handler.rmf
    srm = handler.srm # This is the full, nonsparse matrix

    # Turn into sparse matrix using the same valid positions from the RMF.
    # This allows us to reuse the other definitions from the RMF file,
    # e.g. the F_CHAN and N_CHAN columns.
    rows = []
    for row in range(len(srm)):
        goodinds = rmf[row,:].value >= low_threshold
        srm_row = srm[row,goodinds].value
        srm_row = np.array(srm_row, dtype=rmf.dtype)
        rows.append(srm_row)
    srm = np.array(rows, dtype=object)

    # Create the SRM file.
    # TODO: Need to determine whether this conforms to OGIP.
    # I don't think it does since there's no formal definition of an SRM file.
    shutil.copyfile(rmf_file, out_file)
    with fits.open(out_file) as hdu:
        orig_header = hdu[2].header
        orig_data = hdu[2].data
        num_fields = int(orig_header['TFIELDS'])
        for k in orig_header:
            if orig_header[k] == 'MATRIX':
                index = int(k[-1])
                break

        # Copy the columns into the new file.
        cols = []
        for col in range(num_fields):
            if col != (index-1):
                name = orig_header[f'TTYPE{col+1}']
                form = orig_header[f'TFORM{col+1}']
                new_col = fits.Column(name=name, format=form, array=orig_data[name])
                cols.append(new_col)
            else:
                new_col = fits.Column(name='MATRIX', format=orig_header[f'TFORM{index}'], array=srm)
                cols.append(new_col)
        
        hdu[2] = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        hdu[2].header = orig_header
        hdu.writeto(out_file, overwrite=True)

    # Point the data file to the newly created SRM.
    if data_file is not None:
        update_header_response_files(
            data_file,
            respfile=out_file,
            ancrfile='none',
            outfile=new_data_file,
            overwrite=True
        )


@dataclass
class ResponseHandler():
    """
    Reads the provided RMF and ARF files, computes the SRM,
    and acts as a container for the response matrices.
    """
    
    rmf_file: str
    arf_file: str


    @property
    def energy_bins(self) -> np.ndarray:

        with fits.open(self.rmf_file) as hdu:
            matrix_data = hdu[2].data
            matrix_header = hdu[2].header

        elow = matrix_data['ENERG_LO'] * u.Unit(matrix_header['TUNIT1'])
        ehigh = matrix_data['ENERG_HI'] * u.Unit(matrix_header['TUNIT2'])
        energy_bins = np.vstack( [elow, ehigh] ).T

        return energy_bins


    @property
    def energy_edges(self) -> np.ndarray:
        
        bins = self.energy_bins
        edges = np.append(
            bins.flatten()[::2],
            bins[-1,-1] # Include the last edge
        )
        
        return edges


    @property
    def rmf(self) -> np.ndarray:
        return self._read_rmf()
    

    @property
    def arf(self) -> np.ndarray:
        return self._read_arf()
    

    @property
    def srm(self) -> np.ndarray:
        return self._compute_srm()
    

    def plot_matrix(self, which: str, **kwargs):
        """
        which specifies which matrix to plot: "RMF" or "SRM".
        """

        which = which.upper()
        match which:
            case 'RMF':
                matrix = self.rmf
            case 'SRM':
                matrix = self.srm
            case _:
                print(f'Matrix \'{which}\' is not a valid selection.')
                return

        fig, ax = plt.subplots(layout='constrained')
        
        im = ax.matshow(matrix.value, **kwargs)
        ax.set(title=f'{which}')

        return fig, ax, im


    def _read_rmf(self) -> u.Quantity:
        return construct_matrix(self.rmf_file) * u.ct / u.ph


    def _read_arf(self) -> u.Quantity:

        if self.arf_file is not None:
            with fits.open(self.arf_file) as hdu:
                data = hdu[1].data
                hdr = hdu[1].header
            arf = data['SPECRESP'] * u.Unit(hdr['TUNIT3'])
        else:
            arf = None

        return arf
    

    def _compute_srm(self) -> np.ndarray:
        """
        Computes the SRM from the RMF and ARF.
        b_return_edges allows returning the SRM energy edges instead of the
        bins. The returned edge array is the low energy edge of each bin
        and the last edge being the high energy edge of the last bin.
        """

        arf = self.arf
        rmf = self.rmf

        if arf is not None:
            srm = arf[:, None] * rmf
        else:
            srm = rmf # TODO: Is this true?

        return srm


def construct_matrix(file: str) -> np.ndarray:
    """
    Constructs the full matrix from the provided FITS file
    containing the sparse data.
    """

    with fits.open(file) as hdu:
        matrix_data = hdu[2].data

    # From: https://github.com/KriSun95/nustarFittingExample/blob/master/nustarFittingExample/NuSTAR%20Spectrum.ipynb
    fchan_array = col2arr_py(matrix_data['F_CHAN'])
    nchan_array = col2arr_py(matrix_data['N_CHAN'])
    
    matrix = vrmf2arr_py(
        data=matrix_data['MATRIX'],  
        n_grp_list=matrix_data['N_GRP'],
        f_chan_array=fchan_array, 
        n_chan_array=nchan_array
    )

    return matrix


# From: https://github.com/ianan/nustar_sac/blob/master/python/ns_tresp.py
def col2arr_py(data, **kwargs):
    ''' Takes a list of parameters for each energy channel from a .rmf file and returns it in the correct format.
    From: https://lost-contact.mit.edu/afs/physics.wisc.edu/home/craigm/lib/idl/util/vcol2arr.pro
    
    Parameters
    ----------
    data : array/list-like object
            One parameter's array/list from the .rmf file.
    kwargs : idl_check=Bool or idl_way=Bool
            If idl_check=True the funciton will throw an error if the Python and IDL methods give different answers (they shouldn't).
            If idl_way=True the IDL method's result with be returned instead of the new Python method described.
            
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

    ## this is the quicker way I have chosen to do in Python (this may be revised later but is ~30x faster than way below in Python)
    max_len = np.max([len(r) for r in data]) # find max row length
    chan_array_py = np.array([[*r, *(max_len-len(r))*[0]] for r in data]) # make each row that length (padding with 0)

    #*************************************************************************************************************************************************
    # if you want to involve the IDL way
    # set defaults to help check how it is done in IDL (second dict rewrites the keys of the first)
    defaults = {**{"idl_check":False, "idl_way":False}, **kwargs}

    if defaults["idl_check"] or defaults["idl_way"]:
        ## this is the way IDL does col2arr.pro
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

        if defaults["idl_check"]:
            assert np.array_equal(chan_array_py, chan_array), \
            "The IDL way and the Python way here do not produce the same result. \nPlease check this but trust the IDL way more (set idl_way=True)!"
        if defaults["idl_way"]:
            return chan_array
    #*************************************************************************************************************************************************

    return chan_array_py


def vrmf2arr_py(data=None, n_grp_list=None, f_chan_array=None, n_chan_array=None, **kwargs):
    ''' Takes redistribution parameters for each energy channel from a .rmf file and returns it in the correct format.
    From: https://lost-contact.mit.edu/afs/physics.wisc.edu/home/craigm/lib/idl/spectral/vrmf2arr.pro
    
    Parameters
    ----------
    data : array/list-like object
            Redistribution matrix parameter array/list from the .rmf file. Units are counts per photon.
            Default : None
            
    no_of_channels : int
            Number of entries is the total number of photon channels, the entries themselves show the total number 
            of count channels to which that photon channel contributes.
            Default : None
            
    f_chan_array : numpy.array
            The index of each sub-set channel from each energy bin from the .rmf file run through col2arr_py().
            Default : None
            
    n_chan_array : numpy.array
            The number of sub-set channels in each index for each energy bin from the .rmf file run through col2arr_py().
            Default : None
    kwargs : idl_check=Bool or idl_way=Bool
            If idl_check=True the funciton will throw an error if the Python and IDL methods give different answers (they shouldn't).
            If idl_way=True the IDL method's result with be returned instead of the new Python method described.
            
    Returns
    -------
    A 2D numpy array of the correctly ordered input data with dimensions of energy in the rows and channels in 
    the columns.
    
    Code Example
    -------
    >>> d_rmf = 'directory/'
    >>> f_rmf = 'file.rmf'
    >>> e_lo, e_hi, ngrp, fchan, nchan, matrix = nu_spec.read_rmf(d_rmf+f_rmf)
    >>> fchan_array = nu_spec.col2arr_py(fchan)
    >>> nchan_array = nu_spec.col2arr_py(nchan)
    >>> rmf = nu_spec.vrmf2arr_py(data=matrix,  
                                n_grp_list=ngrp,
                                f_chan_array=fchan_array, 
                                n_chan_array=nchan_array)
    >>> rmf
    array([[0.00033627, 0.0007369 , 0.00113175, ..., 0.        , 0.        , 0.        ],
        [0.00039195, 0.00079259, 0.00138341, ..., 0.        , 0.        , 0.        ],
        [0.00042811, 0.00083381, 0.00157794, ..., 0.        , 0.        , 0.        ],
                                                ...,
        [0.        , 0.        , 0.        , ..., 0.00408081, 0.00409889, 0.00403308],
        [0.        , 0.        , 0.        , ..., 0.00405333, 0.00413722, 0.00413216],
        [0.        , 0.        , 0.        , ..., 0.        , 0.        , 0.        ]])
    ## rows = photon/energy channels, columns = counts channels 
    What's Going On?
    ----------------
    The RMF file has the photon-to-counts conversion information in it. 
    The martix has the photon-to-count conversion value for each count channel (columns) that is involved with theach photon channel (rows). 
            E.g., matrix = [ [a, b, c, d, e, f, ...] , 
                            [        ...          ] , 
                            [        ...          ] , 
                                    ...             ]
    F_chan is the starting index of contiguous counts channels that are involved with the photon channel. 
            E.g., f_chan = [ [0, 5, 0, 0, 0, ...] , 
                            [       ...        ] , 
                            [       ...        ] , 
                                    ...           ] 
                            For the first photon channel, there are rows of counts channels starting at index 0 and 5
    N_chan is the corresponding number of counts channels from each index in the f_chan array.
            E.g., n_chan = [ [2, 3, 0, 0, 0, ...] , 
                            [        ...        ] , 
                            [        ...        ] , 
                                    ...           ]
                            Starting at index 0 for the first photon channel we have the first 2 matrix values, then at index 5 we have the next 3.
                            The total of each row is the same as the n_grp_list and the number of entries in each row of the matrix entry.
    Putting all this together, the rmf matrix is:
            rmf_matrix = [ [a, b, 0, 0, 0, c , d , e, 0 , 0 , ...] ,   #<-- index 0 (f_chan) with 2 entries (n_chan) with photon-to-counts conversion (matrix)
                        [                 ...                   ] , 
                        [                 ...                   ] , 
                                        ...                      ] 
    '''
    
    # this was is about >6x quicker in than the IDL code written in Python
    
    # find the non-zero entries in Nchan, this is the number to counts channels 
    #  in a row that contribute so will have a value if it is useful
    b = np.nonzero(n_chan_array)
    
    # now only want the useful entries from the pre-formatted Nchan and Fchan arrays
    c = f_chan_array[b]
    d = n_chan_array[b]
    
    # to help with indexing, this provides a running sum of the number of counts 
    #  channels that a single photon channel contributes to
    e = np.cumsum(n_chan_array, axis=1)

    # these entries will give the final indices in the row on counts channels
    final_inds = e[b]

    # need to find the starting index so -1, but that means any entry that is 
    #  -1 will be where a zero is needed
    starting_inds = b[1]-1

    # get the  starting indices but the ones that should be 0 are replaced with 
    #  the final on in the list at the minute (-1 in starting_inds)
    start_inds = np.cumsum(n_chan_array, axis=1)[(b[0], starting_inds)] 

    # where starting_inds==-1 that value should be 0, i.e. starting from the first 
    #  value in the rmf matrix
    new_e = np.where(starting_inds!=-1, start_inds, 0)

    # initialise the rmf matrix
    mat_array_py = np.zeros((len(data),len(n_grp_list)))
    
    # now go through row by row (this is the slowest part and needs to be made faster).
    #  Here we go through each photon channel's number of discrete rows of counts channels.
    for r in range(len(c)):
        mat_array_py[b[0][r], c[r]:c[r]+d[r]] = data[b[0][r]][new_e[r]:final_inds[r]]


    #*************************************************************************************************************************************************
    # if you want to involve the IDL way
    # set defaults to help check how it is done in IDL (second dict rewrites the keys of the first)
    defaults = {**{"idl_check":False, "idl_way":False}, **kwargs}

    if defaults["idl_check"] or defaults["idl_way"]:
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

        if defaults["idl_check"]:
            assert np.array_equal(mat_array_py, mat_array), \
            "The IDL way and the Python way here do not produce the same result. \nPlease check this but trust the IDL way more (set idl_way=True)!"
        if defaults["idl_way"]:
            return mat_array
    #*************************************************************************************************************************************************
                    
    return mat_array_py


def _test_ResponseHandler(in_dir: str):
    """
    in_dir contains the rmf and arf files.
    """

    from time_process import time_process

    rmf_fileA = f'{in_dir}/nu20801024001A06_cl_g0_sr.rmf'
    arf_fileA = f'{in_dir}/nu20801024001A06_cl_g0_sr.arf'
    
    handler = ResponseHandler(rmf_fileA, arf_fileA)
    rmf, ebins = time_process(handler.get_rmf)
    arf, ebins = time_process(handler.get_arf)
    srm, ebins = time_process(handler.get_srm)
    print('rmf byte usage:', rmf.nbytes)
    print('arf byte usage:', arf.nbytes)
    print('srm byte usage:', srm.nbytes)
    # Since the RMF and SRM are over 100 megabytes in size,
    # we should probably not store them in memory.
    # Generating the SRM only takes a little over one second.

    fig, ax, im = handler.plot_matrix('RMF')
    fig.colorbar(im, ax=ax)
    plt.show()

    fig, ax, im = handler.plot_matrix('SRM')
    fig.colorbar(im, ax=ax)
    plt.show()