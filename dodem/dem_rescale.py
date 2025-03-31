
import numpy as np


def do_rescale(data, emd_int=True, rscl_factor=[]):

    import copy
    sclf=1E15

    #define the rescaling factor
    mnrat=np.mean(data['dn_in']/data['dn_reg'])
    #or set it to the input value
    if rscl_factor:
        mnrat = rscl_factor

    print('Rescaled by: ', mnrat)
    
    #rescale DEM and edem
    dem=data['DEM']*mnrat
    edem=[ed*mnrat for ed in data['edem']]


    #Make the rmatrix
    temps = data['mts']
    dn_in = np.array(data['dn_in'])
    tresp = data['trmatrix']
    tresp_logt = data['ts_']

    rmatrix, sclf = make_rmatrix(temps, tresp, tresp_logt, dn_in, emd_int=True)

    dn_reg=(rmatrix.T @ dem).squeeze()/sclf

    sze=dn_in.shape
    nf=sze[0]

    chisq=np.sum(((data['dn_in']-dn_reg)/data['edn_in'])**2)/nf

    newdata = copy.deepcopy(data)
    newdata['chisq'] = chisq
    newdata['dn_reg'] = dn_reg
    newdata['DEM'] = dem
    newdata['edem'] = edem

    return newdata, mnrat


def make_rmatrix(temps, tresp, tresp_logt, dn_in, emd_int=True):
    
    #create our bin averages:
    logt=([np.mean([(np.log10(temps[i])),np.log10((temps[i+1]))]) for i in np.arange(0,len(temps)-1)])
    #and widths
    dlogt=(np.log10(temps[1:])-np.log10(temps[:-1]))
    nt=len(dlogt)
    logt=(np.array([np.log10(temps[0])+(dlogt[i]*(float(i)+0.5)) for i in np.arange(nt)]))
    #number of DEM entries
    
    sze=dn_in.shape
    nf=sze[0]
    
    truse=np.zeros([tresp[:,0].shape[0],nf])
    #check the tresp has no elements <0
    #replace any it finds with the mimimum tresp from the same filter
    for i in np.arange(0,nf):
        #keep good TR data
        truse[tresp[:,i] > 0]=tresp[tresp[:,i] > 0]
        #set bad data to the minimum
        truse[tresp[:,i] <= 0,i]=np.min(tresp[tresp[:,i] > 0],axis=0)[i]
    
    
        tr=np.zeros([nt,nf])
        for i in np.arange(nf):
        # Ideally should be interp in log-space, so changed
        # Not as big an issue for purely AIA filters, but more of an issue for steeper X-ray ones
            tr[:,i]=10**np.interp(logt,tresp_logt,np.log10(truse[:,i]))
    
        rmatrix=np.zeros([nt,nf])
        #Put in the 1/K factor (remember doing it in logT not T hence the extra terms)
        dlogTfac=10.0**logt*np.log(10.0**dlogt)
        # Do regularization of EMD or DEM 
        if emd_int:
            l_emd=True
            for i in np.arange(nf):
                rmatrix[:,i]=tr[:,i]
        else:
            for i in np.arange(nf):
                rmatrix[:,i]=tr[:,i]*dlogTfac
                
        #Just scale so not dealing with tiny numbers
        sclf=1E15
        rmatrix=rmatrix*sclf

    return rmatrix, sclf






    